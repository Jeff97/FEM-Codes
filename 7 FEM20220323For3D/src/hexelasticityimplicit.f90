!
! Program for 2D explicit dynamics using linear Quadrilateral elements
! Only linear elasticity is assumed in this implementation
!
! This implementation is for multiple processors
!
! This implementation uses Petsc Vectors
!
! In this implementation, several entries are inserted at a time,
! into PETSc matrices/vectors, as recommended by PETSc.
! This is efficient for large-scale models.
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 10-Jan-2018
! Place : Swansea, UK
!
!
!

      PROGRAM QuadElasticityImplicit

      USE ElementUtilitiesElasticity2D
      USE ElementUtilitiesElasticity3D
      USE WriterVTK
      ! USE Module_SolverPetsc

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
      
!#if defined(PETSC_USE_FORTRAN_MODULES)
      USE petscvec
      USE petscmat
      USE petscksp
      USE petscpc
!#endif
      IMPLICIT NONE

! #include "metis.h"

      ! declare variables

      PetscErrorCode errpetsc;

      DOUBLE PRECISION :: tstart, tend

      LOGICAL :: FILEEXISTS, ISOPEN

      INTEGER :: ndim=3
      INTEGER :: ndof=3
      INTEGER :: npElem=8 ! number of nodes per element
      INTEGER :: io, nn, iter
      DOUBLE PRECISION :: PI=ACOS(-1.0)
      DOUBLE PRECISION :: xc, yc, area, fact, val
      DOUBLE PRECISION :: xNode(8), yNode(8), zNode(8)
      DOUBLE PRECISION :: Klocal(24,24), Flocal(24)
      DOUBLE PRECISION :: dispElem(24), veloElem(24), acceElem(24)
      DOUBLE PRECISION :: elemData(50), timeData(50)


      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DirichletBCs
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ForceBCs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnApplied
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnVTK

      INTEGER, DIMENSION(:), ALLOCATABLE :: DirichletBC_row_nums
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: elemNodeConn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeDofArrayOld
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeDofArrayNew
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemDofArray
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: forAssyMat
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeTypeOld
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeTypeNew
      INTEGER, DIMENSION(:), ALLOCATABLE :: assyForSoln
      INTEGER, DIMENSION(:), ALLOCATABLE :: diag_nnz, offdiag_nnz
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnzVec

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_proc_id

      INTEGER :: forAssyVec(24), elemDofGlobal(24)

      CHARACTER(len=32) :: arg

      Vec  tempVec
      Vec  rhsVec      ! RHS vector
      Vec  reacVec     ! Reaction vector
      Vec  solnVec     ! Solution vector for the matrix system
      ! Vec  solnVecFull ! Solution vector for the mesh
      Mat  matA        ! Linear system matrix
      KSP  ksp         ! Linear solver context
      PC   pc          ! Preconditioner context
      PetscReal rNorm   !norm of solution error
      PetscLogStage stage;

      INTEGER :: stepsCompleted, stepsMax, fileCount
      DOUBLE PRECISION :: timeNow, timeFinal, dt
      DOUBLE PRECISION :: af, am, rhoInf, beta, gamm


      DOUBLE PRECISION :: DTT, IDTT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  disp, dispPrev, dispPrev2, dispCur
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  velo, veloPrev, veloCur
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  acce, accePrev, acceCur

      ! Vec solnTemp
      ! Vec globalM, rhsVec


      PetscInt  nNode  ! number of nodes in the whole mesh
      PetscInt  nElem  ! number of global elements (in the whole model)
      PetscInt  nDBC   ! number of Dirichlet BCs
      PetscInt  nFBC   ! number of force BCs - specified nodal forces

      PetscInt totalDOF  ! total number of global DOF in the whole model

      PetscInt  ee, ii, jj, kk, ind
      PetscInt  count, row, col
      PetscInt  n1, n2, n3, n4, nsize, loadstep, nsteps

      PetscScalar, pointer :: xx_v(:)
      PetscOffset xx_i

      MatInfo info;
      KSPConvergedReason reason;
      PetscInt  its;

      CHARACTER (LEN=100) :: inputFileName
      CHARACTER (LEN=100) :: outFileName
      CHARACTER (LEN=100) :: chartemp
      CHARACTER (LEN=30)  :: string1,string2

      !Set file names
      !The file names are specified as inputs from the command line
      IF( iargc() < 1 ) THEN
        WRITE(*,*) "Number of input files is not sufficient "
        WRITE(*,*) "You must enter the name of the input files"
        STOP "Aborting..."
      END IF

      CALL getarg(1, arg)
      inputFileName = arg
      WRITE(*,*) inputFileName

      ! intialise PETSc environment
      !
      call PetscInitialize(PETSC_NULL_CHARACTER, errpetsc)
      ! CHKERRQ(errpetsc)

      ! tstart = MPI_Wtime()

      !
      ! Read nodal data files
      !

      ! check if the file exists

      write(*,*) " FILEEXISTS = ", FILEEXISTS

      INQUIRE(FILE=inputFileName, EXIST=FILEEXISTS)
      write(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", inputFileName, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and read the contents
      OPEN(1, file=inputFileName,STATUS="OLD",ACTION="READ")

      ! Read number of nodes
      READ(1,*, iostat=io) chartemp, nNode
      write(*,*) chartemp, nNode
      write(*,*) "Number of nodes = ", nNode

      ! Read nodal coordinates
      ALLOCATE(coords(nNode,ndim))
      DO ii=1,nNode
        READ(1,*, iostat=io) nn, coords(ii,1), coords(ii,2), coords(ii,3)
        !WRITE(*,*) nn, coords(ii,1), coords(ii,2)
      END DO 

      ! Initialize NodeType and NodeDofArray arrays
      !
      ALLOCATE(NodeTypeOld(nNode, ndof))
      ALLOCATE(NodeTypeNew(nNode, ndof))
      ALLOCATE(NodeDofArrayOld(nNode, ndof))
      ALLOCATE(NodeDofArrayNew(nNode, ndof))
      ind = nNode*ndof
      ALLOCATE(solnApplied(ind))

      kk=1
      DO ii=1,nNode
        DO jj=1,ndof
          NodeTypeOld(ii,jj)   = 0
          NodeTypeNew(ii,jj)   = 0
          NodeDofArrayOld(ii,jj) = 0
          NodeDofArrayNew(ii,jj) = 0
          solnApplied(kk) = 0.0
          kk = kk+1
        END DO
      END DO


      ! Read number of elements
      READ(1,*, iostat=io) chartemp, nElem
      write(*,*) chartemp, nElem
      write(*,*) "Number of Elements = ", nElem

      ! This program is hardcoded for quadrilateral elements
      npElem = 8
      ALLOCATE(elemNodeConn(nElem,npElem))

      DO ii=1,nElem
        READ(1,*, iostat=io) nn, elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3), elemNodeConn(ii,4), &
        elemNodeConn(ii,5), elemNodeConn(ii,6), elemNodeConn(ii,7), elemNodeConn(ii,8)
        !write(*,*) nn, n1, n2, n3, n4
        !WRITE(*,*) nn, elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3), elemNodeConn(ii,4)
      END DO 

      ! Read Dirichlet BC data
      READ(1,*, iostat=io) chartemp, nDBC
      write(*,*) chartemp, nDBC
      write(*,*) "Number of Dirichlet BCs = ", nDBC

      ALLOCATE(DirichletBCs(nDBC,3))
      DO ii=1,nDBC
        READ(1,*, iostat=io) n1, n2, fact
        DirichletBCs(ii,1) = n1
        DirichletBCs(ii,2) = n2
        DirichletBCs(ii,3) = fact
        ! write(*,*) n1, n2, fact

        NodeTypeOld(n1, n2) = 1
        ind = (n1-1)*ndof + n2
        solnApplied(ind) = fact
      END DO 


      ! Read nodal forces
      READ(1,*, iostat=io) chartemp, nFBC
      write(*,*) chartemp, nFBC
      write(*,*) "Number of nodal forces = ", nFBC

      ALLOCATE(ForceBCs(nFBC,3))
      DO ii=1,nFBC
        READ(1,*, iostat=io) n1, n2, fact

        ForceBCs(ii,1) = n1
        ForceBCs(ii,2) = n2
        ForceBCs(ii,3) = fact

        ! write(*,*) n1, n2, fact
      END DO 

      CLOSE(1)


      WRITE(*,*) " Input files have been read successfully \n\n"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      totalDOF = 1
      DO ii=1,nNode
        DO jj=1,ndof
          IF(NodeTypeOld(ii,jj) == 0) THEN
            NodeDofArrayOld(ii,jj) = totalDOF
            totalDOF = totalDOF + 1
          END IF
        END DO
      END DO

      totalDOF = totalDOF-1


      WRITE(charTemp,*) " Mesh statistics .....\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " nElem   = ", nElem, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " npElem         = ", npElem, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " ndof           = ", ndof, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " Total DOF      = ", totalDOF, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! Partition the mesh. Here METIS is used.
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ind = nNode*ndof
      ALLOCATE( disp(ind) )
      ALLOCATE( velo(ind) )
      ALLOCATE( acce(ind) )
      ALLOCATE( dispPrev(ind) )
      ALLOCATE( dispPrev2(ind) )
      ALLOCATE( veloPrev(ind) )
      ALLOCATE( accePrev(ind) )
      ALLOCATE( dispCur(ind) )
      ALLOCATE( veloCur(ind) )
      ALLOCATE( acceCur(ind) )


      ALLOCATE(elem_proc_id(nElem))
      elem_proc_id = 0

      NodeTypeNew     = NodeTypeOld
      NodeDofArrayNew = NodeDofArrayOld

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !prepare the global matrix pattern
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD, " Preparing matrix pattern \n", errpetsc)

      ! write(*,*) " Creating arrays 1 "

      ! ElemDofArray is used for element matrix/vector assembly
      ! As the starting index in PETSc is ZERO (unlike ONE in Fortran),
      ! ONE is  substracted from each entry of ElemDofArray

      ind = npElem*ndof
      ALLOCATE(ElemDofArray(nElem, ind))

      DO ee=1,nElem
        DO ii=1,npElem
          n1 = ndof*(ii-1)
          n2 = elemNodeConn(ee,ii)

          DO jj=1,ndof
            ElemDofArray(ee, n1+jj) = NodeDofArrayNew(n2,jj) - 1
          END DO
        END DO
        ! IF(this_mpi_proc == 0 ) THEN
          ! write(*,*) ee, elemNodeConn(ee,:), ElemDofArray(ee,:)
        ! END IF
      END DO

      ! write(*,*) " ElemDofArray "
      nsize=npElem*ndof
      ! DO ee=1,nElem
      !   write(*,*)  ElemDofArray(ee,1), ElemDofArray(ee,2), ElemDofArray(ee,3)
      ! END DO

      ! write(*,*) " Creating arrays 2 "
      ALLOCATE(assyForSoln(totalDOF))
      ! write(*,*) " Creating arrays 2 "

      count = 1
      DO ii=1,nNode
        DO jj=1,ndof
          ! write(*,*) ii, jj, NodeDofArray(ii,jj)
          IF(NodeDofArrayNew(ii,jj) /= 0) THEN
            assyForSoln(count) = (ii-1)*ndof + jj;
            count = count + 1
          END IF
        END DO
      END DO

      call PetscPrintf(PETSC_COMM_WORLD, " Preparing matrix pattern DONE \n", errpetsc)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! create PETSc variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " Creating PETSc vectors 1"

      call PetscPrintf(PETSC_COMM_WORLD, " Initialising petsc solver \n", errpetsc)



     call VecCreate(PETSC_COMM_WORLD, solnVec, errpetsc)
     ! call VecCreateSeq(PETSC_COMM_WORLD, totalDOF, solnVec, errpetsc)

     ! write(*,*) " Creating PETSc vectors 2"

     call VecSetSizes(solnVec, PETSC_DECIDE, totalDOF, errpetsc)
     !CHKERRQ(errpetsc)

     ! write(*,*) " Creating PETSc vectors 3"

     call VecSetFromOptions(solnVec, errpetsc)
     !CHKERRQ(errpetsc)
     call VecDuplicate(solnVec, rhsVec, errpetsc)
     !CHKERRQ(errpetsc)

     ! write(*,*) " Creating PETSc vectors 4"

     ind = nNode*ndof
     call VecCreate(PETSC_COMM_WORLD, reacVec, errpetsc)
     !CHKERRQ(errpetsc)
     call VecSetSizes(reacVec, PETSC_DECIDE, ind, errpetsc)
     !CHKERRQ(errpetsc)

     call VecSetFromOptions(reacVec, errpetsc)
     !CHKERRQ(errpetsc)
     ! call VecDuplicate(reacVec, solnVecFull, errpetsc)
     !CHKERRQ(errpetsc)

     write(*,*) " Creating PETSc matrices 1"
     call MatCreate(PETSC_COMM_WORLD, matA, errpetsc)
     !CHKERRQ(errpetsc)

     ! write(*,*) " Creating PETSc matrices 2"

     call MatSetSizes(matA, PETSC_DECIDE, PETSC_DECIDE,  totalDOF, totalDOF, errpetsc)
     !CHKERRQ(errpetsc)

     ! write(*,*) " Creating PETSc matrices 3"

     call MatSetFromOptions(matA, errpetsc)
     !CHKERRQ(errpetsc)

     ! write(*,*) " Creating PETSc matrices 4"

     ALLOCATE(nnzVec(totalDOF))

     DO ii=1,totalDOF
       nnzVec(ii) = 500
     END DO

   !   call MatMPIAIJSetPreallocation(matA,
   !  1  50, PETSC_NULL_CHARACTER,
   !  2   0, PETSC_NULL_CHARACTER, errpetsc); !CHKERRQ(errpetsc);

   !   call MatSeqAIJSetPreallocation(matA, 
   !  1  PETSC_DEFAULT, PETSC_NULL_CHARACTER, errpetsc); !CHKERRQ(errpetsc);

     ind=50
     call MatSeqAIJSetPreallocation(matA, ind, nnzVec, errpetsc) !CHKERRQ(errpetsc)

   !   call MatSeqAIJSetPreallocation(matA, 
   !  1  ind, PETSC_NULL_CHARACTER, errpetsc) !CHKERRQ(errpetsc)

     ! call MatCreateSeqAIJ(MPI_COMM_WORLD,
   !  1  totalDOF, totalDOF, 5, PETSC_NULL_CHARACTER, matA, errpetsc)

     call MatSetOption(matA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, errpetsc)

     call MatSetOption(matA, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, errpetsc)

     write(*,*) " Initialise the Matrix pattern "

     nsize = npElem*ndof
     fact = 0.0
     ! Initialise the Matrix pattern
     ElemLoop: DO ee=1,nElem
       LoopII: DO ii=1,nsize
         row = ElemDofArray(ee,ii)
         IF( row /= -1) THEN
           LoopJJ: DO jj=1,nsize
             col = ElemDofArray(ee,jj)
             IF( col /= -1) THEN
               ! write(*,*) ee, ii, jj, row, col
               call MatSetValue(matA, row, col, fact, INSERT_VALUES, errpetsc)
             END IF
           END DO LoopJJ
         END IF
       END DO LoopII
     END DO ElemLoop

   !   ! call MatGetInfo(matA, MAT_LOCAL, &info);

   !   ! pcout << " mallocs      = " <<  info.mallocs << endl;
   !   ! pcout << " nz_allocated = " <<  info.nz_allocated << endl;
   !   ! pcout << " nz_used      = " <<  info.nz_used << endl;
   !   ! pcout << " nz_unneeded  = " <<  info.nz_unneeded << endl;


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! compute element matrices and vectors and assemble them
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " Matrix pattern initialised. "
      write(*,*) " Computing element matrices now. "

      call MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY, errpetsc)
      !CHKERRQ(errpetsc)
      call MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY, errpetsc)
      !CHKERRQ(errpetsc)
      call MatZeroEntries(matA, errpetsc)
      !CHKERRQ(errpetsc)

      call VecAssemblyBegin(rhsVec, errpetsc)
      call VecAssemblyEnd(rhsVec, errpetsc)
      call VecZeroEntries(rhsVec, errpetsc)

      call VecAssemblyBegin(reacVec, errpetsc)
      call VecAssemblyEnd(reacVec, errpetsc);
      call VecZeroEntries(reacVec, errpetsc)


      disp = 0.0; dispPrev = 0.0; dispPrev2 = 0.0; dispCur = 0.0
      velo = 0.0; veloPrev = 0.0; veloCur = 0.0
      acce = 0.0; accePrev = 0.0; acceCur = 0.0

      write(*,*) " Generating element matrices and vectors "

      ! time integration parameters
      timeData(2) = 1.0;   timeData(3) = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      af = 1.0;
      am = 1.0;
      stepsMax = 20000
      stepsCompleted=0
      dt = 0.00001
      timeNow=0.0
      timeFinal = stepsMax*dt

      DTT  = dt*dt
      IDTT = 1.0/DTT

      ! Young's modulus and Poisson's ratio
      elemData(1) = 240.565;   elemData(2) = 0.3
      !elemData(1) = 10.0**9;   elemData(2) = 0.0

      ! Thickness
      elemData(3) = 1.0;
      ! Density
      elemData(4) = 0.0;
      ! Body force in X-, Y- and Z- direction
      elemData(5) = 0.0;   elemData(6) = 0.0; elemData(7) = 0.0


      ! call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      call PetscPrintf(PETSC_COMM_WORLD, "Computing the solution \n", errpetsc)


      nsteps = 1
      fileCount = 1

      ! Time loop
      DO loadstep=1, nsteps

        WRITE(*,*) " loadstep = ", loadstep

        ! compute time step
        ! dt = 0.01;
        ! dispCur = af*disp + (1.0-af)*dispPrev
        ! veloCur = af*velo + (1.0-af)*veloPrev
        ! acceCur = am*acce + (1.0-am)*accePrev
      
      ! Loop iteration check 2nd result, it should be zero
      Do iter=1, 5

        ! Loop over elements and compute the RHS

        call MatZeroEntries(matA, errpetsc)
        call VecZeroEntries(rhsVec, errpetsc)
        call VecZeroEntries(reacVec, errpetsc)
        
        LoopElem: DO ee=1, nElem
          write(*,*) " ee = ", ee

          dispElem = 0.0; veloElem = 0.0; acceElem = 0.0
          DO ii=1, npElem
            n1 = elemNodeConn(ee,ii)

            xNode(ii) = coords(n1,1)
            yNode(ii) = coords(n1,2)
            zNode(ii) = coords(n1,3)

            jj = (ii-1)*ndof
            kk = (n1-1)*ndof

            dispElem(jj+1) = disp(kk+1)
            dispElem(jj+2) = disp(kk+2)
            dispElem(jj+3) = disp(kk+3)

            veloElem(jj+1) = velo(kk+1)
            veloElem(jj+2) = velo(kk+2)
            veloElem(jj+3) = velo(kk+3)

            acceElem(jj+1) = acce(kk+1)
            acceElem(jj+2) = acce(kk+2)
            acceElem(jj+3) = acce(kk+3)
          END DO

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute striffness matrix and assemble it
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WRITE(*,*) "Computing stiffness"
          ! Compute the element force vector, including residual force
          call StiffnessResidualElasticityLinearHex(xNode, yNode, zNode, elemData, timeData, &
          dispElem, veloElem, acceElem, Klocal, Flocal)

          ! write(*,*) "dispElem: \n"
          ! write(*,*) dispElem

          ! Print the stiffness matrix, if you want
          write(*,*) "Klocal: \n"
          DO ii=1,nsize
            write(*,*) Klocal(ii,1), Klocal(ii,2), Klocal(ii,3), Klocal(ii,4),Klocal(ii,5), Klocal(ii,6), Klocal(ii,7), &
            Klocal(ii,8), Klocal(ii,9), Klocal(ii,10), Klocal(ii,11), Klocal(ii,12), Klocal(ii,13), Klocal(ii,14), &
            Klocal(ii,15), Klocal(ii,16), Klocal(ii,17), Klocal(ii,18), Klocal(ii,19), Klocal(ii,20), &
            Klocal(ii,21), Klocal(ii,22), Klocal(ii,23), Klocal(ii,24)
          END DO

          write(*,*) "Flocal: \n"
          write(*,*) Flocal

          ! Assemble the element matrix
  
          forAssyVec = ElemDofArray(ee,:)

          WRITE(*,*) "nsize = ", nsize
          !WRITE(*,*) forAssyVec

          LoopI1: DO ii=1,nsize
            row = forAssyVec(ii)
            IF( row /= -1) THEN
              LoopJ1: DO jj=1,nsize
                col = forAssyVec(jj)
                IF( col /= -1) THEN
                  call MatSetValue(matA, row, col, Klocal(ii,jj), ADD_VALUES, errpetsc)
                END IF
              END DO LoopJ1
            END IF
          END DO LoopI1
  
          write(*,*) "applying BCs for element ", ee

          ! Also, apply Dirichlet BCs while assembling
          ! Commented for now 
          ! TODO fix it for non-zero BC value
          ! LoopI: DO ii=1,nsize
          !   row = forAssyVec(ii)
          !   IF( row == -1) THEN
          !     fact = solnApplied(elemNodeConn(ee,ii))
          !     ! write(*,) ii, row, fact
          !     LoopJ: DO jj=1,nsize
          !       col = forAssyVec(jj)
          !       IF( col /= -1) THEN
          !         Flocal(jj) = Flocal(jj) - Klocal(jj,ii)*fact
          !       END IF
          !     END DO LoopJ
          !   END IF
          ! END DO LoopI
  
          ! Print the effective force vector, if you want
          !write(*,*) ee, Flocal(1), Flocal(2), Flocal(3), Flocal(4)
          !write(*,*) "\n"
  
          ! Assemble the element vector
          write(*,*) "assembling vector for element ", ee
          DO ii=1,nsize
            row = forAssyVec(ii)
            IF( row /= -1) THEN
              call VecSetValue(rhsVec, row, Flocal(ii), ADD_VALUES, errpetsc)
            END IF
          END DO

          ! END IF
        END DO LoopElem

        call PetscPrintf(PETSC_COMM_WORLD, "Assembly done. Fianlising the solver. \n", errpetsc)

        ! Add specified nodal force 
        DO ii=1,nFBC
          n1 = ForceBCs(ii,1)
          n2 = ForceBCs(ii,2)

          !row  = (n1-1)*ndof+n2
          ! subtract 1 since PETSc's indices start at zero
          row = NodeDofArrayNew(n1,n2)-1
          fact = ForceBCs(ii,3)

          !write(*,*) n1, n2, row, fact

          call VecSetValue(rhsVec, row, fact, ADD_VALUES, errpetsc)
        END DO 

        WRITE(*,*) " Computing velocity and acceleration"

        ! compute velocity and acceleration
        !velo = (disp - dispPrev2)/(2.0*dt);
        !acce = (disp - 2.0*dispPrev + dispPrev2)/DTT;

        write(*,*) "Assembly done. Fianlising the matrices and vectors"

        call MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY, errpetsc)
        !CHKERRQ(errpetsc)
        call MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY, errpetsc)
        !CHKERRQ(errpetsc)
  
        call VecAssemblyBegin(rhsVec, errpetsc)
        !CHKERRQ(errpetsc)
        call VecAssemblyEnd(rhsVec, errpetsc)
        !CHKERRQ(errpetsc)

        call VecAssemblyBegin(reacVec, errpetsc)
        !CHKERRQ(errpetsc)
        call VecAssemblyEnd(reacVec, errpetsc)
        !CHKERRQ(errpetsc)

        !call MatView(matA, PETSC_VIEWER_STDOUT_WORLD, errpetsc)
        !call VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

        call VecNorm(rhsVec, NORM_2, rNorm, errpetsc)

        write(*,*) iter, "... Norm = ", rNorm
        write(*,*) ""

        write(*,*) "Creating KSP context"
  
        !///////////////////////
        ! Create the linear solver and set various options
        !///////////////////////
  
        call KSPCreate(PETSC_COMM_WORLD, ksp, errpetsc)
        !CHKERRQ(errpetsc)
        call KSPSetOperators(ksp, matA, matA, errpetsc)
        !CHKERRQ(errpetsc)
        !call KSPSetType(ksp, KSPCG, errpetsc)
        call KSPSetType(ksp, KSPBCGS, errpetsc)
        !CHKERRQ(errpetsc)

        fact=0.0000001
        ind=5000
        !call KSPSetTolerances(ksp, fact, fact, fact, ind, errpetsc)

        ! Set KSP options from the input file
        call KSPSetFromOptions(ksp, errpetsc) !CHKERRQ(errpetsc);

        write(*,*) "Creating PC context"

        call KSPGetPC(ksp, pc, errpetsc) !CHKERRQ(errpetsc);
        !call PCSetType(pc, PCJACOBI, errpetsc); !CHKERRQ(errpetsc);
        !call PCSetType(pc, PCILU, errpetsc) !CHKERRQ(errpetsc);
        call PCSetType(pc, PCLU, errpetsc) !CHKERRQ(errpetsc);
        ! Set PC options from the input file
        call PCSetFromOptions(pc, errpetsc)
  
        call KSPSetUp(ksp, errpetsc); !CHKERRQ(errpetsc);
  
        write(*,*) "Solving the matrix system"
  
        call KSPSolve(ksp, rhsVec, solnVec, errpetsc)
        CHKERRQ(errpetsc);
  
        call KSPGetIterationNumber(ksp, its, errpetsc) !CHKERRQ(errpetsc);
        call KSPGetConvergedReason(ksp, reason, errpetsc)
  
        IF(reason<0) THEN
          write(*,*) "Divergence."
        ELSE
          call KSPGetIterationNumber(ksp, its, errpetsc);
          write(*,*) "Convergence in", its, " iterations."
          !PetscSynchronizedPrintf(MPI_COMM_WORLD, "Convergence in %d iterations.\n", (int)its);
        END IF
  
        call VecView(solnVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)
  
        write(*,*) "Writing VTK file"
  
        ! Applied BC values
        ! DO ii=1,nNode
          ! disp(ii) = solnApplied(ii)
          !write(*,*) ii, ind, solnApplied(ii)
        ! END DO
  
        ! Add solution for the free DOF
        call VecGetArrayF90(solnVec, xx_v, errpetsc)
  
        write(*,*) "totalDOF = ", totalDOF
  
        DO ii=1,totalDOF
          ind = assyForSoln(ii)
          fact = xx_v(ii)
          disp(ind) = disp(ind) + fact 
        END DO
        CLOSE(1)
        
        write(*,*) "disp at the end of iteration loop: \n"
        write(*,*) disp
      
      END DO  ! End of iteration loop 
        ind = SCAN(inputFileName,".")
        !string1 = inputFileName(1:ind-1)
        !string2 = inputFileName(ind+1:)
        outFileName = trim(inputFileName(1:ind-1)) // ".vtk"
        !WRITE(*,*) outFileName
  
        call PetscPrintf(PETSC_COMM_WORLD, "Writing VTK file \n", errpetsc)

        WRITE(outFileName,'(A,I5.5,A)') "Elasticity-soln-", fileCount, ".vtk"

        ! write the solution to the VTK file
        call writeoutputvtk(ndim, nElem, nNode, npElem, ndof, coords, elemNodeConn, elem_proc_id, disp, outFileName)

        fileCount = fileCount+1

        ! store the variables
        dispPrev2 = dispPrev
        dispPrev  = disp
        veloPrev  = velo
        accePrev  = acce

        stepsCompleted = stepsCompleted + 1
        timeNow = timeNow + dt

      END DO ! Time loop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! End of computations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD, "Deallocating the memory \n", errpetsc)

      call VecRestoreArrayF90(solnVec, xx_v, errpetsc)

      IF( ALLOCATED(coords) )        DEALLOCATE(coords)
      IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
      IF( ALLOCATED(DirichletBCs) )  DEALLOCATE(DirichletBCs)
      IF( ALLOCATED(ForceBCs) )      DEALLOCATE(ForceBCs)
      IF( ALLOCATED(NodeTypeOld) )      DEALLOCATE(NodeTypeOld)
      IF( ALLOCATED(NodeTypeNew) )      DEALLOCATE(NodeTypeNew)
      IF( ALLOCATED(NodeDofArrayOld) )  DEALLOCATE(NodeDofArrayOld)
      IF( ALLOCATED(NodeDofArrayNew) )  DEALLOCATE(NodeDofArrayNew)
      IF( ALLOCATED(ElemDofArray) )  DEALLOCATE(ElemDofArray)
      IF( ALLOCATED(assyForSoln) )   DEALLOCATE(assyForSoln)
      IF( ALLOCATED(solnApplied) )   DEALLOCATE(solnApplied)
      IF( ALLOCATED(solnVTK) )       DEALLOCATE(solnVTK)

      IF( ALLOCATED(elem_proc_id) ) THEN
          DEALLOCATE(elem_proc_id)
      END IF

      IF( ALLOCATED( disp ) ) DEALLOCATE(disp)
      IF( ALLOCATED( velo ) ) DEALLOCATE(velo)
      IF( ALLOCATED( acce ) ) DEALLOCATE(acce)
      IF( ALLOCATED( dispPrev ) ) DEALLOCATE(dispPrev)
      IF( ALLOCATED( dispPrev2 ) ) DEALLOCATE(dispPrev2)
      IF( ALLOCATED( veloPrev ) ) DEALLOCATE(veloPrev)
      IF( ALLOCATED( accePrev ) ) DEALLOCATE(accePrev)
      IF( ALLOCATED( dispCur ) ) DEALLOCATE(dispCur)
      IF( ALLOCATED( veloCur ) ) DEALLOCATE(veloCur)
      IF( ALLOCATED( acceCur ) ) DEALLOCATE(acceCur)

      call VecDestroy(rhsVec, errpetsc)
      call VecDestroy(solnVec, errpetsc)
      call VecDestroy(reacVec, errpetsc)
      call MatDestroy(matA, errpetsc)
      call KSPDestroy(ksp, errpetsc)


      ! tend = MPI_Wtime()

      ! write(*,*) "That took ", (tend-tstart), "seconds"

      call PetscPrintf(PETSC_COMM_WORLD, " Program is successful \n", errpetsc)

      call PetscFinalize(errpetsc)
      ! CHKERRQ(errpetsc)

      END PROGRAM QuadElasticityImplicit
