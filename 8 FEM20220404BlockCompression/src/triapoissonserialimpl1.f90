!
! Program for 2D Poisson equation using linear triangular elements
!
! This is the basic implementation without using OOPs concept. However, 
! it serves as the reference for further implementations
!
! This implementation is for single processor
!
! Insertions to PETSc matrices/vectors are done entry by entry.
! This is not efficient for large-scale models. So, in further implementations 
! several entries are inserted at a time, as recommended by PETSc
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-Oct-2017
!  Original file
! Date  : 22-Jan-2022
!  Merged input files and updated the code
! Place : Swansea, UK
!
!
!

      PROGRAM TriaMeshPoissonEquation

      USE WriterVTK

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"


!#define PETSC_USE_FORTRAN_MODULES


!#if defined(PETSC_USE_FORTRAN_MODULES)
      USE petscvec
      USE petscmat
      USE petscksp
      USE petscpc
!#endif
      IMPLICIT NONE

      ! declare variables
      ! IMPLICIT NONE
      PetscErrorCode errpetsc;

      DOUBLE PRECISION :: tstart, tend

      LOGICAL :: FILEEXISTS=.FALSE., ISOPEN

      INTEGER :: ndim=2
      INTEGER :: ndof=1
      INTEGER :: io, nn
      DOUBLE PRECISION :: PI=ACOS(-1.0)
      DOUBLE PRECISION :: xc, yc, area, fact, val
      DOUBLE PRECISION :: x1, x2, x3, y1, y2, y3
      DOUBLE PRECISION :: x13, x21, x32, y31, y12, y23
      DOUBLE PRECISION :: Bmat(3,2), dispC(3)
      DOUBLE PRECISION :: BmatTrans(2,3)
      DOUBLE PRECISION :: Klocal(3,3), Flocal(3)


      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DirichletBCs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnApplied
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnVTK

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: elemNodeConn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeDofArray
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemDofArray
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: forAssyMat
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeType
      INTEGER, DIMENSION(:), ALLOCATABLE :: assyForSoln
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnzVec

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_procid

      INTEGER :: forAssyVec(3)

      CHARACTER(len=32) :: arg

      Vec  tempVec
      Vec  rhsVec      ! RHS vector
      Vec  reacVec     ! Reaction vector
      Vec  solnVec     ! Solution vector for the matrix system
      ! Vec  solnVecFull ! Solution vector for the mesh
      Mat  matA        ! Linear system matrix
      KSP  ksp         ! Linear solver context
      PC   pc          ! Preconditioner context
      PetscReal norm   !norm of solution error
      PetscLogStage stage;


      PetscInt  nNode  ! number of nodes
      PetscInt  nElem  ! number of global elements (in the whole model)
      ! PetscInt  nElem  ! number of local  elements (owned by the local processor)
      PetscInt  npElem ! number of nodes per element
      PetscInt  nDBC   ! number of Dirichlet BCs
      PetscInt  nFBC   ! number of force BCs - specified nodal forces
      PetscInt  totalDOF ! total number of global DOF in the whole model
      PetscInt  ndofs_local ! total number of local DOF owned by the current processor
      PetscInt  ee, ii, jj, kk, ind
      PetscInt  count, row, col
      PetscInt  n1, n2, n3, nsize

      ! PetscScalar *xx_v
      ! PetscScalar xx_v(1)
      PetscScalar, pointer :: xx_v(:)

      PetscOffset xx_i


      PetscInt  n_mpi_procs     ! total number of processors in the group
      PetscInt  this_mpi_proc   ! rank of the current processor

      MatInfo info;
      KSPConvergedReason reason;
      PetscInt  its;


      CHARACTER (LEN=100) :: inputFileName
      CHARACTER (LEN=100) :: outFileName
      CHARACTER (LEN=20)  :: chartemp

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

      !tstart = MPI_Wtime()

      this_mpi_proc = -1
      !call MPI_Comm_size(PETSC_COMM_WORLD, n_mpi_procs, errpetsc);
      !call MPI_Comm_rank(PETSC_COMM_WORLD, this_mpi_proc, errpetsc);

      write(*,*) " this_mpi_proc = ", this_mpi_proc


      ! check if the file exists

      write(*,*) " FILEEXISTS = ", FILEEXISTS

      INQUIRE(FILE=inputFileName, EXIST=FILEEXISTS)
      write(*,*) " this_mpi_proc = ", this_mpi_proc
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
        READ(1,*, iostat=io) nn, coords(ii,1), coords(ii,2)
        !WRITE(*,*) nn, coords(ii,1), coords(ii,2)
      END DO 


      ! Read number of elements
      READ(1,*, iostat=io) chartemp, nElem
      write(*,*) chartemp, nNode
      write(*,*) "Number of Elements = ", nElem

      ! This program is hardcoded for triangular elements
      npElem = 3
      ALLOCATE(elemNodeConn(nElem,npElem))

      DO ii=1,nElem
        READ(1,*, iostat=io) nn, n1, n2, n3
        elemNodeConn(ii,1) = n1
        elemNodeConn(ii,2) = n2
        elemNodeConn(ii,3) = n3
        !write(*,*) nn, n1, n2, n3, n4
        !WRITE(*,*) nn, elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3)
      END DO 

      ! Read Dirichlet BC data
      READ(1,*, iostat=io) chartemp, nDBC
      write(*,*) chartemp, nDBC
      write(*,*) "Number of Dirichlet BCs = ", nDBC

      ALLOCATE(DirichletBCs(nDBC,3))
      DO ii=1,nDBC
        READ(1,*, iostat=io) nn, xc, yc
        DirichletBCs(ii,1) = nn
        DirichletBCs(ii,2) = xc
        DirichletBCs(ii,3) = yc
        ! write(*,*) nn, xc, yc
      END DO 
      CLOSE(1)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !prepare the global matrix pattern
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " "
      write(*,*) " "
      write(*,*) " Input files have been read successfully "
      write(*,*) " "
      write(*,*) " "
      write(*,*) " Mesh statistics ..... "
      write(*,*) " nElem    ", nElem
      write(*,*) " nNode    ", nNode
      write(*,*) " npElem   ", npElem
      write(*,*) " ndof     ", ndof

      ALLOCATE(elem_procid(nElem))

      elem_procid = 0

      write(*,*) " Preparing matrix pattern "

      ! Initialize NodeType and NodeDofArray arrays
      !
      ALLOCATE(NodeType(nNode, ndof))
      ALLOCATE(NodeDofArray(nNode, ndof))
      n1 = nNode*ndof
      ALLOCATE(solnApplied(n1))

      DO ii=1,nNode
        DO jj=1,ndof
          NodeType(ii,jj)   = 0
          NodeDofArray(ii,jj) = 0
        END DO
      END DO

      DO ii=1,nDBC
        NodeType(DirichletBCs(ii,1), DirichletBCs(ii,2)) = 1
        n1 = (DirichletBCs(ii,1)-1)*ndof + DirichletBCs(ii,2)
        solnApplied(n1) = DirichletBCs(ii,3)
      END DO

      totalDOF = 1;
      DO ii=1,nNode
        DO jj=1,ndof
          IF(NodeType(ii,jj) == 0) THEN
            NodeDofArray(ii,jj) = totalDOF
            totalDOF = totalDOF + 1
          END IF
        END DO
      END DO

      totalDOF = totalDOF-1

      write(*,*) " Total number of DOF = ", totalDOF

      write(*,*) " Creating arrays 1 "

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
            ElemDofArray(ee, n1+jj) = NodeDofArray(n2,jj)-1
          END DO
        END DO
      END DO

      ! write(*,*) " LM array "
      ! nsize=npElem*ndof
      ! DO ee=1,nElem
      !   write(*,*)  ElemDofArray(ee,1), ElemDofArray(ee,2), ElemDofArray(ee,3)
      ! END DO

      write(*,*) " Creating arrays 2 "

      ALLOCATE(assyForSoln(totalDOF))
      write(*,*) " Creating arrays 2 "

      count = 1
      DO ii=1,nNode
        DO jj=1,ndof
          ! write(*,*) ii, jj, NodeDofArray(ii,jj)
          IF(NodeDofArray(ii,jj) /= 0) THEN
            assyForSoln(count) = (ii-1)*ndof + jj;
            count = count + 1
          END IF
        END DO
      END DO

      write(*,*) " assyForSoln array "

      ! DO ii=1,totalDOF
      !   write(*,*) ii, assyForSoln(ii)
      ! END DO

    !   ! ALLOCATE(forAssyMat(totalDOF, totalDOF))
    
    !   ! nsize = npElem*ndof
    !   ! DO ee=1,nElem
    !   !   DO ii=1,nsize
    !   !     kk = ElemDofArray(ee,ii)

    !   !     IF(kk /= 0) THEN
    !   !       DO jj=1,nsize
    !   !         IF(kk /= 0) THEN
    !   !       !     forAssyMat[count1].push_back(tt[jj]);
    !   !         END IF
    !   !       END DO
    !   !     END IF
    !   !   END DO
    !   ! END DO

    !   write(*,*) "Preparing matrix pattern DONE"


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! create PETSc variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " totalDOF = ", totalDOF

      ndofs_local=totalDOF

      write(*,*) " Creating PETSc vectors 1"

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
        nnzVec(ii) = 50
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
      LoopElem: DO ee=1,nElem
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
      END DO LoopElem

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

      write(*,*) " Generating element matrices and vectors "

      DO ee=1,nElem
          Klocal = 0.0
          Flocal = 0.0

          n1 = elemNodeConn(ee,1)
          n2 = elemNodeConn(ee,2)
          n3 = elemNodeConn(ee,3)
          
          x1 = coords(n1,1);  y1 = coords(n1,2)
          x2 = coords(n2,1);  y2 = coords(n2,2)
          x3 = coords(n3,1);  y3 = coords(n3,2)

          ! Compute the element stiffness matrix and force vector
  
          area = 0.5*(x2*y3 - x3*y2 + x3*y1 - x1*y3 + x1*y2 - x2*y1)
          ! write(*,*) " area = ", ee, area

          x13=x1-x3;  x21=x2-x1;  x32=x3-x2;
          y31=y3-y1;  y12=y1-y2;  y23=y2-y3;

          Bmat(1,1) = y23; Bmat(2,1) = y31; Bmat(3,1) = y12;
          Bmat(1,2) = x32; Bmat(2,2) = x13; Bmat(3,2) = x21;
  
          Bmat = Bmat/(2.0*area)

          BmatTrans = TRANSPOSE(Bmat)
  
          Klocal = MATMUL(Bmat, BmatTrans)
          Klocal = area*Klocal

          ! Print the stiffness matrix, if you want
          !DO ii=1,nsize
            !write(*,*) Klocal(ii,1), Klocal(ii,2), Klocal(ii,3)
          !END DO
          !write(*,*) "\n"

          ! Assemble the element matrix

          forAssyVec = ElemDofArray(ee,:)

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


          ! write(*,*) "applying BCs for element ", ee

          ! Also, apply Dirichlet BCs while assembling
          LoopI: DO ii=1,nsize
            row = forAssyVec(ii)
            IF( row == -1) THEN
              fact = solnApplied(elemNodeConn(ee,ii))
              ! write(*,) ii, row, fact
              LoopJ: DO jj=1,nsize
                col = forAssyVec(jj)
                IF( col /= -1) THEN
                  Flocal(jj) = Flocal(jj) - Klocal(jj,ii)*fact
                END IF
              END DO LoopJ
            END IF
          END DO LoopI

          ! Print the effective force vector, if you want
          !write(*,*) Flocal(1), Flocal(2), Flocal(3)
          !write(*,*) "\n"

          ! Assemble the element vector
          ! write(*,*) "assembling vector for element ", ee
          DO ii=1,nsize
            row = forAssyVec(ii)
            IF( row /= -1) THEN
              call VecSetValue(rhsVec, row, Flocal(ii), ADD_VALUES, errpetsc)
            END IF
          END DO

      END DO

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

      write(*,*) "Creating KSP context"

      !///////////////////////
      ! Create the linear solver and set various options
      !///////////////////////

      call KSPCreate(PETSC_COMM_WORLD, ksp, errpetsc)
      !CHKERRQ(errpetsc)
      call KSPSetOperators(ksp, matA, matA, errpetsc)
      !CHKERRQ(errpetsc)
      call KSPSetType(ksp, KSPCG, errpetsc)
      !CHKERRQ(errpetsc)

      fact=0.0000001
      ind=5000
      !call KSPSetTolerances(ksp, fact, fact, fact, ind, errpetsc)

      ! Set KSP options from the input file
      call KSPSetFromOptions(ksp, errpetsc) !CHKERRQ(errpetsc);

      write(*,*) "Creating PC context"

      call KSPGetPC(ksp, pc, errpetsc) !CHKERRQ(errpetsc);
      !call PCSetType(pc, PCJACOBI, errpetsc); !CHKERRQ(errpetsc);
      call PCSetType(pc, PCILU, errpetsc) !CHKERRQ(errpetsc);
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

      !call VecView(solnVec, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

      write(*,*) "Writing VTK file"

      ind = nNode*ndof
      ALLOCATE(solnVTK(ind))

      ! Applied BC values
      DO ii=1,ind
        solnVTK(ii) = solnApplied(ii)
        !write(*,*) ii, ind, solnApplied(ii)
      END DO

      ! Add solution for the free DOF
      !call VecGetArray(solnVec, &xx_v, errpetsc)
      !call VecGetArray(solnVec, xx_v, xx_i, errpetsc)
      call VecGetArrayF90(solnVec, xx_v, errpetsc)

      OPEN(1, file="temp.dat", STATUS="UNKNOWN", ACTION="WRITE")

      write(*,*) "totalDOF = ", totalDOF, xx_i

      DO ii=1,totalDOF
        ind = assyForSoln(ii)
        !fact = xx_v(xx_i+ii)
        fact = xx_v(ii)
        ! x1 = coords(ind,1); y1 = coords(ind,2);
        ! val = (cosh(PI*y1)-sinh(PI*y1)/tanh(PI))*sin(PI*x1)
        !write(*,*) ii, ind, fact
        write(1,*) ii, ind, fact
        solnVTK(ind) = fact
      END DO
      CLOSE(1)

      ind = SCAN(inputFileName,".")
      !string1 = inputFileName(1:ind-1)
      !string2 = inputFileName(ind+1:)
      outFileName = trim(inputFileName(1:ind-1)) // ".vtk"
      !WRITE(*,*) outFileName

      call writeoutputvtk(ndim, nElem, nNode, npElem, ndof, coords, elemNodeConn, elem_procid, solnVTK, outFileName)


      write(*,*) "Deallocating the memory"

      !call VecRestoreArray(solnVec, &xx_v, errpetsc)
      !call VecRestoreArray(solnVec, xx_v, xx_i, errpetsc)
      call VecRestoreArrayF90(solnVec, xx_v, errpetsc)

      IF( ALLOCATED(coords) )        DEALLOCATE(coords)
      IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
      IF( ALLOCATED(DirichletBCs) )  DEALLOCATE(DirichletBCs)
      IF( ALLOCATED(NodeType) )      DEALLOCATE(NodeType)
      IF( ALLOCATED(NodeDofArray) )  DEALLOCATE(NodeDofArray)
      IF( ALLOCATED(ElemDofArray) )  DEALLOCATE(ElemDofArray)
      IF( ALLOCATED(assyForSoln) )   DEALLOCATE(assyForSoln)
      IF( ALLOCATED(solnApplied) )   DEALLOCATE(solnApplied)
      IF( ALLOCATED(solnVTK) )       DEALLOCATE(solnVTK)
      IF( ALLOCATED(nnzVec) )        DEALLOCATE(nnzVec)

      IF( ALLOCATED(elem_procid) ) THEN
          DEALLOCATE(elem_procid)
      END IF


      call VecDestroy(rhsVec, errpetsc)
      call VecDestroy(solnVec, errpetsc)
      call VecDestroy(reacVec, errpetsc)
      call MatDestroy(matA, errpetsc)
      call KSPDestroy(ksp, errpetsc)


      !tend = MPI_Wtime()

      !write(*,*) "That took ", (tend-tstart), "seconds"

      call PetscFinalize(errpetsc)
      ! CHKERRQ(errpetsc)

      write(*,*) "Program is successful"

      END PROGRAM TriaMeshPoissonEquation
