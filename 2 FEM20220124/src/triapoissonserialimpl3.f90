!
! Program for 2D Poisson equation using linear triangular elements
!
! This implementation is for single processor
!
! This implementation follow the second one in 
! "triapoissonserialimpl2.F"
!
! In this implementation, several entries are inserted at a time,
! into PETSc matrices/vectors, as recommended by PETSc.
! This is efficient for large-scale models.
!
! The change in this implementation is that PETSc is used a module
! with a class object
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-Oct-2017
! Place : Swansea, UK
!
!
!

      PROGRAM TriaMeshPoissonEquation

      USE WriterVTK
      USE Module_SolverPetsc

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

      PetscErrorCode errpetsc;

      DOUBLE PRECISION :: tstart, tend

      LOGICAL :: FILEEXISTS, ISOPEN

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
      INTEGER, DIMENSION(:), ALLOCATABLE :: DirichletBC_row_nums
      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_procid

      INTEGER :: forAssyVec(3)

      CHARACTER(len=32) :: arg

      TYPE(PetscSolver) :: solverpetsc


      Vec solnTemp

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

      PetscScalar xx_v(1)
      ! PetscScalar, pointer :: xx_v(:)

      PetscOffset xx_i

      PetscInt  n_mpi_procs     ! total number of processors in the group
      PetscInt  this_mpi_proc   ! rank of the current processor

      PetscInt  its;


      CHARACTER (LEN=100) :: infileNodes
      CHARACTER (LEN=100) :: infileElems
      CHARACTER (LEN=100) :: infileDBCs
      CHARACTER (LEN=100) :: outFileName

      !Set file names
      !The file names are specified as inputs from the command line
      IF( iargc() < 3 ) THEN
        WRITE(*,*) "Number of input files is not sufficient "
        WRITE(*,*) "You must enter names of THREE files"
        ! WRITE(*,*) "a.) Node file, b.) Element file, and c.) Dirichlet BC file"
        STOP "Aborting..."
      END IF

      CALL getarg(1, arg)
      infileNodes = arg
      WRITE(*,*) infileNodes

      CALL getarg(2, arg)
      infileElems = arg
      WRITE(*,*) infileElems

      CALL getarg(3, arg)
      infileDBCs = arg
      WRITE(*,*) infileDBCs


      ! intialise PETSc environment
      !
      call PetscInitialize(PETSC_NULL_CHARACTER, errpetsc)
      ! CHKERRQ(errpetsc)

      tstart = MPI_Wtime()

      call MPI_Comm_size(PETSC_COMM_WORLD, n_mpi_procs, errpetsc);
      call MPI_Comm_rank(PETSC_COMM_WORLD, this_mpi_proc, errpetsc);

      write(*,*) " this_mpi_proc = ", this_mpi_proc


      !
      ! Read nodal data files
      !

      ! check if the file exists

      INQUIRE(file=infileNodes, EXIST=FILEEXISTS)
      write(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileNodes, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nNode = 0
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")
      DO
        READ(1,*, iostat=io)
        IF (io/=0) EXIT
        nNode = nNode + 1
      END DO 
      write(*,*) "Number of nodes = ", nNode
      CLOSE(1)

      ALLOCATE(coords(nNode,ndim))

      ! Open the file and store nodal coordinates
       
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")

      DO ii=1,nNode
        READ(1,*, iostat=io) nn, coords(ii,1), coords(ii,2)
        ! WRITE(*,*) nn, coords(ii,1), coords(ii,2)
      END DO 
      CLOSE(1)


      !
      ! Read element data files
      !

      ! check if the file exists

      INQUIRE(file=infileElems, EXIST=FILEEXISTS)
      write(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileElems, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nElem = 0
      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO
        READ(2,*, iostat=io)
        IF (io/=0) EXIT
        nElem = nElem + 1
      END DO 
      write(*,*) "Number of Elements = ", nElem
      CLOSE(2)

      ! This program is hardcoded for triangular elements
      npElem = 3
      ALLOCATE(elemNodeConn(nElem,3))

      ! Open the file and store nodal coordinates
       
      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO ii=1,nElem
        READ(2,*, iostat=io) nn, n1, n2, n3
        elemNodeConn(ii,1) = n1
        elemNodeConn(ii,2) = n2
        elemNodeConn(ii,3) = n3
        ! write(*,*) nn, n1, n2, n3
      !   WRITE(*,*) nn, elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3)
      END DO 
      CLOSE(2)


      !
      ! Read Dirichlet BC data
      !

      ! check if the file exists

      INQUIRE(file=infileDBCs, EXIST=FILEEXISTS)
      write(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileDBCs, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nDBC = 0
      OPEN(3, file=infileDBCs,STATUS="OLD",ACTION="READ")
      DO
        READ(3,*, iostat=io)
        IF (io/=0) EXIT
        nDBC = nDBC + 1
      END DO 
      write(*,*) "Number of Dirichlet BCs = ", nDBC
      CLOSE(3)

      ! This program is hardcoded for triangular elements

      ALLOCATE(DirichletBCs(nDBC,3))

      call VecCreate(PETSC_COMM_WORLD, solnTemp, errpetsc)
      !CHKERRQ(errpetsc)
      ind = nNode*ndof
      call VecSetSizes(solnTemp, ind, ind, errpetsc)
      !CHKERRQ(errpetsc)
      call VecSetFromOptions(solnTemp, errpetsc)
      !CHKERRQ(errpetsc)

      call VecZeroEntries(solnTemp, errpetsc)

      ALLOCATE( DirichletBC_row_nums(nDBC) )

      ! Open the file and store nodal coordinates
       
      OPEN(3, file=infileDBCs,STATUS="OLD",ACTION="READ")
      DO ii=1,nDBC
        READ(3,*, iostat=io) n1, n2, fact
        DirichletBCs(ii,1) = n1
        DirichletBCs(ii,2) = n2
        DirichletBCs(ii,3) = fact
        ! write(*,*) nn, xc, yc

        ind = (n1-1)*n2
        DirichletBC_row_nums(ii) = ind

        call VecSetValue(solnTemp, ind, fact,
     1    INSERT_VALUES, errpetsc)
      END DO 
      CLOSE(3)



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
          ! IF(NodeType(ii,jj) == 0) THEN
            NodeDofArray(ii,jj) = totalDOF
            totalDOF = totalDOF + 1
          ! END IF
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
            ElemDofArray(ee, n1+jj) = NodeDofArray(n2,jj) - 1
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

      write(*,*) "Preparing matrix pattern DONE"


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! create PETSc variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " totalDOF = ", totalDOF

      ndofs_local=totalDOF

      write(*,*) " Creating PETSc vectors 1"

      ALLOCATE(nnzVec(totalDOF))

      ind = 50
      IF(totalDOF < 50) THEN
        ind = totalDOF
      END IF
      DO ii=1,totalDOF
        nnzVec(ii) = ind
      END DO

      write(*,*) " Initialising petsc solver"

      !Initialize the petsc solver
      call solverpetsc%initialise(totalDOF, totalDOF, nnzVec, nnzVec)

      write(*,*) " Initialise the Matrix pattern "

      nsize = npElem*ndof
      fact = 0.0
      Klocal = 0.0
      WRITE(*,*) Klocal
      ! Initialise the Matrix pattern
      LoopElem: DO ee=1,nElem
        forAssyVec = ElemDofArray(ee,:)

        ! WRITE(*,*) forAssyVec
        call MatSetValues(solverpetsc%mtx, nsize, forAssyVec,
     1                          nsize, forAssyVec,
     2                          Klocal, INSERT_VALUES, errpetsc)
      END DO LoopElem


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! compute element matrices and vectors and assemble them
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) " Matrix pattern initialised. "
      write(*,*) " Computing element matrices now. "

      call solverpetsc%setZero()

      write(*,*) " Generating element matrices and vectors "

      DO ee=1,nElem
        ! IF(elems[ee]->get_subdomain_id() == this_mpi_proc) THEN
          ! write(*,*) " ee = ", ee

          Klocal = 0.0; Flocal = 0.0

          n1 = elemNodeConn(ee,1)
          n2 = elemNodeConn(ee,2)
          n3 = elemNodeConn(ee,3)
          
          ! xcoord(1) = coords(n1,1);  ycoord(1) = coords(n1,2)
          ! xcoord(2) = coords(n2,1);  ycoord(2) = coords(n2,2)
          ! xcoord(3) = coords(n3,1);  ycoord(3) = coords(n3,2)

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

          ! IF(ee < 10) THEN
          !   write(*,*) Klocal
          ! END IF

          ! Assemble the element matrix
          ! write(*,*) "assembling matrix for element ", ee
          forAssyVec = ElemDofArray(ee,:)
          call MatSetValues(solverpetsc%mtx, nsize, forAssyVec,
     1                            nsize, forAssyVec,
     2                            Klocal, ADD_VALUES, errpetsc)

          ! write(*,*) "applying BCs for element ", ee

          ! Assemble the element vector
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

          ! write(*,*) "assembling vector for element ", ee
          call VecSetValues(solverpetsc%rhsVec, nsize, forAssyVec,
     1                              Flocal, ADD_VALUES, errpetsc)

        ! END IF
      END DO

      ! write(*,*) "Assembly done. Fianlising the matrices and vectors"

      write(*,*) "Assembly done. Applying boundary conditions "

      WRITE(*,*) " nDBC = ", nDBC
      ! WRITE(*,*) DirichletBC_row_nums

      ! Apply Dirichlet BCs

      call MatAssemblyBegin(solverpetsc%mtx,
     1  MAT_FINAL_ASSEMBLY, errpetsc)
      ! CHKERRQ(errpetsc)
      call MatAssemblyEnd(solverpetsc%mtx, 
     1 MAT_FINAL_ASSEMBLY, errpetsc)

    !   call MatMult(solverpetsc%mtx, solnTemp, 
    !  1 solverpetsc%rhsVec, errpetsc)

    !   fact = -1.0
    !   call  VecScale(solverpetsc%rhsVec, fact, errpetsc)

      fact = 1.0
      call MatZeroRows(solverpetsc%mtx, nDBC, DirichletBC_row_nums, 
     1  fact, solnTemp, solverpetsc%rhsVec, errpetsc)


      call solverpetsc%factoriseAndSolve()

      write(*,*) "Writing VTK file"

      ind = nNode*ndof
      ALLOCATE(solnVTK(ind))

      ! Applied BC values
      DO ii=1,ind
        solnVTK(ii) = solnApplied(ii)
      END DO

      ! Add solution for the free DOF
      call VecGetArray(solverpetsc%solnVec, xx_v, xx_i, errpetsc)
      ! call VecGetArrayF90(solnVec, xx_v, errpetsc)

      OPEN(1, file="temp.dat", STATUS="UNKNOWN", ACTION="WRITE")

      DO ii=1,totalDOF
        ind = assyForSoln(ii)
        fact = xx_v(xx_i+ii)
        ! fact = xx_v(ii)
        ! x1 = coords(ind,1); y1 = coords(ind,2);
        ! val = (cosh(PI*y1)-sinh(PI*y1)/tanh(PI))*sin(PI*x1)
        write(1,*) ii, ind, fact
        solnVTK(ind) = fact
      END DO
      CLOSE(1)

      WRITE(outFileName,'(A,I5.5,A)') "PoissonTria-soln.vtk"

      ! write the solution to the VTK file
        call writeoutputvtk(
     1     ndim,
     2     nElem,
     3     nNode,
     4     npElem,
     5     ndof, 
     6     coords,
     7     elemNodeConn,
     8     elem_procid,
     9     solnVTK,
     +     outFileName)

      write(*,*) "Deallocating the memory"

      ! call VecRestoreArray(solnVec, xx_v, xx_i, errpetsc)
      ! call VecRestoreArrayF90(solnVec, xx_v, errpetsc)

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
      
      call solverpetsc%free()
      call VecDestroy(solnTemp, errpetsc)
      ! call VecDestroy(rhsVec, errpetsc)
      ! call VecDestroy(solnVec, errpetsc)
      ! call VecDestroy(reacVec, errpetsc)
      ! call MatDestroy(solverpetsc%mtx, errpetsc)
      ! call KSPDestroy(ksp, errpetsc)


      tend = MPI_Wtime()

      write(*,*) "That took ", (tend-tstart), "seconds"

      call PetscFinalize(errpetsc)
      ! CHKERRQ(errpetsc)

      write(*,*) "Program is successful"

      END PROGRAM TriaMeshPoissonEquation
