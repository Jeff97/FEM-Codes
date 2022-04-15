! Subroutines for calculating the stiffness
! and residual (or force) for the simple
! Lagrange family of elements for
!
! 2D Linear Elasticity using displacement formulation
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!

      MODULE ElementUtilitiesElasticity3D

      USE QuadratureUtilities
      USE ElementUtilitiesBasisFuncs

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityLinearPrism(
     1   xNode, yNode, zNode,    ! X,Y and Z coordinates of the 6 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 3 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 6 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)   :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(6,6) :: Dmat
      DOUBLE PRECISION, DIMENSION(18)   :: Flocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(3)
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6), dNdz(6)

      INTEGER, PARAMETER :: nGP=1, nlbf=6
      INTEGER :: DEGREE=1, ETYPE=6, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E    = elemData(1)
      nu   = elemData(2)
      dens = elemData(3)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)


      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = Dmat(1,1);   Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = Dmat(1,1)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2


      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsPrism(nGP, gpts1, gpts2, gpts3, gwts)


      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., ETYPE, DEGREE, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        !WRITE(*,*) " VOLUME = ", gp, Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3

          b1 = dispC(jj+1)
          b2 = dispC(jj+2)
          b3 = dispC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dNdz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dNdz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dNdz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1) !sxx
        strain(2) = grad(2,2) !syy
        strain(3) = grad(3,3) !szz
        strain(4) = 0.5*(grad(1,2)+grad(2,1))  !sxy
        strain(5) = 0.5*(grad(2,3)+grad(3,2))  !syz
        strain(6) = 0.5*(grad(1,3)+grad(3,1))  !szx

        stress = MATMUL(Dmat, strain)

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b3 = dvol*dNdz(ii)
          b4 = (dens*dvol)*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(4) - b3*stress(6)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(4) - b2*stress(2) - b3*stress(5)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3) - b1*stress(6) - b2*stress(5) - b3*stress(3)
        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearPrism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearPrism(
     1   xNode, yNode, zNode,    ! X,Y and Z coordinates of the 6 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(18)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(18,18) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6), dNdz(6)

      INTEGER, PARAMETER :: nGP=2, nlbf=6
      INTEGER :: DEGREE=1, ETYPE=6, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=18

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsPrism(nGP, gpts1, gpts2, gpts3, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., ETYPE, DEGREE, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = (dens*dvol)*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
          END DO
        END DO
      END DO

      ! WRITE(*,*) Klocal

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearPrism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualElasticityLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)     :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(12,12) :: Klocal
      DOUBLE PRECISION, DIMENSION(6,6)   :: Dmat
      DOUBLE PRECISION, DIMENSION(12)    :: Flocal

      DOUBLE PRECISION, DIMENSION(6,12) :: Bmat
      DOUBLE PRECISION, DIMENSION(12,6) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: timefact, af, E, nu, bforce(3)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dNdz(4)

      INTEGER, PARAMETER :: nGP=8, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)

      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = b1*(1.0-nu); Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = b1*(1.0-nu)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts)
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25;      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3
          b1 = valC(jj+1)
          b2 = valC(jj+2)
          b3 = valC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dNdz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dNdz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dNdz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = grad(3,3)
        strain(4) = 0.5*(grad(1,2)+grad(2,1))
        strain(5) = 0.5*(grad(2,3)+grad(3,2))
        strain(6) = 0.5*(grad(1,3)+grad(3,1))

        stress = MATMUL(Dmat, strain)

        Bmat = 0.0
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2
        
          Bmat(1,TI) = dNdx(ii); Bmat(1,TIp1) = 0.0;       Bmat(1,TIp2) = 0.0;
          Bmat(2,TI) = 0.0;       Bmat(2,TIp1) = dNdy(ii); Bmat(2,TIp2) = 0.0;
          Bmat(3,TI) = 0.0;       Bmat(3,TIp1) = 0.0;       Bmat(3,TIp2) = dNdz(ii);
          Bmat(4,TI) = dNdy(ii); Bmat(4,TIp1) = dNdx(ii); Bmat(4,TIp2) = 0.0;
          Bmat(5,TI) = 0.0;       Bmat(5,TIp1) = dNdz(ii); Bmat(5,TIp2) = dNdy(ii)
          Bmat(6,TI) = dNdz(ii); Bmat(6,TIp1) = 0.0;       Bmat(6,TIp2) = dNdx(ii)
        END DO

        BmatTrans = TRANSPOSE(Bmat)
        Bmat = MATMUL(Dmat, Bmat)
  
        Klocal = MATMUL(BmatTrans, Bmat)
        Klocal = dvol*Klocal

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3)
        END DO
      END DO

      END SUBROUTINE StiffnessResidualElasticityLinearTetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(12)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(12,12) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dN_dz(4)

      INTEGER, PARAMETER :: nGP=8, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=12

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts)
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25;      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*dens)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
          END DO
        END DO
      END DO

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearTetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearHexa(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(8)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(24)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(24,24) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(8), dNdx(8), dNdy(8), dN_dz(8)

      INTEGER, PARAMETER :: nGP=8, nlbf=8
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=24

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! call getQuadraturePointsHex(nGP, gpts1, gpts2, gpts3, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

    !     call computeBasisFunctions3D(.FALSE., 1, degree, param,
    !  1    xNode, yNode, zNode,
    !  2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*dens)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
          END DO
        END DO
      END DO

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearHexa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)     :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(6,6)   :: Dmat
      DOUBLE PRECISION, DIMENSION(12)    :: Flocal

      DOUBLE PRECISION, DIMENSION(6,12) :: Bmat
      DOUBLE PRECISION, DIMENSION(12,6) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: timefact, af, E, nu, bforce(3)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dN_dz(4)

      INTEGER, PARAMETER :: nGP=1, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)

      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = b1*(1.0-nu); Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = b1*(1.0-nu)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25
      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3
          b1 = valC(jj+1)
          b2 = valC(jj+2)
          b3 = valC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dN_dz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dN_dz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dN_dz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = grad(3,3)
        strain(4) = 0.5*(grad(1,2)+grad(2,1))
        strain(5) = 0.5*(grad(2,3)+grad(3,2))
        strain(6) = 0.5*(grad(1,3)+grad(3,1))

        stress = MATMUL(Dmat, strain)

        Bmat = 0.0
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2
        
          Bmat(1,TI) = dNdx(ii); Bmat(1,TIp1) = 0.0;       Bmat(1,TIp2) = 0.0;
          Bmat(2,TI) = 0.0;       Bmat(2,TIp1) = dNdy(ii); Bmat(2,TIp2) = 0.0;
          Bmat(3,TI) = 0.0;       Bmat(3,TIp1) = 0.0;       Bmat(3,TIp2) = dN_dz(ii);
          Bmat(4,TI) = dNdy(ii); Bmat(4,TIp1) = dNdx(ii); Bmat(4,TIp2) = 0.0;
          Bmat(5,TI) = 0.0;       Bmat(5,TIp1) = dN_dz(ii); Bmat(5,TIp2) = dNdy(ii)
          Bmat(6,TI) = dN_dz(ii); Bmat(6,TIp1) = 0.0;       Bmat(6,TIp2) = dNdx(ii)
        END DO

        BmatTrans = TRANSPOSE(Bmat)
  
        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b3 = dvol*dN_dz(ii)
          b4 = dvol*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3)

          Flocal(TI)   = Flocal(TI)   - ( b1*stress(1) + b2*stress(4) + b3*stress(6) )
          Flocal(TIp1) = Flocal(TIp1) - ( b1*stress(4) + b2*stress(2) + b3*stress(5) )
          Flocal(TIp2) = Flocal(TIp2) - ( b1*stress(6) + b2*stress(5) + b3*stress(3) )

        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearTetra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getElasticityMooneyRivlin3D(
     1   Elasticity,      ! Elasticity tensor d
     2   Mu, Lambda, kappa, c1, c2,   ! material parameters
     3   F, b, C,         ! left and right Cauchy-Green tensor 
     4   DetF, IC)            ! det of F

         IMPLICIT NONE

         DOUBLE PRECISION :: Mu, Lambda, kappa, c1, c2, DetF
         DOUBLE PRECISION :: Elasticity(9,9), b(3,3), F(3,3), C(3,3)
         DOUBLE PRECISION :: IC 

         INTEGER ::  ii, jj, kk, ll
         INTEGER ::  PP, QQ, RR, SS

         Elasticity = 0.0 

         ! right Cauchy-Green tensor C
         C = 0.0
         DO ii=1, 3
           DO jj=1, 3
             DO PP=1, 3
               C(ii,jj) = C(ii,jj) + F(PP,ii)*F(PP,jj)
             END DO
           END DO
         END DO

         ! the first variant of C  2D case
         IC = C(1,1) + C(2,2)+ C(3, 3)

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                IF (ii==jj .AND. ll==kk) THEN
                 Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1             Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2             + kappa*DetF
                END IF
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                IF (kk==jj .AND. ll==ii) THEN
                 Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1             Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2             - kappa*0.5*( DetF*DetF-1.0 )/DetF
                END IF
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                 DO PP=1, 3
                   DO QQ=1, 3
                     Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1               Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2               + 4.0*c2*(F(ii,PP)*F(jj,PP)*F(kk,QQ)*F(ll,QQ))/DetF
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                 DO PP=1, 3
                   IF (ii==kk) THEN
                     Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1               Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2               + 2.0*(c1+c2*IC)*(F(jj,PP)*F(ll,PP))/DetF
                   END IF
                 END DO
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                 DO PP=1, 3
                   DO QQ=1, 3
                     IF (ii==kk) THEN
                       Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1                  Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2                 - 2.0*c2*C(QQ,PP)*(F(jj,PP)*F(ll,QQ))/DetF
                     END IF
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                 DO PP=1, 3
                   DO QQ=1, 3
                     Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1               Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2               - 2.0*c2*(F(ii,QQ)*F(ll,QQ)*F(kk,PP)*F(jj,PP))/DetF
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

         DO ii=1, 3
           DO jj=1, 3
             DO kk=1, 3
               DO ll=1, 3
                 DO PP=1, 3
                   DO QQ=1, 3
                     Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1               Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2               - 2.0*c2*(F(ii,QQ)*F(kk,QQ)*F(jj,PP)*F(ll,PP))/DetF
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

      END SUBROUTINE getElasticityMooneyRivlin3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getElasticityNeoHookean3D(
     1   Elasticity,      ! Elasticity tensor d
     2   Mu, Lambda,      ! Lamé parameters
     3   F,b,             ! left Cauchy-Green tensor F F^T
     4   DetF)            ! det of F

         IMPLICIT NONE

         DOUBLE PRECISION :: Mu, Lambda, DetF
         DOUBLE PRECISION :: Elasticity(9,9), b(3,3), F(3,3)

         INTEGER ::  ii, jj

         Elasticity = 0.0 
            DO ii=1, 3
              DO jj=1, 3
                Elasticity(3*ii-2,3*jj-2) = Mu*b(ii,jj)/DetF
                Elasticity(3*ii-1,3*jj-1) = Mu*b(ii,jj)/DetF
                Elasticity(3*ii,3*jj) = Mu*b(ii,jj)/DetF
              END DO
            END DO

            DO ii=1, 3
              Elasticity(3*ii-2,ii) = Elasticity(3*ii-2,ii) + ( Mu-Lambda*Log(DetF) )/DetF
              Elasticity(3*ii-1,ii+3) = Elasticity(3*ii-1,ii+3) + ( Mu-Lambda*Log(DetF) )/DetF
              Elasticity(3*ii,ii+6) = Elasticity(3*ii,ii+6) + ( Mu-Lambda*Log(DetF) )/DetF
            END DO

            DO ii=1, 3
              DO jj=1, 3
                Elasticity(4*ii-3,4*jj-3) = Elasticity(4*ii-3,4*jj-3) + Lambda/DetF
              END DO
            END DO
      END SUBROUTINE getElasticityNeoHookean3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getElasticityLinear3D(
     1   Elasticity,      ! Elasticity tensor c
     2   Mu, Lambda)       

         IMPLICIT NONE

         DOUBLE PRECISION :: Mu, Lambda, DetF
         DOUBLE PRECISION :: Elasticity(9,9)

         INTEGER ::  ii, jj

         ! Elasticity tensor c
         Elasticity = 0.0
         Elasticity(2,4) = Mu
         Elasticity(3,7) = Mu
         Elasticity(4,2) = Mu
         Elasticity(6,8) = Mu
         Elasticity(7,3) = Mu
         Elasticity(8,6) = Mu
         DO ii=1, 3
             Elasticity(1,4*ii-3) = Lambda
             Elasticity(5,4*ii-3) = Lambda
             Elasticity(9,4*ii-3) = Lambda
             Elasticity(4*ii-3,4*ii-3) = Elasticity(4*ii-3,4*ii-3) + Mu
         END DO
         DO ii=1, 9
             Elasticity(ii,ii) = Elasticity(ii,ii) + Mu
         END DO
      END SUBROUTINE getElasticityLinear3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE StiffnessResidualElasticityLinearHex(
     1   xNode, yNode, zNode,  ! X and Y coordinates of the 8 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispElem,        ! displacement for the 4 nodes
     5   veloElem,        ! velocity for the 4 nodes
     6   acceElem,        ! acceleration for the 4 nodes
     7   Klocal, Flocal)  ! element stiffness and residual

         IMPLICIT NONE

         DOUBLE PRECISION, DIMENSION(:) :: elemData, timeData
         DOUBLE PRECISION :: xNode(8), yNode(8), zNode(8)
         DOUBLE PRECISION :: dispElem(24), veloElem(24), acceElem(24)
         DOUBLE PRECISION :: Klocal(24,24), Flocal(24), Elasticity(9,9)
         DOUBLE PRECISION :: Gc(3,9), dphi_dx(8), dphi_dy(8), dphi_dz(8)
         DOUBLE PRECISION :: Em, Nu, Mu, Lambda, kappa, thElem
         DOUBLE PRECISION :: c1, c2
         DOUBLE PRECISION :: b1, b2, b3, b4, F(3,3), Stress(9), Strain(9)
         DOUBLE PRECISION :: b(3,3), C(3,3), IC
         DOUBLE PRECISION :: FInv(3,3), DetF, dvol0, dvol
         DOUBLE PRECISION :: area, timefact, af, Jac, force
         DOUBLE PRECISION :: N(8), dN_dX(8), dN_dY(8), dN_dZ(8)

         INTEGER :: degree=1, nlbf=8, nGP=8, gp, ii, jj, pp, KK
        ! I think degree means order of element, like linear, quadratic interpolation
        ! I think nlbf is number of element nodes
   
         DOUBLE PRECISION :: gpts1(8), gpts2(8), gpts3(8), gwts(8)
         DOUBLE PRECISION :: param(3)
         logical :: finite
   
         finite = .TRUE. 
        !  finite = .FALSE. 
         ! material parameters, 
         Em = elemData(1);       Nu = elemData(2)
         kappa = 1000;
         ! Lamé parameters
         Mu = Em/(2+2*Nu);       Lambda = Em*Nu/(1+Nu)/(1-2*Nu)
         c1 = 0.5*Mu;      c2 = -0.5*c1
         ! density and body force are zero

         ! thickness
         thElem = elemData(3)
         ! time integration parameters
         af   = timeData(2)
         timefact = timeData(3)
   
         ! Gauss point coordinates and weights
         ! This function is defined in elemutilitiesquadrature.F 
         call getGaussPointsHexa(ngp, gpts1, gpts2, gpts3, gwts)

         ! loop over Gauss points and compute element stiffness and force vectors
         Klocal = 0.0;      Flocal = 0.0
   
         DO gp=1,nGP
          !  WRITE(*,*) "Gauss point number = ", gp
           param(1) = gpts1(gp)
           param(2) = gpts2(gp)
           param(3) = gpts3(gp)

           ! This function is defined in elementutilitiesbasisfuncs.F 
           ! ETYPE == 8 for hexa elements
           call computeBasisFunctions3D(.FALSE., 8, degree, param, xNode, yNode, zNode, 
     1                                   N, dN_dX, dN_dY, dN_dZ, Jac)
           
          !  WRITE(*,*) "Jac = ", Jac

           IF(Jac < 0.0) THEN
             STOP " Negative Jacobian for the hexa element" 
           END IF
           
          ! Calculate the deformation gradient tensor F
           F = 0
           F(1,1) = 1.0;      F(2,2) = 1.0;      F(3,3) = 1.0

           IF(finite) THEN ! 3D case
             DO ii=1, nlbf
               F(1,1) = F(1,1) + dispElem(3*ii-2) * dN_dX(ii)
               F(1,2) = F(1,2) + dispElem(3*ii-2) * dN_dY(ii)
               F(1,3) = F(1,3) + dispElem(3*ii-2) * dN_dZ(ii)
               F(2,1) = F(2,1) + dispElem(3*ii-1) * dN_dX(ii)
               F(2,2) = F(2,2) + dispElem(3*ii-1) * dN_dY(ii)
               F(2,3) = F(2,3) + dispElem(3*ii-1) * dN_dZ(ii)
               F(3,1) = F(3,1) + dispElem(3*ii)   * dN_dX(ii)
               F(3,2) = F(3,2) + dispElem(3*ii)   * dN_dY(ii)
               F(3,3) = F(3,3) + dispElem(3*ii)   * dN_dZ(ii)
             END DO
           END IF

          ! Inverse of F  3D case
           DetF = -F(1,3)*F(2,2)*F(3,1)+F(1,2)*F(2,3)*F(3,1)+F(1,3)*F(2,1)*F(3,2) 
     1            -F(1,1)*F(2,3)*F(3,2)-F(1,2)*F(2,1)*F(3,3)+F(1,1)*F(2,2)*F(3,3)

           FInv(1,1) =  (F(2,2)*F(3,3)-F(2,3)*F(3,2))/DetF
           FInv(1,2) =  (F(1,3)*F(3,2)-F(1,2)*F(3,3))/DetF
           FInv(1,3) =  (F(1,2)*F(2,3)-F(1,3)*F(2,2))/DetF
           FInv(2,1) =  (F(2,3)*F(3,1)-F(2,1)*F(3,3))/DetF
           FInv(2,2) =  (F(1,1)*F(3,3)-F(1,3)*F(3,1))/DetF
           FInv(2,3) =  (F(1,3)*F(2,1)-F(1,1)*F(2,3))/DetF
           FInv(3,1) =  (F(2,1)*F(3,2)-F(2,2)*F(3,1))/DetF
           FInv(3,2) =  (F(1,2)*F(3,1)-F(1,1)*F(3,2))/DetF
           FInv(3,3) =  (F(1,1)*F(2,2)-F(1,2)*F(2,1))/DetF

          ! left Cauchy-Green tensor b=F F^T
           b = 0.0
           DO ii=1, 3
             DO jj=1, 3
               DO KK=1, 3
                b(ii,jj) =b(ii,jj) + F(ii,KK)*F(jj,KK)
               END DO
             END DO
           END DO

          ! Elasticity tensor dd
           IF(finite) THEN
             call getElasticityMooneyRivlin3D(Elasticity, Mu, Lambda, kappa, c1, c2, F, b, C, DetF, IC)
             ! Calculate stress {11,12,13,21,22,23,31,32,33}
             Stress = 0.0
             Stress(1) = 0.5*kappa*(DetF - 1.0/DetF) + 2.0*(c1+c2*IC)*b(1,1)/DetF 
     1                 - 2.0*c2*(b(1,1)*b(1,1)+b(1,2)*b(2,1)+b(1,3)*b(3,1))/DetF
             Stress(2) =                               2.0*(c1+c2*IC)*b(1,2)/DetF
     1                 - 2.0*c2*(b(1,1)*b(1,2)+b(1,2)*b(2,2)+b(1,3)*b(3,2))/DetF
             Stress(3) =                               2.0*(c1+c2*IC)*b(1,3)/DetF
     1                 - 2.0*c2*(b(1,1)*b(1,3)+b(1,2)*b(2,3)+b(1,3)*b(3,3))/DetF
             Stress(4) =                               2.0*(c1+c2*IC)*b(2,1)/DetF
     1                 - 2.0*c2*(b(2,1)*b(1,1)+b(2,2)*b(2,1)+b(2,3)*b(3,1))/DetF
             Stress(5) = 0.5*kappa*(DetF - 1.0/DetF) + 2.0*(c1+c2*IC)*b(2,2)/DetF 
     1                 - 2.0*c2*(b(2,1)*b(1,2)+b(2,2)*b(2,2)+b(2,3)*b(3,2))/DetF
             Stress(6) =                               2.0*(c1+c2*IC)*b(2,3)/DetF
     1                 - 2.0*c2*(b(2,1)*b(1,3)+b(2,2)*b(2,3)+b(2,3)*b(3,3))/DetF
             Stress(7) =                               2.0*(c1+c2*IC)*b(3,1)/DetF
     1                 - 2.0*c2*(b(3,1)*b(1,1)+b(3,2)*b(2,1)+b(3,3)*b(3,1))/DetF
             Stress(8) =                               2.0*(c1+c2*IC)*b(3,2)/DetF
     1                 - 2.0*c2*(b(3,1)*b(1,2)+b(3,2)*b(2,2)+b(3,3)*b(3,2))/DetF
             Stress(9) = 0.5*kappa*(DetF - 1.0/DetF) + 2.0*(c1+c2*IC)*b(3,3)/DetF 
     1                 - 2.0*c2*(b(3,1)*b(1,3)+b(3,2)*b(2,3)+b(3,3)*b(3,3))/DetF
           else
             call getElasticityLinear3D(Elasticity, Mu, Lambda)
            ! Calculate linear strain {11,12,13,21,22,23,31,32,33}, I multiply 2.0 on each component
             Strain = 0.0
   
             Strain(1) = F(1,1) - 1.0 
             Strain(2) = 0.5*( F(1,2)+ F(2,1) ) 
             Strain(3) = 0.5*( F(1,3)+ F(3,1) ) 
             Strain(4) = Strain(2) 
             Strain(5) = F(2,2) - 1.0 
             Strain(6) = 0.5*( F(2,3)+ F(3,2) ) 
             Strain(7) = Strain(3) 
             Strain(8) = Strain(6)
             Strain(9) = F(3,3) - 1.0 
   
            ! Calculate stress {11,12,13,21,22,23,31,32,33}
             Stress = 0.0
             DO ii=1, 9
              Stress(ii) = 2*Mu*Strain(ii)
             END DO
             Stress(1) = Stress(1) + Lambda*(Strain(1)+Strain(5)+Strain(9))
             Stress(5) = Stress(5) + Lambda*(Strain(1)+Strain(5)+Strain(9))
             Stress(9) = Stress(9) + Lambda*(Strain(1)+Strain(5)+Strain(9))
           END IF
           
          ! Calculate derivative to the current coordinate denoted by dphi_dx, dphi_dy and dphi_dz
           dphi_dx = 0.0;     dphi_dy = 0.0;     dphi_dz = 0.0
           DO ii=1, nlbf
           dphi_dx(ii) = dN_dX(ii)*FInv(1,1) + dN_dY(ii)*FInv(2,1) + dN_dZ(ii)*FInv(3,1)
           dphi_dy(ii) = dN_dX(ii)*FInv(1,2) + dN_dY(ii)*FInv(2,2) + dN_dZ(ii)*FInv(3,2)
           dphi_dz(ii) = dN_dX(ii)*FInv(1,3) + dN_dY(ii)*FInv(2,3) + dN_dZ(ii)*FInv(3,3)
           END DO

          ! volume in the current configuration
           dvol0 = gwts(gp) * Jac;
           dvol = dvol0 * DetF
           !WRITE(*,*) "dvol = ", Jac, gwts(gp), dvol

           Gc = 0.0
           Do ii=1, nlbf
             b1 = dphi_dx(ii)*dvol;
             b2 = dphi_dy(ii)*dvol;
             b3 = dphi_dz(ii)*dvol;
             b4 = N(ii)*dvol;
             Do pp=1, 9
              Gc(1,pp) = b1*Elasticity(1,pp) + b2*Elasticity(4,pp) + b3*Elasticity(7,pp)
              Gc(2,pp) = b1*Elasticity(2,pp) + b2*Elasticity(5,pp) + b3*Elasticity(8,pp)
              Gc(3,pp) = b1*Elasticity(3,pp) + b2*Elasticity(6,pp) + b3*Elasticity(9,pp)
             END DO

             ! rewrite it with G^T \sigma to calculate F_int in eq.(1.31)
             Flocal(3*ii-2) = Flocal(3*ii-2) - b1*Stress(1) - b2*Stress(2) - b3*Stress(3)
             Flocal(3*ii-1) = Flocal(3*ii-1) - b1*Stress(4) - b2*Stress(5) - b3*Stress(6)
             Flocal(3*ii)   = Flocal(3*ii)   - b1*Stress(7) - b2*Stress(8) - b3*Stress(9)
   
             DO jj=1, nlbf
               Klocal(3*ii-2, 3*jj-2) = Klocal(3*ii-2, 3*jj-2) 
     1                  + Gc(1,1)*dphi_dx(jj) + Gc(1,4)*dphi_dy(jj) + Gc(1,7)*dphi_dz(jj)
               Klocal(3*ii-2, 3*jj-1) = Klocal(3*ii-2, 3*jj-1)
     1                  + Gc(1,2)*dphi_dx(jj) + Gc(1,5)*dphi_dy(jj) + Gc(1,8)*dphi_dz(jj)
               Klocal(3*ii-2, 3*jj)   = Klocal(3*ii-2, 3*jj)  
     1                  + Gc(1,3)*dphi_dx(jj) + Gc(1,6)*dphi_dy(jj) + Gc(1,9)*dphi_dz(jj)

               Klocal(3*ii-1, 3*jj-2) = Klocal(3*ii-1, 3*jj-2)
     1                  + Gc(2,1)*dphi_dx(jj) + Gc(2,4)*dphi_dy(jj) + Gc(2,7)*dphi_dz(jj)
               Klocal(3*ii-1, 3*jj-1) = Klocal(3*ii-1, 3*jj-1)
     1                  + Gc(2,2)*dphi_dx(jj) + Gc(2,5)*dphi_dy(jj) + Gc(2,8)*dphi_dz(jj)
               Klocal(3*ii-1, 3*jj)   = Klocal(3*ii-1, 3*jj)  
     1                  + Gc(2,3)*dphi_dx(jj) + Gc(2,6)*dphi_dy(jj) + Gc(2,9)*dphi_dz(jj)

               Klocal(3*ii, 3*jj-2) = Klocal(3*ii, 3*jj-2)
     1                  + Gc(3,1)*dphi_dx(jj) + Gc(3,4)*dphi_dy(jj) + Gc(3,7)*dphi_dz(jj)
               Klocal(3*ii, 3*jj-1) = Klocal(3*ii, 3*jj-1)
     1                  + Gc(3,2)*dphi_dx(jj) + Gc(3,5)*dphi_dy(jj) + Gc(3,8)*dphi_dz(jj)
               Klocal(3*ii, 3*jj)   = Klocal(3*ii, 3*jj)  
     1                  + Gc(3,3)*dphi_dx(jj) + Gc(3,6)*dphi_dy(jj) + Gc(3,9)*dphi_dz(jj)
             END DO
           END DO
         END DO

       !   WRITE(*,*) Klocal

      END SUBROUTINE StiffnessResidualElasticityLinearHex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE ElementUtilitiesElasticity3D
