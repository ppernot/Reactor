MODULE define_IVP

  IMPLICIT NONE
  ! DX assigned a value in program that Uses define_IVP

  ! Physical constants
  DOUBLE PRECISION,parameter :: R = 8.3144621d0, &  ! Universal gas constant
                                                    ! (J.mol^-1.K^-1 =
                                                    !  kg.m^2.s^-2.mol^-1.K^-1)
                                Avogadro = 6.02214129d23, &  ! Avogadro constant
                                pi = 3.14159265359, &  ! pi
                                conc_min = 0.D0,    &  ! Treshold concentration
                                infinity = huge(1.d0)

  integer :: neqn !, indxCG

  DOUBLE PRECISION                            :: advectionRate, dx, scaleFactor
  DOUBLE PRECISION, allocatable, dimension(:) :: gasFlux, photonFlux !,molm,cd,if0

  ! Sparse stoechiometry matrices and reaction rates
  character*10, allocatable      :: speciesList(:)
  integer                        :: nbSpecies, nbReac, nbPhotoReac,&
                                    spectrumRange(1:2)=(/ 50, 200 /),&
                                    sp1, sp2
  integer, allocatable           :: D(:,:), L(:,:), LL(:,:)
  integer, allocatable           :: Dphoto(:,:), Lphoto(:,:)
  DOUBLE PRECISION , allocatable :: reactionRates(:), Dk(:),&
                                    crossSections(:,:), photoRates(:), yloc(:)
  DOUBLE PRECISION, allocatable  :: v(:), ly(:)
  DOUBLE PRECISION, allocatable  :: absorb(:,:), sumabs(:), intabs(:)

CONTAINS

  subroutine THRESH(Y)
    integer           :: i
    DOUBLE PRECISION  :: Y(:)
    do i = 1, size(Y)
       Y(i) = max(Y(i),conc_min)
    enddo
  END subroutine THRESH

  SUBROUTINE F_E(NEQN,T,Y,DY)
    INTEGER, intent(in)            :: NEQN
    DOUBLE PRECISION, intent(in)   :: T, Y(NEQN)
    DOUBLE PRECISION, intent(out)  :: DY(NEQN)

    ! Restore real concentrations
    call thresh(Y)
    yloc = y * scaleFactor

    DY = 0D0

    ! Photolysis
    if (nbPhotoReac /= 0) call  PHOTO_KINET_SPARSE(yloc,DY)

    ! Transport
    DY = DY + ( gasFlux - advectionRate * yloc) / dx ! molec.cm^-3.s^-1
    !+ cd*(conc0-conc(1:))/dx**2/2
 
    ! Back to scaled concentrations
    dy = dy / scaleFactor

  END SUBROUTINE F_E

  SUBROUTINE F_I(GRID_POINT,NPDES,T,Y,DY,WANT_JAC,JAC)
    INTEGER, intent(in)           :: GRID_POINT,NPDES
    LOGICAL, intent(in)           :: WANT_JAC
    DOUBLE PRECISION, intent(in)  :: T, Y(NPDES)
    DOUBLE PRECISION, intent(out) :: DY(NPDES), JAC(NPDES,NPDES)

    if (nbReac == 0) then
      DY = 0D0
      if (WANT_JAC) JAC = 0D0
    else
      call thresh(Y)
      yloc = Y * scaleFactor
      call  KINET_SPARSE(yloc,DY,WANT_JAC,JAC)
      dy  =  dy / scaleFactor
    endif

  END SUBROUTINE F_I

  !DOUBLE PRECISION FUNCTION SR(NEQN,T,Y)
  !  INTEGER, intent(in) :: NEQN
  !  DOUBLE PRECISION, intent(in)  :: T, Y(NEQN)
  !  SR = 4d0 / DX**2
  !END FUNCTION SR

  SUBROUTINE PHOTO_KINET_SPARSE(y,dy)
    DOUBLE PRECISION, intent(in)    :: y(:)
    DOUBLE PRECISION, intent(inout) :: dy(:)
    integer                         :: imaxD, imaxL, ii, jj, i, j, k

    imaxD  = size(Dphoto,1)
    imaxL  = size(Lphoto,1)

    absorb = 0D0
    do i = 1, imaxL
      ii = Lphoto(i,1)
      absorb(ii,sp1:sp2) = crossSections(ii,sp1:sp2) * y(Lphoto(i,2))
    enddo
    do k = sp1, sp2
      sumabs(k) = sum(absorb(1:nbPhotoReac,k))
    end do
    intabs = photonFlux * (1d0-exp(-sumabs*dx))/dx ! ph.cm^-3.s^-1
    do k = sp1, sp2
        if(sumabs(k) == 0D0) cycle
        intabs(k) = intabs(k) / sumabs(k)
    end do

    photoRates = 0D0
    do i = 1, imaxL
      ii = Lphoto(i,1)
      photoRates(ii) = photoRates(ii)+dot_product(intabs,absorb(ii,:))
    enddo

    do i = 1, imaxD
       ii = Dphoto(i,2)
       dy(ii) = dy(ii) + Dphoto(i,3) * photoRates(Dphoto(i,1))
    enddo

  END SUBROUTINE PHOTO_KINET_SPARSE

  SUBROUTINE KINET_SPARSE(y,dy,want_jac,jac)
  ! Adapted from Cangiani2012 : Biochemical Pathways Simulation
  ! Dk = t(D*k) (done in calling program)
  ! dC = Dk %*% apply(parms$L,1,function(x) prod(y^x))
  ! NB: Dense matrices dimensions : reactions x species
  ! NB: Description of sparse matrices D & L as triplets (reac:i,species:j,val)

    DOUBLE PRECISION, intent(in)           :: y(:)
    DOUBLE PRECISION, intent(out)          :: dy(:)
    logical, intent(in)                    :: want_jac
    DOUBLE PRECISION, optional,intent(out) :: jac(:,:)
    integer                                :: imaxD, imaxL, imaxLL, &
                                              ii, jj, i, j
    DOUBLE PRECISION                       :: tmp

    imaxD  = size(D,1)
    imaxL  = size(L,1)

    ! v_Dense = L_Sparse %*% ly_Dense
    ly = log(y)
    v = 0D0
    do i = 1, imaxL
      ii = L(i,1)
      v(ii) = v(ii) + L(i,3)*ly(L(i,2))
    enddo
    v = exp(v)

    ! dy_Dense = Dk_Sparse %*% v_Dense
    dy = 0D0
    do i = 1, imaxD
       ii = D(i,2)
       dy(ii) = dy(ii) + Dk(i)*v(D(i,1))
    enddo

    if(want_jac) then ! Compute jacobian
      jac = 0D0
      do j = 1, size(y)
        ! Contract L matrix

        ! 1- select reactions involving species j
        v = 0D0
        do i = 1, imaxL
          if (L(i,2) == j) then
            ii = L(i,1)
            v(ii) = v(ii) + 1D0
          endif
        enddo

        ! 2- reduce L matrix according to v/=0
        LL = 0
        jj = 0
        do i = 1, imaxL
           ii = L(i,1)
           if (v(ii) == 0D0) cycle
           jj = jj + 1
           LL(jj,:) = L(i,:)
        enddo
        imaxLL = jj

        ! v_Dense = LL_Sparse %*% ly_Dense
        v = 0D0
        do i = 1, imaxLL
          ii  = LL(i,1)
          jj  = LL(i,2)
          tmp = LL(i,3)
          if(v(ii) == 0D0) v(ii) = 1D0 ! Leave unconcerned reactions to 0 
          if(jj == j) then
             v(ii) = v(ii) * tmp * y(jj)**(tmp-1)
          else
             v(ii) = v(ii) * y(jj)**tmp
          endif
        enddo

        ! jac_Dense = Dk_Sparse %*% v_Dense
        do i = 1, imaxD
           ii = D(i,2)
           jac(ii,j) = jac(ii,j) + Dk(i)*v(D(i,1))
        enddo

      enddo

    endif

  END SUBROUTINE KINET_SPARSE

END MODULE define_IVP
