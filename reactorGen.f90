MODULE define_IVP

  IMPLICIT NONE
  ! DX assigned a value in program that USEs define_IVP

  ! Physical constants
  real(8),parameter :: R = 8.3144621d0,          &  ! Universal gas constant (J.mol^-1.K^-1 = kg.m^2.s^-2.mol^-1.K^-1)
                       Avogadro = 6.02214129d23, &  ! Avogadro constant
                       pi = 3.14159265359,       &  ! pi
                       conc_min = 1.d-30,        &  ! Treshold concentration
                       infinity = huge(1.d0)

  integer :: neqn !, indxCG

  DOUBLE PRECISION :: advectionRate, dx, scaleFactor
  real(8), allocatable, dimension(:) :: gasFlux, photonFlux !,molm,cd,if0

  ! Sparse stoechiometry matrices and reaction rates
  character*10, allocatable :: speciesList(:)
  integer              :: nbSpecies, nbReac, nbPhotoReac,&
                          spectrumRange(1:2)=(/ 50, 200 /),&
                          sp1, sp2
  integer, allocatable :: D(:,:), L(:,:), LL(:,:)
  integer, allocatable :: Dphoto(:,:), Lphoto(:,:)
  real*8 , allocatable :: reactionRates(:), Dk(:),&
                          crossSections(:,:), photoRates(:), yloc(:)
  real*8, allocatable  :: v(:), ly(:)
  real*8, allocatable :: absorb(:,:), sumabs(:), intabs(:)

CONTAINS

  subroutine THRESH(Y)
    integer :: i
    real*8  :: Y(:)
    do i = 1, size(Y)
       Y(i) = max(Y(i),conc_min)
    enddo
  END subroutine THRESH

  SUBROUTINE F_E(NEQN,T,Y,DY)
    INTEGER, intent(in)  :: NEQN
    real*8, intent(in)   :: T, Y(NEQN)
    real*8, intent(out)  :: DY(NEQN)

    call thresh(Y)

    ! Restore real concentrations
    yloc = y * scaleFactor
    DY = 0

    ! Photolysis
    if (nbPhotoReac /= 0) call  PHOTO_KINET_SPARSE(yloc,DY)

    ! Transport
    DY = DY + ( gasFlux - advectionRate * yloc) / dx ! molec.cm^-3.s^-1
    !+ cd*(conc0-conc(1:))/dx**2/2

    ! Back to scaled concentrations
    dy = dy / scaleFactor

    !print *, 'F_E:',real(DY(1:10))
   END SUBROUTINE F_E

  DOUBLE PRECISION FUNCTION SR(NEQN,T,Y)
    INTEGER, intent(in) :: NEQN
    real*8, intent(in)  :: T, Y(NEQN)

    SR = 4d0 / DX**2

  END FUNCTION SR

  SUBROUTINE PHOTO_KINET_SPARSE(y,dy)
    real*8, intent(in)    :: y(:)
    real*8, intent(inout) :: dy(:)
    integer               :: imaxD, imaxL, ii, jj, i, j, k
    real*8                :: vsum

    imaxD  = size(Dphoto,1)
    imaxL  = size(Lphoto,1)

    absorb = 0d0
    do i = 1, imaxL
      ii = Lphoto(i,1)
      absorb(ii,sp1:sp2) = crossSections(ii,sp1:sp2) * y(Lphoto(i,2))
    enddo
    do k = sp1, sp2
      sumabs(k)=sum(absorb(1:nbPhotoReac,k))
    end do
    intabs = photonFlux * (1d0-exp(-sumabs*dx))/dx ! ph.cm^-3.s^-1
    do k = sp1, sp2
        if(sumabs(k) == 0) cycle
        intabs(k) = intabs(k)/sumabs(k)
    end do

    photoRates = 0
    do i = 1, imaxL
      ii = Lphoto(i,1)
      vsum=dot_product(intabs,absorb(ii,:))
      photoRates(ii) = photoRates(ii) + vsum
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

    real*8, intent(in)  :: y(:)
    real*8, intent(out) :: dy(:)
    logical, intent(in) :: want_jac
    real*8, optional,&
            intent(out) :: jac(:,:)

    integer             :: imaxD, imaxL, imaxLL, ii, jj, i, j
    real*8              :: tmp

    imaxD  = size(D,1)
    imaxL  = size(L,1)

    ! v_Dense = L_Sparse %*% ly_Dense
    ly = log(y)
    v = 0
    do i = 1, imaxL
      ii = L(i,1)
      v(ii) = v(ii) + L(i,3)*ly(L(i,2))
    enddo
    v = exp(v)

    ! dy_Dense = Dk_Sparse %*% v_Dense
    dy = 0
    do i = 1, imaxD
       ii = D(i,2)
       dy(ii) = dy(ii) + Dk(i)*v(D(i,1))
    enddo

    if(want_jac) then ! Compute jacobian
      jac = 0
      do j = 1, size(y)
        ! Contract L matrix

        ! 1- select reactions involving species j
        v = 0
        do i = 1, imaxL
          if (L(i,2) == j) then
            ii = L(i,1)
            v(ii) = v(ii) + 1
          endif
        enddo

        ! 2- reduce L matrix according to v/=0
        LL = 0
        jj = 0
        do i = 1, imaxL
           ii = L(i,1)
           if (v(ii) == 0) cycle
           jj = jj + 1
           LL(jj,:) = L(i,:)
        enddo
        imaxLL = jj

        ! v_Dense = LL_Sparse %*% ly_Dense
        v = 0
        do i = 1, imaxLL
          ii  = LL(i,1)
          jj  = LL(i,2)
          tmp = LL(i,3)
          if(v(ii) == 0) v(ii) = 1 ! Leave to 0 unconcerned reactions
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

  SUBROUTINE F_I(GRID_POINT,NPDES,T,YG,DYG,WANT_JAC,JAC)
    INTEGER, intent(in) :: GRID_POINT,NPDES
    LOGICAL, intent(in) :: WANT_JAC
    real*8, intent(in)  :: T, YG(NPDES)
    real*8, intent(out) :: DYG(NPDES), JAC(NPDES,NPDES)

    if (nbReac == 0) then
      DYG=0
      if(WANT_JAC) JAC = 0

    else
      call thresh(YG)
      yloc = YG * scaleFactor
      call  KINET_SPARSE(yloc,DYG,WANT_JAC,JAC)
      dyg =  dyg / scaleFactor
    endif
!     print *, 'F_I:',real(DYG); stop

  END SUBROUTINE F_I

END MODULE define_IVP
PROGRAM REACTOR

  USE define_IVP
  USE IRKC_M

!  IMPLICIT NONE
  ! Control parameters: Reactor geometry; T,P conditions; gas mixture...
  character*20                :: runId, beamSpectrumFile
  character*10, dimension(10) :: reactantsSpecies
  real*8,       dimension(10) :: reactantsComposition
  logical                     :: ifRestart=.FALSE., debug=.FALSE.
  integer                     :: nbSnapshots
  real*8                      :: reactorLength, reactorSection,&
                                 beamSection, reactantsFlux,&
                                 gasTemperature, electronsTemperature,&
                                 totalPressure, reactantsPressure,&
                                 reactionTime, beamIntensity

  namelist /REAC_DATA/ runId, debug, ifRestart, nbSnapshots,&
                       beamSpectrumFile,&
                       beamIntensity,&                        ! ph.cm^-2.s^-1
                       spectrumRange,&                        ! nm
                       reactorLength,&                        ! cm
                       reactorSection, beamSection,&          ! cm^2
                       reactantsFlux,&                        ! sccm
                       gasTemperature, electronsTemperature,& ! K
                       totalPressure, reactantsPressure,&     ! Pa
                       reactionTime,&                         ! s
                       reactantsSpecies,reactantsComposition


  ! Integration parameters
  TYPE(IRKC_SOL) :: SOL
  DOUBLE PRECISION,allocatable :: TOUT(:)
  DOUBLE PRECISION    :: T0, TEND, DTOUT, TEPS, tau

  ! Misc variables
  integer :: ngrid, i, next
  real*8  :: pToConc
  real*8, allocatable       :: Y0(:), initialConcentrations(:)
  real*8, allocatable       :: speciesMass(:)
  integer, allocatable      :: speciesCharge(:)


  ! Get control data ==================================================
  open(10,file = 'control.dat', status = 'old')
  read(10, nml = REAC_DATA)
  close(10)
  write(*,REAC_DATA)

  ! Initializations ==================================================
  call initGeometry
  if (debug) print *, 'initGeometry done...'
  call initGasMixture
  if (debug) print *, 'initGasMixture done...'
  call initFlow
  if (debug) print *, 'initFlow done...'
  call initPhotoChemistry
  if (debug) print *, 'initPhotoChemistry done...'
  call initChemistry
  if (debug) print *, 'initChemistry done...'

  ! Integration time and log-scaled snapshots time ===================
  allocate(TOUT(nbSnapshots))
  T0   = 0D0
  TEND = reactionTime
  TEPS = 1e-14*TEND   ! Shortest output time
  DTOUT = log(TEND/TEPS) / (nbSnapshots-2)
  TOUT(1)=T0
  DO I= 2, nbSnapshots
     TOUT(I)= TEPS*exp((I-2)*DTOUT)
  ENDDO

  ! Open results file
  open(10,file='fracmol_out.dat')
  write(10,'(500a15)') 'Time', speciesList

  ! ==================================================================
  ! Integration
  neqn = ngrid * nbSpecies
  allocate( Y0(neqn), yloc(neqn) )
  allocate( v(size(reactionRates)), ly(neqn) )

  ! Scale concentrations for numerical accurarcy
  Y0 = initialConcentrations/scaleFactor

  NEXT = 1
  WRITE(10,'(500E15.6e3)') TOUT(NEXT), Y0*scaleFactor

  DO NEXT = 2, nbSnapShots
    if(debug) print *,'Step #',NEXT,'/',nbSnapShots,'++++++++++++++++'
    SOL = IRKC_SET(TOUT(NEXT-1),Y0,TOUT(NEXT),&
                   AE=1D-12, RE=1D-8,&
                   NPDES=nbSpecies,&
                   ONE_STEP=.FALSE.)
    CALL IRKC(SOL,F_E,F_I,SR)
    Y0 = SOL%Y
    call thresh(Y0)
    WRITE(10,'(500E15.6e3)') TOUT(NEXT), Y0*scaleFactor
  !     if(debug) CALL IRKC_STATS(SOL)
  END DO

  CLOSE(10)

  CONTAINS

    SUBROUTINE initGeometry

      ! Grid geometry : box model
      dx = reactorLength
      print *,'Reacteur *************************'
      print *,'Longueur  / cm   :',real(reactorLength)
      print *,'Section   / cm^2 :',real(reactorSection)
      print *,'-> Volume / cm^3 :',real(reactorSection*reactorLength)
      print *,'**********************************'

      ngrid = 1

    END SUBROUTINE initGeometry

    SUBROUTINE initGasMixture
      character*10 :: species
      real*8       :: dummy
      integer      :: indx 

      ! Info about chemical species
      open(20, file = 'species_aux.dat', status = 'old')
      read(20,*) nbSpecies
      allocate(speciesList(nbSpecies), &
               speciesMass(nbSpecies), &
               speciesCharge(nbSpecies))
      read(20,*) speciesList
      read(20,*) speciesMass
      read(20,*) speciesCharge
      close(20)
      ! Molar masses
!       speciesMass = speciesMass*1d-3 ! ?????????????????????????????????????????????????

      ! Initial concentrations
      allocate(initialConcentrations(nbSpecies))

      ! From pressure to concentration : n/V = P/(RT)
      pToConc = reactantsPressure / (R*gasTemperature) & ! P/RT : mol.m^-3
                  * Avogadro / 1d6                         ! -> molec.cm^-3

      if (ifRestart) then
        ! Restart from previous state
        ! First number on line is time
        open(20, file="init_conc.dat", status="old")
        read(20, *) dummy, initialConcentrations
        close(20)

      else
        print *,'Gas Mixture **********************'
        ! Use initial composition and pressure to define conc. in molec.cm^-3
        initialConcentrations = 0
        do i= 1, size(reactantsSpecies)
          species = reactantsSpecies(i)
          if (species == '') exit
          indx = speciesIndex(species)
          initialConcentrations(indx) = reactantsComposition(i) * pToConc
          print *,trim(species),' / molec.cm^-3 :', real(initialConcentrations(indx))
          print *,trim(species),' / molec       :', real(initialConcentrations(indx)* &
                                                      reactorSection*reactorLength)
        end do
        print *,'**********************************'
!         initialConcentrations(indxCG) = totalPressure - reactantsPressure

      end if
      scaleFactor = maxval(initialConcentrations)

    END SUBROUTINE initGasMixture

    SUBROUTINE initFlow
      character*10 :: species
      real*8       :: reacTime, cFlux
      integer      :: indx

      reactantsFlux = 1.689d-3 * reactantsFlux ! Convert sccm to Pa.m^3/s
      cFlux = reactantsFlux *Avogadro / (R * gasTemperature * reactorSection) ! molec.cm^-2.s^-1

      ! Advection rate / cm.s^-1
      advectionRate = cFlux / pToConc
      tau = reactorLength/advectionRate

      ! Max. integration time
      if( reactionTime < 0 ) then
        ! Use characteristic time, if not infinite
        if (tau > infinity) stop 'Please set a reaction time'
        reactionTime = 5*tau
      endif

      print *,'Advection ************************'
      print *,'v_ad / cm.s^-1     =',real(advectionRate)
      print *,'tau  / s           =',real(tau)
      print *,'Integration time / s=',real(reactionTime)
      print *,'**********************************'

      ! Molecular flux
      allocate(gasFlux(nbSpecies))
      gasFlux = 0
      print *,'Gas flow *************************'
      print *,'qre / mol.s^-1  =',real(cFlux/Avogadro*reactorSection)
      do i= 1, size(reactantsSpecies)
        species = reactantsSpecies(i)
        if (species == '') exit
        indx = speciesIndex(species)
        gasFlux(indx) = cFlux * reactantsComposition(i)  ! molec.cm^-2.s^-1
        print *,trim(species),' / molec.cm^-2.s^-1 :', real(gasFlux(indx))
        print *,trim(species),' / molec.s^-1       :', real(gasFlux(indx)*reactorSection)
        print *,trim(species),' / molec over reactionTime :', real(gasFlux(indx)*&
                                                           reactorSection*reactionTime)
      end do
      print *,'**********************************'
!       gasFlux(indxCG) = cFlux * (totalPressure/reactantsPressure - 1d0)     ! Carrier Gas

    END SUBROUTINE initFlow

    SUBROUTINE initPhotoChemistry
      integer             :: imaxL, ii, jj, k, ios, nnz, nw
      integer             :: i, j, irate, imax, il
      real*8              :: lambda, s
      character(len=100)  :: words(100), fname
      character(len=3000) :: line

      nbPhotoReac = 0

       ! D and L sparse matrices'irradiation_flux.dat'
      open(20, file = "photo_DL.dat", status = "old", iostat=ios)

      if (ios /= 0) RETURN ! No photodissociations; exit

      read(20,*) nnz
      allocate(Dphoto(nnz,3))
      read(20,*) ((Dphoto(i,j),i=1,nnz),j=1,3)
      read(20,*) nnz
      allocate(Lphoto(nnz,3))
      read(20,*) ((Lphoto(i,j),i=1,nnz),j=1,3)
      close(20)
      nbPhotoReac = maxval(Dphoto(:,1))

      ! Build partial cross sections
      open(20, file = "photo_params.dat", status = "old")
      allocate(crossSections(nbPhotoReac,spectrumRange(1):spectrumRange(2)),&
               photoRates(nbPhotoReac) )
      crossSections = 0
      do irate = 1, nbPhotoReac
        read( 20, '(a)', iostat = ios) line

        if (ios /= 0) stop "Pb reading photo_params.dat"

        call split(adjustl(line), nw, words)

        read(words(1),*) fname ! File containing the cross-sections
        open(21,file='Photo/'//fname, status="old")
        do
          read(21,*,iostat=ios) lambda, s ! nm, cm^2.molec-1.nm^-1
          if( ios/=0 ) exit
          il = int(lambda+0.5)
          if( il < spectrumRange(1) ) cycle
          if( il > spectrumRange(2) ) exit
          crossSections(irate,il) = s
        enddo
        close(21)

        if( words(2) /= '') then ! Use branching ratios
          read(words(2),*) fname ! File containing the BRs
          open(21,file='Photo/'//fname, status="old")
          do
            read(21,*,iostat=ios) lambda, s ! nm, no unit
            if( ios/=0 ) exit
            il = int(lambda+0.5)
            if( il < spectrumRange(1) ) cycle
            if( il > spectrumRange(2) ) exit
           crossSections(irate,il) = crossSections(irate,il) * s
          enddo
          close(21)
        endif

      end do
      close(20)

      ! Irradiation flux
      allocate(photonFlux(spectrumRange(1):spectrumRange(2)))
      open(20,file=beamSpectrumFile,status="old")
      photonFlux = 0d0
      do
        read(20,*,iostat=ios) lambda, s ! nm, ph.cm^-2.s^-1.nm^-1
        if( ios/=0 ) exit
        il = int(lambda+0.5)
        if( il < spectrumRange(1) ) cycle
        if( il > spectrumRange(2) ) exit
        photonFlux(il) = s
      enddo
      close(20)
      open(55,file='photonFlux.out')
      write(55,*) real(photonFlux)
      close(55)

      ! Renormalization of spectrum
      if(beamIntensity >= 0) photonFlux = photonFlux/sum(photonFlux) * beamIntensity

      ! Redistribute beam intensity to reactorSection
      photonFlux = photonFlux * (beamSection/reactorSection)


      print *,'Photons **************************'
      print *,'flux  / ph.cm^-2.s^-1 =',real(sum(photonFlux))
      print *,'flux  / ph.s^-1       =',real(sum(photonFlux)*reactorSection)
      print *,'integ / ph            =',real(sum(photonFlux)*reactorSection*reactionTime)
      print *,'**********************************'
!       photonFlux=photonFlux/9.537**2 !*20.2 ! Titan
!       photonFlux=0d0;photonFlux(82)=6.7d13*1d4/Avogadro    ! irradiation at 82 nm in Imanaka and Smith's reactor

      print *,'Absorption ***********************'
      imaxL  = size(Lphoto,1)
      sp1 = spectrumRange(1)
      sp2 = spectrumRange(2)
      ! Allcocate tables also used in PHOTO_KINET_SPARSE
      allocate(absorb(nbPhotoReac,sp1:sp2), &
               sumabs(sp1:sp2), &
               intabs(sp1:sp2))
      absorb = 0d0
      do i = 1, imaxL
        ii = Lphoto(i,1)
        absorb(ii,sp1:sp2) = crossSections(ii,sp1:sp2) * &
                             initialConcentrations(Lphoto(i,2))
      enddo
      do k = sp1, sp2
        sumabs(k)=sum(absorb(1:nbPhotoReac,k))
      end do
      print *,'% absorbed ', (sum(1-exp(-sumabs*dx)))* 100 / (sp2-sp1+1)
      print *,'**********************************'

    END SUBROUTINE initPhotoChemistry

    SUBROUTINE initChemistry
      integer :: i, j, nbSpec, irate, imax, ios, nnz, nw
      real*8  :: tmp, k_para(20)
      character*10 :: type = 'kooij'
      character(len=100)  :: words(100), fname
      character(len=3000) :: line

      nbReac = 0

       ! D and L sparse matrices
      open(20, file = "reac_DL.dat", status = "old", iostat=ios)

      if (ios /= 0) RETURN ! No reactions; exit

      read(20,*) nnz

      if (nnz == 0) RETURN ! No reactions; exit

      allocate(D(nnz,3),Dk(nnz))
      read(20,*) ((D(i,j),i=1,nnz),j=1,3)
      read(20,*) nnz
      allocate(L(nnz,3),LL(nnz,3))
      read(20,*) ((L(i,j),i=1,nnz),j=1,3)
      close(20)
      nbReac = maxval(D(:,1))

      ! Reactions list
      open(20, file = "reac_params.dat", status = "old")
      allocate(reactionRates(nbReac))
      do irate = 1, nbReac
        read( 20, '(a)', iostat = ios ) line
        if (ios /= 0) stop "Pb reading reac_params.dat"
        call split(adjustl(line), nw, words)

        ! Get formula type
        type = words(nw)

        ! Get parameters
        select case (type)

          case ('3-body')
            imax = 10
            tmp  = gasTemperature

          case('dr')
            imax = 5
            tmp  = electronsTemperature

          case ('kooij','ionpol1','ionpol2')
            imax = 5
            tmp  = gasTemperature

          case default
            print *,'Unknown reaction type :'//type
            stop 'initChemistry'

        end select

        do i = 1, imax
          read(words(i),*) k_para(i)
        end do

        ! Calculate rate constant
        reactionRates(irate) = krate(type,k_para,tmp,pToConc)

      end do
      close(20)
      print *, real(pToConc)
      print *, real(minval(reactionRates)),real(maxval(reactionRates))
      open(55,file='rrates.out')
      write(55,*) real(reactionRates)
      close(55)

      ! Compute Dk to accelerate KINET_SPARSE
      ! Dk_Sparse = D_Sparse * rate_Dense
      do i = 1, size(D,1)
        Dk(i) = D(i,3) * reactionRates(D(i,1))
      enddo

    END SUBROUTINE initChemistry

    integer function speciesIndex(name_spec) result(no)
      character(len=*) :: name_spec

      no = 1
      do
        if (no > size(speciesList) ) then
          no = 0
          exit
        end if
        if (name_spec(1:len_trim(name_spec))== &
            speciesList(no)(1:len_trim(speciesList(no)))) exit
        no=no+1
      end do
      if (no == 0) then
        print *,"Species not in list: "//name_spec
        stop
      endif
    end function speciesIndex

    real*8 function krate(type,c,T,P) result(k)
      character(len=*)     :: type
      real*8, dimension(:) :: c
      real*8               :: T, F, k0, kInf, P
      real*8, parameter    :: T0 = 300   ! Reference temp. for kooij expression

      select case (type)
        case ('kooij', 'dr')
          k = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k = k * F

        case ('3-body')
          k0   = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k0 = k0 * F
          kInf = c(6) * (T/T0)**c(7) * exp(-c(8)/T)
          F = c(9) * exp( c(10)*abs(1/T-1/T0) )
          kInf = kInf * F
          k = troe(k0,kInf,P)
       
        case ('ionpol1')
          k = c(1) * c(2) * (0.62 + 0.4767*c(3)*(T0/T)**0.5)
          F = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k = k * F

        case ('ionpol2')
          k = c(1) * c(2) * (1 + 0.0967*c(3)*(T0/T)**0.5 &
                                  + c(3)**2 / 10.526 * T0/T   )
          F = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k = k * F

      end select

    end function krate

    real*8 function troe(k0,kInf,P) result(k)
	real*8               :: k0, kInf, P
	real*8               :: cE, nE, dE, fE, Pr, Fc, lFc, c1

        ! Simplified Troe formula
        Fc = 0.64

        lFc= log10(Fc)
        Pr = k0*P/kInf

        cE = -0.4 -0.67*lFc
        nE = 0.75 -1.27*lFc
        dE = 0.14
        c1 = log10(Pr)+cE
        fE = 1+(c1/(nE-dE*c1))**2

        k  = kInf * (Pr/(1+Pr)) * Fc**(1/fE)  

    end function troe


END PROGRAM REACTOR

subroutine split(chain,nw,words)
  IMPLICIT NONE
  integer,parameter :: lmax=5000
  integer :: nw,lw
  character(len=lmax) line
  character(len=*   ) chain, words(*)
  logical continue

  if(len(chain).gt.lmax) stop 'Split: Problem of string length!'
  line=chain
  nw=0
  continue=.true.
  if(len_trim(line).le.0) continue=.false.
  do while(continue)
    do while(line(1:1).eq.' ')
      line=line(2:)
    enddo
    lw=index(line,' ')-1
    nw=nw+1
    words(nw)=line(:lw)
    line=line(lw+1:)
    if(len_trim(line).le.0) continue=.false.
  enddo

end subroutine split





!   ! Diffusion coefficients in N2 - ref: PhD Eric H. Wilson
!   allocate(cd(nbSpecies))
!   no_n2=name_no("N2",speciesList(1:nspec),nspec)   ! serial number of N2
!   if (no_n2==0) stop "Data of N2 do not exist in spec_data.dat"
!   A_N2 = 5.09d18
!   S_N2 = 0.81d00
!   cDiff_N2=A_N2*gasTemperature**S_N2
!   mass_N2 = molm(no_n2)
!   do i=1,nspec
!     if (molm(i) < mass_N2) then
!       cd(i)=sqrt(5d-1*mass_N2/molm(i)+5d-1)
!     else
!       cd(i)=sqrt(mass_N2/molm(i))
!     end if
!   end do
!   cd=cd*cDiff_N2/Avogadro*1d5/totalPressure ! ?????????????
