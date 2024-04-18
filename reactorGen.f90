! v2.0
PROGRAM REACTOR

  USE define_IVP
  USE IRKC_M

  IMPLICIT NONE

  ! Declare the interface for POSIX fsync function
            interface
              function fsync (fd) bind(c,name="fsync")
              use iso_c_binding, only: c_int
                integer(c_int), value :: fd
                integer(c_int) :: fsync
              end function fsync
            end interface

  ! Control parameters: Reactor geometry; T,P conditions; gas mixture...
  character*20                    :: runId, beamSpectrumFile
  character*10, dimension(10)     :: reactantsSpecies
  DOUBLE PRECISION, dimension(10) :: reactantsComposition
  logical                         :: ifRestart=.FALSE., debug=.FALSE.
  integer                         :: nbSnapshots
  DOUBLE PRECISION                :: reactorLength,  reactorSection,       &
                                     beamSection,    reactantsFlux,        &
                                     gasTemperature, electronsTemperature, &
                                     totalPressure,  reactantsPressure,    &
                                     reactionTime,   beamIntensity,        &
                                     relativeError = 1D-4,                 &
                                     absoluteError = 1D-4,                 &
                                     TEPS = 1D-10,                         &
                                     wallFactor = 1D0

  namelist /REAC_DATA/ runId, debug, ifRestart, nbSnapshots,&
                       beamSpectrumFile,&
                       beamIntensity,&                        ! ph.cm^-2.s^-1
                       spectrumRange,&                        ! nm
                       spectralResolution,&                   ! nm
                       reactorLength,&                        ! cm
                       reactorSection, beamSection,&          ! cm^2
                       reactantsFlux,&                        ! sccm
                       gasTemperature, electronsTemperature,& ! K
                       totalPressure, reactantsPressure,&     ! Pa
                       reactionTime,&                         ! s
                       reactantsSpecies, reactantsComposition, &
                       relativeError, absoluteError, &
                       useSR, SRmax, &
                       wallFactor


  ! Integration parameters
  TYPE(IRKC_SOL) :: SOL
  DOUBLE PRECISION,allocatable  :: TOUT(:)
  DOUBLE PRECISION              :: T0, TEND, DTOUT, tau
  ! Misc variables
  integer                       :: ngrid, i, next, ret
  DOUBLE PRECISION              :: pToConc
  DOUBLE PRECISION, allocatable :: Y0(:), initialConcentrations(:)
  DOUBLE PRECISION, allocatable :: speciesMass(:)
  integer, allocatable          :: speciesCharge(:)


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
  T0      = 0D0
  TEND    = reactionTime
  DTOUT   = log(TEND/TEPS) / (nbSnapshots-2)
  TOUT(1) = T0
  DO I= 2, nbSnapshots
     TOUT(I)= TEPS * exp((I-2)*DTOUT)
  ENDDO

  ! Open results file
  open(10,file='fracmol_out.dat')
  write(10,'(500a15)') 'Time', speciesList

  ! Statistics file
  open(11,file='integ_stats.out')

  ! ==================================================================
  ! Integration
  neqn = ngrid * nbSpecies
  allocate( Y0(neqn), yloc(neqn) )
  allocate( v(size(reactionRates)), ly(neqn) )

  ! Scale concentrations for numerical accurarcy
  ! Defined in initGasMixture
  Y0 = initialConcentrations / scaleFactor

  NEXT = 1
  ! Output unscaled results
  WRITE(10,*) TOUT(NEXT), Y0 * scaleFactor

  DO NEXT = 2, nbSnapshots

    if(debug) then
      print *,'Step #',NEXT,'/',nbSnapShots
      FLUSH(6) ! For graphical interface
    end if

    SOL = IRKC_SET(TOUT(NEXT-1), Y0, TOUT(NEXT) , &
                   AE            = absoluteError, &
                   RE            = relativeError, &
                   NPDES         = nbSpecies    , &
                   ONE_STEP      = .FALSE.      , &
                   STOP_ON_ERROR = .TRUE. )

    if(useSR) then   
      CALL IRKC(SOL,F_E,F_I,SR)
    else
      CALL IRKC(SOL,F_E,F_I)
    endif 

    Y0 = SOL%Y

    WRITE(10,*) TOUT(NEXT), Y0 * scaleFactor

    FLUSH(10) ! For graphical interface
    ret = fsync(fnum(10))
    if (ret /= 0) stop "Error calling FSYNC on file 10"

    ! Integration stats
    WRITE(11,*) SOL%T, SOL%NFE, SOL%NFI, SOL%NSTEPS, SOL%NACCPT, &
                SOL%NREJCT, SOL%NFESIG, SOL%MAXM, SOL%SPRAD
    FLUSH(11) ! For graphical interface
    ret = fsync(fnum(11))
    if (ret /= 0) stop "Error calling FSYNC on file 11"

  END DO

  if(debug) then
      print *,'*** Done ***'
      FLUSH(6) ! For graphical interface
  end if

  CLOSE(10)
  CLOSE(11)

  ! Output Photolysis final rates
  open(55,file='phrates.out')
  write(55,*) photoRates
  close(55)

  CONTAINS

    SUBROUTINE initGeometry
      IMPLICIT NONE

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
      IMPLICIT NONE
      character*10     :: species
      DOUBLE PRECISION :: dummy
      integer          :: indx 

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
!       speciesMass = speciesMass*1d-3 ! ???????????????????????????

      ! Initial concentrations
      allocate(initialConcentrations(nbSpecies))

      ! From pressure to concentration : n/V = P/(RT)
      pToConc = reactantsPressure / (R*gasTemperature) & ! P/RT : mol.m^-3
                  * Avogadro / 1D6                         ! -> molec.cm^-3

      if (ifRestart) then
        ! Restart from previous state
        ! First number on line is time
        open(20, file="init_conc.dat", status="old")
        read(20, *) dummy, initialConcentrations
        close(20)

      else
        print *,'Gas Mixture **********************'
        ! Use initial composition and pressure to define conc. in molec.cm^-3
        initialConcentrations = 0D0
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
      IMPLICIT NONE
      character*10     :: species
      DOUBLE PRECISION :: reacTime, cFlux
      integer          :: indx

      reactantsFlux = 1.689d-3 * reactantsFlux ! Convert sccm to Pa.m^3/s
      cFlux = reactantsFlux * Avogadro / (R * gasTemperature * reactorSection) ! molec.cm^-2.s^-1

      ! Advection rate / cm.s^-1
      advectionRate = cFlux / pToConc
      tau = reactorLength / advectionRate

      ! Max. integration time
      if( reactionTime < 0D0 ) then
        ! Use characteristic time, if not infinite
        if (tau > infinity) stop 'Please set a reaction time'
        reactionTime = 5D0 * tau
      endif

      print *,'Advection ************************'
      print *,'v_ad / cm.s^-1       :',real(advectionRate)
      print *,'tau  / s             :',real(tau)
      print *,'Integration time / s :',real(reactionTime)
      print *,'**********************************'

      ! Molecular flux
      allocate(gasFlux(nbSpecies))
      gasFlux = 0D0
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
      IMPLICIT NONE
      integer             :: imaxL, ii, jj, k, ios, nnz, nw
      integer             :: i, j, irate, imax, il
      DOUBLE PRECISION    :: lambda, s, qy, spF
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

      ! Vector limits
      sp1 = int(spectrumRange(1)/spectralResolution)
      sp2 = int(spectrumRange(2)/spectralResolution)

      ! Build partial cross sections
      open(20, file = "photo_params.dat", status = "old")
      allocate(crossSections(nbPhotoReac,sp1:sp2),&
               photoRates(nbPhotoReac) )
      crossSections = 0D0
      do irate = 1, nbPhotoReac
        read( 20, '(a)', iostat = ios) line

        if (ios /= 0) stop "Pb reading photo_params.dat"

        call split(adjustl(line), nw, words)

        read(words(1),*) fname ! File containing the cross-sections
        open(21,file='Photo/'//fname, status="old")
        do
          read(21,*,iostat=ios) lambda, s ! nm, cm^2
          if( ios/=0 ) exit
          il = int( (lambda+0.5*spectralResolution)/ &
                       spectralResolution             )
          if( il < sp1 ) cycle
          if( il > sp2 ) exit
          crossSections(irate,il) = s
        enddo
        close(21)

        if( words(2) /= '') then ! Use branching ratios / quantum yields
          read(words(2),*) fname ! File containing the BRs
          open(21,file='Photo/'//fname, status="old")
          do
            read(21,*,iostat=ios) lambda, qy ! nm, no unit
            if( ios/=0 ) exit
            il = int( (lambda+0.5*spectralResolution)/ &
                       spectralResolution               )
            if( il < sp1 ) cycle
            if( il > sp2 ) exit
            crossSections(irate,il) = crossSections(irate,il) * qy
          enddo
          close(21)
        endif

      end do
      close(20)

      ! Irradiation flux
      allocate(photonFlux(sp1:sp2))
      open(20,file=beamSpectrumFile,status="old")
      photonFlux = 0D0
      do
        read(20,*,iostat=ios) lambda, s ! nm, ph.cm^-2.s^-1.nm^-1
        if( ios/=0 ) exit
        il = int( (lambda+0.5*spectralResolution)/ &
                   spectralResolution                )
        if( il < sp1 ) cycle
        if( il > sp2 ) exit
        photonFlux(il) = s
      enddo
      close(20)
      open(55,file='photonFlux.out')
      do il = sp1, sp2
        write(55,*) real(il*spectralResolution),real(photonFlux(il))
      enddo
      close(55)

      ! Renormalization of spectrum
      if(beamIntensity >= 0) photonFlux = photonFlux /&
                                          (sum(photonFlux) * spectralResolution) *&
                                          beamIntensity

      ! Redistribute beam intensity to reactorSection
      photonFlux = photonFlux * (beamSection/reactorSection)

      spF = sum(photonFlux) * spectralResolution

      print *,'Photons **************************'
      print *,'flux  / ph.cm^-2.s^-1 =',real(spF)
      print *,'flux  / ph.s^-1       =',real(spF * reactorSection)
      print *,'integ / ph            =',real(spF * reactorSection * reactionTime)
      print *,'**********************************'
!       photonFlux=photonFlux/9.537**2 !*20.2 ! Titan
!       photonFlux=0d0;photonFlux(82)=6.7d13*1d4/Avogadro    ! irradiation at 82 nm in Imanaka and Smith's reactor

      print *,'Absorption ***********************'
      
      imaxL  = size(Lphoto,1)    
      ! Allcocate tables also used in PHOTO_KINET_SPARSE
      allocate(absorb(nbPhotoReac,sp1:sp2), &
               sumabs(sp1:sp2), &
               intabs(sp1:sp2))
      absorb = 0d0
      do i = 1, imaxL
        ii = Lphoto(i,1)
        absorb(ii,sp1:sp2) = crossSections(ii,sp1:sp2) * &
                             initialConcentrations(Lphoto(i,2)) ! cm^-1
      enddo
      do k = sp1, sp2
        sumabs(k)=sum(absorb(1:nbPhotoReac,k))
      end do

      !!! WORK A LITTLE MORE ON THAT !!! Mean of the absorption per wl ???
      print *,'% absorbed (WIP!!!) ', (sum(1-exp(-sumabs*dx)))* 100 / & 
                              (spectrumRange(2)-spectrumRange(1)) ! ????
      print *,'**********************************'
      
      open(55, file='init_absorption.out')
      do il = sp1, sp2
        write(55,*) real(il*spectralResolution), real(1-exp(-sumabs(il)*dx))
      enddo
      close(55)
      
    END SUBROUTINE initPhotoChemistry

    SUBROUTINE initChemistry
      IMPLICIT NONE
      integer             :: i, j, nbSpec, irate, imax, ios, nnz, nw
      DOUBLE PRECISION    :: tmp, k_para(20)
      character*10        :: type = 'kooij'
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
      open(55,file='check_rrates.out')
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

          case ('3body')
            imax = 10
            tmp  = gasTemperature

          case ('assocVV','assocMD')
            imax = 16
            tmp  = gasTemperature

          case('dr')
            imax = 5
            tmp  = electronsTemperature

          case ('kooij','assoc0','ionpol1','ionpol2')
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
        write(55,*) type, reactionRates(irate)
 
      end do
      close(55)
      close(20)
      ! print *, real(pToConc)
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
      IMPLICIT NONE
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

    DOUBLE PRECISION function krate(type,c,T,P) result(k)
      IMPLICIT NONE
      character(len=*)     :: type
      DOUBLE PRECISION, dimension(:) :: c
      DOUBLE PRECISION               :: T, F, k0, k0P, kInf, P, kR, Fc, F1
      DOUBLE PRECISION               :: Ci, Ni, kInf1, lPr, fExp, lF1 
      DOUBLE PRECISION, parameter    :: T0 = 300D0 ! Reference temp. for kooij expression

      select case (type)

        case ('kooij', 'dr')
          k = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k = k * F

        case ('assocMD')
          k0   = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F    = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k0P  = k0 * F * P
          kInf = c(6) * (T/T0)**c(7) * exp(-c(8)/T)
          F    = c(9) * exp( c(10)*abs(1/T-1/T0) )
          kInf = kInf * F
          kR   = c(11) * (T/T0)**c(12) * exp(-c(13)/T)
          F    = c(14) * exp( c(15)*abs(1/T-1/T0) )
          kR   = kR * F
          Fc   = c(16)

          Ni   = 0.75d0 - 1.27d0 * log10(Fc) ! approx. 1 for Fc=0.6 only...
          fExp = 1.d0 + (log10(k0P / kInf) / Ni)**2
          lF1  = log10(Fc) / fExp 
          k    = kInf * (k0P * 10.d0**lF1 + kR) / (k0P + kInf) 
          k    = k * wallFactor
      
        case ('assocVV') ! Different Kooij and Kassoc expressions
          kInf = c(1) * T**c(2) * exp(-c(3)/T)
          F    = c(4) * exp( c(5)*abs(1.d0/T-1.d0/T0) )
          kInf = kInf * F
          k0   = c(6) * T**c(7) * exp(-c(8)/T)
          F    = c(9) * exp( c(10)*abs(1.d0/T-1.d0/T0) )
          k0P  = k0 * F * P
          kR   = c(11) * T**c(12) * exp(-c(13)/T)
          F    = c(14) * exp( c(15)*abs(1.d0/T-1.d0/T0) )
          kR   = kR * F
          Fc   = c(16)

          if (kR >= 0.99d0 * kInf) then
            k = kInf * wallFactor
          else
            Ci    = -0.4d0 - 0.67d0 * log10(Fc)
            Ni    = 0.75d0 - 1.27d0 * log10(Fc)
            lPr   = log10(k0P / kInf)
            fExp  = 1.d0 + ((lPr + Ci) / (Ni - 0.14d0 * (lPr + Ci)))**2
            lF1   = log10(Fc) / fExp
            kInf1 = kInf - kR
            k     = kR + (10.d0**lF1 * kInf1 * k0P) / (kInf1 + k0P)
            k     = k * wallFactor
          end if
          
        case ('assoctroe')
          k0   = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F    = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k0   = k0 * F
          kInf = c(6) * (T/T0)**c(7) * exp(-c(8)/T)
          F    = c(9) * exp( c(10)*abs(1/T-1/T0) )
          kInf = kInf * F
          k    = troe(k0,kInf,P) * wallFactor
          
        case ('assoc0')
          k0   = c(1) * (T/T0)**c(2) * exp(-c(3)/T)
          F    = c(4) * exp( c(5)*abs(1/T-1/T0) )
          k0P  = k0 * F * P
          k    = k0P * wallFactor
          
        case ('ionpol1')
          k = c(1) * c(2) * (0.62d0 + 0.4767d0*c(3)*sqrt(T0/T))
          F = c(4) * exp( c(5) * abs(1.d0/T-1.d0/T0) )
          k = k * F

        case ('ionpol2')
          k = c(1) * c(2) * (1.d0 + 0.0967d0*c(3)*sqrt(T0/T) &
                                  + c(3)*c(3) / 10.526d0 * T0/T   )
          F = c(4) * exp( c(5) * abs(1.d0/T-1.d0/T0) )
          k = k * F

      end select

    end function krate

    DOUBLE PRECISION function troe(k0,kInf,P) result(k)
      IMPLICIT NONE
	    DOUBLE PRECISION :: k0, kInf, P
	    DOUBLE PRECISION :: cE, nE, dE, fE, Pr, Fc, lFc, c1

        ! Simplified Troe formula
        Fc = 0.64D0

        lFc= log10(Fc)
        Pr = k0 * P / kInf

        cE = -0.4D0 -0.67D0*lFc
        nE = 0.75D0 -1.27D0*lFc
        dE = 0.14D0
        c1 = log10(Pr) + cE
        fE = 1.d0 + (c1/(nE-dE*c1))**2

        k  = kInf * (Pr/(1.d0+Pr)) * Fc**(1.d0/fE)  

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
