      SUBROUTINE Nongrey2(P,Teq0,Tint,grav,Output,mu,f,&
      & ROSS,COEFF,COMP,STAR,ALBEDO,CONV,T,TEFF,Ab,TAU,gradad,gradrad,&
      & verbose,PRC,N)
      IMPLICIT NONE
!!!!!!INPUT VARIABLES!!!!!!!!!!!!!!!!
      !Model parameters 
!      REAL*8, intent(in) :: Pmin,Pmax !minimal and maximal pressure of the model. 
      INTEGER :: N ! Number of levels in the model
      REAL*8 :: P(N)
      CHARACTER*150 :: Output ! Name of the output file
      CHARACTER*150 :: OutputFIle, OutputFIleinfo ! Name of the output file
      !Planet parameters
      REAL*8, intent(in) :: Tint !Temperature for the internal flux
      REAL*8, intent(in) :: Teq0 !Temperature of the incoming flux
      REAL*8 :: mu,f
      CHARACTER*10 :: ROSS,COEFF,COMP,STAR,ALBEDO
      logical :: CONV,verbose
!f2py intent(in,out,copy) :: p(n)
!f2py integer intent(hide),depend(p) :: n=size(p)
!f2py CHARACTER*150, intent(in), optional :: Output=''
!f2py REAL*8, intent(in), optional :: mu=1/sqrt(3.0) !Cosine of the stellar angle mu=1/sqrt(3) for averaged profiles
!f2py REAL*8, intent(in), optional :: f=0.25 !Parameter for averaged profiles f=0.5 for dayside average and f=0.25 for planet-average
      REAL*8, intent(in) :: grav !Gravity of the planet (m/s**2)
      !Options
!f2py CHARACTER*10, intent(in), optional :: ROSS='auto' ! Can be either 'user' when Kappa is defined by the user or 'auto' to use the opacities of Freedman et al. 2008 as fitted by Valencia et al. 2013
!f2py CHARACTER*10, intent(in), optional :: COEFF='auto' ! Can be either 'user' when the input coefficients for Gv1,Gv2,Betav,Beta,Gp or 'auto' to use the fit provided by Parmentier et al. 2014
!f2py CHARACTER*10, intent(in), optional :: COMP='solar' ! Can be 'solar' for a solar-composition atmosphere or 'NOTIO' for an atmosphere without TiO/VO.
!f2py CHARACTER*10, intent(in), optional :: STAR='sun' ! Can only be 'sun' for now. 
!f2py CHARACTER*10, intent(in), optional :: ALBEDO='auto' !If 'user' uses the value of Ab, if 'auto', uses the fit of Parmentier at al. 2014.
!f2py logical, intent(in), optional :: CONV=1 !If YES calculates a radiative/convective profile, if NO, calculates only the radiative solution
!f2py logical, intent(in), optional :: verbose=1
      !Atmosphere parameters
      REAL*8 :: Gv1 !Gamma_v1=Kv1/Kr is the ratio of the first visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gv2 !Gamma_v2=Kv2/Kr is the ratio of the second visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gv3 !Gamma_v2=Kv3/Kr is the ratio of the third visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gp  !Gamma_P=Kp/Kr ratio of the Planck to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Beta ! Relative spectral width of the first thermal band
!      REAL*8 :: Ab ! Bond albedo of the planet

!!!!!!INPUTS/OUTPUTS!!!!!!!!!!
      REAL*8, DIMENSION(N) :: Kappa ! Rosseland mean opacity at each atmospheric level, either calculated or set by the user

!!!!!!!OUTPUTS!!!!!!!!!!!!!
      REAL*8, intent(out) :: Teff ! Effective temperature of the model 
      REAL*8, intent(out) :: Ab
      REAL*8, intent(out) :: TAU(N) ! Rosseland optical depth
!      REAL*8, DIMENSION(N) :: P   ! Pressure
      REAL*8, intent(out) :: T(N)   ! Temperature
      REAL*8, intent(out) :: gradad(N) ! Adiabatic gradient
      REAL*8, intent(out) :: gradrad(N) ! Radiative gradients
      REAL*8, intent(out) :: PRC ! Pressure at the radiative/convective boudary

!!!!!!!INTERNAL VARIABLES
      INTEGER :: i,K
      CHARACTER*4 :: xtemp
!      print *, ROSS,STAR,COEFF,ALBEDO
!      print *, ALBEDO==''
!       print *, size(P)
!       print *, transfer( conv , 1 )
!       OutputFile='te1500ti500g5.tp'
       if (Output/='') then
       OutputFIleinfo=trim(Output)//'.inf'
       OutputFile = trim(Output)//'.tp'
       endif
       

!!!!!!Model parameters
!      Pmin=1e1 !Minimal pressure (Pa)
!      Pmax=1e6 ! Maximal pressure (Pa)
!!!!!!!Planet parameters
!      Tint=500!Internal temperature (K)
!      Teq0=1500!!Equilibrium temperature for zero albedo (K)
!      f=0.25 !f=0.25 for a planet-average profile, f=0.5 for a dayside average profile
!      mu=1/sqrt(3.0) ! Cosine of the irradiation angle. mu=1/sqrt(3) is the mean mu to use for a dayside or a planet average.
!!      Teff0=1000 !Effective temperature for zero albedo
!      grav=5! gravity (m/s**2)


!!!!!Options
!      ROSS="AUTO" !This can be AUTO or USER. If USER, the Rosseland mean opacity must be provided below for each level. If AUTO, the Rosseland opacities are calculated from the fit of Valencia et al. 2013
!      COEFF="AUTO" !This can be AUTO or USER whether you want to use the fit of the coefficients provided in Parmentier et al. 2013 or set your own coefficients.
!      STAR="SUN"
!      COMP="SOLAR"    !This can be SOLAR for a solar-composition atmosphere or NOTIO for an atmosphere without TiO/VO
!      ALBEDO="AUTO" !This can be AUTO or USER. If USER, the Bond albedo must be provided below. If AUTO, the Bond Albedo is calculated from the fit of Parmentier et al. 2014
!      CONV="YES" ! If "YES", a convective zone is calculated at the bottom of the model using the convective gradient of Parmentier et al. 2014. If "NO" only the radiative solution is provided

!!!!!Input variables, used only if COEFF="user"
      if(COEFF.eq."user") then
      Gv1=10 
      Gv2=0.1
      Gv3=0.01 ! 0.5 is two bands of equal width for the absorption of the stellar flux
      Gp=10! Gp=10 Planck opacities are ten times higher than Rosseland ones.
      beta=0.5! The two thermal bands have the same width
      endif

!!!!!Input variables, used only if ALBEDO="user"
      if(ALBEDO.eq."user") then
      Ab=0!!Bond albedo set by the user if Albedo="user".
      endif

!!!!!Input variables, used only if ROSS="user"
      if(ROSS.eq."user") then
      DO i=1,N
      Kappa(i)=1e-1 ! kg/m^2 Here I set the Rosseland mean opacity, if ROSS="CONSTANT". Kappa(i) can be any functional form of the pressure P(i)
      ENDDO
      endif

!!!!!!!!!!!!!!END OF INPUTS!!!!!!!!!!!!!!!!!!!!!

!      !Set the pressure levels
!      DO i=1,N
!       P(i)=10**(log10(Pmin)+(log10(Pmax)-log10(Pmin))*                 &
!     &1.0*(i-1)/(N-1))!The pressure varies from Pmin to Pmax in logspace with N points!     
!      ENDDO

      !Calculate the Temperature 
      call tprofile(Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Kappa,      &
     &Ab,Teff,ROSS,COEFF,COMP,STAR,CONV,ALBEDO,TAU,P,T,N,PRC,gradad,    &
     &gradrad,verbose)
    
      !Output formats
      if (Output/='') then
 1400    format (f10.2,1X,f10.2,1X,f10.2,1X,f10.2,1X,f10.2,1X,f6.2,1X,  &
      &e9.3,1X,e9.3,1X,f5.4,1X,e9.3,1X,f5.4,1X,f5.4,1X                  & 
      &,A10,A10,A10,A10,A10,A10' Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp&  
      &,Beta,Ab,ROSS,COEFF,COMP,STAR,CONV,ALBEDO')

 1500   format(                                                         &
     &      'c',t1,'(1)N',t6,'(2)P(Pa)',t17,'(3)T(K)'                   &
     &      ,t27,'(4)Tau', t37,'(5)Kr(m^2/kg)', t50, '(6)Gradad'        &
     &      , t60, '(7)Gradrad', t72,'(8)PRC (Pa)')

 9900   format (1X,I3,",",E9.3,",",F10.3,",",E9.3,",",E9.3,",",F8.5,"," &
     &,F8.5,",",E9.3,",",E9.3,",",E9.3,",",F5.3,",",E9.3,",",F5.3,",",  &
     &F6.3,",",F6.1,",",F6.1,",",F6.1,",",F5.3,",",F5.3,",",F8.3,",",   &
     &A10,",",A10,",",A10,",",A10,",",L,",",A10)

    
      write(xtemp,'(I4.4)') Int((1-Ab)**0.25*(4*f)**0.25*Teq0)
!       OutputFile='./profilesNoTiO/PTprofile-Tmu-'//xtemp//'.dat'
      !Open and write output file
       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFIleinfo&
     &)
       write(40,'(A)')'N,P(Pa),T(K),Tau,Kappa(m^2/kg),gradad,gradrad,PRC&
     &,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff(K),Tint(K),Teq0(K),mu,f,grav(m/s^2), &
     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO'
        DO i=1,N
         write(40,9900)i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i),&
     &PRC,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff,Tint,Teq0,mu,f,grav,ROSS,COEFF,   &
     &COMP,STAR,CONV,ALBEDO        
        ENDDO
       close(40)
       
!9800   format (1X,I3,1X,E9.3,1X,F10.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,      &
!     &E9.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,E9.3)
!       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFile)
!       write(40,1400) Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Ab,&
!     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO
!       write(40,*),' '
!       write(40,*),' '
!       write(40,1500)
!       DO i=1,N
!         write(40,9800) i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i)&
!     &,PRC
!       ENDDO
!       close(40)
9800   format (4X,I3,4X,E9.3,4X,F10.3)
       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFile)
       write(40,'(4X,A3,4X,A9,4X,A10)') 'I','P','T'
!       write(40,9800) 'i','P','T'
!       write(40,*),' '
!       write(40,*),' '
!       write(40,1500)
       DO i=1,N
         write(40,9800) i,P(i),T(i)
       ENDDO
      close(40)
!      ENDDO
      write(*,*) 'Done, PTprofile is available in files: "',            &
     &TRIM(OutputFile),'" or "',TRIM(OutputFIleinfo),'"'
     endif
      END subroutine
