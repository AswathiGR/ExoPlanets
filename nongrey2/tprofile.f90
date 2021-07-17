!This subroutine computes the PT profile for an irradiated plane-parallel atmosphere using the analytical solution of Parmentier & Guillot 2013 and the coefficients of Parmentier et al. 2013
! Units in the code are SI. (P in Pa, Kappa in kg/m^2, T in K, grav in m²/s etc...)
! This code makes use of the fit to the Rosseland mean of Freedman et al. 2008 given by Valencia et al. 2013
! This code is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License (see http://creativecommons.org/licenses/by-nc-sa/3.0/ for more details).

      SUBROUTINE tprofile(Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,     &
     &Kappa,Ab,Teff,ROSS,COEFF,COMP,STAR,CONV,ALBEDO,TAU,P,T,N,PRC,     &
     &gradad,gradrad,verbose)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: Tint !Temperature for the internal flux
      REAL*8, INTENT(IN) :: Teq0 !Temperature of the incoming flux
      REAL*8, INTENT(IN) :: mu !Stellar angle
      REAL*8, INTENT(IN) :: f !Parameter for averaged profiles
      REAL*8, INTENT(IN) :: grav      !Gravity of the planet
      LOGICAL, INTENT(IN) :: verbose
      REAL*8, INTENT(INOUT) :: Gv1 !Gamma_v1=Kv1/Kr is the ratio of the first visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8, INTENT(INOUT) :: Gv2 !Gamma_v2=Kv2/Kr is the ratio of the second visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8, INTENT(INOUT) :: Gv3 !Gamma_v2=Kv2/Kr is the ratio of the second visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8, INTENT(INOUT) :: Gp  !Gamma_P=Kp/Kr ratio of the Planck to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8, INTENT(INOUT) :: Beta ! Relative spectral width of the first thermal band
      REAL*8, INTENT(INOUT) :: Ab ! Bond albedo of the planet
      REAL*8, INTENT(INOUT) :: Teff !Effective temperature of the model
      INTEGER, INTENT(IN) :: N ! Number of levels in the atmosphere
      REAL*8, DIMENSION(N), INTENT(IN) :: P !Pressure
      CHARACTER*10, INTENT(IN) :: ROSS !Can be either "user" when Kappa is defined by the user or 'auto' to use the opacities of Freedman et al. 2008 as fitted by Valencia et al. 2013
      CHARACTER*10, INTENT(IN) :: COEFF !Can be either "user" when the input coefficients for Gv1,Gv2,Betav,Beta,Gp or 'auto' to use the fit provided by Parmentier et al. 2014
      CHARACTER*10, INTENT(IN) :: COMP !Can be 'solar' for a solar-composition atmosphere or 'NOTIO' for an atmosphere without TiO/VO.
      CHARACTER*3, INTENT(IN) :: STAR ! Can only be 'sun' for now. 
      logical, INTENT(IN) :: CONV !If YES calculates a radiative/convective profile, if NO, calculates only the radiative solution
      CHARACTER*10, INTENT(IN) :: ALBEDO !If USER uses the value of Ab, if AUTO, uses the fit of Parmentier at al. 2014.
      REAL*8, DIMENSION(N), INTENT(INOUT) :: Kappa !Rosseland mean opacity used in the model.

      REAL*8, DIMENSION(N), INTENT(OUT) :: TAU ! Optical depth
      REAL*8, DIMENSION(N), INTENT(OUT) :: T ! Temperature
      REAL*8, INTENT(OUT) :: PRC !Pressure of the radiative/convective boudary
      REAL*8,PARAMETER :: met=0 !Metalicity use for the calculation of the Rosseland mean opacities when using the fit by Valencia et al. 2013 by default it is zero. Could be used as an input parameter later.
      REAL*8, DIMENSION(N), INTENT(OUT) :: gradad !Convective gradient d(ln(T))/d(ln(P))
      REAL*8, DIMENSION(N), INTENT(OUT) :: gradrad !Raiative gradient
! Now all the internal variables
      REAL*8 :: Betav1,Betav2,Betav3
      REAL*8 :: aGp,bGp,cGp,aBa,aGp2,bGp2,cGp2
      REAL*8 :: l1Gv1,l2Gv1,l3Gv1,l4Gv1,l5Gv1,l6Gv1,k1Gv1,k2Gv1,k3Gv1,  &
     &k4Gv1,k5Gv1,k6Gv1
      REAL*8 :: l1Gv2,l2Gv2,l3Gv2,l4Gv2,l5Gv2,l6Gv2,k1Gv2,k2Gv2,k3Gv2,  &
     &k4Gv2,k5Gv2,k6Gv2
      REAL*8 :: l1Gv3,l2Gv3,l3Gv3,l4Gv3,l5Gv3,l6Gv3,k1Gv3,k2Gv3,k3Gv3,  &
     &k4Gv3,k5Gv3,k6Gv3
      REAL*8 :: l1B,l2B,l3B,l4B,l5B,l6B,k1B,k2B,k3B,k4B,k5B,k6B,l7B,k7B
      REAL*8 :: At1,At2,AA1v1,AA1v2,AA2v1,AA2v2,AA1v3,AA2v3
      REAL*8 :: a0,a1,a2v1,a2v2,a3v1,a3v2,b0,b1v1,b1v2,b2v1,b2v2,b3v1,  &
     &b3v2,b1v3,b2v3,b3v3,a3v3,a2v3
      REAL*8 :: A,B,Cv1,Cv2,Dv1,Dv2,Ev1,Ev2,Cv3,Dv3,Ev3
      REAL*8 :: R,Gv1s,Gv2s,Gv3s,TauLim,G1,G2 
      REAL*8 :: T1,T2,T3,T4,T5,T6
      REAL*8 :: X
      INTEGER :: i,j,iRC1
      REAL*8 :: aAlb,bAlb
      REAL*8 :: Teff0
      REAL*8 :: TAURC
      REAL*8 :: Tmu,Tmu0
      INTEGER :: iRC
!Error message if the input parameters are out of bounds
      if (Tint<0) then
       Print*, "Error : Tint must be positive"
       STOP
      endif

      if (Teq0<0) then 
       Print*, "Error : Teq0 must be positive"
       STOP
      endif

      if (f>1) then
       PRINT*, "Error : f must be smaller than one"
       STOP
      elseif (f<0) then
       PRINT*, "Error : f must be positive"
       STOP
      endif

      if (mu>1) then
       PRINT*, "Error : mu=cos(theta) must be smaller than one"
       STOP
      elseif (mu<0) then
       PRINT*, "Error : mu=cos(theta) must be positive"
       STOP
      endif

      if (grav<0) then
       Print*, "Error, gravity must be positive"
       STOP
      endif
      
      Tmu0=(4*f)**0.25*Teq0
      Teff0=(Tint**4+Tmu0**4)**0.25
!Calculate the Albedo
     if(ALBEDO.eq."user") then
      if(Ab<0) then
        Print*, "Error: Bond albedo Ab should be positive"
        STOP
      elseif(Ab>1) then
        Print*, "Error: Bond albedo Ab should be smaller than 1"
        STOP
      endif
     elseif(ALBEDO.eq.'auto') then
       if(Teff0<250) then
        aAlb=-0.335*grav**0.070
        bAlb=0
       elseif(Teff0<750) then
        aAlb=-0.335*grav**0.070+2.149*grav**0.135
        bAlb=-0.896*grav**0.135
       elseif(Teff0<1250) then
        aAlb=-0.335*grav**0.070-0.428*grav**0.135
        bAlb=0
       else
        aAlb=16.947-3.174*grav**0.070-4.051*grav**0.135
        bAlb=-5.472+0.917*grav**0.070+1.170*grav**0.135
       endif
      
       Ab=aAlb+bAlb*log10(Teff0)
       Ab=10**Ab
      else
       Print*,'ERROR: ALBEDO must be either "AUTO" or "user"'
       STOP
     endif
       
      Tmu=(4*f*(1-Ab))**0.25*Teq0
      Teff=(Tint**4+Tmu**4)**0.25

! If the coefficients do not come from the fit, checking the values of the coefficients
      if(COEFF.eq."user") then
       if (Gv1<0) then
        Print*, "Error : Gv1 must be positive"
        STOP
       endif

       if (Gv2<0) then
        Print*, "Error : Gv2 must be positive"
        STOP
       endif

       if (Gp<1) then
        Print*, "Error : Gp must be bigger than 1"
        STOP
       endif

       if (Beta<0) then
        Print*, "Error : Beta must be positive"
        STOP
       elseif (Beta>1) then 
        Print*, "Error : Beta must be lower than 1"
        STOP
       endif 
      elseif(COEFF.eq.'auto') then
        if (Tmu.lt.Tint) then
        if(verbose) write(*,*)'WARNIG : Tmu<Tint, the coefficients of &
     &Parmentier et al. 2014 might not be valid in this regime.'
        endif
       if(COMP.eq.'solar') then      
        if(STAR.eq.'sun') then
         l1Gv1=-5.51
         l2Gv1=1.23
         l3Gv1=8.65
         l4Gv1=-12.96
         l5Gv1=-23.75
         l6Gv1=12.65
         
         k1Gv1=2.48
         k2Gv1=-0.45
         k3Gv1=-3.45
         k4Gv1=4.33
         k5Gv1=7.76
         k6Gv1=-3.27
         
         l1Gv2=-7.37
         l2Gv2=13.99
         l3Gv2=-15.18
         l4Gv2=-10.41
         l5Gv2=-19.95
         l6Gv2=13.56
         
         k1Gv2=2.53
         k2Gv2=-6.75
         k3Gv2=5.02
         k4Gv2=3.31
         k5Gv2=6.34
         k6Gv2=-3.81
         
         l1Gv3=-3.03
         l2Gv3=-13.87
         l3Gv3=-11.95
         l4Gv3=-6.97
         l5Gv3=-3.65
         l6Gv3=-6.02
         
         k1Gv3=-0.2
         k2Gv3=4.51
         k3Gv3=3.74
         k4Gv3=1.94
         k5Gv3=0.89
         k6Gv3=1.61
         
         l1B=0.84
         l2B=0.84
         l3B=0.84
         l4B=0.84
         l5B=0.84
         l6B=6.21
         
         k1B=0.0
         k2B=0.0
         k3B=0.0
         k4B=0.0
         k5B=0.0
         k6B=-1.63
         aGp=-2.36
         bGp=13.92
         cGp=-19.38
         aGp2=aGp
         bGp2=bGp
         cGp2=cGp
        elseif (STAR.eq.'GJ1214') then
         Print*, "Error : STAR 'GJ1214' not implemented yet, please use &
     &'sun'"
         STOP
        else 
         Print*, "Error : STAR must be either 'sun' or 'GJ1214"
         STOP
        endif !End the condition on STAR

       elseif(COMP.eq.'NOTIO') then
        if(STAR.eq.'sun') then 
         l1Gv1=-5.51
         l2Gv1=1.23
         l3Gv1=8.65
         l4Gv1=-12.96
         l5Gv1=-1.68
         l6Gv1=10.37
         
         k1Gv1=2.48
         k2Gv1=-0.45
         k3Gv1=-3.45
         k4Gv1=4.33
         k5Gv1=0.75
         k6Gv1=-2.91
         
         l1Gv2=-7.37
         l2Gv2=13.99
         l3Gv2=-15.18
         l4Gv2=-10.41
         l5Gv2=6.96
         l6Gv2=-2.4
         
         k1Gv2=2.53
         k2Gv2=-6.75
         k3Gv2=5.02
         k4Gv2=3.31
         k5Gv2=-2.21
         k6Gv2=0.62
         
         l1Gv3=-3.03
         l2Gv3=-13.87
         l3Gv3=-11.95
         l4Gv3=-6.97
         l5Gv3=0.02
         l6Gv3=-16.54
         
         k1Gv3=-0.2
         k2Gv3=4.51
         k3Gv3=3.74
         k4Gv3=1.94
         k5Gv3=-0.28
         k6Gv3=4.74
         
         l1B=0.84
         l2B=0.84
         l3B=0.84
         l4B=0.84
         l5B=3.0
         l6B=3.0
         
         k1B=0.0
         k2B=0.0
         k3B=0.0
         k4B=0.0
         k5B=-0.69
         k6B=-0.69
         aGp=-2.36
         bGp=13.92
         cGp=-19.38
         aGp2=-12.45
         bGp2=82.25
         cGp2=-134.42


        elseif (STAR.eq.'GJ1214') then
         Print*, "Error : STAR 'GJ1214' not implemented yet, please use &
     &'sun'"
         STOP
        else 

        Print*, "Error : STAR must be either 'sun' or 'GJ1214"
        STOP
        endif !End the condition on STAR

       else
       Print*, "Error : COMP must be either 'solar' or 'NOTIO'"
       STOP
       endif ! End the condition on TiO

!In the case I used Fit I recalculate the coefficients using in function of the irradiation temperature

        X=log10(Teff)
           T1=200.000000
           T2=300.000000
           T3=600.000000
           T4=1400.000000
           T5=2000.000000
 
         if(Teff.le.T1) then        
          Gv1=l1Gv1+k1Gv1*X
          Gv2=l1Gv2+k1Gv2*X
          Gv3=l1Gv3+k1Gv3*X
          Beta=l1B+k1B*X
          Gp=aGp*X**2+bGp*X+cGp
        elseif(Teff.le.T2) then
          Gv1=l2Gv1+k2Gv1*X
          Gv2=l2Gv2+k2Gv2*X
          Gv3=l2Gv3+k2Gv3*X
          !Gp=l2Gp+k2Gp*X
          Beta=l2B+k2B*X
          Gp=aGp*X**2+bGp*X+cGp
        elseif(Teff.le.T3) then
          Gv1=l3Gv1+k3Gv1*X
          Gv2=l3Gv2+k3Gv2*X
          Gv3=l3Gv3+k3Gv3*X
          !Gp=l3Gp+k3Gp*X
          Beta=l3B+k3B*X
          Gp=aGp*X**2+bGp*X+cGp
        elseif(Teff.le.T4) then
          Gv1=l4Gv1+k4Gv1*X
          Gv2=l4Gv2+k4Gv2*X
          Gv3=l4Gv3+k4Gv3*X
          !Gp=l4Gp+k4Gp*X
          Beta=l4B+k4B*X
          Gp=aGp*X**2+bGp*X+cGp
        elseif(Teff.le.T5) then
          Gv1=l5Gv1+k5Gv1*X
          Gv2=l5Gv2+k5Gv2*X
          Gv3=l5Gv3+k5Gv3*X
          !Gp=l5Gp+k5Gp*X
          Beta=l5B+k5B*X
          Gp=aGp2*X**2+bGp2*X+cGp2
        elseif(Teff.ge.T5) then
          Gv1=l6Gv1+k6Gv1*X
          Gv2=l6Gv2+k6Gv2*X
          Gv3=l6Gv3+k6Gv3*X
          !Gp=l6Gp+k6Gp*X
          Beta=l6B+k6B*X
          Gp=aGp2*X**2+bGp2*X+cGp2
      endif !End of the condition on temperature
          Gv1=10**Gv1
          Gv2=10**Gv2
          Gv3=10**Gv3
          Gp=10**Gp
          Betav1=1.0/3
          Betav2=1.0/3
          Betav3=1.0/3 
        if(Gp.lt.1) then
         Print*, "WARNING : the fit gives Gp<1, so Gp is set to 1.      &
     &This correspond to a grey model for the thermal opacities"
       Gp=1+1e-6
      endif
      else
       Print*, "Error : COEFF must be either 'user' or 'auto'"
       STOP
      endif ! End the condition on Fit

      !We calculate different useful quantities
      Gv1s=Gv1/mu
      Gv2s=Gv2/mu 
      Gv3s=Gv3/mu 
      R=1+(Gp-1)/(2*Beta*(1-Beta))+sqrt(((Gp-1)/(2*Beta*(1-Beta)))**2+  &
     &(Gp-1)/(2*Beta*(1-Beta)))
      G1=Beta+R-Beta*R
      G2=G1/R
      TauLim=1/(G1*G2)*sqrt(Gp/3.0)

      !Now we calculate the coefficients 
      At1=G1**2*log(1+1/(TauLim*G1))
      At2=G2**2*log(1+1/(TauLim*G2))
      AA1v1=G1**2*log(1+Gv1s/G1)
      AA2v1=G2**2*log(1+Gv1s/G2)
      AA1v2=G1**2*log(1+Gv2s/G1)
      AA2v2=G2**2*log(1+Gv2s/G2)
      AA1v3=G1**2*log(1+Gv3s/G1)
      AA2v3=G2**2*log(1+Gv3s/G2)
      a0=1/G1+1/G2

      a1=-1.0/(3*TauLim**2)*(Gp/(1-Gp)*(G1+G2-2)/(G1+G2)+(G1+G2)*TauLim-&
     &(At1+At2)*TauLim**2)!a1->infinity for Gp->1 but th    e product a1*b0 -> 0. Check that the computer does the calculation correctly.

      a2v1=TauLim**2/(Gp*Gv1s**2)*((3*G1**2-Gv1s**2)*(3*G2**2-Gv1s**2)* &
     &(G1+G2)-3*Gv1s*(6*G1**2*G2**2-Gv1s**2*(G1**2+G2**2)))/(1-Gv1s**2* &
     &TauLim**2)

      a2v2=TauLim**2/(Gp*Gv2s**2)*((3*G1**2-Gv2s**2)*(3*G2**2-Gv2s**2)* &
     &(G1+G2)-3*Gv2s*(6*G1**2*G2**2-Gv2s**2*(G1**2+G2**2)))/(1-Gv2s**2* &
     &TauLim**2)

      a2v3=TauLim**2/(Gp*Gv3s**2)*((3*G1**2-Gv3s**2)*(3*G2**2-Gv3s**2)* &
     &(G1+G2)-3*Gv3s*(6*G1**2*G2**2-Gv3s**2*(G1**2+G2**2)))/(1-Gv3s**2* &
     &TauLim**2)

      a3v1=-TauLim**2*(3*G1**2-Gv1s**2)*(3*G2**2-Gv1s**2)*(AA2v1+AA1v1) &
     &/(Gp*Gv1s**3*(1-Gv1s**2*TauLim**2))

      a3v2=-TauLim**2*(3*G1**2-Gv2s**2)*(3*G2**2-Gv2s**2)*(AA2v2+AA1v2) &
     &/(Gp*Gv2s**3*(1-Gv2s**2*TauLim**2))

      a3v3=-TauLim**2*(3*G1**2-Gv3s**2)*(3*G2**2-Gv3s**2)*(AA2v3+AA1v3) &
     &/(Gp*Gv3s**3*(1-Gv3s**2*TauLim**2))

      b0=1.0/(G1*G2/(G1-G2)*(At1-At2)/3-(G1*G2)**2/sqrt(3*Gp)-(G1*G2)**3&
     &/((1-G1)*(1-G2)*(G1+G2)))

      b1v1=G1*G2*(3*G1**2-Gv1s**2)*(3*G2**2-Gv1s**2)*TauLim**2/         &
     &(Gp*Gv1s**2*(Gv1s**2*TauLim**2-1))

      b1v2=G1*G2*(3*G1**2-Gv2s**2)*(3*G2**2-Gv2s**2)*TauLim**2/         &
     &(Gp*Gv2s**2*(Gv2s**2*TauLim**2-1))

      b1v3=G1*G2*(3*G1**2-Gv3s**2)*(3*G2**2-Gv3s**2)*TauLim**2/         &
     &(Gp*Gv3s**2*(Gv3s**2*TauLim**2-1))

      b2v1=3*(G1+G2)*Gv1s**3/((3*G1**2-Gv1s**2)*(3*G2**2-Gv1s**2))
      b2v2=3*(G1+G2)*Gv2s**3/((3*G1**2-Gv2s**2)*(3*G2**2-Gv2s**2))
      b2v3=3*(G1+G2)*Gv3s**3/((3*G1**2-Gv3s**2)*(3*G2**2-Gv3s**2))
      b3v1=(AA2v1-AA1v1)/(Gv1s*(G1-G2))
      b3v2=(AA2v2-AA1v2)/(Gv2s*(G1-G2))
      b3v3=(AA2v3-AA1v3)/(Gv3s*(G1-G2))

      A=1.0/3*(a0+a1*b0)
      B=-1.0/3*(G1*G2)**2/Gp*b0
      Cv1=-1.0/3*(b0*b1v1*(1+b2v1+b3v1)*a1+a2v1+a3v1)
      Cv2=-1.0/3*(b0*b1v2*(1+b2v2+b3v2)*a1+a2v2+a3v2)
      Cv3=-1.0/3*(b0*b1v3*(1+b2v3+b3v3)*a1+a2v3+a3v3)
      Dv1=1.0/3*(G1*G2)**2/Gp*b0*b1v1*(1+b2v1+b3v1)
      Dv2=1.0/3*(G1*G2)**2/Gp*b0*b1v2*(1+b2v2+b3v2)
      Dv3=1.0/3*(G1*G2)**2/Gp*b0*b1v3*(1+b2v3+b3v3)
      Ev1=(3-(Gv1s/G1)**2)*(3-(Gv1s/G2)**2)/(9*Gv1s*((Gv1s*TauLim)**2-1)&
     &)
      Ev2=(3-(Gv2s/G1)**2)*(3-(Gv2s/G2)**2)/(9*Gv2s*((Gv2s*TauLim)**2-1)&
     &)
      Ev3=(3-(Gv3s/G1)**2)*(3-(Gv3s/G2)**2)/(9*Gv3s*((Gv3s*TauLim)**2-1)&
     &)
!Now I calculate the temperature in function of Tau
      if (ROSS.eq."user") then !If there is no opacity file provided, we use constant opacity Kappa(1)
       DO i=1,N
        if(i.eq.1) then
         TAU(1)=Kappa(1)/grav*P(1)
        else
         TAU(i)=TAU(i-1)+sqrt(Kappa(i-1)*Kappa(i))/grav*(P(i)-P(i-1))
        endif
      T(i)=(3.0*Tint**4/4*(TAU(i)+A+B*exp(-TAU(i)/TauLim))+3.0*Betav1   &
     &*Tmu**4/4*(Cv1+Dv1*exp(-TAU(i)/TauLim)+Ev1*exp(-Gv1s*TAU(i)))     &
     &+3.0*(Betav2)*Tmu**4/4*(Cv2+Dv2*exp(-TAU(i)/TauLim)+Ev2*          &
     &exp(-Gv2s*TAU(i)))                                                &
     &+3.0*(Betav3)*Tmu**4/4*(Cv3+Dv3*exp(-TAU(i)/TauLim)+Ev3*          &
     &exp(-Gv3s*TAU(i))))**0.25

       ENDDO

      elseif (ROSS.eq.'auto') then !Here I use the fit of the Rosseland mean opacities from Valencia at al. 2013
      !I first Estimate the skin temperature setting tau(1)=0 to have an estimate of kappa(1)
       TAU(1)=0
      T(1)=(3.0*Tint**4/4*(TAU(1)+A+B*exp(-TAU(1)/TauLim))+3.0*Betav1   &
     &*Tmu**4/4*(Cv1+Dv1*exp(-TAU(1)/TauLim)+Ev1*exp(-Gv1s*TAU(1)))     &
     &+3.0*Betav2*Tmu**4/4*(Cv2+Dv2*exp(-TAU(1)/TauLim)+Ev2*            &
     &exp(-Gv2s*TAU(1)))                                                &
     &+3.0*Betav3*Tmu**4/4*(Cv3+Dv3*exp(-TAU(1)/TauLim)+Ev3*            &
     &exp(-Gv3s*TAU(1))))**0.25

      call valencia(P(1),T(1),met,Kappa(1)) !I estimate Kappa(1)

      TAU(1)=Kappa(1)/grav*P(1) ! I can recalculate tau(1) and so T(1) with the new value of kappa(1)
      T(1)=(3.0*Tint**4/4*(TAU(1)+A+B*exp(-TAU(1)/TauLim))+3.0*Betav1   &
     &*Tmu**4/4*(Cv1+Dv1*exp(-TAU(1)/TauLim)+Ev1*exp(-Gv1s*TAU(1)))     &
     &+3.0*Betav2*Tmu**4/4*(Cv2+Dv2*exp(-TAU(1)/TauLim)+Ev2*            &
     &exp(-Gv2s*TAU(1)))                                                & 
     &+3.0*Betav3*Tmu**4/4*(Cv3+Dv3*exp(-TAU(1)/TauLim)+Ev3*            &
     &exp(-Gv3s*TAU(1))))**0.25 

       DO i=2,N
        call valencia(sqrt(P(i-1)*P(i)),T(i-1),met,Kappa(i))
        TAU(i)=TAU(i-1)+Kappa(i)/grav*(P(i)-P(i-1))
      T(i)=(3.0*Tint**4/4*(TAU(i)+A+B*exp(-TAU(i)/TauLim))+3.0*Betav1   &
     &*Tmu**4/4*(Cv1+Dv1*exp(-TAU(i)/TauLim)+Ev1*exp(-Gv1s*TAU(i)))     &
     &+3.0*Betav2*Tmu**4/4*(Cv2+Dv2*exp(-TAU(i)/TauLim)+Ev2*            &
     &exp(-Gv2s*TAU(i)))                                                &
     &+3.0*Betav3*Tmu**4/4*(Cv3+Dv3*exp(-TAU(i)/TauLim)+Ev3*            &
     &exp(-Gv3s*TAU(i))))**0.25
        DO j=1,5
        call valencia(sqrt(P(i-1)*P(i)),sqrt(T(i-1)*T(i)),met,      &
     &Kappa(i))
        TAU(i)=TAU(i-1)+Kappa(i)/grav*(P(i)-P(i-1))
      T(i)=(3.0*Tint**4/4*(TAU(i)+A+B*exp(-TAU(i)/TauLim))+3.0*Betav1   &
     &*Tmu**4/4*(Cv1+Dv1*exp(-TAU(i)/TauLim)+Ev1*exp(-Gv1s*TAU(i)))     &
     &+3.0*Betav2*Tmu**4/4*(Cv2+Dv2*exp(-TAU(i)/TauLim)+Ev2*            &
     &exp(-Gv2s*TAU(i)))                                                &
     &+3.0*Betav3*Tmu**4/4*(Cv3+Dv3*exp(-TAU(i)/TauLim)+Ev3*            &
     &exp(-Gv3s*TAU(i))))**0.25
        ENDDO
       ENDDO
      else 
       Print*, "ROSS must be either 'user' or 'auto'"
       STOP
      endif
!We now provide a convective adjustment based on the Schwarchild criterion. 
!We only handle a convective zone from the bottom of the model to the first radiative/convective boundary
!Convective zone above this RC boundary (also called "detached" convective zones) are not handled
!We use a simple equation of state for the convective gradient T/T0=(P/P0)**gradad  with gradad=0.33-0.1*(T/3000K) for molecular hydrogen

      gradrad(N)=0
      gradad(N)=0
      DO i=1,N-1
         gradrad(i)=(log10(T(i))-log10(T(i+1)))/(log10(P(i))-           &
     &log10(P(i+1)))
         gradad(i)=0.32-0.10*T(i)/3000
      ENDDO
      iRC=N-1
      iRC1=N-1
      if(CONV) then
       DO i=N-1,1,-1
        if(IRC1.le.i+1) then
!Here we allow to search for the radiative/convective boundary in the zone where gradrad>0.7*gradad
!The RC boundary is still set where gradrad>gradad
         if(gradrad(i).gt.0.7*gradad(i)) then
         iRC1=i
         endif
         if(gradrad(i).gt.0.98*gradad(i)) then
          iRC=i
          PRC=P(iRC)
          TAURC=TAU(iRC)
         endif
        endif 
       ENDDO
       if(iRC.lt.N) then
        DO i=iRC,N-1
         gradad(i)=0.32-0.10*T(i)/3000
         if(gradad(i).lt.0) then
          write(*,*)'The adiabatic gradient is negative, probably becaus&
     &e temperatures are too high (T=',T(i),', P=',P(i),'. It is set to &
     &zero'
          gradad(i)=0
         endif
         T(i+1)=T(i)*(P(i+1)/P(i))**gradad(i)
        ENDDO
       endif
      else
      !ALL NaN
       PRC=0
       PRC=PRC/PRC
       iRC=PRC
       TAURC=PRC
!      else 
!       print*,"CONV must be either 'YES' or 'NO'"
!       STOP
      endif   
      END SUBROUTINE tprofile
