       SUBROUTINE valencia(P,T,met,Kappa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine implement the fit of the Rosseland mean opacities from Valencia et al. 2013
!We added the condition that Kappa cannot be bigger than 10^4 m^2/kg in order to avoid divergence in the integration of the Tau(P) function. In any case, for kappa>10^4m^2/kg, the atmosphere should be convective and the radiative solution becomes irrelevant.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      REAL*8,INTENT(IN)::P !Pressure in Pa
      REAL*8,INTENT(IN)::T ! Temperature in K
      REAL*8,INTENT(IN)::met !metalicity

      REAL*8,INTENT(OUT)::Kappa  ! Rosseland mean Opacity (m2/kg)

!Now the internal variables
      REAL*8 :: logP
      REAL*8 :: logT
      REAL*8 :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11
      REAL*8 :: KlowP,KhighP
      REAL*8,PARAMETER :: Pi=3.1415926535
      logP=log10(P)+1 !Convert Pascals into Dynes/cm2
      logT=log10(T)
!      logT=5

      c1=-37.5
      c2=0.00105
      c3=3.2610
      c4=0.84315
      c5=-2.339

       if(T.lt.800) then
       c6=-14.051
       c7=3.055
       c8=0.024
       c9=1.877
       c10=-0.445
       c11=0.8321
      else
       c6=82.241
       c7=-55.456
       c8=8.754
       c9=0.7048
       c10=-0.0414
       c11=0.8321
      endif
      KlowP=c1*(logT-c2*logP-c3)**2+(c4*met+c5)
      KhighP=(c6+c7*logT+c8*logT**2)+logP*(c9+c10*logT)+met*c11*        &
     &(0.5+1/Pi*atan((logT-2.5)/0.2))
      KlowP=10**KlowP
      KhighP=10**KhighP
      Kappa=KlowP+KhighP
      Kappa=Kappa/10 !(convert from cm2/g to m2/kg)

      ! In order to avoid divergence when Kappa becomes two big, we modified the fit such that Kappa can't be bigger than 10000 m^2/kg
      if(Kappa.gt.1e10) then 
      Kappa=1e10
      endif
      RETURN
      END SUBROUTINE valencia
