  module przesuniecia
!-------------------------------------------
! geoproc.f95
! module geoproc
!
!-------------------------------------------

  implicit none
  
  public   ::  efektydp
  contains   
  
    subroutine efektydp(uemisji,l,q2,gfaktor,thetaemisji,phiemisji)
      use const
      implicit none
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
      double precision  ::  uemisji,l,q2,gfaktor,thetaemisji,phiemisji!,gfac
!-------------------------------zmienne lokalne------------------------------------!
      double precision  ::  domega,omega,da,dsig,ddel,dvel,q,miem,re
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
      
      re = one/uemisji
      ddel = re**two-two*re+a**two
      dsig = re**two
      da = (re**two+a**two)**two-a**two*ddel
      omega = two*a*re/da
      domega = one/(a+re**ohf)
      dvel = (domega-omega)*da/(dsig*dsqrt(ddel))
      
      gfaktor = dsqrt(dsig*ddel/da)*(one-dvel**two)**hf/(one-domega*l)
!      gfac = (re**hf*(re**two-try*re+two*a*re**hf)**hf)/(re**ohf+a-l)
      q = dsqrt(q2)
      thetaemisji = dacos(q*gfaktor/re)
      miem = q*gfaktor/re
!      phiem = dasin((l*dsig*dsqrt(ddel)*da**(-one)+dvel*(l*domega-one))/(dsin(thetem)*(one-l*domega))) 
      phiemisji = dasin((l*dsig*ddel**(one/two)*da**(-one)+l*omega*dvel-dvel)/(dsqrt(one-miem**two)*(one-l*domega)))
     
      return
    end subroutine

  end module przesuniecia
