  module redshift
!-------------------------------------------
! redshift.f95
! module redshift
!
!-------------------------------------------

  implicit none
  
  public   ::  redshif,backshif
  contains   
  

    subroutine redshif(a,ue,l,q2,gfac,thetem,phiem)
      use const
      implicit none
      double precision  ::  a,ue,re,l,q2,gfac,thetem,phiem
      double precision  ::  domega,omega,da,dsig,ddel,dvel,q,miem
      double precision, parameter  ::  one=1.0, two=2.0
      
      re = one/ue
      ddel = re**two-two*re+a**two
      dsig = re**two
      da = (re**two+a**two)**two-a**two*ddel
      omega = two*a*re/da
      domega = one/(a+re**(3.0/two))
      dvel = (domega-omega)*da/(dsig*dsqrt(ddel))
      
      gfac = dsqrt(dsig*ddel/da)*(one-dvel**two)/(one-domega*l)
      q = dsqrt(q2)
      thetem = dacos(q*gfac/re)
      miem = q*gfac/re
!      phiem = dasin((l*dsig*dsqrt(ddel)*da**(-one)+dvel*(l*domega-one))/(dsin(thetem)*(one-l*domega))) 
      phiem = dasin((l*dsig*ddel**(one/two)*da**(-one)+l*omega*dvel-dvel)/(dsqrt(one-miem**two)*(one-l*domega)))
     
      return
    end subroutine

    subroutine backshif(a,ue,nuthem,nuphem,l,q2,gfac)
      use const
      implicit none
      double precision  ::  a,ue,re,l,q2,gfac,nuthem,nuphem
      double precision  ::  domega,omega,da,dsig,ddel,dvel,q,miem ! to miem to w ogole co innego niz miem wyzej
      double precision, parameter  ::  one=1.0, two=2.0
      
      re = one/ue
      ddel = re**two-two*re+a**two
      dsig = re**two
      da = (re**two+a**two)**two-a**two*ddel
      omega = two*a*re/da
      domega = one/(a+re**(3.0/two))
      dvel = (domega-omega)*da/(dsig*dsqrt(ddel))
      print*,"nuth:",nuthem,"nuph:",nuphem
      miem = nuphem*dsqrt(one-nuthem**two)
      l = (miem+dvel)/(domega*miem+dsig*dsqrt(ddel)*da**(-one)+omega*dvel)
!      l = 0.0
      gfac = dsqrt(dsig*ddel/da)*(one-dvel**two)/(one-domega*l)
      q = (re*nuthem/gfac)**two      
      q2 = q
      print*,"gfac:",gfac,"l:",l,"q:",q,"q2:",q2
      pause
      return
    end subroutine
    
  end module redshift
