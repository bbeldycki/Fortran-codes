  module funkcjeCarlsona
!-------------------------------------------
! functions.f95
! module general_functions
! Module contains all necessary functions used in other modules
! Potem to opisze bardziej
!-------------------------------------------

  implicit none
  
  public   ::  rf,rd,rj,rc,InvHipSec
  contains   


    double precision function rf(x,y,z)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosc calki eliptycznej (Carlsona) pierwszego rodzaju R_{F}
! R_{F}(x,y,z) = 1/2 \int from 0 do infinity dt (t+x)^{-1/2}(t+y)^{-1/2}(t+z)^{-1/2}
! Wartosci x, y i z musza byc nieujemne i conajwyzej jedna z nich moze byc rowna 0
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisana
! petla do while dajaca ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  x,y,z
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  ERRTOL=0.08, TINY=1.5d-38, BIG=3.0d37
    double precision, parameter  ::  C1=1.0/24.0, C2=0.1, C3=3.0/44.0, C4=1.0/14.0
    double precision  ::  alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    
    if (min(x,y,z)<zero .or. min(x+y,x+z,y+z)<TINY .or. max(x,y,z)>BIG) then
      print*,"Invalid arguments in rf, ERROR !!!"
    endif
    xt = x
    yt = y
    zt = z

    sqrtx = dsqrt(xt)
    sqrty = dsqrt(yt)
    sqrtz = dsqrt(zt)
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    xt = 0.25*(xt+alamb)
    yt = 0.25*(yt+alamb)
    zt = 0.25*(zt+alamb)
    ave = trd*(xt+yt+zt)
    if (ave==zero) then
      delx = zero
      dely = zero
      delz = zero
    else
      delx = (ave-xt)/ave
      dely = (ave-yt)/ave
      delz = (ave-zt)/ave
    endif

    do while (dmax1(dabs(delx),dabs(dely),dabs(delz))>ERRTOL)
      sqrtx = dsqrt(xt)
      sqrty = dsqrt(yt)
      sqrtz = dsqrt(zt)
      alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt = 0.25*(xt+alamb)
      yt = 0.25*(yt+alamb)
      zt = 0.25*(zt+alamb)
      ave = trd*(xt+yt+zt)
      if (ave==zero) then
        delx = zero
        dely = zero
        delz = zero
      else
        delx = (ave-xt)/ave
        dely = (ave-yt)/ave
        delz = (ave-zt)/ave
      endif
    enddo
    e2 = delx*dely-delz**two
    e3 = delx*dely*delz
    rf = (one+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)
    return
    end function


    double precision function rd(x,y,z)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosc calki eliptycznej (Carlsona) drugiego rodzaju R_{D}
! R_{D}(x,y,z) = 3/2 \int from 0 do infinity dt (t+x)^{-1/2}(t+y)^{-1/2}(t+z)^{-3/2}
! Wartosci x i y musza byc nieujemne i conajwyzej jedna z nich moze byc rowna 0
! Wartosc z zawsze musi byc dodatnia
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisana
! petla do while dajaca ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  x,y,z
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  ERRTOL=0.05, TINY=1.e-25, BIG=4.5e21, C1=3.0/14.0, C2=1.0/6.0
    double precision, parameter  ::  C3=9.0/22.0, C4=3.0/26.0, C5=0.25*C3, C6=1.5*C4
    double precision  ::  alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sums,xt,yt,zt
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!

    if (min(x,y)<0.d0 .or. min(x+y,z)<TINY .or. max(x,y,z)>BIG) then
      print*,"Invalid arguments in rd, ERROR !!!"
    endif
    xt = x
    yt = y
    zt = z
    sums = zero
    fac = one

    sqrtx = dsqrt(xt)
    sqrty = dsqrt(yt)
    sqrtz = dsqrt(zt)
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    sums = sums+fac/(sqrtz*(zt+alamb))
    fac = 0.25*fac
    xt = 0.25*(xt+alamb)
    yt = 0.25*(yt+alamb)
    zt = 0.25*(zt+alamb)
    ave = 0.2*(xt+yt+try*zt)
    delx = (ave-xt)/ave
    dely = (ave-yt)/ave
    delz = (ave-zt)/ave
 
    do while (dmax1(dabs(delx),dabs(dely),dabs(delz))>ERRTOL)
      sqrtx = dsqrt(xt)
      sqrty = dsqrt(yt)
      sqrtz = dsqrt(zt)
      alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      sums = sums+fac/(sqrtz*(zt+alamb))
      fac = 0.25*fac
      xt = 0.25*(xt+alamb)
      yt = 0.25*(yt+alamb)
      zt = 0.25*(zt+alamb)
      ave = 0.2*(xt+yt+try*zt)
      delx = (ave-xt)/ave
      dely = (ave-yt)/ave
      delz = (ave-zt)/ave
    enddo

    ea = delx*dely
    eb = delz*delz
    ec = ea-eb
    ed = ea-6.0*eb
    ee = ed+ec+ec
    rd = try*sums+fac*(one+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/&
   &(ave*dsqrt(ave))
    return
    end function


    double precision function rc(x,y)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosc calki eliptycznej (Carlsona) trzeciego rodzaju R_{C}
! R_{C}(x,y) = 1/2 \int from 0 do infinity dt (t+x)^{-1/2}(t+y)^{-1}
! Wartosc x musi byc nieujemna, wartosc y musi byc rozna od 0
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisana
! petla do while dajaca ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  x,y
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  ERRTOL=0.04, TINY=1.69e-38, BIG=3.0e38, SQRTNY=1.3e-19
    double precision, parameter  ::  TNBG=TINY*BIG, COMP1=2.236/SQRTNY, COMP2=TNBG*TNBG/25.0
    double precision, parameter  ::  C1=0.3, C2=1.0/7.0, C3=0.375, C4=9.0/22.0
    double precision  ::  alamb,ave,s,w,xt,yt
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    if (x<zero .or. y==zero .or. (x+dabs(y))<TINY .or. (x+dabs(y))>BIG .or. &
   &(y<(-COMP1) .and. x>zero .and. x<COMP2)) then    
      print*,"Invalid arguments in rc, ERROR !!!"
    endif
    if (y>zero) then
      xt = x
      yt = y
      w = one
    else
      xt = x-y
      yt = -y
      w = dsqrt(x)/dsqrt(xt)
    endif
    alamb = two*dsqrt(xt)*dsqrt(yt)+yt
    xt = 0.25*(xt+alamb)
    yt = 0.25*(yt+alamb)
    ave = trd*(xt+yt+yt)
    s = (yt-ave)/ave
    do while (dabs(s)>ERRTOL)
      alamb = two*dsqrt(xt)*dsqrt(yt)+yt
      xt = 0.25*(xt+alamb)
      yt = 0.25*(yt+alamb)
      ave = trd*(xt+yt+yt)
      s = (yt-ave)/ave
    enddo      
    rc = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/dsqrt(ave)
    return
    end function


    double precision function rj(x,y,z,p)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosc calki eliptycznej (Carlsona) trzeciego rodzaju R_{J}
! R_{J}(x,y,z,p) = 3/2 \int from 0 do infinity dt (t+x)^{-1/2}(t+y)^{-1/2}(t+z)^{-1/2}(t+p)^{-1}
! Wartosci x, y i z musza byc nieujemne i conajwyzej jedna moze byc rowna 0
! p musi byc niezerowe, jezeli jest mniejsze od 0 to zwracana jest wartosc glowna "Cauchy'ego"
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisana
! petla do while dajaca ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  x,y,z,p
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  ERRTOL=0.05, TINY=2.5e-13, BIG=9.0e11, C1=3.0/14.0, C2=1.0/3.0
    double precision, parameter  ::  C3=3.0/22.0, C4=3.0/26.0, C5=0.75*C3, C6=1.5*C4
    double precision, parameter  ::  C7=0.5*C2, C8=C3+C3
    double precision  ::  a1,alamb,alpha,ave,b1,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,fac
    double precision  ::  pt,rcx,rho,sqrtx,sqrty,sqrtz,sums,tau,xt,yt,zt
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    if (dmin1(x,y,z)<zero .or. dmin1(x+y,x+z,y+z,dabs(p))<TINY .or. dmax1(x,y,z,dabs(p))>BIG) then
      print*,"Invalid arguments in rj, ERROR !!!"
    endif
    sums = zero
    fac = one
    if (p>zero) then
      xt = x
      yt = y
      zt = z
      pt = p
    else
      xt = dmin1(x,y,z)
      zt = dmax1(x,y,z)
      yt = x+y+z-xt-zt
      a1 = one/(yt-p)
      b1 = a1*(zt-yt)*(yt-xt)
      pt = yt+b1
      rho = xt*zt/yt
      tau = p*pt/yt
      rcx = rc(rho,tau)
    endif
    sqrtx = dsqrt(xt)
    sqrty = dsqrt(yt)
    sqrtz = dsqrt(zt)
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    alpha = (pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**two
    beta = pt*(pt+alamb)**two
    sums = sums+fac*rc(alpha,beta)
    fac = 0.25*fac
    xt = 0.25*(xt+alamb)
    yt = 0.25*(yt+alamb)
    zt = 0.25*(zt+alamb)
    pt = 0.25*(pt+alamb)
    ave = 0.2*(xt+yt+zt+pt+pt)
    delx = (ave-xt)/ave
    dely = (ave-yt)/ave
    delz = (ave-zt)/ave
    delp = (ave-pt)/ave
    do while (max(dabs(delx),dabs(dely),dabs(delz),dabs(delp))>ERRTOL)
      sqrtx = dsqrt(xt)
      sqrty = dsqrt(yt)
      sqrtz = dsqrt(zt)
      alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      alpha = (pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**two
      beta = pt*(pt+alamb)**two
      sums = sums+fac*rc(alpha,beta)
      fac = 0.25*fac
      xt = 0.25*(xt+alamb)
      yt = 0.25*(yt+alamb)
      zt = 0.25*(zt+alamb)
      pt = 0.25*(pt+alamb)
      ave = 0.2*(xt+yt+zt+pt+pt)
      delx = (ave-xt)/ave
      dely = (ave-yt)/ave
      delz = (ave-zt)/ave
      delp = (ave-pt)/ave
    enddo
      
    ea = delx*(dely+delz)+dely*delz
    eb = delx*dely*delz
    ec = delp**two
    ed = ea-try*ec
    ee = eb+two*delp*(ea-ec)
    rj = try*sums+fac*(one+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+delp*ea*&
   &(C2-delp*C3)-C2*delp*ec)/(ave*dsqrt(ave))
    if (p<zero) then
      rj = a1*(b1*rj+try*(rcx-rf(xt,yt,zt)))
    endif
    return
    end function


    double precision function InvHipSec(x)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosc Inverse hyperbolic secant
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  x     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    if (x>zero) then
      InvHipSec = dlog((one+dsqrt(one-x*x))/x)
    else 
      InvHipSec = zero
    endif
    return      
    end function

  end module funkcjeCarlsona
