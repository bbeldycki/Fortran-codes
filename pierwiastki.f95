  module pierwiastki
!-------------------------------------------
! roots.f95
! module roots
! Opis zrobie pozniej
!-------------------------------------------

  implicit none
  
  public   ::  gauleg,laguer,zroots
  contains   
  

    subroutine gauleg(x1,x2,x,w,n)
    use const
    implicit none
! Cel procedury:
! Procedura zwraca tablice wag i "odcietych" dla n punktowej kwadratury Gaussa-Legendrea
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisana
! petla do while dajaca ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    integer n
    double precision  ::  x1,x2
    double precision, dimension(n)  ::  x,w
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  EPS=3.0d-14, quater=0.25
    integer  ::  i,j,m
    double precision  ::  p1,p2,p3,pp,xl,xm,z,z1
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    m = (n+1)/2
    xm = hf*(x2+x1)
    xl = hf*(x2-x1)
    do i=1,m
      z = dcos(pi*(i-quater)/(n+hf))
      do while (dabs(z-z1)>EPS)
        p1 = one
        p2 = zero
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((two*j-one)*z*p2-(j-one)*p3)/j
        enddo
        pp = n*(z*p1-p2)/(z*z-one)
        z1 = z
        z = z1-p1/pp
      enddo
      x(i) = xm-xl*z
      x(n+1-i) = xm+xl*z
      w(i) = two*xl/((one-z*z)*pp*pp)
      w(n+1-i) = w(i)
    enddo
    return
    end subroutine


    subroutine laguer(ax,m,x,its)
    use const
    implicit none
! Cel procedury:
! Opis zastosowania jest w ksiazce (patrz linijke nizej)
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisany
! warunek if dajacy ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    integer  ::  m,its
    double complex  ::  x
    double complex, dimension(m+1)  ::  ax
!-------------------------------zmienne lokalne------------------------------------!
    integer, parameter  ::  MR=8, MT=10, MAXIT=MT*MR
    double precision, parameter  ::  EPSS=1.0d-10
    integer  ::  iter,j
    double precision  ::  abx,abp,abm,err
    double precision, dimension(MR) ::  frac
    double complex  ::  dx,x1,b,d,f,g,h,sq,gp,gm,g2
    data frac /hf, 0.25d0, 0.75d0, 0.13d0, 0.38d0, 0.62d0, 0.88d0, one/
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!

    do iter=1,MAXIT
      its = iter
      b = ax(m+1)
      err = dble(cdabs(b))
      d = dcmplx(zero,zero)
      f = dcmplx(zero,zero)
      abx = dble(cdabs(x))
      do j=m,1,-1
        f = x*f+d
        d = x*d+b
        b = x*b+ax(j)
        err = dble(cdabs(b))+abx*err
      enddo
      err = EPSS*err
      if (dble(cdabs(b))<=err) then
        return
      else
        g = d/b
        g2 = g*g
        h = g2-two*f/b
        sq = cdsqrt((m-1)*(m*h-g2))
        gp = g+sq
        gm = g-sq
        abp = dble(cdabs(gp))
        abm = dble(cdabs(gm))
        if (abp<abm) then
          gp = gm
        endif
        if (max(abp,abm)>zero) then
          dx = m/gp
        else
          dx = cdexp(dcmplx(dlog(one+abx),dble(iter)))
        endif
      endif
      x1 = x-dx
      if (x==x1) then
        return
      endif
      if (mod(iter,MT)/=0) then
        x = x1
      else
        x = x-dx*frac(iter/MT)
      endif
    enddo
    print*,"too many iterations in laguer"
    return
    end subroutine


    subroutine zroots(ax,m,roots,polish)
    use const
    implicit none
! Cel procedury:
! Opis zastosowania jest w ksiazce (patrz linijke nizej)
! Funkcja zaczerpnieta z Numerical Recipies in Fortran Press et al (1992)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    integer  ::  m
    double complex, dimension(m)  ::  roots
    double complex, dimension(m+1)  ::  ax
    logical  ::  polish
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  EPS=1.e-6
    integer, parameter  ::  MAXM=10
    integer  ::  i,j,jj,its,nc1
    double complex  ::  x,b,c
    double complex, dimension(MAXM)  ::  ad
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!

    if (m>MAXM-1) then
      print*,"M too large in ZROOTS"
    endif
    do j=1,m+1
      ad(j) = ax(j)
    enddo
    do j=m,1,-1
      x = dcmplx(zero,zero)
      call laguer(ad,j,x,its)
      if (dabs(aimag(x))<=(two*EPS**two*dabs(dble(x)))) then
        x = dcmplx(dble(x),zero)
      endif
      roots(j) = x
      b = ad(j+1)
      do jj=j,1,-1
        c = ad(jj)
        ad(jj) = b
        b = x*b+c
      enddo
    enddo
    if (polish) then
      do j=1,m
        call laguer(ax,m,roots(j),its)
      enddo
    endif
    do j=2,m
      x = roots(j)
      do i=j-1,1,-1
        if (dble(roots(i))<=dble(x)) then
          nc1 = 1
          exit
        endif
        roots(i+1) = roots(i)
      enddo
      if (nc1/=1) then
        i = 0
      endif
      roots(i+1) = x
    enddo
    return
    end subroutine
      
    
  end module pierwiastki
