  module funkcjeeliptycznejacobiego
!-------------------------------------------
! jacobifunct.f95
! module jacobifunct
! Potem to opisze
!-------------------------------------------

  implicit none
  
  public   ::  sncndn
  contains   

    subroutine sncndn(uu,emmc,sn,cn,dn)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza wartosci funkcji eliptycznych Jacobiego (sn,cn,dn) dla podanych wartosci
! u i k_{c}^{2} oznaczonych tutaj przez uu i emmc.
! Procedura zaczerpnieta z Numerical Recipies in Fortran Press et al (1992) i przepisana w 
! lepszej formie (tzn. nie ma przeskokow (go to) w procedurze tylko jest rozsadnie napisany
! warunek if dajacy ten sam rezultat)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  uu,cn,dn,sn,emmc
!-------------------------------zmienne lokalne------------------------------------!
    double precision, parameter  ::  CA=0.0003
    integer  ::  i,ii,l
    double precision  ::  a1,b1,c1,d1,emc,u
    double precision, dimension(13)  ::  em,en
    logical  ::  bo
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
      
    emc = emmc
    u = uu
    if (emc/=zero) then
      bo = (emc<zero)
      if (bo) then
        d1 = one-emc
        emc = -emc/d1
        d1 = dsqrt(d1)
        u = u*d1
      endif
      a1 = one
      dn = one
      do i=1,13
        l = i
        em(i) = a1
        emc = dsqrt(emc)
        en(i) = emc
        c1 = hf*(a1+emc)
        if (dabs(a1-emc)<(CA*a1)) then
          exit
        endif
        emc = a1*emc
        a1 = c1
      enddo
      u = c1*u
      sn = dsin(u)
      cn = dcos(u)
      if (sn==zero) then
        if (bo) then
          a1 = dn
          dn = cn
          cn = a1
          sn = sn/d1
        endif
      else
        a1 = cn/sn
        c1 = a1*c1
        do ii=l,1,-1
          b1 = em(ii)
          a1 = c1*a1
          c1 = dn*c1
          dn = (en(ii)+a1)/(a1+b1)
          a1 = c1/b1
        enddo
        a1 = one/dsqrt(c1**two+one)
        if (sn<zero) then
          sn = -a1
        else
          sn = a1
        endif
        cn = c1*sn
        if (bo) then
          a1 = dn
          dn = cn
          cn = a1
          sn = sn/d1
        endif
      endif
    else
      cn = one/dcosh(u)
      dn = cn
      sn = dtanh(u)
    endif
    return       
  
    end subroutine sncndn

  end module funkcjeeliptycznejacobiego
