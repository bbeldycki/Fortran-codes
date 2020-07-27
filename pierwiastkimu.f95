  module pierwiastkimu
!-------------------------------------------
! mmuroots.f95
! module mmuroots
! Opis zrobie pozniej
!-------------------------------------------

  implicit none
  
  public   ::  rozwiazaniamu
  contains   
  

    subroutine rozwiazaniamu(q2,l,mup,muminus,muplus)
    use const
    implicit none
! Cel procedury:
! Procedura ta oblicza pierwiastki rownania M(mu) w przypadku, gdy a i q2 sa rozne od zera
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------Opis wszystkich zmiennych znajduje sie w pliku zmienne_i_parametry---------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  q2,l,mup,muminus,muplus
    double precision  ::  mplus,mminus,rozw,ql2     ! zmienne pomocnicze lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!  
    ql2 = q2+l**two
    rozw = -hf*(a**two-ql2+dsign(one,a**two-ql2)*dsqrt((a**two-ql2)**two+4.0*q2*a**two))
    if ((a**two-ql2)<zero) then
      mminus = -one*rozw/a**two
      mplus = q2/rozw
    else
      mminus = q2/rozw
      mplus = -one*rozw/a**two
    endif
    mplus = dmin1(mplus,one)
    if (mminus>zero) then
! Przypadek, gdzie mamy asymetryczne pierwiastki
      muplus = dsign(one,mup)*dsqrt(mplus)
      muminus = dsign(one,mup)*dsqrt(mminus)
      if (dabs(muminus)>dabs(mup)) then
        muminus = mup
      endif
      if (dabs(mup)>dabs(muplus)) then
        muplus = mup
      endif
    else
! Przypadek, gdzie mamy symetryczne pierwiastki. W tym przypadku orbita moze przecinac plaszczyzne dysku
      muplus = dsqrt(mplus)
      if (muplus<mup) then
        muplus = mup
      endif
      muminus = -muplus
    endif
    return
    end subroutine

  end module pierwiastkimu


