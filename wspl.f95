  module zmianazmiennej
!-------------------------------------------
! cormuf.f95
! module cormuf
!
!-------------------------------------------

  implicit none
  
  public   ::  zmiennaniezalezna
  contains   
  

    subroutine zmiennaniezalezna(muminus,muplus,mup,muk,mutp0,mutp,zcmu,k,lpptp,ofset,muv,mutp1)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza lpt punktow na trajektorii okreslona przez wartosci q2 i l w metryce Kerra o spinie a
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------Opis wszystkich zmiennych znajduje sie w pliku zmienne_i_parametry---------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  muminus,muplus,mup,muk,zcmu,ofset,muv
    integer  ::  mutp0,mutp,k,lpptp,mutp1     ! mutp0 jest zmienna lokalna, nie ma jej w pliku zmienne_i_parametry
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  dmu,dmu1,dmu2,dmuk,mu1,mu2,mub
    double precision  ::  kmax1,kmax2,keff
    integer  ::  a1,a2,a3,a4,dtpm,dmukdmu1
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
! Najpierw obliczam calkowita dlugosc trajektorii      
      dtpm = mutp-mutp0
      a1 = int(zcmu*(-one)**mutp0)
      a2 = int(zcmu*(-one)**mutp)
! SPrawdzam czy kolejny turning point napotkany to muminus lub muplus
      mu1 = (a1+one)*muplus/two+(one-a1)*muminus/two
      mu2 = (one-a1)*muplus/two+(a1+one)*muminus/two
      dmu1 = dabs(mu1-mup)
      dmu2 = dabs(mu2-mu1)
      a3 = 2*int((two*dble(dtpm)+try-a1)/4.0)-1
      dmu = dble(a1)*(muplus-mup)+dble(a2)*(muk-muminus)+dble(a3)*(muplus-muminus)
! Teraz obliczam liczbe turning pointow pomiedzy mup, a aktualnym indeksem (k)
      dmuk = dble(k)*dmu/dble(lpptp)
      dmukdmu1 = int(dmuk/dmu1)
      if (int(dmukdmu1)<0) then
        dmukdmu1 = 50
      endif
      mutp1 = min(dmukdmu1,1)+max(int((dmuk-dmu1)/dmu2),0)
      a4 = int(zcmu*(-one)**(mutp1+mutp0))
      mub = mu1*mod(mutp1,2)+mup*((sign(one,dble(hf-mutp1))+one)/two)+mu2*mod(mutp1-1,2)*(sign(one,dble(mutp1-hf))+one)/two
      kmax1 = int(dmu1*lpptp/dmu)+1
      kmax2 = int(dmu2*lpptp/dmu)+1
      keff = k-(sign(one,k-kmax1)+one)*int(dmu1*lpptp/dmu)/two-(int(dmu2*lpptp/dmu)+1)*(int((k-kmax1)/dble(kmax2))&
     &*(sign(one,k-kmax1-kmax2)+one)/two)
      muv = mub+a4*(keff-ofset)*dmu/dble(lpptp)
      mutp1 = mutp1+mutp0
      return

    end subroutine

  end module zmianazmiennej


