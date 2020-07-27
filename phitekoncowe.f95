  module obliczphite
!-------------------------------------------
! calphite.f95
! module caluf
! Potem opisze to
!-------------------------------------------

  implicit none
  
  public   ::  geodezyjnaphite
  contains   

    subroutine geodezyjnaphite(up,uk,mup,muk,l,l2,q2,mutp,rtp,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,dphiu,&
   &dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,ctmupplus,ctmutp,&
   &cphimupplus,cphimutp,elcaltp)
    use const
    use calki
    implicit none
! Cel procedury:
! Procedura oblicza lpt punktow na trajektorii okreslona przez wartosci q2 i l w metryce Kerra o spinie a

! Parametry wejsciowe procedury:
! up - poczatkowa wartosc u. Jezeli up=0, to wartosci tki i lambdaki beda nieskonczone.
! uk - koncowa wartosc u. Jezeli czrad=.TRUE. to policzona, w przeciwnym przypadku musi byc zadana
! upom - policzone wartosci u sa pomiedzy upom i uk. Dla pelnej geodezyjnej ustawiamy up = upom
! mup - poczatkowa wartosc mu
! muk - koncowa wartosc mu, jezeli czrad=.TRUE. to musi byc podana, w przeciwnym razie jest policzona
! a - spin czarnej dziury
! l - bezwymiarowa wartosc momentu pedu wzdluz osi z
! q2 - bezwymiarowa wartosc starej Cartera
! alfa - parametr zderzenia w nieskonczonosci (prostopadly do osi obrotu czarnej dziury). Dla alfa=0 mamy l=0
! beta - parametr zderzenia w nieskonczonosci (rownolegly do osi obrotu czarnej dziury). Dla beta=0 i mup=0
! mamy orbite w plaszczyznie dyskowej
! mutp - liczba mu turning pointow, policzona w przypadku gdy czrad=.FALSE., musi byc podana w przeciwnym przypadku
! zcu - znak calki u (du/dlmabda), = 1/-1 dla przychodzacych/odchodzacych promieni
! zcmu - znac calki mu (dmu/dlambda), moze byc policzony jesli znamy mup, beta i wiemy, ze wartosc 1/up to duza liczba

! Parametry wyjsciowe procedury:
! dtu - wartosc calki (wzor 47 Praca Dextera)
! dtmu - wartosc calki (wzor 45 Praca Dextera)
! dphiu - wartosc calki (wzor 46 Praca Dextera)
! dphimu - wartosc calki (wzor 48 Praca Dextera)
! rd2r - wartosc rd dla zupelnej calki eliptycznej drugiego rodzaju
! rj3r - wartosc rj dla zupelnej calki eliptycznej trzeciego rodzaju
! lambda - wartosc parametru afinicznego
! ctuw1 - wartosc calki od up (czwarty wyraz z wzoru (47) Praca Dextera). Obliczona gdy elcaltp=.TRUE., jezeli 
! nie ma fizycznych turning pointow to zmienna zostanie pominieta. W przeciwnym przypadku dana jako input.
! ctuw2 - wartosc calki od up (trzeci wyraz z wzoru (47) Praca Dextera). Obliczona gdy elcaltp=.TRUE., jezeli 
! nie ma fizycznych turning pointow to zmienna zostanie pominieta. W przeciwnym przypadku dana jako input.
! ctuw3 - wartosc calki od up (pierwszy wyraz z wzoru (47) Praca Dextera). Obliczona gdy elcaltp=.TRUE., jezeli 
! nie ma fizycznych turning pointow to zmienna zostanie pominieta. W przeciwnym przypadku dana jako input.
! ctuw4 - wartosc calki od up (drugi wyraz z wzoru (47) Praca Dextera). Obliczona gdy elcaltp=.TRUE., jezeli 
! nie ma fizycznych turning pointow to zmienna zostanie pominieta. W przeciwnym przypadku dana jako input.
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  up,uk,mup,muk,l,l2,q2,zcu,zcmu,cu,h1,u1,u2,u3,u4,dphimu,dtmu,dphiu,dtu,lambda,rfdcup
    double precision  ::  rfdcuk,ctmupplus,ctmutp
    double precision  ::  rfdcmup,rfdcmuk,rfdcmutp
    double precision  ::  rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,cphimupplus,cphimutp
    integer  ::  mutp,rtp,np
    logical  ::  elcaltp
    integer, dimension(5)  ::  p
    double precision  ::  mne,mpo,upl,yy,pi2,ql2,s1,iumin,ee,dd,f,g,h,f1,f2,g1,g2,h2
    double precision  ::  phiu0,phiu01,phiu02,cphimukplus,ctmukplus!,rfr1,rfr1pplus,rfr1kplus
    double precision  ::  phiu11,phiu12,phiu1,phiu2
    double precision  ::  tu0,tu11,tu12,tu13,tu14,tu1,tu2,tu3,tu4
    double precision  ::  muplus,vfm,vfp,vsm,vsp,qs,ur,a1,a2,a3
    double precision  ::  lambdau
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    pi2 = two*pi
    p(1) = -1
    p(2) = -1
    p(3) = -1
    upl = one/(one+dsqrt(one-a**two))
! Policze teraz 1/u_{-}, bo u_{-} jest rozbiezne gdy a=0
    iumin = one-dsqrt(one-a**two)
    ur = -one/(two*dsqrt(one-a**two))
    qs = dsign(one,q2)
    dd = two*((a-l)**two+q2)
    ee = -a**two*q2
    ql2 = q2+l2
    a1 = zcmu
    a2 = zcmu*(-one)**mutp
    a3 = two*int((two*dble(mutp)+try-zcmu)/4.0)-one
    if (np==0) then
      dtu = zero
      dtmu = zero
      dphiu = zero
      dphimu = zero
      lambdau = zero
    elseif (np<3) then
!----------------------------------------------------------------------------------!
! Tutaj mamy przypadek gdy U(u) ma 3 pierwiastki rzeczywiste
!----------------------------------------------------------------------------------! 
      p(5) = 0
      if (dabs(u3-u2)<(1.d-12)) then
! Przypadek rownych pierwiastkow gdy zachodzi u1<0<u2<=u3       
        dtu = zcu*tuszescian(up,uk,u1,u2,l)/dsqrt(dd)
        dphiu = zcmu*(phiuszescian(uk,u1,u2,l)-phiuszescian(up,u1,u2,l))/dsqrt(dd)
        if (up>=u3) then
          dphiu = -dphiu
          dtu = -dtu
        endif
      elseif (up<=u2) then
!----------------------------------------------------------------------------------!
! Przypadek 1 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
        if (elcaltp.and.(up/=u2)) then
! Pierwsze trzy calki z rownania z wzoru (47) Praca Dextera
          phiu01 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,one,-one/upl,rfdcup,up,u2)
          phiu02 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,one,-iumin,rfdcup,up,u2)
          ctuw2 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,one,rfdcup,up,u2)
          ctuw3 = phiu01
          ctuw4 = phiu02
          p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          ctuw1 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,one,rfdcup,up,u2)  
        elseif (up==u2) then
          ctuw1 = zero
          ctuw2 = zero 
          ctuw3 = zero
          ctuw4 = zero 
          phiu01 = zero
          phiu02 = zero
        else
          phiu01 = ctuw3
          phiu02 = ctuw4
        endif
        if (uk/=u2) then
! Pierwsze trzy calki z rownania z wzoru (47) Praca Dextera        
          phiu11 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,one,-one/upl,rfdcuk,uk,u2)
          phiu12 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,one,-iumin,rfdcuk,uk,u2)
          tu12 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,one,rfdcuk,uk,u2)
          tu13 = phiu11
          tu14 = phiu12
          p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          tu11 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,one,rfdcuk,uk,u2)
        else
          tu11 = zero
          tu12 = zero
          tu13 = zero
          tu14 = zero
          phiu11 = zero
          phiu12 = zero
        endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
! Najpierw zapisze dla czesci zwiazanej z up
        tu0 = (iumin-one/upl)*ctuw1+(iumin**two-one/upl**two)*ctuw2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &ctuw3+(two*a*(a-l)+a**two*iumin+iumin**try)*ctuw4
! A nastepnie zapisze dla czesci zwiazanej z uk
        tu1 = (iumin-one/upl)*tu11+(iumin**two-one/upl**two)*tu12-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu13+(two*a*(a-l)+a**two*iumin+iumin**try)*tu14
! Lacze powyzsze wartosci
        dtu = zcu*ur*(tu0-(-one)**rtp*tu1)/dsqrt(dd)
! Analogicznie robie z rownaniem (48) Praca Dextera
        phiu0 = -(l/upl+two*(a-l))*phiu01+(l*iumin+two*(a-l))*phiu02
        phiu1 = -(l/upl+two*(a-l))*phiu11+(l*iumin+two*(a-l))*phiu12
        dphiu = zcu*ur*(phiu0-(-one)**rtp*phiu1)/dsqrt(dd)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*(ctuw1-(-one)**rtp*tu11)/dsqrt(dd)
      elseif (up>=u3) then
!----------------------------------------------------------------------------------!
! Przypadek 2 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
        if (elcaltp.and.(up/=u3)) then
          p(4) = -2
! Pierwsze trzy calki z rownania z wzoru (47) Praca Dextera
          phiu01 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,one,-one/upl,rfdcup,u3,up)
          phiu02 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,one,-iumin,rfdcup,u3,up)
          ctuw2 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,one,rfdcup,u3,up)
          ctuw3 = phiu01
          ctuw4 = phiu02
          p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          ctuw1 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,one,rfdcup,u3,up)  
        elseif (uk==u3) then
          ctuw1 = zero
          ctuw2 = zero
          ctuw3 = zero
          ctuw4 = zero
          phiu11 = zero
          phiu12 = zero
        else
          phiu01 = ctuw3
          phiu02 = ctuw4
        endif
        if (uk/=u3) then 
          p(4) = -2
! Pierwsze trzy calki z rownania z wzoru (47) Praca Dextera
          phiu11 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,one,-one/upl,rfdcuk,u3,uk)
          phiu12 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,one,-iumin,rfdcuk,u3,uk)
          tu12 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,one,rfdcuk,u3,uk)
          tu13 = phiu01
          tu14 = phiu02
          p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          tu11 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,one,rfdcuk,u3,uk)
        else
          tu11 = zero
          tu12 = zero
          tu13 = zero
          tu14 = zero
          phiu11 = zero
          phiu12 = zero
        endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
! Najpierw zapisze dla czesci zwiazanej z up
        tu0 = (iumin-one/upl)*ctuw1+(iumin**two-one/upl**two)*ctuw2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &ctuw3+(two*a*(a-l)+a**two*iumin+iumin**try)*ctuw4
! A nastepnie zapisze dla czesci zwiazanej z uk
        tu1 = (iumin-one/upl)*tu11+(iumin**two-one/upl**two)*tu12-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu13+(two*a*(a-l)+a**two*iumin+iumin**try)*tu14
! Lacze powyzsze wartosci
        dtu = zcu*ur*(tu0-(-one)**rtp*tu1)/dsqrt(dd)
! Analogicznie robie z rownaniem (48) Praca Dextera
        phiu0 = -(l/upl+two*(a-l))*phiu01+(l*iumin+two*(a-l))*phiu02
        phiu1 = -(l/upl+two*(a-l))*phiu11+(l*iumin+two*(a-l))*phiu12
        dphiu = zcu*ur*(phiu0-(-one)**rtp*phiu1)/dsqrt(dd)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*(ctuw1-(-one)**rtp*tu11)/dsqrt(dd)
      endif      
    elseif (np==3) then
!----------------------------------------------------------------------------------!
! Tutaj mamy przypadek gdy U(u) ma 3 pierwiastki z czego tylko jeden rzeczywisty
! Przypadek 3 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
      np = 3
      f = -one/(dd*u1)
      g = f/u1
      h = one
      if (up<uk) then
        p(4) = -2
! Pierwsze trzy calki z rownania z wzoru (47) Praca Dextera
        tu2 = dwazjedrzecz(p,-u1,one,zero,one,f,g,h,rfdcup,up,uk)
        tu3 = dwazjedrzecz(p,-u1,one,one,-one/upl,f,g,h,rfdcup,up,uk)
        tu4 = dwazjedrzecz(p,-u1,one,one,-iumin,f,g,h,rfdcup,up,uk)
        p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = dwazjedrzecz(p,-u1,one,zero,one,f,g,h,rfdcup,up,uk)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(dd)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(dd)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dd)
      elseif (up>uk) then 
        p(4) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
        tu2 = dwazjedrzecz(p,-u1,one,zero,one,f,g,h,rfdcup,uk,up)
        tu3 = dwazjedrzecz(p,-u1,one,one,-one/upl,f,g,h,rfdcup,uk,up)
        tu4 = dwazjedrzecz(p,-u1,one,one,-iumin,f,g,h,rfdcup,uk,up)
        p(4) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = dwazjedrzecz(p,-u1,one,zero,one,f,g,h,rfdcup,uk,up)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(dd)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(dd)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dd)
      else
        dtu = zero
        dphiu = zero
        lambdau = zero
      endif
    endif
!----------------------------------------------------------------------------------!
! Teraz bedziemy rozwazac przypadeki gdy q2=0 lub gdy a=0
!----------------------------------------------------------------------------------!
    if (q2==zero) then
!----------------------------------------------------------------------------------!
! W tym przypadku t i phi wyrazaja sie przez funkcje elementarne
!----------------------------------------------------------------------------------!
      s1 = dsign(one,mup)
      a1 = s1*zcmu
      a2 = s1*zcmu*(-one)**(mutp+1)
      if (dabs(l)<dabs(a)) then
        muplus = s1*dsqrt(one-l/a**two)
        cphimupplus = dsign(one,muplus*l/a)*(datan((muplus-dtan(dasin(mup/muplus)/two))/dsqrt(one-muplus**two))&
       &+datan((muplus+dtan(dasin(mup/muplus)/two))/dsqrt(one-muplus**two)))
        cphimukplus = dsign(one,muplus*l/a)*(datan((muplus-dtan(dasin(muk/muplus)/two))/dsqrt(one-muplus**two))&
       &+datan((muplus+dtan(dasin(muk/muplus)/two))/dsqrt(one-muplus**two)))
        dtmu = dabs(muplus*a)*(a2*dsqrt(one-muk**two/muplus**two)+a1*dsqrt(one-mup**two/muplus**two))
        dphimu =  a1*cphimupplus+a2*cphimukplus
      else
        dtmu = zero
        dphimu = zero        
      endif
    elseif (a==zero) then
!----------------------------------------------------------------------------------!
! W tym przypadku t i phi wyrazaja sie przez funkcje elementarne rowniez
!----------------------------------------------------------------------------------!
      dtmu = zero
      muplus = dsqrt(q2/ql2)
      vfm = (muk-muplus**two)/((one-muk)*muplus)
      vfp = -(muk+muplus**two)/((one+muk)*muplus)
! komentarz potem
      if (dabs(vfm)>one) then
        vfm = dsign(one,vfm)*one
      endif
      if (dabs(vfp)>one) then
        vfp = dsign(one,vfp)*one
      endif
      if ((mup-muplus**two)==zero) then
        vsm = -one
      else
        vsm = (mup-muplus**two)/((one-mup)*muplus)
      endif
      if ((mup+muplus**two)==zero) then
        vsp = -one
      else
        vsp = -(mup+muplus**two)/((one+mup)*muplus)
      endif
      if (dabs(vsp)>one) then
        vsp = dsign(one,vsp)*one
      endif
      if (dabs(vsm)>one) then
        vsm = dsign(one,vsm)*one
      endif
      cphimupplus = pi-dasin(vsm)+dasin(vsp)
      cphimukplus = pi+dasin(vfm)-dasin(vfp)
      cphimutp = two*pi
      dphimu = -l*cu+dsign(one,l)*hf*(a1*cphimupplus+a2*cphimukplus+a3*cphimutp)
    endif
    if (np==4) then
!----------------------------------------------------------------------------------!
! To tez jest specjalny przypadek gdy mamy q2=0 i a=l, wowczas U(u)=1
! Elementy u calek t i phi wyrazaja sie przez funkcje elementarne
!----------------------------------------------------------------------------------!
! Calka t_u wzor (47) praca Dextera
      dtu = zcu*ur*(iumin-one/upl)*(one/up-one/uk)+(iumin**two-one/upl**two)*dlog(uk/up)+(a**two/upl+&
     &one/upl**try)*upl*dlog((uk/upl-one)/(up/upl-one))-(a**two*iumin+iumin**try)*iumin*dlog((uk*iumin-one)/&
     &(up*iumin-one))
! Calki phi_u wzor (48) praca Dextera
      dphiu = zcu*ur*(l*dlog((uk/upl-one)/(up/upl-one))-l*iumin**two*dlog((uk*iumin-one)/(up*iumin-one)))
! Lambda wzor (49) praca Dextera
      lambdau = zcu*(one/up-one/uk)
    elseif (np==5) then
!----------------------------------------------------------------------------------!
! Tutaj mamy przypadek gdy U(u) ma 4 pierwiastki z czego jedna para pierwiastkow jest zespolona
! Przypadek 5 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
      p(4) = -1
      f = -qs*one/(dabs(ee)*u1*u4)
      g = (u1+u4)*f/(u1*u4)
      h = one
      if (up<uk) then
        p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
        tu2 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,zero,one,f,g,h,rfdcup,up,uk)
        tu3 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,one,-one/upl,f,g,h,rfdcup,up,uk)
        tu4 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,one,-iumin,f,g,h,rfdcup,up,uk)
        p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,zero,one,f,g,h,rfdcup,up,uk)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dabs(ee))      
      elseif (up>uk) then
        p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
        tu2 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,zero,one,f,g,h,rfdcup,uk,up)
        tu3 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,one,-one/upl,f,g,h,rfdcup,uk,up)
        tu4 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,one,-iumin,f,g,h,rfdcup,uk,up)
        p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = dwazdwarzecz(p,-u1,one,u4*qs,-one*qs,zero,one,f,g,h,rfdcup,uk,up)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dabs(ee))    
      else
        dtu = zero
        dphiu = zero
        lambdau = zero
      endif
    elseif (np==6) then
!----------------------------------------------------------------------------------!
! Tutaj mamy przypadek gdy U(u) ma 4 pierwiastki zespolone
! Przypadek 5 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
      p(4) = -1
      h2 = one/h1
      g1 = dd/(ee*(h2-h1))
      g2 = -g1
      f1 = one/dsqrt(ee)
      f2 = f1
      if (up<uk) then
        p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
        tu2 = czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,one,rfdcup,up,uk)
        tu3 = czteryzespolone(p,f1,g1,h1,f2,g2,h2,one,-one/upl,rfdcup,up,uk)
        tu4 = czteryzespolone(p,f1,g1,h1,f2,g2,h2,one,-iumin,rfdcup,up,uk)
        p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,one,rfdcup,up,uk)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dabs(ee))   
      elseif (up>uk) then
        p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
        tu2 = -czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,one,rfdcup,uk,up)
        tu3 = -czteryzespolone(p,f1,g1,h1,f2,g2,h2,one,-one/upl,rfdcup,uk,up)
        tu4 = -czteryzespolone(p,f1,g1,h1,f2,g2,h2,one,-iumin,rfdcup,uk,up)
        p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
        tu1 = -czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,one,rfdcup,uk,up)
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
        dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Minus przez wyrazach phi bierze sie ze zmiany znaku wyrazania (u/u_{\pm} -1), zeby zachowac wszystkie 
! argumenty dodatnie
        phiu0 = -tu3
        phiu1 = -tu4
! Analogicznie robie z rownaniem (48) Praca Dextera
        dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*tu1/dsqrt(dabs(ee))   
      else
        dtu = zero
        dphiu = zero
      endif
    elseif (np>6) then
!----------------------------------------------------------------------------------!
! Tutaj mamy przypadek gdy U(u) ma 4 pierwiastki rzeczywiste
!----------------------------------------------------------------------------------!
      p(4) = -1
      if (dabs(u3-u2)<1.d-12) then
!----------------------------------------------------------------------------------!
! Przypadek gdy mamy dwa rowne pierwiastki rzeczywiste
!----------------------------------------------------------------------------------!
        if (up<uk) then
          if (up<u2) then
!----------------------------------------------------------------------------------!
! Przypadek 7 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
            p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
            phiu1 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/upl,rfdcup,up,uk)
            phiu2 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-iumin,rfdcup,up,uk)
            tu2 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,up,uk)
            tu3 = phiu1
            tu4 = phiu2
            p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
            tu1 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,up,uk)
          elseif (up>u3) then
!----------------------------------------------------------------------------------!
! Przypadek 8 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
            p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
            phiu1 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/upl,rfdcup,up,uk)
            phiu2 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-iumin,rfdcup,up,uk)
            tu2 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,up,uk)
            tu3 = phiu1
            tu4 = phiu2
            p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
            tu1 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,up,uk)
          else
            tu11 = zero
            tu12 = zero
            tu13 = zero
            tu14 = zero
            phiu11 = zero
            phiu12 = zero
          endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
          dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
         &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Analogicznie robie z rownaniem (48) Praca Dextera
          dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
          lambdau = zcu*tu1/dsqrt(-qs*ee)   
        elseif (up>uk) then
          if (up<u2) then
!----------------------------------------------------------------------------------!
! Przypadek 7 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
            p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
            phiu1 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/upl,rfdcup,uk,up)
            phiu2 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-iumin,rfdcup,uk,up)
            tu2 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,uk,up)
            tu3 = -phiu1
            tu4 = -phiu2
            p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
            tu1 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,uk,up)
          elseif (up>u3) then
            p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera           
            phiu1 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/upl,rfdcup,uk,up)
            phiu2 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-iumin,rfdcup,uk,up)
            tu2 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,uk,up)
            tu3 = -phiu1
            tu4 = -phiu2
            p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
            tu1 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,uk,up)
          else
            tu11 = zero
            tu12 = zero
            tu13 = zero
            tu14 = zero
            phiu11 = zero
            phiu12 = zero
          endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
          dtu = ((iumin-one/upl)*tu1+(iumin**two-one/upl**two)*tu2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
         &tu3+(two*a*(a-l)+a**two*iumin+iumin**try)*tu4)*zcu*ur/dsqrt(-qs*ee)
! Analogicznie robie z rownaniem (48) Praca Dextera
          dphiu = zcu*ur*((l/upl+two*(a-l))*phiu0-(l*iumin+two*(a-l))*phiu1)/dsqrt(-qs*ee)
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
          lambdau = zcu*tu1/dsqrt(-qs*ee)  
        else
          dtu = zero
          dphiu = zero
          lambdau = zero
        endif
      elseif (up<=u2) then
!----------------------------------------------------------------------------------!
! Przypadek gdy mamy cztery rozne pierwiastki rzeczywiste
! Przypadek 7 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
        if (elcaltp.and.(up/=u2)) then
! Policze calki w ktorych wystepuje up
! Raz na kazda geodezyjna
          p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
          phiu01 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/upl,rfdcup,up,u2)
          phiu02 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-iumin,rfdcup,up,u2)
          ctuw2 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,up,u2)
          ctuw3 = phiu01
          ctuw4 = phiu02
          p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          ctuw1 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcup,up,u2)
        elseif (up==u2) then
          ctuw1 = zero
          ctuw2 = zero
          ctuw3 = zero
          ctuw4 = zero
          phiu01 = zero
          phiu02 = zero
        else
          phiu01 = ctuw3
          phiu02 = ctuw4
        endif
        if (uk/=u2) then
          p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
          phiu11 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/upl,rfdcuk,uk,u2)
          phiu12 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-iumin,rfdcuk,uk,u2)
          tu12 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcuk,uk,u2)
          tu13 = phiu11
          tu14 = phiu12
          p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          tu11 = -czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,one,rfdcuk,uk,u2)
        else
          tu11 = zero
          tu12 = zero
          tu13 = zero
          tu14 = zero
          phiu11 = zero
          phiu12 = zero
        endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
! Najpierw zapisze dla czesci zwiazanej z up
        tu0 = (iumin-one/upl)*ctuw1+(iumin**two-one/upl**two)*ctuw2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &ctuw3+(two*a*(a-l)+a**two*iumin+iumin**try)*ctuw4
! A nastepnie zapisze dla czesci zwiazanej z uk
        tu1 = (iumin-one/upl)*tu11+(iumin**two-one/upl**two)*tu12-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu13+(two*a*(a-l)+a**two*iumin+iumin**try)*tu14
! Lacze powyzsze wartosci
        dtu = zcu*ur*(tu0-(-one)**rtp*tu1)/dsqrt(dabs(ee))
! Analogicznie robie z rownaniem (48) Praca Dextera
        phiu0 = -(l/upl+two*(a-l))*phiu01+(l*iumin+two*(a-l))*phiu02
        phiu1 = -(l/upl+two*(a-l))*phiu11+(l*iumin+two*(a-l))*phiu12
        dphiu = zcu*ur*(phiu0-(-one)**rtp*phiu1)/dsqrt(dabs(ee))
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*(ctuw1-(-one)**rtp*tu11)/dsqrt(dabs(ee))
      elseif (up>u3) then
!----------------------------------------------------------------------------------!
! Przypadek gdy mamy cztery rozne pierwiastki rzeczywiste
! Przypadek 8 Patrz tablica (uzupelnic numer tablicy potem) po wiecej informacji
!----------------------------------------------------------------------------------!
        if (elcaltp.and.(up/=u3)) then
          p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
          phiu01 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/upl,rfdcup,u3,up)
          phiu02 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-iumin,rfdcup,u3,up)
          ctuw2 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,u3,up)
          ctuw3 = phiu01
          ctuw4 = phiu02
          p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          ctuw1 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcup,u3,up)
        elseif (up==u3) then
          ctuw1 = zero
          ctuw2 = zero
          ctuw3 = zero
          ctuw4 = zero
          phiu01 = zero
          phiu02 = zero
        else
          phiu01 = ctuw3
          phiu02 = ctuw4
        endif
        if (uk/=u3) then
          p(5) = -2
! Pierwsze trzy calki z wzoru (47) Praca Dextera
          phiu11 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/upl,rfdcuk,u3,uk)
          phiu12 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-iumin,rfdcuk,u3,uk)
          tu12 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcuk,u3,uk)
          tu13 = phiu11
          tu14 = phiu12
          p(5) = -4
! Czwarta calka z wzoru (47) Praca Dextera
          tu11 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,one,rfdcuk,u3,uk)
        else
          tu11 = zero
          tu12 = zero
          tu13 = zero
          tu14 = zero
          phiu11 = zero
          phiu12 = zero
        endif
! Zapisuje wartosc z rownania (47) Praca Dextera jako sume tych calek
! Najpierw zapisze dla czesci zwiazanej z up
        tu0 = (iumin-one/upl)*ctuw1+(iumin**two-one/upl**two)*ctuw2-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &ctuw3+(two*a*(a-l)+a**two*iumin+iumin**try)*ctuw4
! A nastepnie zapisze dla czesci zwiazanej z uk
        tu1 = (iumin-one/upl)*tu11+(iumin**two-one/upl**two)*tu12-(two*a*(a-l)+a**two/upl+one/upl**try)*&
       &tu13+(two*a*(a-l)+a**two*iumin+iumin**try)*tu14
! Lacze powyzsze wartosci
        dtu = zcu*ur*(tu0-(-one)**rtp*tu1)/dsqrt(dabs(ee))
! Analogicznie robie z rownaniem (48) Praca Dextera
        phiu0 = -(l/upl+two*(a-l))*phiu01+(l*iumin+two*(a-l))*phiu02
        phiu1 = -(l/upl+two*(a-l))*phiu11+(l*iumin+two*(a-l))*phiu12
        dphiu = zcu*ur*(phiu0-(-one)**rtp*phiu1)/dsqrt(dabs(ee))
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
        lambdau = zcu*(ctuw1-(-one)**rtp*tu11)/dsqrt(dabs(ee))
      endif
    endif
! OPISIK TUTAJ
    if (np>4) then
! Obliczam pierwiastki dwukwadratowe M(mu)     
      yy = -hf*(a**two-ql2+dsign(one,a**two-ql2)*dsqrt((a**two-ql2)**two+4.0*q2*a**two))
      if ((a**two-ql2)<zero) then
        mne = -yy/a**two
        mpo = q2/yy
      else
        mne = q2/yy
        mpo = -yy/a**two
      endif
! Zabezpieczam obliczenia przed zakragleniem bledow w mpo
      if (mpo>one) then
        mpo = one
      endif
      muplus = dsqrt(mpo)
! Teraz uzyje troche innej wzoru na rozdzielenie calku mu, wszystkie calki beda policzone wzgledem
! muplus, ale procedury rozkladu i obliczania wspolczynnikow jest taka sama
      a1 = zcmu
      a2 = zcmu*(-one)**(mutp+1)
      a3 = two*int((two*mutp-zcmu+1)/4.0)
      if (mne<zero) then
! Zabezpieczam sie przed zaokraglaniem bledow przy liczeniu pierwiastkow
        if (muplus<mup) then
          muplus = mup
          mpo = muplus**two
        endif
!        print*,"rffmu1:",rfdcmup,"rffmu2:",rfdcmuk,"rffmu3:",rfdcmutp
!        print*,"muplus:",muplus,"mpo:",mpo
!        pause
! Tutaj rozpatrzy przypadek z symetrycznymi pierwiastkami, orbita moze przecinak plaszczyzne dyskowa
! Obliczam teraz skladowe calek t i phi _mu zgodnie z rownaniami 45 i 46 z pracy Dextera
        call calkiphiitsymetryczne(mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,&
       &ctmukplus,ctmutp,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,elcaltp)
! Zgodnie z rownanie 45
        dtmu = a**two*mne*cu+(a1*ctmupplus+a2*ctmukplus+a3*ctmutp)
! Zgodnie z rownaniem 46
        dphimu = -l*cu+l*(a1*cphimupplus+a2*cphimukplus+a3*cphimutp)
      else
! Tutaz rozpatrze niesymetryczne pierwiastki
        if (dsign(one,mup)==-1) then
          muplus = -muplus
        endif
! Zabezpieczam sie przed zaokraglaniem bledow w przypadku gdy mup jest turning pointem
        if (dabs(muplus)<dabs(mup)) then
          muplus = mup
          mpo = muplus**two
        endif
        if (dabs(mne)>mup**two) then
          mne = mup**two
        endif
! Obliczam teraz skladowe calek t i phi _mu zgodnie z rownaniami 45 i 46 z pracy Dextera
        call calkiphiitniesymetryczne(mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,&
       &ctmukplus,ctmutp,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,elcaltp)
! Zgodnie z rownanie 45
        dtmu = (a1*ctmupplus+a2*ctmukplus+a3*ctmutp)
! Zgodnie z rownaniem 46
        dphimu = -l*cu+l*(a1*cphimupplus+a2*cphimukplus+a3*cphimutp)
      endif
    endif
!    print*,"Phiu:",dphiu,"phimu:",dphimu,"phiu+phimu:",dphiu+dphimu
! I teraz obliczam wartosc lambda wzor (49) praca Dextera
    lambda = lambdau+dtmu
    if (l==zero) then
      dphiu = dphiu-dsign(one,a)*pi*mutp
    endif
    return
    end subroutine

  end module obliczphite
