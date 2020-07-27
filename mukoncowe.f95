  module obliczmuk
!-------------------------------------------
! calmuf.f95
! module calmuf
!
!-------------------------------------------

  implicit none
  
  public   ::  geodezyjnamuk
  contains   

    subroutine geodezyjnamuk(up,uk,mup,muk,l,l2,q2,cu,mutp,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
   &rfdcmup,rfdcmuk,rfdcmut,cu0,cmu1,cmu3,pehate,elcaltp)
    use const
    use funkcjeCarlsona
    use pierwiastki
    use funkcjeeliptycznejacobiego
    use calki
    implicit none
! Cel procedury:
! Procedura oblicza koncowa wartosc wspolrzednej katowej odpowiadajaca koncowej wartosci wspolrzednej radialnej
! Zadane mamy up, uk i mup oraz stale ruchu

! Parametry wejsciowe procedury:
! up - poczatkowa wartosc u. Jezeli up=0, to wartosci tki i lambdaki beda nieskonczone
! uk - koncowa wartosc u
! mup - poczatkowa wartosc mu
! a - spin czarnej dziury
! l - bezwymiarowa wartosc momentu pedu wzdluz osi z
! q2 - bezwymiarowa wartosc starej Cartera
! mutp - liczba mu turning pointow, policzona w przypadku gdy czrad=.FALSE., musi byc podana w przeciwnym przypadku
! zcu - znak calki u (du/dlmabda), = 1/-1 dla przychodzacych/odchodzacych promieni
! zcmu - znac calki mu (dmu/dlambda), moze byc policzony jesli znamy mup, beta i wiemy, ze wartosc 1/up to duza liczba

! Parametry wyjsciowe procedury:
! muk - koncowa wartosc mu
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  up,uk,mup,muk,l,l2,q2,zcu,zcmu,cu,h1,u1,u2,u3,u4,rfdcup,rfdcuk
    double precision  ::  rfdcmup,rfdcmuk,rfdcmut,cu0,cmu1,cmu2,cmu3
    integer  ::  mutp,rtp,np
    logical  ::  elcaltp,pehate
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  aa,a1,a2,a3,sn,cn,dn,cc,dd,ee
    double precision  ::  f,g,f1,f2,g1,g2,h2,ql2,pi2,mne,mpo,muplus,muminus
    double precision  ::  upl,yy,dis,rr,qs,m1,s1,fi,bb,qq,iu1
    double precision  ::  sarg,iarg,uarg,farg,czr,czz
    integer  ::  i,nrzecz
    double complex, dimension(4)  ::  pierwrzeczywiste
    double complex, dimension(5)  ::  wspolprzecz
    double complex, dimension(6)  ::  pierwiastkih
    double complex, dimension(7)  ::  wspolpzesp
    integer, dimension(5)  ::  p
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    pi2 = two*pi
    upl = one/(one+dsqrt(one-a**two))
    ql2 = q2+l2
! Obliczam wspolrzynniki rownania U(u)
    cc = a**two-ql2
    dd = two*((a-l)**two+q2)
    ee = -a**two*q2
! Rozpatruje rownanie U(u)
! Istnieja dwa przypadki:
! 1) U(u) bedzie rownaniem trzciego stopnia 
! Szukam pierwiastkow rownania U(u) = 0
    if ((ee==zero) .and. (dd/=zero)) then
      p(1) = -1
      p(2) = -1
      p(3) = -1
      p(4) = 0
      p(5) = 0
      qq = cc**two/(9.0*dd**two)
      rr = (two*cc**try/dd**try+27.0/dd)/54.0
      dis = rr**two-qq**try
      if (dis<-1.d-16) then
! Tutaj mam przypadek trzech rzeczywistych rozwiazan, gdzie zachodzi u1<0<u2<=u3
! Pierwiastki sa zapisane zawsze w kolejnosci rosnacej
        fi = dacos(rr/qq**ohf)
        u1 = -two*dsqrt(qq)*dcos(fi/try)-cc/(try*dd)
        u2 = -two*dsqrt(qq)*dcos((fi-two*pi)/try)-cc/(try*dd)
        u3 = -two*dsqrt(qq)*dcos((fi+two*pi)/try)-cc/(try*dd) 
        iu1 = zero
! Przeskok, zeby kontrolowac wartosci up
150     continue
        if (up<=u2) then
          if (uk>u2) then
            uk = u2
          endif
          np = 1
! Tablica rozwiazan r patrz wiersz 1
          if (elcaltp.and.(up/=u2)) then
            cu0 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,zero,rfdcup,up,u2)
          elseif (up==u2) then
            cu0 = zero
          endif
          if (uk/=u2) then
            iu1 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,zero,rfdcuk,uk,u2)
          endif
          cu = zcu*(cu0-(-one)**rtp*iu1)/dsqrt(dd)
        elseif (up>=u3) then
          if (uk<u3) then
            uk = u3
          endif
          np = 2
! Tablica rozwiazan r patrz wiersz 2
          if (elcaltp.and.(up/=u3)) then
            cu0 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,zero,rfdcup,u3,up)
          elseif (up==u3) then
            cu0 = zero
          endif
          if (uk/=u3) then
            iu1 = -trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,zero,rfdcuk,u3,uk)
          endif
          cu = zcu*(cu0-(-one)**rtp*iu1)/dsqrt(dd)          
        else
          print*," Niefizyczny przypadek, modyfikuje dane"
          if (zcu==one) then
            up = u3
          else
            up = u2
          endif
          goto 150
        endif
      elseif (dabs(dis)<1.d-16) then
! Tutaj mamy przypadek trzech rozwiazan z dwoma jednakowymi pierwiastkami
        u1 = -two*dsqrt(qq)-cc/(try*dd)
        u2 = -two*dsqrt(qq)*dcos(two*pi/try)-cc/(try*dd)
        u3 = u2
        rtp = 0
        if (up<=u2) then
          np = 1
          if (uk>u2) then
            uk = u2
            cu = zcu*1.d300
          else
            farg = (dsqrt(uk-u1)+dsqrt(u2-u1))/dabs(dsqrt(uk-u1)-dsqrt(u2-u1))
            sarg = (dsqrt(up-u1)+dsqrt(u2-u1))/dabs(dsqrt(up-u1)-dsqrt(u2-u1))
            cu = zcu*(dlog(farg)-dlog(sarg))/dsqrt((u2-u1)*dd)
          endif
        elseif (up>=u2) then
          np = 2
          if (uk<u2) then
            uk = u2
            cu = zcu*1.d300
          else
            farg = (dsqrt(uk-u1)+dsqrt(u2-u1))/dabs(dsqrt(uk-u1)-dsqrt(u2-u1))
            sarg = (dsqrt(up-u1)+dsqrt(u2-u1))/dabs(dsqrt(up-u1)-dsqrt(u2-u1))
            cu = zcu*(dlog(farg)-dlog(sarg))/dsqrt((u2-u1)*dd)
          endif
        endif
      else
! Tutaj mamy przypadek z tylko jednym rzeczywistym pierwiastkiem i dwoma zespolonymi
        np = 3
        aa = -dsign(one,rr)*(dabs(rr)+dsqrt(dis))**trd
        if (aa/=zero) then
          bb = qq/aa
        else
          bb = zero
        endif
! Tablica rozwiazan r patrz wiersz 3
        u1 = (aa+bb)-cc/(try*dd)
        f = -one/(dd*u1)
        g = f/u1
        if (uk>up) then
          cu = zcu*dwazjedrzecz(p,-u1,one,zero,zero,f,g,one,rfdcup,up,uk)/dsqrt(dd)
        elseif (uk<up) then
          cu = -zcu*dwazjedrzecz(p,-u1,one,zero,zero,g,f,one,rfdcup,uk,up)/dsqrt(dd)
        else
          cu = zero
        endif
      endif
      if (q2==zero) then
! Obliczam muk dla przypadku szczegolnego (q2=0)
! W tym przypadku moze byc co najwyzej jeden turning point
        s1 = dsign(one,mup)
        if (l2<a**two .and. mup/=zero) then
          muplus = s1*dsqrt(one-l2/a**two)
! Rozwiazanie na muk wzor (43) z pracy Dextera
          muk = muplus/dcosh(dabs(a*muplus)*cu-s1*zcmu*InvHipSec(mup/muplus))
        else
          muk = zero
        endif
      elseif (a==zero) then
! Obliczam muk dla przypadku szczegolnego (a=0)
        a1 = zcmu
        muplus = dsqrt(q2/ql2)
        if (mup>muplus) then
          muplus = mup
        endif
        cmu1 = dacos(mup/muplus)/dsqrt(ql2)
        cmu3 = pi/dsqrt(ql2)
        if (zcmu==one) then
! Wzor (30) z pracy Dextera
          mutp = int((cu-cmu1)/cmu3)+int((one+dsign(one,(cu-cmu1)/cmu3))/two)
        else
! Wzor (33) z pracy Dextera
          mutp = int((cu+cmu1)/cmu3)
        endif
        a2 = zcmu*(-one)**mutp
! Wzor (35) z pracy Dextera
        a3 = two*int((two*dble(mutp)+try-zcmu)/4.0)-one
! Rozwiazanie na muk wzor (44) z pracy Dextera
! Korzystam z wslanosci muminus = -muplus
        muk = -muplus*dcos(dsqrt(ql2)*(cu-a1*cmu1-a3*cmu3)/a2)          
      endif
    elseif ((ee==zero) .and. (dd==zero)) then
! Znowu mamy specjalny przypadek gdzue q2=0 i l=a, wowczas U(u) = 1
      np = 4
      rtp = 0
      cu = zcu*(uk-up)
      muk = zero
    else
!      print*,"wszedlem tutaj?"
!      pause
! 2) U(u) bedzie rownaniem czwartego stopnia 
! Szukam pierwiastkow rownania U(u) = 0
      p(1) = -1
      p(2) = -1
      p(3) = -1
      p(4) = -1
      p(5) = 0
      if (elcaltp) then
        wspolprzecz(1) = dcmplx(one,zero)
        wspolprzecz(2) = dcmplx(zero,zero)
        wspolprzecz(3) = dcmplx(cc,zero)
        wspolprzecz(4) = dcmplx(dd,zero)
        wspolprzecz(5) = dcmplx(ee,zero)
        call zroots(wspolprzecz,4,pierwrzeczywiste,.TRUE.)
        nrzecz = 0
!        print*,pierwrzeczywiste
!        pause
! Upewniam sie, ze policzone pierwiastki nie maja jakis glupich czesci urojonych w stylu 1.0d-24
        do i=1,4
          czz = dimag(pierwrzeczywiste(i))
          czr = dreal(pierwrzeczywiste(i))
          if (dabs(czz)<1.0d-8) then
              czz = zero
          endif
          pierwrzeczywiste(i) = dcmplx(czr,czz)
        enddo
!        print*,pierwrzeczywiste
!        pause
        do i=1,4
          if (dimag(pierwrzeczywiste(i))==zero) then
            nrzecz = nrzecz+1
          endif
        enddo
        if (nrzecz==2) then
! Przypadek gdzie mamy dwa pierwiastki rzeczywiste i dwa zespolone
          np = 5
          u1 = dble(pierwrzeczywiste(1))
          if (dimag(pierwrzeczywiste(2))==zero) then
            u4 = dble(pierwrzeczywiste(2))
          else
            u4 = dble(pierwrzeczywiste(4))
          endif
        elseif (nrzecz==0) then
          np = 6
        else
          u1 = dble(pierwrzeczywiste(1))
          u2 = dble(pierwrzeczywiste(2))
          u3 = dble(pierwrzeczywiste(3))
          u4 = dble(pierwrzeczywiste(4))
149       continue
! Istnieja przypadki gdzie wszystkie pierwiastki sa rzeczywiste, ale sa niefizyczne
          if (u2>upl .and. u3>upl) then
            np = 5
          elseif (up<=u2) then
            np = 7
          elseif (up>=u3) then
            np = 8
          else
! Jezeli zajdzie przypadek u2<up<u3 to musimy modyfikowac dane
            print*,"Niefizyczny przypadek, modyfikuje dane"
            if (zcu==one) then
              up = zero
            else
              up = upl
            endif
            goto 149
          endif
        endif
      endif
      if (np==5) then
! Przypadek gdzie mamy dwa pierwiastki rzeczywiste i dwa zespolone
! Tablica rozwiazan r patrz wiersz 5
        qs = dsign(one,q2)
        f = -qs/(dabs(ee)*u1*u4)
        g = (u1+u4)*f/(u1*u4)
        if (uk>up) then
          cu = zcu*dwazdwarzecz(p,-u1,one,qs*u4,-qs*one,zero,zero,f,g,one,rfdcup,up,uk)/dsqrt(dabs(ee))
        elseif (uk<up) then
          cu = -zcu*dwazdwarzecz(p,-u1,one,qs*u4,-qs*one,zero,zero,f,g,one,rfdcup,uk,up)/dsqrt(dabs(ee))
        else
          cu = zero
        endif
! Koniec warunku (np==5)
      elseif (np==6) then
! Przypadek gdzie mamy nie mamy pierwiastkow rzeczywistych
! Tablica rozwiazan r patrz wiersz 6
        np = 6
        wspolpzesp(1) = dcmplx(one,zero)
        wspolpzesp(2) = dcmplx(-cc/dsqrt(ee),zero)
        wspolpzesp(3) = dcmplx(-one,zero)
        wspolpzesp(4) = dcmplx(dsqrt(ee)*(two*cc/ee-dd**two/ee**two),zero)
        wspolpzesp(5) = dcmplx(-one,zero)
        wspolpzesp(6) = dcmplx(-cc/dsqrt(ee),zero)
        wspolpzesp(7) = dcmplx(one,zero)
        call zroots(wspolpzesp,6,pierwiastkih,.TRUE.)
        i = 0
        h1 = zero
148     continue
        i = i+1
        if (dimag(pierwiastkih(i))==zero) then
          h1 = dble(pierwiastkih(i))
        endif
        if (h1==zero) then
          goto 148
        endif
        h2 = one/h1
        g1 = dd/((h2-h1)*ee)
        g2 = -g1
        f1 = one/dsqrt(ee)
        f2 = f1
        if (uk>up) then
          cu = zcu*czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,zero,rfdcup,up,uk)/dsqrt(dabs(ee))
        elseif (uk<up) then
          cu = -zcu*czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,zero,rfdcup,uk,up)/dsqrt(dabs(ee))
        else
          cu = zero
        endif
! Koniec warunku (np==6)
      else
! Teraz beda przypadki gdzie sa cztery pierwiastki rzeczywiste
        if (dabs(u3-u2)>1.d-12) then
          iu1 = zero
          if (np==7) then
! Tablica rozwiazan r patrz wiersz 7
            if (uk>u2) then
              uk = u2
            endif
            if (elcaltp.and.(up/=u2)) then
              cu0 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,zero,rfdcup,up,u2)
            elseif (up==u2) then
              cu0 = zero
            endif
            if (uk/=u2) then
              iu1 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,zero,rfdcuk,uk,u2)
            endif
            cu = zcu*(cu0-(-one)**rtp*iu1)/dsqrt(dabs(ee))
!            print*,"cu:",cu,"cu0:",cu0,"zcu:",zcu,"iu1:",iu1
          elseif (np==8) then
! Tablica rozwiazan r patrz wiersz 7
            if (uk<u3) then
              uk = u3
            endif
            if (elcaltp.and.(up/=u3)) then
              cu0 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,zero,rfdcup,u3,up)
            elseif (up==u3) then
              cu0 = zero
            endif
            if (uk/=u3) then
              iu1 = -czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,zero,rfdcuk,u3,uk)
            endif
            cu = zcu*(cu0-(-one)**rtp*iu1)/dsqrt(dabs(ee))
          endif
        endif
! Koniec warunku na (np==7 lub np==8)
      endif 
!      print*,"cu:",cu,"cu0:",cu0,"iu1:",iu1
!      pause
! Obliczam pierwiastki Mu(mu) w przypadku gdy sa dwukwadratowe
      yy = -hf*(a**two-ql2+dsign(one,a**two-ql2)*dsqrt((a**two-ql2)**two+4.0*q2*a**two))
      if ((a**two-ql2)<zero) then
        mne = -yy/a**two
        mpo = q2/yy
      else
        mne = q2/yy
        mpo = -yy/a**two
      endif
      if (mpo>one) then
        mpo = one
      endif
      muplus = dsqrt(mpo)
      if (mne<zero) then
! Teraz bedzie przypadek rozwiazan symetrycznych. W tym przypadku orbita moze przecinac plaszczyzne dysku
        if (muplus<mup) then
          muplus = mup
          mpo = muplus**two
        endif
        muminus = -muplus
        if (elcaltp) then
          call calkaimusymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfdcmup,rfdcmut)
        endif
        a1 = zcmu
        if (zcmu==one) then
! Wzor (30) z pracy Dextera
          mutp = int((cu-cmu1)/cmu3)+int((one+dsign(one,(cu-cmu1)/cmu3))/two)
        else
! Wzor (33) z pracy Dextera
          mutp = int((cu+cmu1)/cmu3)
        endif
        a2 = zcmu*(-one)**mutp
! Wzor (35) z pracy Dextera
        a3 = two*int((two*dble(mutp)+try-zcmu)/4.0)-one
        iarg = dabs(a)*(cu-a1*cmu1-a3*cmu3)/a2
        uarg = dsqrt(mpo-mne)*iarg
        m1 = -mne/(mpo-mne)
        call sncndn(uarg,m1,sn,cn,dn)
        muk = muminus*cn
        if (pehate) then
          call calkaimusymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfdcmuk)
        endif
      else
! Teraz bedzie przypadek rozwiazan niesymetrycznych
        muplus = dsign(one,mup)*muplus
        muminus = dsign(one,mup)*dsqrt(mne)
        if (dabs(muminus)>dabs(mup)) then
          muminus = mup
          mne = muminus**two
        endif
        if (dabs(mup)>dabs(muplus)) then
          muplus = mup
          mpo = muplus**two
        endif
        if (elcaltp) then
          call calkaimuniesymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfdcmup,rfdcmut)
        endif
        a1 = zcmu
        if (zcmu==one) then
! Wzor (30) z pracy Dextera
          mutp = int((cu-cmu1)/cmu3)+int(dabs((one+dsign(one,(cu-cmu1)/cmu3))/two))
        else
! Wzor (33) z pracy Dextera
          mutp = int((cu+cmu1)/cmu3)
        endif
        a2 = zcmu*(-one)**mutp
! Wzor (35) z pracy Dextera
        a3 = two*int((two*dble(mutp)+try-zcmu)/4.0)-one
        iarg = dabs(a)*(cu-a1*cmu1-a3*cmu3)/a2
        uarg = dabs(muplus)*iarg
        m1 = mne/mpo
        call sncndn(uarg,m1,sn,cn,dn)
        muk = muminus/dn
        if (pehate) then
          call calkaimuniesymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfdcmuk)
        endif
      endif
    endif
!    print*,"muk:",muk,"cu:",cu,"cmu:",(a1*cmu1+a2*cmu2+a3*cmu3)
    return
    end subroutine
    
  end module obliczmuk

