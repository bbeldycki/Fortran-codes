  module obliczuk
!-------------------------------------------
! caluf.f95
! module caluf
!
!-------------------------------------------

  implicit none
  
  public   ::  geodezyjnauk
  contains   

    subroutine geodezyjnauk(up,uk,mup,muk,l,l2,q2,cmu,mutp,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,rfdcmup,&
   &rfdcmuk,rfdcmut,cu0,cmu1,cmu3,pehate,elcaltp)
    use const
    use funkcjeCarlsona
    use pierwiastki
    use funkcjeeliptycznejacobiego
    use calki
    implicit none
! Cel procedury:
! Procedura oblicza koncowa wartosc wspolrzednej radialnej odpowiadajaca koncowej wartosci wspolrzednej katawej
! Zadane mamy up, mup i muk oraz stale ruchu

! Parametry wejsciowe procedury:
! up - poczatkowa wartosc u. Jezeli up=0, to wartosci tki i lambdaki beda nieskonczone.
! mup - poczatkowa wartosc mu
! muk - koncowa wartosc mu
! a - spin czarnej dziury
! l - bezwymiarowa wartosc momentu pedu wzdluz osi z
! q2 - bezwymiarowa wartosc starej Cartera
! mutp - liczba mu turning pointow, policzona w przypadku gdy czrad=.FALSE., musi byc podana w przeciwnym przypadku
! zcu - znak calki u (du/dlmabda), = 1/-1 dla przychodzacych/odchodzacych promieni
! zcmu - znac calki mu (dmu/dlambda), moze byc policzony jesli znamy mup, beta i wiemy, ze wartosc 1/up to duza liczba

! Parametry wyjsciowe procedury:
! uk - koncowa wartosc u. Jezeli czrad=.TRUE. to policzona, w przeciwnym przypadku musi byc zadana

!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  up,uk,mup,muk,l,l2,q2,zcu,zcmu,cmu,h1,u1,u2,u3,u4,rfdcup,rfdcuk
    double precision  ::  rfdcmup,rfdcmuk,rfdcmut,cu0,cmu1,cmu2,cmu3
    integer  ::  mutp,rtp,np
    logical  ::  elcaltp,pehate
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  aa,a1,a2,a3,c1,c2,c3,c4,c5,sn,cn,dn,cc,dd,ee
    double precision  ::  f,g,f1,f2,g1,g2,h2,ql2,pi2,cu,mne,mpo,muplus!,muminus
    double precision  ::  upl,yy,dis,rr,dummy,qs,m1,s1,fi,bb,qq,n2,pp,r,ua,ub,mn2,pr2,chw,iut,iu1
    double precision  ::  jarg,sarg,np56,mp56,czz,czr
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
    if (q2==zero) then
! Obliczam wartosc calki mu (przypadej specjalny), gdy q2 jest rowne zero. Wowczas mamy 0 lub co najwyzej 1 mu turning point
      if (l2>=a**two .or. mup==zero) then
        cmu = zero
        uk = -one
        np = 0
        return
      else
        s1 = dsign(one,mup)
        a1 = s1*zcmu
        a2 = s1*zcmu*(-one)**(mutp+1)
        muplus = s1*dsqrt(one-l2/a**two)
        cmu1 = InvHipSec(mup/muplus)/dabs(a*muplus)
        cmu2 = InvHipSec(muk/muplus)/dabs(a*muplus)
        cmu = a1*cmu1+a2*cmu2
      endif
    elseif (a==zero) then
! Obliczam wartosc calki mu (przypadej specjalny), gdy a jest rowne zero
      a1 = zcmu
      muplus = dsqrt(q2/ql2)
      if (mup>muplus) then
        muplus = mup
      endif
      cmu1 = (hf*pi-dasin(mup/muplus))/dsqrt(ql2)
      cmu2 = (hf*pi+dasin(muk/muplus))/dsqrt(ql2)
      cmu3 = pi/dsqrt(ql2)
      a2 = zcmu*(-one)**mutp
      a3 = two*int((two*dble(mutp)+try-zcmu)/4.0)-one
      cmu = a1*cmu1+a2*cmu2+a3*cmu3
    endif
! Po szczegolnych przypadkach teraz bede rozpatrywac rownanie U(u)
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
! Przeskok, zeby kontrolowac wartosci up
151     continue
        if (up<=u2) then
          np = 1
          if (elcaltp.and.(up/=u2)) then
            cu0 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,zero,rfdcup,up,u2)
          elseif (up==u2) then
            cu0 = zero
          endif
! Tablica rozwiazan r patrz wiersz 1
          m1 = (u3-u2)/(u3-u1)
          jarg = dsqrt((u3-u1)*dd*hf*(cmu-zcu*cu0/dsqrt(dd)))
          call sncndn(jarg,m1,sn,cn,dn)
          uk = u1+(u2-u1)*cn**two/dn**two
          rtp = int((dsign(one,zcu*(cmu-zcu*cu0/dsqrt(dd)))+one)/two)
          if (pehate) then
            iu1 = trzyrzeczywiste(p,-u1,one,u2,-one,u3,-one,zero,zero,rfdcuk,uk,u2)/dsqrt(dd)
          endif
        elseif (up>=u3) then
          np = 2
          if (elcaltp.and.(up/=u3)) then
            cu0 = trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,zero,rfdcup,u3,up)
          elseif (up/=u3) then
            cu0 = zero
          endif
! Tablica rozwiazan r patrz wiersz 2
          m1 = (u3-u2)/(u3-u1)
          jarg = dsqrt((u3-u1)*dd*hf*(cmu+zcu*cu0/dsqrt(dd)))
          call sncndn(jarg,m1,sn,cn,dn)
          dummy = uk
          uk = u1+(u3-u1)*dn**two/cn**two
          rtp = int((-dsign(one,zcu*(cmu+cu0/dsqrt(dd)))+one)/two)
          if (pehate) then
            iu1 = trzyrzeczywiste(p,-u1,one,-u2,one,-u3,one,zero,zero,rfdcuk,u3,uk)/dsqrt(dd)
          endif
        else
          print*," Niefizyczny przypadek, modyfikuje dane"
          if (zcu==one) then
            up = u3
          else
            up = u2
          endif
          goto 151
        endif
      elseif (dabs(dis)<1d-16) then
        rtp = 0
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
            sarg = dsqrt((up-u1)/(u2-u1))
            jarg = dsqrt((u2-u1)*dd)*hf*zcu*cmu+hf*dlog((one+sarg)/(one-sarg))
            uk = u1+(u2-u1)*dtanh(jarg)**two
          endif
        elseif (up>=u2) then
          np = 2
          if (uk<u2) then
            uk = u2
            cu = zcu*1.d300
          else
            sarg = dsqrt((u2-u1)/(up-u1))
            jarg = -dsqrt((u2-u1)*dd)*hf*zcu*cmu+hf*dlog((one+sarg)/(one-sarg))
            uk = u1+(u2-u1)/dtanh(jarg)**two
          endif
        endif
      else
! Tutaj mamy przypadek z tylko jednym rzeczywistym pierwiastkiem i dwoma zespolonymi
        np = 3
! Tablica rozwiazan r patrz wiersz 3
        rtp = 0
        aa = -dsign(one,rr)*(dabs(rr)+dsqrt(dis))**trd
        if (aa/=zero) then
          bb = qq/aa
        else
          bb = zero
        endif
        u1 = (aa+bb)-cc/(try*dd)
        f = -one/(dd*u1)
        g = f/u1
! Upewniam sie, ze istnieje prawidlowe rozwiazania
        if (zcu>zero) then
          iut = dwazjedrzecz(p,-u1,one,zero,zero,f,g,one,dummy,up,upl)/dsqrt(dd)
        else
          iut = dwazjedrzecz(p,-u1,one,zero,zero,f,g,one,dummy,zero,up)/dsqrt(dd)
        endif
        if (cmu>iut) then
          uk = -one
          np = 0
          rtp = 0
          return
        endif
        np56 = -g/two
        if (elcaltp) then
          cu0 = zcu*dwazjedrzecz(p,-u1,one,zero,zero,f,g,one,dummy,u1,up)/dsqrt(dd)
        endif
        c3 = -one
        if (a/=zero.or.l/=zero) then
          c3 = (a+l)/(a-l)
        endif
        c2 = dsqrt(u1*(try*u1+c3))
        c1 = dsqrt(c2*dd)
        m1 = hf+(6.0*u1+c3)/(8.0*c2)
        jarg = c1*(cmu+cu0)
        call sncndn(jarg,m1,sn,cn,dn)
        uk = (c2+u1-(c2-u1)*cn)/(one+cn)
        if (pehate) then
          cu = dwazjedrzecz(p,-u1,one,zero,zero,f,g,one,rfdcup,up,uk)/dsqrt(dd)
        endif
      endif
    elseif (ee==zero .and. dd==zero) then
! Znowu mamy specjalny przypadek gdzue q2=0 i l=a i mup = 0
! W tym przypadku nie da sie obrocic calek
      np = 4
      uk = -one
      rtp = 0
    else
! 2) U(u) bedzie rownaniem czwartego stopnia 
! Ale najpierw znajde pierwiastki Mu(mu) w przypadku gdy sa dwukwadratowe
      yy = -hf*(a**two-ql2+dsign(one,a**two-ql2)*dsqrt((a**two-ql2)**two+4.0*q2*a**two))
      if ((a**two-ql2)<zero) then
        mne = -yy/a**two
        mpo = q2/yy
      else
        mne = q2/yy
        mpo = -yy/a**two
      endif
      muplus = dsqrt(mpo)
! Ponizsze wzory na rozdzial calki mu beda sie troche roznic. Wszystkie calki sa liczone wzgledem muplus
      a1 = zcmu
      a2 = zcmu*(-one)**(mutp+1)
      a3 = two*int((two*mutp-zcmu+one)/4.0)
      if (mne<zero) then
! Teraz bedzie przypadek rozwiazan symetrycznych. W tym przypadku orbita moze przecinac plaszczyzne dysku
        if (mup>muplus) then
          muplus = mup
          mpo = muplus**two
        endif
! Obliczam calki ktore uwzgledniaja up i turning point
! Robie to tylko raz na kazda geodezyjna
        if (elcaltp) then
          call calkaimusymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfdcmup,rfdcmut)
        endif
        call calkaimusymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfdcmuk)
        cmu = a1*cmu1+a2*cmu2+a3*cmu3
      else
! Teraz bedzie przypadek rozwiazan niesymetrycznych
        if (dabs(muk)<dsqrt(mne)) then
          uk = -one
          np = 0
          rtp = 0
          return
        else
          if (dsign(one,mup)==-1) then
            muplus = -muplus
          endif
          if (dabs(muplus)<dabs(mup)) then
            muplus = mup
            mpo = muplus**two
          endif
          mne = min(mup**two,mne)
          if (elcaltp) then
            call calkaimuniesymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfdcmup,rfdcmut)
          endif
          call calkaimuniesymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfdcmuk)
        endif
        cmu = a1*cmu1+a2*cmu2+a3*cmu3
      endif
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
! Upewniam sie, ze policzone pierwiastki nie maja jakis glupich czesci urojonych w stylu 1.0d-24
        do i=1,4
          czz = dimag(pierwrzeczywiste(i))
          czr = dreal(pierwrzeczywiste(i))
          if (dabs(czz)<1.0d-3) then
              czz = zero
          endif
          pierwrzeczywiste(i) = dcmplx(czr,czz)
        enddo
        do i=1,4
          if (dimag(pierwrzeczywiste(i))==zero) then
            nrzecz = nrzecz+1
          endif
        enddo
        if (nrzecz==2) then
! Przypadek gdzie mamy dwa pierwiastki rzeczywiste i dwa zespolone
          np = 5
          rtp = 0
          u1 = dble(pierwrzeczywiste(1))
          if (dimag(pierwrzeczywiste(2))==zero) then
            u4 = dble(pierwrzeczywiste(2))
          else
            u4 = dble(pierwrzeczywiste(4))
          endif
        elseif (nrzecz==0) then
          np = 6
          rtp = 0
        else
          u1 = dble(pierwrzeczywiste(1))
          u2 = dble(pierwrzeczywiste(2))
          u3 = dble(pierwrzeczywiste(3))
          u4 = dble(pierwrzeczywiste(4))
152       continue
          if (u2>upl .and. u3>upl) then
            np = 5
          elseif (up<=u2) then
            np = 7
          elseif (up>u3) then
            np = 8
          else
            print*,"Niefizyczny przypadek, modyfikuje dane"
            if (zcu==1) then
              up = u3
            else 
              up = u2
            endif
            goto 152
          endif
        endif
      endif ! Koniec warunku (elcaltp)
      if (np==5) then
! Tablica rozwiazan r patrz wiersz 5
        qs = dsign(one,q2)
        f= -qs/(u1*u4*dabs(ee))
        g = (u1+u4)*f/(u1*u4)
! Upewniam sie, ze istnieje prawidlowe rozwiazania
        if (zcu>zero) then
          iut = dwazdwarzecz(p,-u1,one,qs*u4,-qs*one,zero,zero,f,g,one,dummy,up,upl)/dsqrt(ee)
        else
          iut = dwazdwarzecz(p,-u1,one,qs*u4,-qs*one,zero,zero,f,g,one,dummy,zero,up)/dsqrt(ee)
        endif
        if (cmu>iut) then
          uk = -one
          np = 0
          rtp = 0
          return
        endif
        if (qs==one) then
          ua = u4
          ub = u1
        else
          ua = u1
          ub = u4
        endif
        mp56 = -g*hf
        n2 = f-g**two/4.0
        c4 = dsqrt((mp56-u4)**two+n2)
        c5 = dsqrt((mp56-u1)**two+n2)
        c1 = dsqrt(dabs(ee*c4*c5))
        if (elcaltp) then
          cu0 = zcu*dwazdwarzecz(p,-u1,one,qs*u4,-qs*u4,zero,zero,f,g,one,dummy,utp,up)/dsqrt(dabs(ee))
        endif
        m1 = qs*((c4+qs*c5)**two-(u4-u1)**two)/(4.0*c4*c5)
        jarg = c1*(cmu+cu0)
        call sncndn(jarg,m1,sn,cn,dn)
        uk = (u4*c5+qs*u1*c4-(qs*u4*c5-u1*c4)*cn)/((c4-qs*c5)*cn+(qs*c4+c5))
        if (pehate) then
          cu = dwazdwarzecz(p,-u1,one,qs*u4,-qs*one,zero,zero,f,g,one,rfdcup,up,uk)/dsqrt(dabs(ee))
        endif 
! Koniec warunku (np==5)
      elseif (np==6) then
! Teraz bedzie przypadek gdzie nie mamy wcale pierwiastkow rzeczywistych
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
153     continue
        i = i+1
        if (dimag(pierwiastkih(i))==zero) then
          h1 = dble(pierwiastkih(i))
        endif
        if (h1==zero) then
          goto 153
        endif
        h2 = one/h1
        g1 = dd/((h2-h1)*ee)
        f1 = one/dsqrt(ee)
        f2 = f1
! Upewniam sie, ze istnieje prawidlowe rozwiazania 
        if (zcu>zero) then
          iut = czteryzespolone(p,f1,g2,h1,f2,g2,h2,zero,zero,dummy,up,upl)/dsqrt(dabs(ee))
        else
          iut = czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,zero,dummy,zero,up)/dsqrt(dabs(ee))
        endif
        if (cmu>iut) then
          uk = -one
          rtp = 0
          return
        endif
! Teraz chce zapisac czesc rzeczywista i urojona pierwiastkow m+/- n i p+/-r za pomoca wielkosci rzeczywistych
        wspolpzesp(1) = ee**(-try)
        wspolpzesp(2) = -cc/ee**try
        wspolpzesp(3) = -ee**(-two)
        wspolpzesp(4) = -ee**(-two)*(dd**two/ee-two*cc)
        wspolpzesp(5) = -one/ee
        wspolpzesp(6) = -cc/ee
        wspolpzesp(7) = one
        call zroots(wspolpzesp,6,pierwiastkih,.TRUE.)
        i = 0
        mn2 = zero
154     continue
        if (dimag(pierwiastkih(i))==zero) then
          mn2 = dble(pierwiastkih(i))
        endif
        if (mn2==zero) then
          goto 154
        endif
        pp = dd/(two*ee**two*(mn2**two-one/ee))
        mp56 = -pp-hf*dd/ee
        if (mp56<pp) then
          chw = pp
          pp = mp56
          mp56 = chw
        endif
        pr2 = one/(mn2*ee)
        np56 = dsqrt(mn2-mp56**two)
        r = dsqrt(pr2-pp**two)
        c4 = dsqrt((mp56-pp)**two+(np56+r)**two)
        c5 = dsqrt((mp56-pp)**two+(np56-r)**two)
        c1 = (c4+c5)*hf*dsqrt(dabs(ee))
        c2 = dsqrt((4.0*np56**two-(c4-c5)**two)/((c4+c5)**two-4.0*np56**two))
        c3 = mp56+c2*np56
        if (elcaltp) then
          cu0 = dsign(one,(up-c3))*czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,zero,dummy,c3,up)/dsqrt(dabs(ee))
        endif
        m1 = ((c4-c5)/(c4+c5))**two
        jarg = c1*(zcu*cmu+cu0)
        call sncndn(jarg,m1,sn,cn,dn)
        uk = c3+(np56*(one+c2**two)*sn/cn)/(one-c2*sn/cn)
        if (pehate) then
          cu = czteryzespolone(p,f1,g1,h1,f2,g2,h2,zero,zero,rfdcmup,up,uk)/dsqrt(-ee)
        endif
! Koniec warunku (np==6)
      else
! Teraz beda przypadki gdzie sa cztery pierwiastki rzeczywiste
        if (np==7) then
! Tablica rozwiazan r patrz wiersz 7
          if (elcaltp.and.(up/=u2)) then
            cu0 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,zero,rfdcup,up,u2)
          endif
          jarg = dsqrt(dabs(ee)*(u3-u1)*(u4-u2))*hf*(cmu-zcu*cu0/dsqrt(-ee))
          m1 = (u4-u1)*(u3-u2)/((u4-u2)*(u3-u1))
          call sncndn(jarg,m1,sn,cn,dn)
          dummy = uk
          uk = ((u2-u1)*u3*sn**two-u2*(u3-u1))/((u2-u1)*sn**two-(u3-u1))
          rtp = int((dsign(one,zcu*(cmu-zcu*cu0/dsqrt(dabs(ee))))+one)/two)
          if (pehate) then
            iu1 = czteryrzeczywiste(p,-u1,one,u2,-one,u3,-one,u4,-one,zero,zero,rfdcuk,uk,u2)/dsqrt(-ee)
          endif
        elseif (np==8) then
! Tablica rozwiazan r patrz wiersz 8
          if (elcaltp.and.(up/=u3)) then
            cu0 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,zero,rfdcup,u3,up)
          endif
          jarg = dsqrt(dabs(ee)*ee*(u3-u1)*(u4-u2))*hf*(cmu+zcu*cu0/dsqrt(-ee))
          m1 = (u4-u1)*(u3-u2)/((u4-u2)*(u3-u1))
          call sncndn(jarg,m1,sn,cn,dn)
          uk = ((u4-u3)*u2*sn**two-u3*(u4-u2))/((u4-u3)*sn**two-(u4-u2))
          rtp = int((-dsign(one,zcu*(cmu+zcu*(cu0/dsqrt(-ee))))+one)/two)
          if (pehate) then
            iu1 = czteryrzeczywiste(p,-u1,one,-u2,one,-u3,one,u4,-one,zero,zero,rfdcuk,u3,uk)/dsqrt(-ee)
          endif
        endif
! Koniec warunku na (np==7 lub np==8)
      endif
    endif
    return
    end subroutine

  end module obliczuk
