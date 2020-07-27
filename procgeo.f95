  module procgeo
!-------------------------------------------
! geoproc.f95
! module geoproc
!
!-------------------------------------------

  implicit none
  
  public   ::  geodezyjna
  contains   
  
    subroutine geodezyjna(up,uk,upom,mup,muk,l,q2,alfa,beta,mutp,rtp,zcu,zcmu,lpt,ofset,phiit,czrad,zmzm,np,&
   &uki,muki,tki,phiki,mutpi,rtpi,lambdaki)
    use const
    use zmianazmiennej
    use pierwiastkimu
    use obliczphite
    use obliczuk
    use obliczmuk
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
! lpt - liczba puntkow na trajektorii 

! Parametry wyjsciowe procedury:
! uki(maxlet) - tablica z wartosciamy uk pomiedzy up i uk lub w przypadku rtp=1 pomiedzy up, utp i uk
! muki(maxlet) - tablica z wartosciamu muk odpowiadajaca wartosciom uki
! tki(maxlet) - tablica z wartosciami delta t, czyli t(uki(i))-t(up)
! phiki(maxlet) - tablica z wartosciami delta phi, czyli phi(uki(i))-phi(up) 
! mutpi(maxlet) - tablica z liczba mu turning pointow napotkanych na trajektorii 
! rtpi(maxlet) - tablica z liczba u turning pointow napotkanych na trajektorii 
! lambdaki(maxlet) - tablica z wartosciami parametru afinicznego wzdluz trajektorii lambda(uki(i))
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    logical  ::  phiit,czrad,zmzm,elcaltp
    integer  ::  lpt,np,mutp,rtp
    double precision  ::  up,uk,mup,muk,upom,l,q2,zcu,zcmu,ofset,alfa,beta
    double precision, dimension(maxlet)  ::  uki,muki,tki,phiki,lambdaki
    integer, dimension(maxlet)  ::  mutpi,rtpi
!-------------------------------zmienne lokalne------------------------------------!
    integer  ::  k,kmax,kex,mutp1,rtp1,lpptp
    double precision  ::  muminus,muplus,muv,uni,l2,upl,du,cu,u1,u2,u3,u4,h1,rfdcup,rfdcuk,rfdcmup
    double precision  ::  rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,lambda,rd2r,rj3r,ctmupplus,ctmutp
    double precision  ::  dtu,dtmu,dphiu,dphimu,ctuw1,ctuw2,ctuw3,ctuw4,cphimupplus,cphimutp
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    alfa = alfa
! Zadaje wartosci zmiennej logicznej oraz obliczam kwadram momentu pedu i zadaje wartosci intow    
    l2 = l**two
    lpptp = 0
    kex = 0
    elcaltp = .TRUE.
    upl = one/(one+dsqrt(one-a**two))
    !print*,"upom:",upom,"up:",up
! Przechodze teraz do warunku na tryb obliczen
    if (czrad) then
! Obliczam pierwiastki rownania M(mu)
      call rozwiazaniamu(q2,l,mup,muminus,muplus)
      k = 1
! Bede teraz rozwiazywac geodezyjna po uk w rownych krokach pomiedzy mup i muk
! Zakladam, ze nie bedzie turning pointow po drodze
! Zanim to robie jeszcze zmiane zmiennych
!      print*,"wartosc cu 1:",cu
      call zmiennaniezalezna(muminus,muplus,mup,muk,0,mutp,zcmu,k,lpt,ofset,muv,mutp1)
      if (mup/=zero .or. beta/=zero) then
! Wolam procedure geodezyjnauk
!        print*,"wartosc cu 2:",cu
        call geodezyjnauk(up,uk,mup,muv,l,l2,q2,cu,mutp1,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,rfdcmup,&
       &rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.TRUE.)
!        print*,"wartosc cu 3:",cu
        if (phiit) then
! Wolam procedure geodezyjnaphite
          call geodezyjnaphite(up,uk,mup,muv,l,l2,q2,mutp1,rtp,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,dphiu,&
         &dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,ctmupplus,&
         &ctmutp,cphimupplus,cphimutp,.TRUE.)
        endif
!        print*,"wartosc cu 4:",cu
        !print*,"up:",up,"uk:",uk,"mup:",mup,"muk:",muv,"k:",k,mutp1,rtp
!        pause
      else
        uk = up
      endif     
      uki(k) = uk
      muki(k) = muv
      tki(k) = dtu+dtmu
      phiki(k) = dphiu+dphimu
      mutpi(k) = mutp1
      rtpi(k) = rtp
      lambdaki(k) = lambda
      do k=2,lpt
! Obliczam pierwiastki rownania M(mu)
        call rozwiazaniamu(q2,l,mup,muminus,muplus)
! Bede teraz rozwiazywac geodezyjna po uk w rownych krokach pomiedzy mup i muk
! Zakladam, ze nie bedzie turning pointow po drodze
! Zanim to robie jeszcze zmiane zmiennych
!        print*,"wartosc cu 1:",cu
        call zmiennaniezalezna(muminus,muplus,mup,muk,0,mutp,zcmu,k,lpt,ofset,muv,mutp1)
        if (mup/=zero .or. beta/=zero) then
! Wolam procedure geodezyjnauk
!          print*,"wartosc cu 2:",cu
          call geodezyjnauk(up,uk,mup,muv,l,l2,q2,cu,mutp1,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
         &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.TRUE.)
!          print*,"wartosc cu 3:",cu
          if (phiit) then
! Wolam procedure geodezyjnaphite
            call geodezyjnaphite(up,uk,mup,muv,l,l2,q2,mutp1,rtp,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,&
           &dphiu,dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,&
           &ctmupplus,ctmutp,cphimupplus,cphimutp,.TRUE.)
          endif
!          print*,"wartosc cu 4:",cu
          !print*,"up:",up,"uk:",uk,"mup:",mup,"muk:",muv,"k:",k,mutp1,rtp
!          pause
        else
          uk = up
        endif     
        uki(k) = uk
        muki(k) = muv
        tki(k) = dtu+dtmu
        phiki(k) = dphiu+dphimu
        mutpi(k) = mutp1
        rtpi(k) = rtp
        lambdaki(k) = lambda
      enddo
    else
!      muk = zero
!      print*,"up:",up,"uk:",uk,"mup:",mup,"muk:",muk,"lpt:",lpt
! Wywolujemy procedura geodezyjnamuk, ktora znajdzie trajektorie, zeby wziac pod uwage odpowiedni przypadek i pierwiastki
      call geodezyjnamuk(up,uk,mup,muk,l,l2,q2,cu,mutp,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,rfdcmup,&
     &rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.TRUE.)
      !print*,"upom:",upom,"up:",up
!      print*,"up:",up,"uk:",uk,"mup:",mup,"muk:",muk,"lpt:",lpt
!      print*,"np:",np
      !pause
! Rozpatrujemy wszystkie fizyczne przypadki
      if (np>2 .and. np<7 .or. (np==1 .and. zcu<zero) .or. (np==7 .and. zcu<zero) .or. (np==2 .and. zcu>zero)&
     & .or. (np==8 .and. zcu>zero)) then
! Nie bedzie turning pointow radialnych
        rtp = 0
      endif
      if (np>2 .and. np<7) then
        if (zcu>zero) then
          utp = upl
        else
          utp = zero
        endif
      elseif (np==1 .or. np==7) then
        utp = u2
        if (rtp==1 .and. uk==utp) then
!          print*,"uk:",uk,"upom:",upom
          uk = upom
        endif
      else
        utp = u3
        if (rtp==1 .and. uk==utp) then
          uk = upom
        endif
      endif
! Ustawiam dodatkowa zmienna od liczby turning pointow na zero
      rtp1 = 0
!      print*,"utp:",utp,"upom:",upom,"zcu:",zcu,"uk:",uk,"rtp:",rtp
! Ustawiam krok do obliczania punktow wzdluz trajekorii
      du = sign(one,utp-upom)*zcu*((utp-upom)+(two*rtp-one)*(utp-uk))
!      print*,"du:",du
      du = du/dble(lpt)
!      print*,"du:",du
      kmax = min(int((utp-upom)/du+ofset),lpt)    
! Zabezpieczam sie przed krokiem du=0
      if (du==zero) then
        kmax = 0
      endif
      !print*,"upom:",upom,"up:",up
      !print*,"znak:",sign(one,utp-upom),"du:",du
      !pause
      if ((sign(one,upom-utp)/=sign(one,upom-up)).or.(upom==up)) then    
        if (kmax/=0) then
! Obliczam pierwszy punkt na trajektorii oddzielnie dla kazdej trajektorii, zeby policzyc calki w subroutynie "nazwa"
          k = 1
          !print*,"upom:",upom,"ofset:",ofset,"du:",du
          uni = upom+(k-ofset)*du
          !print*,"up:",up,"uni:",uni,"mup:",mup,"muv:",muv
          !pause
! Wolam procedure geodezyjnamuk
          call geodezyjnamuk(up,uni,mup,muv,l,l2,q2,cu,mutp,rtp1,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
         &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.FALSE.)
!          print*,"k:",k,"rfdcup:",rfdcup,"rfdcuk:",rfdcuk,"rffmu1:",rfdcmup,"rffmu2:",rfdcmuk,"rffmu3:",rfdcmutp
          if (phiit) then
! Wolam procedure geodezyjnaphite
            call geodezyjnaphite(up,uni,mup,muv,l,l2,q2,mutp,rtp1,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,&
           &dphiu,dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,&
           &ctmupplus,ctmutp,cphimupplus,cphimutp,.TRUE.)
          endif
          uki(k) = uni
          muki(k) = muv
          tki(k) = dtu+dtmu
          phiki(k) = dphiu+dphimu
          mutpi(k) = mutp
          rtpi(k) = 0
          lambdaki(k) = lambda
!          print*,"up:",up,"uk:",uk,"mup:",mup,"muk:",muk,"lpt:",lpt
          print*,"k:",k,"uki(1):",uki(1),"rki(1):",one/uki(1),"muki(1):",muki(1),"phi(1):",phiki(1)
!          print*,"upom:",upom,"up:",up
        endif 
! Lecimy naszym fotonem z up do uk lub do turning pointu
        do k=2,kmax
          uni = upom+(k-ofset)*du
! Wolam procedure geodezyjnamuk
          call geodezyjnamuk(up,uni,mup,muv,l,l2,q2,cu,mutp,rtp1,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
         &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.FALSE.)
          if (phiit) then
! Wolam procedure geodezyjnaphite
            call geodezyjnaphite(up,uni,mup,muv,l,l2,q2,mutp,rtp1,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,&
           &dphiu,dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,&
           &ctmupplus,ctmutp,cphimupplus,cphimutp,.FALSE.)
          endif   
          uki(k) = uni
          muki(k) = muv
          tki(k) = dtu+dtmu
          phiki(k) = dphiu+dphimu
          mutpi(k) = mutp
          rtpi(k) = 0
          lambdaki(k) = lambda
          print*,"k:",k,"uki(k):",uki(k),"rki(k):",one/uki(k),"muki(k):",muki(k),"phi(k):",phiki(k)
!          print*,"upom:",upom,"up:",up
        enddo
        !pause
        if (rtp==1) then
          if (lpt==1) then
            elcaltp = .TRUE.
          endif
          if (zmzm) then
! Ta opcja uzywa procedury geodezyjnauk, zeby zagescic punkty w poblizu (r/u) turning pointu, gdzie male kroki w (r/u)
! moga prowadzic do duzych zmian kroku w mu.
! Dlatego najpierw obliczam muk i mutp1 po przeciwnej stronie (r/u) turning pointu (za nim tak jakby)
            kex = 8
            k = kmax+kex+1
            lpptp = 120
! Wolam procedure geodezyjnamuk
            call geodezyjnamuk(uk,uni,mup,muv,l,l2,q2,cu,mutp1,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
           &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.FALSE.)
            muk = muv
! Obliczam pierwiastki rownania M(mu)
            call rozwiazaniamu(q2,l,mup,muminus,muplus)
            do k=1,lpptp
! Wolam procedure do zmiany zmiennych
              call zmiennaniezalezna(muminus,muplus,muki(kmax-kex),muk,mutpi(kmax-kex),mutp1,zcmu,k,lpptp,&
             &ofset,muv,mutp)
! Teraz wolam procedure geodezyjnauk
              call geodezyjnauk(up,uni,mup,muv,l,l2,q2,cu,mutp,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
             &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,.FALSE.)
! Wolam procedure geodezyjnaphite
              call geodezyjnaphite(up,uni,mup,muv,l,l2,q2,mutp,rtp,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,&
             &dphiu,dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,&
             &ctmupplus,ctmutp,cphimupplus,cphimutp,.FALSE.)
              uki(k+kmax-kex) = uni
              muki(k+kmax-kex) = muv
              tki(k+kmax-kex) = dtu+dtmu
              phiki(k+kmax-kex) = dphiu+dphimu
              mutpi(k+kmax-kex) = mutp
              rtpi(k+kmax-kex) = rtp
              lambdaki(k+kmax-kex) = lambda
            enddo
            lpptp = lpptp-2*kex
          endif
! Teraz bede sledzil foton od turning pointu (jezeli obecny) do uk
          do k=1+kmax+kex,lpt
            uni = two*utp-(upom+(k-ofset)*du)
! Wolam procedure geodezyjnamuk
            call geodezyjnamuk(uk,uni,mup,muv,l,l2,q2,cu,mutp,rtp,zcu,zcmu,np,h1,u1,u2,u3,u4,rfdcup,rfdcuk,&
           &rfdcmup,rfdcmuk,rfdcmutp,cu0,cmu1,cmu3,phiit,elcaltp)
            if (phiit) then
! Wolam procedure geodezyjnaphite
              call geodezyjnaphite(uk,uni,mup,muv,l,l2,q2,mutp,rtp,zcu,zcmu,cu,h1,dphimu,dtmu,np,u1,u2,u3,u4,&
             &dphiu,dtu,lambda,rfdcup,rfdcuk,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,ctuw1,ctuw2,ctuw3,ctuw4,&
             &ctmupplus,ctmutp,cphimupplus,cphimutp,elcaltp)
            endif
            uki(k+lpptp) = uni
            muki(k+lpptp) = muv
            tki(k+lpptp) = dtu+dtmu
            phiki(k+lpptp) = dphiu+dphimu
            mutpi(k+lpptp) = mutp
            rtpi(k+lpptp) = 1
            lambdaki(k+lpptp) = lambda
            print*,"k:",k,"uki(k):",uki(k),"rki(k):",one/uki(k),"muki(k):",muki(k),"phi(k):",phiki(k)
!            print*,"upom:",upom,"up:",up
          enddo
          lpt = lpt+lpptp
        endif
        !pause
      else
        lpt = 0
      endif
    endif

    end subroutine 
    
  end module procgeo
