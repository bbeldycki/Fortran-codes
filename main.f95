    program main
! tutaj potem dodam moduly
    use const
    use procgeo
    use przesuniecia
    implicit none
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
!    double precision  ::  kp,rmin,rmax,rmaxi     ! chwilowo bez uzycia
    double precision  ::  z1,z2,ums,thetem,phiem
    double precision  ::  alfa,beta,l,l2,q2,almin,almax,bmin,bmax,maxab     ! zmienne zwiazane z warunkami brzegowymi/poczatkowymi
    double precision  ::  up,uk,upom,mup,muk,upl,iumin    ! zmienne zwiazane z warunkami brzegowymi/poczatkowymi
    integer  ::  lptz,lpt,ltp,tsd,lpro,lpphi,ql     ! zmienne zwiazane z warunkami brzegowymi/poczatkowymi
    double precision  ::  zcu,zcmu,ofset,dist,gfac
    integer  ::  np,mutp,rtp
    integer  ::  zst1     ! zmienna sterujaca
    integer  ::  i,j,ii,k,ij    ! zmienne petlowe (i-w petli po alfach, j-w petli po betach, ii-w petli po trajektoriach, k-w petli zapisu)  
    character(LEN=20)  ::  outplik     ! zmienne stringi do zapisu lub czytania danych
    integer  ::  nplikt    ! "wskaznik" na plik wyjsciowy
    double precision, dimension(maxlet)  ::  uki,muki,tki,phiki,lambdaki    ! tablice z odpowiednimi wynikami
    integer, dimension(maxlet)  ::  mutpi,rtpi    ! tablice z odpowiednimi wynikami

!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
! Ustawiam liczbe punktow na trajektorie rowna 1 (wartosc ta moze ale nie musi sie zmienic dalej)
    lpt = 1
! Ustawiam wartosc ofsetu 
    ofset = hf
! Na ten moment wszystkie parametry beda ustawiane i zmienane w petli lub recznie.
! Czytanie danych wejsciowych z pliku zaimplementuje w swoim czasie. 
    nplikt = 15
    outplik = 'data.dat'
! Otwieram plik do zapisu
    open(unit=nplikt,file=outplik)

!----------------------------------------------------------------------------------!
! Na ten moment zakladam wartosci roznych parametrow.
! Czesc wartosci znajduje sie w pliku const.f95, a pozostale zostana zadeklarowane ponizej.
!----------------------------------------------------------------------------------!
 !   print*,int(dble(1528883)/dble(1500))+1,mod(1528883-1,1500)+1
    almin = -10.0
    almax = 10.0
    bmin = -10.0
    bmax = 10.0
    i = 1020
    j = 383
    alfa = almin+(almax-almin)*(dble(i-1)+hf)/dble(1500)    
    beta = bmin+(bmax-bmin)*(dble(j-1)+hf)/dble(1500)
!    print*,alfa,beta
!    stop
    mup = dcos(inc*pi/180.0)     ! mu poczatkowe (musi zawierac sie w przedziale -1,1)
!    mup = 0.5
    ql = 1     ! zmienna ta bedzie w przyszlosci czytana z pliku, bedzie ona mowic programowi jak zmienne sie przypisza
! Narazie pracujemy na kwadratowej siatce na detektorze
    tsd = 2
! Siatka rozciaga sie od almin do almax i bmin do bmax
    almin = alfa
    almax = alfa
    bmin = beta
    bmax = beta
    maxab = DMAX1(almin**two,almax**two)**two+DMAX1(bmin**two,bmax**two)**two
! Liczba punktow w kierunku (x/r) i (y/phi)
    lpro = 1
    lpphi = 1 
! Liczba punktow na poszczegolnych trajektoriach
    lpt = 500
! Liczba trajektorii do policzenia
    ltp = lpro*lpphi
! Uplus, to rozmiar event horizonu
    upl = one/(one+dsqrt(one-a**two))
! Odwrotnosc Umin, w programie bedzie wygodniej korzystac z odrotnosci Umin niz z samego Umin
    iumin = one-dsqrt(one-a**two)
! Testowa wielkosc
    dist = 1.d4*parsek*ccc**two/(graw*10.0*eslon)
! Jesli nasz obserwator bedzie znajdowal sie w nieskonczonosci to nasz policzony czas bedzie sie rozbiegal
! Dlatego zrobimy tak, zeby 1/r_obs bylo skonczone
    up = DMIN1(1.d-4,one/(zr*maxab))
    up = 5.d-11
!    print*,"up:",up,"dist:",dist,"rp:",one/up,"rp/dist:",dist/up
!    stop
! Jezeli chcemy policzyc tylko koncowa wartosc u lub mu, to wtedy liczba punktow na trajketorii wynosi lpt = 1
! Musimy wtedy zmienic wartosc ofsetu
    if (lpt==1) then
      ofset = 1.d-8
    endif
! Nadaje wartosci dla zmiennej sterujacej
! W przypadku gdy zst1 = 2, bedziemy obliczac uk, gdy znamy up, mup i muk
! W pozostalych przypadkach bedziemy obliczac muk, gdy znamy up, uk i mup
! Standardowo bedziemy korzystac z czesci, ktora oblicza muk
    zst1 = 1
! determining isco in R_g=G*M/c**2
    z1 = one+(one-a**two)**trd*((one+a)**trd+(one-a)**trd)
    z2 = (try*a**two+z1**two)**hf
    rms = try+z2-((try-z1)*(try+z1+two*z2))**hf
    ums = one/rms              ! more useful unit 
!----------------------------------------------------------------------------------!
!------------Przechodzimy do samego obliczania trajektorii teraz-------------------!
!----------------------------------------------------------------------------------!

    if (zst1==2) then
! Rozpatrujemy przypadek w ktorym obliczamy uk, gdy znamy up, mup i muk
! Petla po wszystkich trajektoriach
      do ii=1,ltp
        i = idint(dble(ii-1)/dble(lpphi))+1
        j = mod(ii-1,lpphi)+1
! Pracuje narazie na kwadratowej siatce, dlatego warunek na siatki zapisze i zostawie zakomentowany
        alfa = almin+(almax-almin)*(dble(i-1)+hf)/dble(lpro)
        beta = bmin+(bmax-bmin)*(dble(j-1)+hf)/dble(lpphi)
        !if (tsd>1) then        
          !alfa = almin+(almax-almin)*(dble(i-1)+hf)/dble(lpro)
          !beta = bmin+(bmax-bmin)*(dble(j-1)+hf)/dble(lpphi)
        !else
          !if (lpphi/=1) then
            !kp = two*pi*dble(j-1)/(dble(lpphi-1))
          !else
            !hp = zero
          !endif
! Obliczam parametry zderzenia dla obserwatora w nieskonczonosci gdy phi0=-pi/2
          ! alfa = pri(i)*dsin(kp)
          ! beta = -pri(i)*dcos(kp)
        !endif
! Oblizam moment pedu i stala Cartera
        l = -alfa*dsqrt(one-mup**two)
        l2 = l**two
        q2 = beta**two-(a**two-alfa**two)*mup**two
! Sprawdzam czy stale l,a,q2 nie sa zbyt male. Jezeli sa to ustawiam ich wartosc na zero
        if (dabs(q2)<dgzw**two) then
          q2 = zero
        endif
        if (dabs(a)<dgzw) then
          a = zero
        endif
        if (dabs(l)<dgzw) then
          l = zero
        endif
! Okreslam znaki calek u i mu
! Znak calki U jest rowny pochodnej u po lambdzie (du/dlambda)
        zcu = one
! Teraz znak calku MU
! Jezeli beta>0 i up jest daleko od czarnej dziury oraz zcu = 1, to mu jest dodatnie
        if (beta>zero .and. mup<one) then
          zcmu = one
        else
          zcmu = -one
        endif
! Teraz oblicze liczba turning pointow w kierunku mu. To dziala tylko gdy muk jest blisko plaszczyzny dysku (muk~0)
        mutp = int((sign(one,mup)*zcmu+one)/two)
        muk = zero
! Wywoluje procedure geodezyjna, ktora obliczy wspolrzedne geodezyjnej i liczbe turning pointow w kazdym z lpt pomiedzy mup i muk
        call geodezyjna(up,uk,upom,mup,muk,l,q2,alfa,beta,mutp,rtp,zcu,zcmu,lpt,ofset,.TRUE.,.TRUE.,.FALSE.,&
       &np,uki,muki,tki,phiki,mutpi,rtpi,lambdaki)
! Zapis geodezyjnej do pliku
        write(nplikt,*) alfa,beta,lpt,np
        write(nplikt,200) uki(1),muki(1),tki(1),phiki(1),lambdaki(1),mutpi(1),rtpi(1)
        do k=2,lpt
          write(nplikt,200) uki(k),muki(k),tki(k),phiki(k),lambdaki(k),mutpi(k),rtpi(k)
        enddo
! Koniec petli po wszystkich trajektoriach
      enddo
    else
! Teraz rozpatruje przypadek standardowy, w ktorym obliczamy muk, gdy znamy up, uk i mup
! Zapamietuje ile powinno byc policzonych punktow dla kazdej poszczegolnej trajektorii
      lptz = lpt
! Petla po wszystkich trajektoriach
      do ii=1,ltp
        if (zst1==1) then
! uk obiera rowne wartsci od upom do horyzontu lub do turning pointu i wraca do upom (dla calych geodezyjnych zachodzi upom=up)   
          upom = up
          uk = upl
! Ustawiamy rtp=1, ale to nie oznacza, ze musimy napotkac na turning point na trajektorii
! Jezeli trafimy na turning point, to policzymy tez punkty za nim
          rtp = 1
          i = idint(dble(ii-1)/dble(lpphi))+1
          j = mod(ii-1,lpphi)+1
! Pracuje narazie na kwadratowej siatce, dlatego warunek na siatki zapisze i zostawie zakomentowany
          alfa = almin+(almax-almin)*(dble(i-1)+hf)/dble(lpro)
          beta = bmin+(bmax-bmin)*(dble(j-1)+hf)/dble(lpphi)
          !if (tsd>1) then        
            !alfa = almin+(almax-almin)*(dble(i-1)+hf)/dble(lpro)
            !beta = bmin+(bmax-bmin)*(dble(j-1)+hf)/dble(lpphi)
          !else
            !if (lpphi/=1) then
              !kp = two*pi*dble(j-1)/(dble(lpphi-1))
            !else
              !hp = zero
            !endif
! Obliczam parametry zderzenia dla obserwatora w nieskonczonosci gdy phi0=-pi/2
            ! alfa = pri(i)*dsin(kp)
            ! beta = -pri(i)*dcos(kp)
          !endif
          alfa = -50.231670647863893
          beta = 0.35987528063294921
! Oblizam moment pedu i stala Cartera
          l = -alfa*dsqrt(one-mup**two)
          l2 = l**two
          q2 = beta**two-(a**two-alfa**two)*mup**two
! Sprawdzam czy stale l,a,q2 nie sa zbyt male. Jezeli sa to ustawiam ich wartosc na zero
! Okreslam znaki calek u i mu
! Znak calki U jest rowny pochodnej u po lambdzie (du/dlambda)
          zcu = one
! Teraz znak calku MU
! Jezeli beta>0 i up jest daleko od czarnej dziury oraz zcu = 1, to mu jest dodatnie
          if (beta>zero .and. mup<one) then
            zcmu = one
          else
            zcmu = -one
          endif
        else
! W tej czesci beda przypisawane wartosci zmiennych w zaleznosci od wczytanej wartosci ql
! Narazie zapisze ql = 1, bo ta czesc kodu uzupelnie pozniej
          ql = 1
        endif
! Sprawdzam czy stale l,a,q2 nie sa zbyt male. Jezeli sa to ustawiam ich wartosc na zero
        if (dabs(q2)<dgzw**two) then
          q2 = zero
        endif
        if (dabs(a)<dgzw) then
          a = zero
        endif
        if (dabs(l)<dgzw) then
          l = zero
        endif 
! Dla testu wyprintuje jeszcze alfa i beta przed obliczeniem trajketorii
        print*,"alfa:",alfa,"beta:",beta
        print*,"l:",l,"q2:",q2
!        l = -49.550285926981090
!        q2 = 4.3720420661632948
        !print*,"l:",l,"q2:",q2
! Ustawiam, ze liczba turning pointow w kierunku mu wyjsciowo jest rowna zero, co moze sie zmienic w trakcie obliczen
        mutp = 0
! Wywoluje procedure geodezyjna, ktora obliczy wspolrzedne geodezyjnej i liczbe turning pointow w kazdym z lpt pomiedzy mup i muk
        call geodezyjna(up,uk,upom,mup,muk,l,q2,alfa,beta,mutp,rtp,zcu,zcmu,lpt,ofset,.TRUE.,.FALSE.,.FALSE.,&
       &np,uki,muki,tki,phiki,mutpi,rtpi,lambdaki)   
!        stop
        ij = 1
        do k=1,lpt              ! lookin for the part of trajectory from observer to the disk
!          print*,"k:",k,"muki:",muki(k),"uki:",uki(k),"ums:",ums,"phi:",phiki(k)
          if (muki(k)<zero) then
!            print*,"k:",k,"muki:",muki(k) 
            if (uki(k-1)>(one/100.0) .and. ums>=(uki(k-1)) .and. muki(k-1)>=zero) then
              ij = k-1
              exit
            endif
          endif             
        enddo
        if (phiki(ij)<zero) then
          phiem = two*pi+phiki(ij)
        else
          phiem = phiki(ij)
        endif
        print*,"Jestem tutaj?"
!        print*,"ij:",ij,"uki:",uki(ij),"rki:",one/uki(ij),"phiki:",phiem
        if (ums>=uki(ij) .and. uki(ij)>(one/100.0) .and. muki(ij)>=zero) then         ! photons hitting the disk 
!          call nazwa(it,pname)   
!          open (56,file=pname,status='new') 
!          write(56,*) alfa,beta,ij
!          write(56,*) ' u        mu      phi'
!          do k=1,ij
!            write(56,*) uki(k),muki(k),phiki(k)
!          enddo
!          close (56) 
          call efektydp(uki(ij),l,q2,gfac,thetem,phiem)
          print*,"gfac:",gfac,"l:",l,"q:",dsqrt(q2)
!          write(nplikt,*) ii,alfa,beta,l,q2,uki(ij),muki(ij),phiki(ij),gfac
        endif    
        stop
! Zapis geodezyjnej do pliku
        write(nplikt,*) alfa,beta,lpt,np
        write(nplikt,200) uki(1),muki(1),tki(1),phiki(1),lambdaki(1),mutpi(1),rtpi(1)
        do k=2,lpt
          write(nplikt,200) uki(k),muki(k),tki(k),phiki(k),lambdaki(k),mutpi(k),rtpi(k)
        enddo  
! Koniec petli po wszystkich trajektoriach
      enddo
    endif
! formatowanie danych wyjsciowych
200 format(5(1x,1d16.8),2(i2))

!----------------------------------------------------------------------------------!
!------------------------Koniec obliczania trajektorii-----------------------------!
!----------------------------------------------------------------------------------!

    end program
