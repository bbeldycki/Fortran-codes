  module const
!------------------------------------
! agata_const.f95:  25.07.2015 
! global declaretions for agata_cc1.f95 
!
!------------------------------------

  implicit none

! wybrane stale fizyczne w jednostkach cgs:

  double precision, parameter  ::  hhh = 6.62620d-27                 ! stala plancka [erg*s]
  double precision, parameter  ::  skkk = 1.38062d-16                ! stala boltzmana [erg/K]
  double precision, parameter  ::  ccc = 2.99792458d+10              ! predkosc swiatla [cm/s]
  double precision, parameter  ::  pmas = 1.6726d-24                 ! masa protonu [g]
  double precision, parameter  ::  nmas = 1.6749d-24                 ! masa neutronu [g]
  double precision, parameter  ::  emas = 9.1093d-28                 ! masa elektronu [g]
  double precision, parameter  ::  sigma = 5.670367d-05              ! stala stefana-boltzmana [erg/s/cm^2/K^4]
  double precision, parameter  ::  eslon = 1.989d+33                 ! masa slonca [g]
  double precision, parameter  ::  graw = 6.67408d-08                ! stala grawitacji [cm^3/g/s^2]   
  double precision, parameter  ::  ekrg = 1.60217d-9                 ! przeliczanie energii z keV na erg
  double precision, parameter  ::  ergf = 1.50916d+26                ! przeliczanie energii z erg na czestosc Hz  
  double precision, parameter  ::  tcros = 6.65245871d-25            ! przekroj czynny na oddzialywanie Thomsona [cm^2]
  double precision, parameter  ::  esk = 0.343                       ! Thomson el. scat.
  double precision, parameter  ::  parsek = 3.0856d+18               ! parsek [cm]

! parametry liczby

  double precision, parameter  ::  one=1.0, two=2.0, hf=0.5, trd=1.0/3.0, zero = 0.0, pi = 3.141592654
  double precision, parameter  ::  try=3.0, ohf=1.5
  
! zmienne globalne,parametry i parametry startowe: 

  integer, parameter  ::  maxlet = 5000                 ! opis parametru patrz plik (zmienne_i_parametry)
  integer, parameter  ::  maxlwr = 5000                 ! opis parametru patrz plik (zmienne_i_parametry)
  integer, parameter  ::  ndr = 501                     ! ndr - liczba punktow w siatce radialnej 
  integer, parameter  ::  nf = 1000                     ! nf - liczba punktow w siatce energii
  integer, parameter  ::  mm = 1000                     ! zmienne do petli gdzie bedzieym chcieli 1000 krokow
  double precision  ::  dgzw = 1.d-5                    ! opis parametru patrz plik (zmienne_i_parametry) 
  double precision, parameter  ::  zr = 100.0           ! opis parametru patrz plik (zmienne_i_parametry)                    
  double precision  ::  etacc = 12.0                    ! visc., 16,0.06 (Hans,Olek) 12,0.08 (Ja) 1 (Andrzej)
!  double precision  ::  em = 10.0                       ! masa BH ! trzeba to zmienic
  double precision  ::  utp                             ! opis parametru patrz plik (zmienne_i_parametry)
  double precision  ::  emin = 0.00001, emax = 100.0    ! energia w keV
!  double precision  ::  rA,rB,rC,rD,rE,rstru            ! parameters for relativistic corrections from Novikov and Thorne 1973  ! trzeba to zmienic
  double precision  ::  rms                             ! orbita marginalnie stabilna
  double precision, parameter  ::  rout = 1.d2          ! promien zewnetrzny dysku 100 R_g 
  double precision  ::  a = 0.998                       ! spin BH
  double precision, parameter  ::  inc = 80.0           ! inklinacja w stopniach
  double precision, parameter  ::  mdot = 3.d-2         ! tempo akrecji (L/L_{edd})


  end module const             


