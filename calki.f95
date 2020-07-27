  module calki
!-------------------------------------------
! roots.f95
! module roots
! Opis zrobie pozniej
!-------------------------------------------

  implicit none
  
  public  ::  trzyrzeczywiste,dwazjedrzecz,czteryrzeczywiste,dwazdwarzecz,czteryzespolone,&
             &calkieliptycznelegendra,calkiphiitniesymetryczne,calkiphiitsymetryczne,calkaimuniesymetryczna,&
             &calkaimuniesymetrycznak,calkaimusymetryczna,calkaimusymetrycznak,&
             &tuszescian,phiuszescian
  contains   
  
    double precision function trzyrzeczywiste(p,a1,b1,a2,b2,a3,b3,a4,b4,ffr,y,x)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki eliptycznej (Carlsona) 
! \int od y do x dt Iloczyn_{i=1,4} (a_{i}+b_{i}t)^{p_{i}/2}
! Wartosci x, y musza byc nieujemne oraz zachodzi y<x
! Funkcja zaczerpnieta z Carlson (1989) tam jest wiecej informacji
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    integer,dimension(5)  ::  p
    double precision  ::  a1,a2,a3,a4,b1,b2,b3,b4,ffr,y,x
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  d12,d13,d14,d24,d34
    double precision  ::  r12,r13,r14i,r24i,r34i
    double precision  ::  x1,x2,x3,x4,y1,y2,y3,y4
    double precision  ::  u1,u2,u3,u12,u22,u32
    double precision  ::  w22,q22,p22
    double precision  ::  i1c,i2c,i3c,k2c,a111m2
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
! Carlson (1989) rownanie 2.1      
    d12 = a1*b2-a2*b1
    d13 = a1*b3-a3*b1
    d14 = a1*b4-a4*b1
    d24 = a2*b4-a4*b2
    d34 = a3*b4-a4*b3
! Carlson (1989) rownanie 2.2
    x1 = dsqrt(a1+b1*x)
    x2 = dsqrt(a2+b2*x)
    x3 = dsqrt(a3+b3*x)
    x4 = dsqrt(a4+b4*x)
    y1 = dsqrt(a1+b1*y)
    y2 = dsqrt(a2+b2*y)
    y3 = dsqrt(a3+b3*y)
    y4 = dsqrt(a4+b4*y)
! Carlson (1989) rownanie 2.3
    u1 = (x1*y2*y3+y1*x2*x3)/(x-y)
    u2 = (x2*y1*y3+y2*x1*x3)/(x-y)
    u3 = (x3*y1*y2+y3*x1*x2)/(x-y)
    u12 = u1**two
    u22 = u2**two
    u32 = u3**two
! Carlson (1989) rownanie 2.4
    w22 = u12-b4*d12*d13/d14
! Carlson (1989) rownanie 2.5
    q22 = w22*(x4*y4/(x1*y1))**two
    p22 = q22 + b4*d24*d34/d14
! Musze policzyc 3 calki dla p=[-1,-1,-1,0], p=[-1,-1,-1,-2], p=[-1,-1,-1,-4]
    if (p(4)==0) then
! Carlson (1989) rownanie 2.12
      ffr = rf(u32,u22,u12)
      trzyrzeczywiste = two*ffr
      return
    else
! Carlson (1989) rownanie 2.12
      i1c = two*ffr
! Carlson (1989) rownanirownanie 2.14
      i3c = -two*b1*d12*d13*rj(u32,u22,u12,w22)/(try*d14)+two*rc(p22,q22)
      if (p(4)==-2) then
! Carlson (1989) rownanie 2.49
        trzyrzeczywiste = (b4*i3c-b1*i1c)/d14
        return
      elseif (p(4)==-4) then
! Carlson (1989) rownanie 2.13
        i2c = two*d12*d13*rd(u32,u22,u12)/try+two*x1*y1/u1
! Carlson (1989) rownanie 2.1
        r12 = (a1*b2-a2*b1)/(b1*b2)
        r13 = (a1*b3-a3*b1)/(b1*b3)
        r14i = b1*b4/(a1*b4-a4*b1)
        r24i = b2*b4/(a2*b4-a4*b2)
        r34i = b3*b4/(a3*b4-a4*b3)
! Carlson (1989) rownanie 2.6
        a111m2 = x1*x2/(x3*x4*x4)-y1*y2/(y3*y4*y4)
! Carlson (1989) rownanie 2.59
        k2c = b2*b3*i2c-two*b4*a111m2
! Carlson (1989) rownanie 2.62
        trzyrzeczywiste = b4*k2c/(two*d14*d24*d34)+(b1/d14)**two*(one-r12*r13*r24i*r34i/two)*i1c
!        trzyrzeczywiste = -one*(r14i+r24i+r34i)*i3c/(two*d14)+b4*k2c/(two*d14*d24*d34)+(b1/d14)**two*&
!       &(one-r12*r13*r24i*r34i/two)*i1c
        return
      endif
    endif
    end function


    double precision function dwazjedrzecz(p,a1,b1,a4,b4,f,g,h,ffr,y,x)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki eliptycznej (Carlsona) 
! \int od y do x dt Iloczyn_{i=1;4} (a_{i}+b_{i}t)^{p_{i}/2}(f+gt+ht^2)^{p_2/2}
! Wartosci x, y musza byc nieujemne oraz zachodzi y<x
! Funkcja zaczerpnieta z Carlson (1991) tam jest wiecej informacji
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  a1,b1,a4,b4,f,g,h,ffr,y,x
    integer,dimension(5)  ::  p
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  x1,y1,x4,y4,d14
    double precision  ::  beta1,beta4
    double precision  ::  c11,c44,c142
    double precision  ::  xi,eta
    double precision  ::  m2
    double precision  ::  lp2,lm2,wp2
    double precision  ::  u,u2,w2
    double precision  ::  q,q2,p2
    double precision  ::  rho
    double precision  ::  r12xr13,r24xr34,r14i,mr24pmr34
    double precision  ::  i1c,k2c,i3c,n2c
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
! Carlson (1991) rownanie 2.1
    x1 = dsqrt(a1+b1*x)
    y1 = dsqrt(a1+b1*y)
    x4 = dsqrt(a4+b4*x)
    y4 = dsqrt(a4+b4*y)
    d14 = a1*b4-a4*b1
! Carlson (1991) rownanie 2.2
    beta1 = g*b1-two*h*a1
    beta4 = g*b4-two*h*a4
! Carlson (1991) rownanie 2.3
    c11 = dsqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
    c44 = dsqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
    c142 = two*f*b1*b4-g*(a1*b4-a4*b1)+two*h*a1*a4
! Carlson (1991) rownanie 2.4
    xi = dsqrt(f+g*x+h*x*x)
    eta = dsqrt(f+g*y+h*y*y)
! Carlson (1991) rownanie 3.1
    m2 = ((x1+y1)*dsqrt(two*xi*eta+two*f+g*(x+y)+two*h*x*y)/(x-y))**two
! Carlson (1991) rownanie 3.2
    lp2 = m2-beta1+dsqrt(two*h)*c11
    lm2 = m2-beta1-dsqrt(two*h)*c11
! Musze policzyc 3 calki dla p=[-1,-1,-1,0], p=[-1,-1,-1,-2], p=[-1,-1,-1,-4]
    if (p(4)==0) then
      ffr = rf(m2,lm2,lp2)
! Carlson (1991) rownanie 3.8
      dwazjedrzecz = 4.0*ffr
      return
    else
! Carlson (1991) rownanie 3.8
      i1c = 4.0*ffr
! Carlson (1991) rownanie 3.2
      wp2 = m2-b1*(c142+c11*c44)/d14
! Carlson (1991) rownanie 3.3
      u = (x1*eta+y1*xi)/(x-y)
      u2 = u**two
      w2 = u2-c11*c11*b4/(two*d14)
! Carlson (1991) rownanie 3.4
      q = x4*y4*dsqrt(w2)/(x1*y1)
      q2 = w2*(x4*y4/(x1*y1))**two
      p2 = q2+c44*c44*b4/(two*d14)
! Carlson (1991) rownanie 3.5
      rho = dsqrt(two*h)*c11-beta1
! Carlson (1991) rownanie 3.10
      i3c = (two*c11/(try*c44))*((-4.0*b1/d14)*(c142+c11*c44)*rj(m2,lm2,lp2,wp2)-&
     &6.0*rf(m2,lm2,lp2)+try*rc(u2,w2))+two*rc(p2,q2)
      if (p(4)==-2) then
! Carlson (1989) rownanie 2.49          
        dwazjedrzecz = (b4*i3c-b1*i1c)/d14
        return
      elseif (p(4)==-4) then
! Carlson (1991) rownanie 3.11
        n2c = dsqrt(8.0*h/(9.0*c11*c11))*(4.0*rho*rd(m2,lm2,lp2)-6.0*rf(m2,lm2,lp2)+&
       &try/u)+two/(x1*y1*u)
! Carlson (1991) rownanie 3.12
        k2c = c11*c11*n2c/two-two*d14*(xi/(x1*x4*x4)-eta/(y1*y4*y4))
! Carlson (1991) rownanie 2.19
        r12xr13 = c11*c11/(two*h*b1*b1)
        r24xr34 = c44*c44/(two*h*b4*b4)
        r14i = b1*b4/d14
        mr24pmr34 = -two*b4*beta4/(c44*c44)
! Carlson (1989) rownanie 2.62
        dwazjedrzecz = k2c/(4.0*d14*b4*h*r24xr34)+(b1/d14)**two*(one-r12xr13/(two*r24xr34))*i1c
!          dwazjedrzecz = -one*i3c*(r14i+mr24pmr34)/(two*d14)+k2c/(4.0*d14*b4*h*r24xr34)+&
!         &(b1/d14)**two*(one-r12xr13/(two*r24xr34))*i1c
        return
      endif
    endif
    end function  


    double precision function czteryrzeczywiste(p,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,ffr,y,x)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki eliptycznej (Carlsona) 
! \int od y do x dt Iloczyn_{i=1,5} (a_{i}+b_{i}t)^{p_{i}/2}
! Wartosci x, y musza byc nieujemne oraz zachodzi y<x
! Funkcja zaczerpnieta z Carlson (1988) tam jest wiecej informacji
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    integer,dimension(5)  ::  p
    double precision  ::  a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,ffr,y,x
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  d12,d13,d14,d15,d24,d25,d34,d35,d45,r12,r13,r15,r25,r35,r45!,r15i,r25i,r35i,r45i
    double precision  ::  x1,x2,x3,x4,x5,y1,y2,y3,y4,y5
    double precision  ::  u122,u132,u142,u14
    double precision  ::  w,w2
    double precision  ::  q2,p2
    double precision  ::  a111m1m2,sums
    double precision  ::  i1,i2,i3
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
!  Carlson (1988) rownanie 2.1
    d12 = a1*b2-a2*b1
    d13 = a1*b3-a3*b1
    d14 = a1*b4-a4*b1
    d15 = a1*b5-a5*b1
    d24 = a2*b4-a4*b2
    d25 = a2*b5-a5*b2
    d34 = a3*b4-a4*b3
    d35 = a3*b5-a5*b3
    d45 = a4*b5-a5*b4
    r12 = (a1*b2-a2*b1)/(b1*b2)
    r13 = (a1*b3-a3*b1)/(b1*b3)
    r15 = (a1*b5-a5*b1)/(b1*b5)
    r25 = (a2*b5-a5*b2)/(b2*b5)
    r35 = (a3*b5-a5*b3)/(b3*b5)
    r45 = (a4*b5-a5*b4)/(b4*b5)
!  Carlson (1988) rownanie 2.2
    x1 = dsqrt(a1+b1*x)
    x2 = dsqrt(a2+b2*x)
    x3 = dsqrt(a3+b3*x)
    x4 = dsqrt(a4+b4*x)
    x5 = dsqrt(a5+b5*x)
    y1 = dsqrt(a1+b1*y)
    y2 = dsqrt(a2+b2*y)
    y3 = dsqrt(a3+b3*y)
    y4 = dsqrt(a4+b4*y)
    y5 = dsqrt(a5+b5*y)
!  Carlson (1988) rownanie 2.3
    u122 = ((x1*x2*y3*y4+y1*y2*x3*x4)/(x-y))**two
    u132 = ((x1*x3*y2*y4+y1*y3*x2*x4)/(x-y))**two
    u14 = (x1*x4*y2*y3+y1*y4*x2*x3)/(x-y)
    u142 = u14**two
!  Carlson (1988) rownanie 2.4
    w2 = u122-d13*d14*d25/d15
    w = dsqrt(w2)
!  Carlson (1988) rownanie 2.5
    q2 = w2*(x5*y5/(x1*y1))**two
    p2 = q2+d25*d35*d45/d15
!  Carlson (1988) rownanie 2.6
    a111m1m2 = x1*x2*x3/(x4*x5*x5)-y1*y2*y3/(y4*y5*y5)
! Musze policzyc 3 calki dla p=[-1,-1,-1,-1], p=[-1,-1,-1,-1,-2], p=[-1,-1,-1,-1,-4]
    if (p(5)==0) then
!  Carlson (1988) rownanie 2.13
      ffr = rf(u122,u132,u142)
      czteryrzeczywiste = two*ffr
      return
    else
!  Carlson (1988) rownanie 2.13
      i1 = two*ffr
!  Carlson (1988) rownanie 2.15
      i3 = two*d12*d13*d14*rj(u122,u132,u142,w2)/(try*d15)+two*rc(p2,q2)
      if (p(5)==-2) then
!  Carlson (1988) rownanie 2.35
        czteryrzeczywiste = (b5*i3-b1*i1)/d15
        return
      elseif (p(5)==-4) then
!  Carlson (1988) rownanie 2.14
        i2 = two*d12*d13*rd(u122,u132,u142)/try+two*x1*y1/(x4*y4*u14)
        sums = (r25*r35*r45+r15*r35*r45+r15*r25*r45+r15*r25*r35)/(r15*r25*r35*r45)
!  Carlson (1988) rownanie 2.49
!          czteryrzeczywiste = -one*i3*sums/(two*d15)+b5*b5*d24*d34*i2/(two*d15*d25*d35*d45)+&
!         &(b1/d15)**two*(one-r12*r13/(two*r25*r35))*i1-b5*b5*a111m1m2/(d15*d25*d35)
        czteryrzeczywiste = b5*b5*d24*d34*i2/(two*d15*d25*d35*d45)+&
       &(b1/d15)**two*(one-r12*r13/(two*r25*r35))*i1-b5*b5*a111m1m2/(d15*d25*d35)
        return
      endif
    endif
    czteryrzeczywiste = zero
    return
    end function


    double precision function dwazdwarzecz(p,a1,b1,a4,b4,a5,b5,f,g,h,ffr,y,x)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki eliptycznej (Carlsona) 
! \int od y do x dt Iloczyn_{i=1;4;5} (a_{i}+b_{i}t)^{p_{i}/2}(f+gt+ht^2)^{p_2/2}
! Wartosci x, y musza byc nieujemne oraz zachodzi y<x
! Funkcja zaczerpnieta z Carlson (1991) tam jest wiecej informacji
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    integer,dimension(5)  ::  p
    double precision  ::  a1,b1,a4,b4,a5,b5,f,g,h,ffr,y,x
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  x1,y1,x4,y4,x5,y5,d14,d15,d24d34,d25d35,d45,beta5
    double precision  ::  c11,c14,c15,c44,c55,c112,c142,c152,c442,c552
    double precision  ::  r15,r45,r12r13,r25r35,r25mr35m
    double precision  ::  xi,eta
    double precision  ::  a111m1m2
    double precision  ::  m2
    double precision  ::  lp2,lm2,wp2
    double precision  ::  u,u2,w2
    double precision  ::  q,q2,p2
    double precision  ::  i1,i2,i3
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
!  Carlson (1991) rownanie 2.1     
    x1 = dsqrt(a1+b1*x)
    x4 = dsqrt(a4+b4*x)
    x5 = dsqrt(a5+b5*x)
    y1 = dsqrt(a1+b1*y)
    y4 = dsqrt(a4+b4*y)
    y5 = dsqrt(a5+b5*y)
    d14 = a1*b4-a4*b1
    d15 = a1*b5-a5*b1
    d45 = a4*b5-a5*b4
!  Carlson (1991) rownanie 2.2
    beta5 = g*b5-two*h*a5
!  Carlson (1991) rownanie 2.3
    c11 = dsqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
    c14 = dsqrt(two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4)
    c15 = dsqrt(two*f*b1*b5-g*(a1*b5+a5*b1)+two*h*a1*a5)
    c44 = dsqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
    c55 = dsqrt(two*f*b5*b5-two*g*a5*b5+two*h*a5*a5)
    c112 = two*f*b1*b1-two*g*a1*b1+two*h*a1*a1
    c142 = two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4
    c152 = two*f*b1*b5-g*(a1*b5+a5*b1)+two*h*a1*a5
    c442 = two*f*b4*b4-two*g*a4*b4+two*h*a4*a4
    c552 = two*f*b5*b5-two*g*a5*b5+two*h*a5*a5
!  Carlson (1991) rownanie 2.19
    d24d34 = c442/two
    d25d35 = c552/two
    r15 = d15/(b1*b5)
    r45 = d45/(b4*b5)
    r12r13 = c112/(two*h*b1*b1)
    r25r35 = c552/(two*h*b5*b5)
    r25mr35m = two*b5*beta5/c552
!  Carlson (1991) rownanie 2.4
    xi = dsqrt(f+g*x+h*x*x)
    eta = dsqrt(f+g*y+h*y*y)
!  Carlson (1991) rownanie 2.5
    a111m1m2 = x1*xi/(x4*x5*x5)-y1*eta/(y4*y5*y5)
!  Carlson (1991) rownanie 2.6
    m2 = ((x1*y4+y1*x4)*dsqrt((xi+eta)**two-h*(x-y)**two)/(x-y))**two
!  Carlson (1991) rownanie 2.7
    lp2 = m2+c142+c11*c44
    lm2 = m2+c142-c11*c44
    wp2 = m2+d14*(c152+c11*c55)/d15
!  Carlson (1991) rownanie 2.8
    u = (x1*x4*eta+y1*y4*xi)/(x-y)
    u2 = u*u
    w2 = u2-c112*d45/(two*d15)
!  Carlson (1991) rownanie 2.9
    q = one
    q2 = w2*(x5*y5/(x1*y1))**two
    p2 = q2+c552*d45/(two*d15)
! Musze policzyc 3 calki dla p=[-1,-1,-1,-1], p=[-1,-1,-1,-1,-2], p=[-1,-1,-1,-1,-4]
    if (p(5)==0) then
      ffr = rf(m2,lm2,lp2)
      dwazdwarzecz = 4.0*ffr
      return
    else
!  Carlson (1991) rownanie 2.14
      i1 = 4.0*ffr
!  Carlson (1991) rownanie 2.16
      i3 = two*c11*(4.0*d14*(c152+c11*c55)*rj(m2,lm2,lp2,wp2)/d15-6.0*rf(m2,lm2,lp2)+&
     &3.0*rc(u2,w2))/(try*c55)+two*rc(p2,q2)
!        print*,"i3:",i3
      if (p(5)==-2) then
!  Carlson (1988) rownanie 2.35
        dwazdwarzecz = (b5*i3-b1*i1)/d15
        return
      elseif (p(5)==-4) then
!  Carlson (1991) rownanie 2.15        
        i2 = two*c11*(4.0*(c142+c11*c44)*rd(m2,lm2,lp2)-6.0*rf(m2,lm2,lp2)+try/u)/(try*c44)+&
       &two*x1*y1/(x4*y4*u)
!         print*,"i2:",i2
!  Carlson (1988) rownanie 2.49 
        dwazdwarzecz = b5*b5*d24d34*i2/(two*d15*d25d35*d45)+b1*b1*(one-r12r13/(two*r25r35))*i1/(d15*d15)&
       &-b5*b5*a111m1m2/(d15*d25d35)
!        dwazdwarzecz = -one*i3*(one/r15+r25mr35m+one/r45)/(two*d15)+b5*b5*d24d34*i2&
!       &/(two*d15*d25d35*d45)+b1*b1*(one-r12r13/(two*r25r35))*i1/(d15*d15)-b5*b5*a111m1m2/(d15*d25d35)
        return
      endif
    endif
    dwazdwarzecz = zero
    return
    end function


    double precision function czteryzespolone(p,f1,g1,h1,f2,g2,h2,a5,b5,ffr,y,x)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki eliptycznej (Carlsona) 
! \int od y do x dt (f_{1}+g_{1}t+h_{1}t^2)^{p_1/2}(f_{2}+g_{2}t+h_{2}t^2)^{p_2/2}(a_{5}+b_{5}t)^{p_{5}/2}
! Wartosci x, y musza byc nieujemne oraz zachodzi y<x
! Funkcja zaczerpnieta z Carlson (1992) tam jest wiecej informacji
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    integer,dimension(5)  ::  p
    double precision  ::  f1,g1,h1,f2,g2,h2,a5,b5,ffr,y,x
!-------------------------------zmienne lokalne------------------------------------!
    double precision  ::  xi1,xi2,eta1,eta2
    double precision  ::  xi1p,eta1p
    double precision  ::  b,e
    double precision  ::  th1,th2
    double precision  ::  hi1,hi2
    double precision  ::  u,u2,m,m2
    double precision  ::  dl11,dl12,dl22,dl112,dl122,dl222,dd
    double precision  ::  ddp,ddm,lp2,lm2
    double precision  ::  g
    double precision  ::  sig
    double precision  ::  alfa15,alfa25,beta15,beta25
    double precision  ::  gm1,gm2
    double precision  ::  dlam,dsig2
    double precision  ::  psi,psi2
    double precision  ::  xi5,eta5
    double precision  ::  am11m1,a1111m2,a1111m4
    double precision  ::  dx
    double precision  ::  ds
    double precision  ::  mu,dt,v2
    double precision  ::  b2
    double precision  ::  a2
    double precision  ::  dh
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
!  Carlson (1992) rownanie 2.1
    xi1 = dsqrt(f1+g1*x+h1*x*x)
    xi2 = dsqrt(f2+g2*x+h2*x*x)
    eta1 = dsqrt(f1+g1*y+h1*y*y)
    eta2 = dsqrt(f2+g2*y+h2*y*y)      
!  Carlson (1992) rownanie 2.2
    xi1p = (g1+two*h1*x)/(two*xi1)
    eta1p = (g1+two*h1*y)/(two*eta1)
!  Carlson (1992) rownanie 2.3
    b = xi1p*xi2-eta1p*eta2
    e = xi1p*xi1**two*xi2-eta1p*eta1**two*eta2
!  Carlson (1992) rownanie 2.4
    th1 = two*f1+g1*(x+y)+two*h1*x*y
    th2 = two*f2+g2*(x+y)+two*h2*x*y
!  Carlson (1992) rownanie 2.5
    hi1 = dsqrt(two*xi1*eta1+th1)
    hi2 = dsqrt(two*xi2*eta2+th2)
!  Carlson (1992) rownanie 2.6
    u = (xi1*eta2+eta1*xi2)/(x-y)
    u2 = u**two
    m = hi1*hi2/(x-y)
    m2 = m**two
!  Carlson (1992) rownanie 2.7
    dl11 = dsqrt(4.0*f1*h1-g1*g1)
    dl12 = dsqrt(two*f1*h2+two*f2*h1-g1*g2)
    dl22 = dsqrt(4.0*f2*h2-g2*g2)
    dl112 = 4.0*f1*h1-g1*g1
    dl122 = two*f1*h2+two*f2*h1-g1*g2
    dl222 = 4.0*f2*h2-g2*g2
    dd = dsqrt(dl122**two-dl112*dl222)      
!  Carlson (1992) rownanie 2.8
    ddp = dl122+dd
    ddm = dl122-dd
    lp2 = m2+ddp
    lm2 = m2+ddm
!  Carlson (1992) rownanie 2.9
    g = two*dd*ddp*rd(m2,lm2,lp2)/try+dd/(two*u)+(dl122*th1-dl112*th2)/(4.0*xi1*eta1*u)
! Musze policzyc 3 calki dla p=[-1,-1,-1,-1], p=[-1,-1,-1,-1,-2], p=[-1,-1,-1,-1,-4]
    if (p(5)==0) then 
!  Carlson (1992) rownanie 2.10
      ffr = rf(m2,lm2,lp2)
!  Carlson (1992) rownanie 2.36
      czteryzespolone = 4.0*ffr
      return
    else   
!  Carlson (1992) rownanie 2.10
      sig = g-ddp*ffr+b 
!  Carlson (1992) rownanie 2.11
      alfa15 = two*f1*b5-g1*a5
      alfa25 = two*f2*b5-g2*a5
      beta15 = g1*b5-two*h1*a5
      beta25 = g2*b5-two*h2*a5
!  Carlson (1992) rownanie 2.12
      gm1 = f1*b5*b5-g1*a5*b5+h1*a5*a5
      gm2 = f2*b5*b5-g2*a5*b5+h2*a5*a5
!  Carlson (1992) rownanie 2.13
      dlam = dl112*gm2/gm1
      dsig2 = m2+dlam
!  Carlson (1992) rownanie 2.14
      psi = (alfa15*beta25-alfa25*beta15)/two
      psi2 = psi**two
!  Carlson (1992) rownanie 2.15
      xi5 = a5+b5*x
      eta5 = a5+b5*y
!  Carlson (1992) rownanie 2.16
      am11m1 = xi2/xi1-eta2/eta1
      a1111m2 = xi1*xi2/xi5-eta1*eta2/eta5
      a1111m4 = xi1*xi2/(xi5**two)-eta1*eta2/(eta5**two)
!  Carlson (1992) rownanie 2.17
      dx = xi5*eta5*(th1*am11m1/two-xi5*eta5*a1111m4)/(x-y)**two
!  Carlson (1992) rownanie 2.18
      ds = (m2+dl122)/two-u2
!  Carlson (1992) rownanie 2.19
      mu = gm1*xi5*eta5/(xi1*eta1)
      dt = mu*ds+two*gm1*gm2
      v2 = mu**two*(ds**two+dlam*u2)
!  Carlson (1992) rownanie 2.20
      b2 = (ds**two/u2+dlam)*dsig2**two
!  Carlson (1992) rownanie 2.21
      a2 = b2 + dlam**two*psi2/(gm1*gm2)
!  Carlson (1992) rownanie 2.22
      dh = dl112*psi*(rj(m2,lm2,lp2,dsig2)/try+rc(a2,b2)/two)/(gm1**two)-dx*rc(dt**two,v2)
      if (p(5)==-2) then
!  Carlson (1992) rownanie 2.39
        czteryzespolone = -two*(b5*dh+beta15*ffr/gm1)
        return
      elseif (p(5)==-4) then
!  Carlson (1992) rownanie 2.41
        czteryzespolone = b5*(beta15/gm1+beta25/gm2)*dh+beta15**two*ffr/(gm1**two)+b5**two*&
      &(sig-b5*a1111m2)/(gm1*gm2)
        return
      endif
    endif
    end function


    subroutine calkieliptycznelegendra(phi0,phif,m1,n,elcaltp,rf0,rff,rfc,rdc,rjc,elle0,ellef,ellec,ellpi0,&
   &ellpif,ellpic)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Procedura oblicza wartosci calek eliptycznych Legendrea (wiecej informacji w ksiazce podanej
! ponizej)
! Procedura jest kombinacja kilki funkcji z Numerical Recipies in Fortran Press et al (1992)

! elle0 - calka eliptyczna Legendrea drugiego rodzaju od mup do muplus
! ellef - calka eliptyczna Legendrea drugiego rodzaju od muk do muplus
! ellec - zupelna calka eliptyczna Legendrea drugiego rodzaju
! ellpi0 - calka eliptyczna Legendrea trzeciego rodzaju od mup do muplus
! ellpif - calka eliptyczna Legendrea trzeciego rodzaju od muk do muplus
! ellpic - zupelna calka eliptyczna Legendea trzeciego rodzaju
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  phi0,phif,m1,n,rf0,rff,rfc,rdc,rjc
    double precision  ::  elle0,ellef,ellec,ellpi0,ellpif,ellpic
    logical  ::  elcaltp
    double precision  ::  s0,sf,ak,q0,qf,rj0,rjf,rd0,rdf     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
! Przechodze na parametr k powiazany z parametrem m1 nastepujaca relacja
    ak = dsqrt(one-m1)
    if (elcaltp) then
! Najpierw obliczam mup, muplus skladniki o ile nie byly policzone wczesniej
      s0 = dsin(phi0)
      q0 = (one-s0*ak)*(one+s0*ak)
! Tu korzystam z calek Carlsona, zeby znalezc calke Legendre'a drugiego rodzaju
      rd0 = rd(one-s0*s0,q0,one)
      rdc = rd(zero,m1,one)
      elle0 = s0*(rf0-((s0*ak)**two)*rd0/try)
      ellec = two*(rfc-ak*ak*rdc/try)
! jezeli bylo ustawione n=1, to w tym przypadku skladniki zwiazane z phi sa rowne 0
      if (n/=one) then
! Tu korzystam z calek Carlsona, zeby znalezc calke Legendre'a trzeciego rodzaju
        rj0 = rj(one-s0*s0,q0,one,one-n*s0*s0)
        rjc = rj(zero,m1,one,one-n)
        ellpi0 = s0*(rf0+n*s0*s0*rj0/try)
        ellpic = two*(rfc+n*rjc/try)
      else
        ellpi0 = zero
        ellpic = zero
      endif
! Procedury z Press et al (1992) dzialaja tylko gdy (0<=phi<=pi/2). W przypadku gdy (phi>pi/2)
! obliczam calka od 0 do pi/2 - calka od pi/2 do pi. Zeby otrzymac poprawna wartosci, robimy tak
! calka od 0 do pi = 2 calka od 0 do pi/2 - (calka od 0 do pi/2 i calka od pi/2 do pi)
      if (phi0>(pi/two)) then
        ellpi0 = ellpic-ellpi0
        elle0 = ellec-elle0
      endif
    endif
! Powtarzam to co wyzej dla czesci zwiazanej z muk
    sf = sin(phif)
    qf = (one-sf*ak)*(one+sf*ak)
    rdf = rd(one-sf*sf,qf,one)
    ellef = sf*(rff-((sf*ak)**two)*rdf/try)
    if (n/=one) then
      rjf = rj(one-sf*sf,qf,one,one-n*sf*sf)
      ellpif = sf*(rff+n*sf*sf*rjf/try)
    else
      ellpif = zero
    endif
    if (phif.gt.(pi/two)) then
      ellpic = two*(rfc+n*rjc/try)
      ellec = two*(rfc-ak*ak*rdc/try)
      ellpif = ellpic-ellpif
      ellef = ellec-ellef
    endif
    return
    end subroutine


    subroutine calkiphiitniesymetryczne(mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,&
   &ctmukplus,ctmutp,rfdcmup,rfdcmuk,rfdcmut,rd2r,rj3r,elcaltp)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza calki z PHI_{mu} i T_{mu} w przypadku antysymetrycznym tzn. muminus jest wieksze od zero
! M_{-} > 0
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,ctmukplus,ctmutp
    double precision  ::  rfdcmup,rfdcmuk,rfdcmut,rd2r,rj3r
    logical  ::  elcaltp
    double precision  ::  fac,fac2,m1,n,elle0,ellef,ellec,ellpi0,ellpif,ellpic,phi0,phif     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!      
    fac = dabs(a)*muplus
    m1 = mne/mpo
! Jezeli orbita przecina biegun, to nie licze calek phi
    if (mpo/=one) then
      n = (mpo-mne)/(one-mpo)
      fac2 = one-mpo
    else
      n = -one
      fac2 = one
    endif
    phi0 = dasin(dsqrt((mpo-mup**two)/(mpo-mne)))
    phif = dasin(dsqrt((mpo-muk**two)/(mpo-mne)))
! Wolam procedure do obliczenia calek Legendra 2 i 3 rodzaju
    call calkieliptycznelegendra(phi0,phif,m1,-n,elcaltp,rfdcmup,rfdcmuk,rfdcmut,rd2r,rj3r,elle0,ellef,ellec&
   &,ellpi0,ellpif,ellpic)
    if (elcaltp) then
      ctmupplus = fac*elle0
      ctmutp = fac*ellec/two
      cphimupplus = ellpi0/(fac*fac2)
      cphimutp = ellpic/(two*fac*fac2)
    endif
    ctmukplus = fac*ellef
    cphimukplus = ellpif/(fac*fac2)
    return
    end subroutine


    subroutine calkiphiitsymetryczne(mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,&
   &ctmukplus,ctmutp,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,elcaltp)
    use const
    implicit none
! Cel procedury:
! Procedura oblicza calki z PHI_{mu} i T_{mu} w przypadku symetrycznym tzn. muminus jest mniejsze od zero
! M_{-} < 0
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,mup,muk,muplus,cphimupplus,cphimukplus,cphimutp,ctmupplus,ctmukplus,ctmutp
    double precision  ::  rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r
    logical  ::  elcaltp
    double precision  ::  fac,fac2,m1,n,elle0,ellef,ellec,ellpi0,ellpif,ellpic,phi0,phif     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!   
    fac = dabs(a)*dsqrt(mpo-mne)
    m1 = -mne/(mpo-mne)
! Jezeli orbita przecina biegun, to nie licze calek phi
    if (mpo/=one) then
      n = mpo/(one-mpo)
! Licze dodatkowy wyraz, zeby miec pewnosc ze calki phi beda zbiezne jak bedziemy przechodzic przez
! biegun
      fac2 = one-mpo
    else
! Ustawiam n takie, zeby moc wykryc zmiane
      n = -one
      fac2 = one
    endif
    phi0 = dacos(mup/muplus)
    phif = dacos(muk/muplus)
    call calkieliptycznelegendra(phi0,phif,m1,-n,elcaltp,rfdcmup,rfdcmuk,rfdcmutp,rd2r,rj3r,elle0,ellef,ellec&
   &,ellpi0,ellpif,ellpic)
    if (elcaltp) then
      ctmupplus = fac*elle0
      ctmutp = fac*ellec
      cphimupplus = ellpi0/(fac*fac2)
      cphimutp = ellpic/(fac*fac2)
    endif
    ctmukplus = fac*ellef
    cphimukplus = ellpif/(fac*fac2)
    return   
    end subroutine


    subroutine calkaimuniesymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfr1pplus,rfr1)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Procedura oblicza czesc calki I_{mu} od mup do mu turning pointu w przypadku niesymetrycznym
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,mup,muplus,cmu1,cmu3,rfr1pplus,rfr1
    double precision  ::  m1,fac,phi0,s0,ak,q0     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------! 
       
    fac = dabs(a)*muplus
    m1 = mne/mpo
    s0 = dsqrt((mpo-mup**two)/(mpo-mne))
    phi0 = dasin(s0)
    ak = dsqrt(one-m1)
    q0 = (one-s0*ak)*(one+s0*ak)
    rfr1pplus = rf(one-s0**two,q0,one)
    rfr1 = rf(zero,m1,one)
    cmu3 = two*rfr1
    cmu1 = s0*rfr1pplus
    if (dtan(phi0)<zero) then
      cmu1 = cmu3-cmu1
    endif
    cmu1 = cmu1/fac
    cmu3 = cmu3/(two*fac)
    return
    end subroutine


    subroutine calkaimuniesymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfr1kplus)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Procedura oblicza czesc calki I_{mu} od mu turning point do muk w przypadku niesymetrycznym
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,muk,muplus,cmu2,cmu3,rfr1kplus
    double precision  ::  m1,fac,ak,sf,qf!,phif   ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------! 
    cmu3 = cmu3
!    phif = phif
    fac = dabs(a)*muplus
    m1 = mne/mpo
    sf = dsqrt((mpo-muk**two)/(mpo-mne))
    ak = dsqrt(one-m1)
    qf = (one-sf*ak)*(one+sf*ak)
    rfr1kplus = rf(one-sf*sf,qf,one)
    cmu2 = sf*rfr1kplus/fac
    cmu2 = cmu2/fac
    return
    end subroutine


    subroutine calkaimusymetryczna(mne,mpo,mup,muplus,cmu1,cmu3,rfr1pplus,rfr1)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Procedura oblicza czesc calki I_{mu} od mup do mu turning point w przypadku symetrycznym
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,mup,muplus,cmu1,cmu3,rfr1pplus,rfr1
    double precision  ::  m1,fac,phi0,s0,ak,q0     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    fac = dabs(a)*dsqrt(mpo-mne)
    m1 = -mne/(mpo-mne)
    phi0 = dacos(mup/muplus)
    s0 = dsin(phi0)
    ak = dsqrt(one-m1)
    q0 = (one-s0*ak)*(one+s0*ak)
    rfr1pplus = rf(one-s0**two,q0,one)
    rfr1 = rf(zero,m1,one)
    cmu3 = two*rfr1
    cmu1 = s0*rfr1pplus
    if (dtan(phi0)<zero) then
      cmu1 = cmu3-cmu1
    endif
    cmu1 = cmu1/fac
    cmu3 = cmu3/fac
    return
    end subroutine


    subroutine calkaimusymetrycznak(mne,mpo,muk,muplus,cmu2,cmu3,rfr1kplus)
    use const
    use funkcjeCarlsona
    implicit none
! Cel procedury:
! Procedura oblicza czesc calki I_{mu} od mu turning point do muk w przypadku symetrycznym
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!---------------Opis zmiennych znajduje sie w pliku zmienne_i_parametry------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  mne,mpo,muk,muplus,cmu2,cmu3,rfr1kplus
    double precision  ::  m1,fac,ak,sf,qf,phif     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    fac = dabs(a)*dsqrt(mpo-mne)
    m1 = -mne/(mpo-mne)
    phif = dacos(muk/muplus)
    sf = dsin(phif)
    ak = dsqrt(one-m1)
    qf = (one-sf*ak)*(one+sf*ak)
    rfr1kplus = rf(one-sf**two,qf,one)
    cmu2 = sf*rfr1kplus
    if (dtan(phif)<zero) then
      cmu2 = cmu3*fac-cmu2
    endif
    cmu2 = cmu2/fac
    return
    end subroutine


    double precision function tuszescian(up,uk,u1,u2,l)
    use const
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki tu w przypadku gdy mamy trzy pierwiastki (2 rowne, u2=u3)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  up,uk,u1,u2,l
    double precision  ::  upl,umi,sf,ss,su21,su1,sup,summ     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    upl = one/(one+dsqrt(one-a**two))
    umi = one/(one-dsqrt(one-a**two))
    sf = dsqrt(uk-u1)
    ss = dsqrt(up-u1)
    su21 = dsqrt(u2-u1)
    su1 = dsqrt(-u1)
    sup = dsqrt(one-u1/upl)
    summ = dsqrt(one-u1/umi)
    tuszescian = sf/(uk*u1*u2)-ss/(up*u1*u2)-((-u2-two*u1*(one+u2*umi+u2/upl))*dlog((sf+su1)*(su1-ss)/&
   &((-sf+su1)*(ss+su1))))/(two*(-u1)**(try/two)*u2*u2)+((one-two*a*l*u2*u2+a*a*u2**two*(one+two*u2))*&
   &dlog((su21+sf)*(su21-ss)/((su21-sf)*(su21+ss))))/(u2*u2*su21*(u2/umi-one)*(u2/upl-one))+&
   &((-two*a*l+a*a*(two+one/umi)+one/umi**try)*dsqrt(one/umi)*dlog((sup+sf*sqrt(one/umi))*(summ-ss*&
   &dsqrt(one/umi))/((summ-sf*dsqrt(one/umi))*(summ+ss*dsqrt(one/umi)))))/(summ*(u2/umi-one)*(one/umi-one/upl))&
   &+((-two*a*l+a*a*(two+one/upl)+one/upl**try)*dsqrt(one/upl)*dlog((sup+sf*dsqrt(one/upl))*&
   &(sup-ss*dsqrt(one/upl))/((sup-sf*dsqrt(one/upl))*(sup+ss*dsqrt(one/upl))))/((one/upl-one/umi)*sup*(u2/upl-one)))
    return
    end function


    double precision function phiuszescian(u,u1,u2,l)
    use const
    implicit none
! Cel procedury:
! Funkcja oblicza wartosc calki phiu w przypadku gdy mamy trzy pierwiastki (2 rowne, u2=u3)
!----------------------------------------------------------------------------------!
!----------------------------Definicje zmiennych-----------------------------------!
!-------------------Opis modulow znajduje sie w pliku moduly-----------------------!
!--------------------Wszystkie ponizsze zmienne sa lokalne-------------------------!
!----------------------------------------------------------------------------------!
    double precision  ::  u,u1,u2,l
    double precision  ::  upl,umi,s,su21,su1,sup,summ     ! zmienne lokalne
!----------------------------------------------------------------------------------!
!------------------------Koniec definicji zmiennych--------------------------------!
!----------------------------------------------------------------------------------!
    upl = one/(one+dsqrt(one-a**two))
    umi = one/(one-dsqrt(one-a**two))
    s = dsqrt(u-u1)
    su21 = dsqrt(u2-u1)
    su1 = dsqrt(-u1)
    sup = dsqrt(one-u1/upl)
    summ = dsqrt(one-u1/umi)
    phiuszescian = (l+two*a*u2-two*l*u2)*dlog((su21+s)/(su21-s))/(su21*(u2/umi-one)*(u2/upl-one))+&
   &(two*(a-l)+l/umi)*dsqrt(one/umi)*dlog((summ+s*dsqrt(one/umi))/(summ-s*dsqrt(one/umi)))/(summ*&
   &(u2/umi-two)*(one/umi-one/upl)+(two*(a-l)+l/upl)*dsqrt(one/upl)*dlog((sup+s*dsqrt(one/upl))/&
   &(sup-s*dsqrt(one/upl)))/(sup*(u2/upl-one)*(one/upl-one/umi)))
    return
    end function
    
  end module calki
