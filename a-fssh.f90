!!Program to simulate the second example of Tully's paper.
  
PROGRAM trajau
  implicit none
  Real*8 p, q, pdot, rdot, timestep, k1p, k2p, k3p, k4p
  real*8 k1q, k2q, k3q, k4q,random,random2
  Real*8, dimension(2,2) :: V,d
  complex*16, dimension(2) :: coef, coefdot
  complex*16, dimension(2,2):: a, adot, delR, delP, delRdot, delPdot
  Real*8 t,m, runtime, gamma_c, gamma_r, temp3
  real*8, dimension(4):: plist
  real*8, parameter:: hbar=1 !1.05457173d-34
  Integer state, counter, flag, flag2, loopcount, run, idex, sh_seed, rnd_size,maxrun
  Real*8 e1, e2, d12, poriginal, b12, b21, probability
  Real*8, dimension(2):: l
  Real*8, dimension(2,2) :: CC, CC1, Vp, Vpd, Vd
  Integer, dimension(6):: stat
  Integer, allocatable :: seedy(:)
  Real*8 ran2
  call random_seed(size=rnd_size)
  !plist= (/4.4,4.6,7.7,7.9,8.1,8.3,8.7,8.9,9.2,9.7,10.3,10.7,11.3/)
  plist= (/7.9,8.1,8.3,8.7/)
  allocate(seedy(rnd_size))
  m=2000
  flag=0
  flag2=0
  seedy=-1234
  sh_seed = -1234
  open(90, FILE='statistics.k10-deco')
  open(8, File = 'gamma_c_1h0_nodeco')
  open(4, File = 'rnd_1')
  open(88, File = 'collapse_')
  open(65, File = 'Force')
  open(69, File = 'print_dRdP')
  timestep=1
  idex=1
  poriginal =10.0
  maxrun = 200
  write(90,*) 'poriginal, counter, lowtrans, uptrans, lowref, upref, coll, reset'
  do while(poriginal.le.10.0)

      write(*,*) 'poriginal', poriginal
      counter = 0
      run=0
      stat = 0
      do while(run<maxrun) 
          p=poriginal
          runtime =1
          !runtime=(40.0/poriginal)*m 
          run=run+1
          t = 0
          state = 1
          coef = (0,0)
          a=(0,0)
          coef(1)=(1,0)
          a(1,1)=(1,0)
          flag = 0
          delP = 0
          delR = 0 
          delPdot = 0
          delRdot = 0

          q=-10.0 ! Initial condition set up, each of them making a start from -7
          loopcount = 0

 
          call Potential(q,V)!The diabatic potential matrix is evaluated
          e1= V(1,1)+V(2,2)
          e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
          l(2) = (e1+Sqrt(e2))/2 !The two eigen values are evaluated 
          l(1) = (e1-Sqrt(e2))/2
          Vd=0.0
          Vd(1,1) = l(1)
          Vd(2,2)=l(2)
          call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
          call inverse2(CC,CC1)
          call PotentialP(q, Vp)
          Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
 
          d12 = Vpd(2,1)/(l(2)-l(1))
       
          d=0.0
          d(2,1) = -d12
          d(1,2) = d12 
          rdot =p/m
          pdot = -Vpd(state,state)
 
          call getcoefDot(coefdot, coef,Vd, rdot, d)
          call getADot(adot, a, Vd, rdot, d)

          do while (abs(q)<10.01) !This is the trajectory loop
              if (t.ge.5000000) exit
              loopcount=loopcount+1
              call propagate(p, q, coef, a, timestep, state, m)
              call Potential(q,V)!The diabatic potential matrix is evaluated

              e1= V(1,1)+V(2,2)
              e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
              l(2) = 0.5*(e1+Sqrt(e2)) !The two eigen values are evaluated 
              l(1) = 0.5*(e1-Sqrt(e2))

              Vd(1,1)=l(1)
              Vd(2,2)=l(2)
                   
              call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
              call inverse2(CC,CC1)
              call PotentialP(q, Vp) !The derivative of the potential is calculated in the adiabatic state
 
              Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
              d12 = Vpd(2,1)/(l(2)-l(1))
 
              d(2,1) = -d12
              d(1,2) = d12 
              call au_propagate(p,q,coef, timestep, state,m, delR,delP, Vd, d, Vpd) 
!
             gamma_r = -timestep*real((-Vpd(state,state)+Vpd(3-state,3-state))*&
                       &(-delR(state,state) + delR(3-state,3-state))/2.0)
             gamma_c = -gamma_r 
             gamma_r = 0.0

             gamma_c = sign(gamma_c, Real( (delR(1,1)-delR(2,2))*(delP(1,1)-delP(2,2)) ))
!
!              gamma_r = -timestep*real((-Vpd(state,state)+Vpd(3-state,3-state))*&
!                       &(-delR(state,state) + delR(3-state,3-state))/2.0)
!              gamma_c = -gamma_r - 2.0* timestep*abs(Vpd(1,2)*(delR(3-state,3-state) - delR(state,state)))
              
              b21 = -2.0*Realpart(coef(2)*conjg(coef(1)))*p*d(2,1)/m!*Real(Real(conjg(a(2,1))*rdot*d(2,1)))+  2*Real(Aimag(conjg(a(2,1))*Vd(2,1)))/hbar!These are the real deal
              b12 = -2.0*Realpart(coef(1)*conjg(coef(2)))*p*d(1,2)/m!+ 2*Real(Aimag(conjg(a(1,2))*Vd(1,2)))/hbar !These are the probablities to be calculated
              a(1,1)=coef(1)*conjg(coef(1))
              a(2,2)=coef(2)*conjg(coef(2))
              call init_random_seed()
              call random_number(random)
              random=ran2(sh_seed)    
              seedy = sh_seed
              call random_seed(put=seedy) 
              call random_number(random)
              if(state.eq.1) then
                  probability = timestep*b21/abs(a(1,1))
                  if(probability>random) then
                      flag=1
                  end if
             end if
             if(state.eq.2) then
                 probability = timestep*b12/abs(a(2,2))
                 if(probability>random) then
                     flag=1
                 end if
             end if
 
             if(flag==1) then 
                 if(l(3-state)-l(state)<0.5*((p**2)/m)) then
                     p = sign(sqrt(p**2 +2*m*(l(state)-l(3-state))),p)
                     state = 3-state
                     pdot=-Vpd(state,state)
                     flag2=1
                     temp3 = Real(delP(state,state))
                     delP(1,1) = delP(1,1) - temp3
                     delP(2,2) = delP(2,2) - temp3
                     temp3 = Real(delR(state,state))
                     delR(1,1) = delR(1,1) - temp3
                     delR(2,2) = delR(2,2) - temp3
                 end if
                 flag=0
             end if
             random=ran2(sh_seed)    
             seedy = sh_seed
             call random_seed(put=seedy) 
             call random_number(random)
                   
                      
             if(gamma_c>random) then
                 coef(state) = 1.0
                 coef(3-state) = 0.0
                 stat(5) = stat(5)+1
                 delP=0
                 delR=0
             end if
             if(gamma_r>random) then
                 stat(6) = stat(6)+1
                 delP=0
                 delR=0
             end if
                   
             t=t+timestep
         end do!The time domain
         !write(*,*) 'run done'
         if(q>0) then
           if(state==1) then
               stat(1)=stat(1)+1
           else
               stat(2)=stat(2)+1
           end if
         end if
         if(q<0) then
           if(state ==1) then
               stat(3) = stat(3) + 1
           else
               stat(4) = stat(4) + 1
           end if
         end if


         !write(7,*) 
   
         if (flag2==1) then

                 counter=counter+1
                 flag2=0
         end if
 end do !run
 
 write(90,*) poriginal, counter, stat(1), stat(2), stat(3), stat(4),stat(5),stat(6)



 poriginal=poriginal+2.0
 idex=idex+1

 counter = 0
            end do! the p step

   close(6)

end program

!****************************************************THE MAIN PROGRAM ENDS HERE***************************************************


!*******************************************ALL THAT FOLLOWS ARE THE NECESSARY FUNCTIONS***************************************














subroutine au_propagate(p,q,coef, h, state,m, delR,delP, V, d, Vpd) 
           
    implicit none
    real*8, intent(inout):: p, q, h, m
    integer, intent(inout):: state
    complex*16, intent(inout), dimension(2):: coef
    complex*16, intent(inout), dimension(2,2):: delP, delR
    complex*16, dimension(2,2):: k1p,k2p,k3p,k4p
    complex*16,dimension(2,2) :: k1r,k2r,k3r,k4r
    real*8, dimension(2,2) :: V, d, Vpd
    complex*16,dimension(2,2) :: delRdot, delPdot
    real*8 resolution, minih
    integer ii
    resolution = 100.0
    minih = h/resolution
    !write(69,*) V
    !write(69,*) Vpd
    !write(69,*) d
    !write(69,*) delR
    !write(69,*) delP
    !write(69,*) coef 
   
    do ii = 1,int(resolution)
        call getdelRdelPdot(delR,delP,delRdot,delPdot, V, d, Vpd, coef, m, p, state)
    !    write(69,*) delRdot
    !    write(69,*) delPdot

        k1p = minih*delPdot
        k1r = minih*delRdot

        call getdelRdelPdot((delR+0.5*k1r),(delP+0.5*k1p),delRdot,delPdot, V, d, Vpd, coef, m, p, state)

        k2p = minih*delPdot
        k2r = minih*delRdot

        call getdelRdelPdot((delR+0.5*k2r),(delP+0.5*k2p),delRdot,delPdot, V, d, Vpd, coef, m, p, state)

        k3p = minih*delPdot
        k3r = minih*delRdot
        call getdelRdelPdot((delR+k3r),(delP+k3p),delRdot,delPdot, V, d, Vpd, coef, m, p, state)

        k4p = minih*delPdot
        k4r = minih*delRdot
        


        delP = delP+(k1p + 2.0*k2p + 2.0*k3p +k4p)/6.0
        delR = delR+(k1r + 2.0*k2r + 2.0*k3r +k4r)/6.0
    end do
end subroutine


subroutine propagate(p, q, coef, a, h, state, m)
    implicit none
    real*8, intent(inout):: p, q, h, m
    integer, intent(inout):: state
    complex*16, intent(inout), dimension(2):: coef
    complex*16, intent(inout), dimension(2,2):: a
    complex*16, dimension(2,2):: k1a,k2a,k3a,k4a
    complex*16,dimension(2) :: k1c,k2c,k3c,k4c
    Real*8 k1q,k2q,k3q,k4q,pdot
    Real*8 k1p,k2p,k3p,k4p,rdot
    Real*8 e1, e2, d12
    Real*8, dimension(2):: l
    complex*16, dimension(2):: coefdot
    Real*8, dimension(2,2) :: CC, CC1, Vp, Vpd, Vd, V, d
    complex*16, dimension(2,2)::adot
    call Potential(q,V)!The diabatic potential matrix is evaluated
    d=0
    e1= V(1,1)+V(2,2)
    e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
    l(2) = 0.5*(e1+Sqrt(e2)) !The two eigen values are evaluated 
    l(1) = 0.5*(e1-Sqrt(e2))

    Vd(1,1)=l(1)
    Vd(2,2)=l(2)
                                        
    call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
    call inverse2(CC,CC1)
    call PotentialP(q, Vp) !The derivative of the potential is calculated in the adiabatic state
  
    Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
    d12 = Vpd(2,1)/(l(2)-l(1))
         
    d(2,1) = -d12
    d(1,2) = d12 
    rdot =p/m
    pdot = -Vpd(state,state)
          
    call getcoefDot(coefdot, coef, Vd, rdot, d)
    call getADot(adot, a, Vd, rdot, d)
          
    k1p = h * pdot
    k1q= h*rdot
    k1a =h*adot
    k1c= h*coefdot

    
    call Potential(q+0.5*k1q,V)!The diabatic potential matrix is evaluated
    e1= V(1,1)+V(2,2)
    e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
    l(2) = 0.5*(e1+Sqrt(e2)) !The two eigen values are evaluated 
    l(1) = 0.5*(e1-Sqrt(e2))

    Vd(1,1)=l(1)
    Vd(2,2)=l(2)
                                        
    call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
    call inverse2(CC,CC1)
    call PotentialP(q+0.5*k1q, Vp) !The derivative of the potential is calculated in the adiabatic state
  
    Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
    d12 = Vpd(2,1)/(l(2)-l(1))
         
    d(2,1) = -d12
    d(1,2) = d12 
    rdot =(p+0.5*k1p)/m
    pdot = -Vpd(state,state)
          
    call getcoefDot(coefdot, coef + 0.5*k1c, Vd, rdot, d)
    call getADot(adot, a+0.5*k1a, Vd, rdot, d)
    
    k2p = h * pdot
    k2q= h*rdot
    k2a =h*adot
    k2c= h*coefdot

    call Potential(q+0.5*k2q,V)!The diabatic potential matrix is evaluated
    e1= V(1,1)+V(2,2)
    e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
    l(2) = 0.5*(e1+Sqrt(e2)) !The two eigen values are evaluated 
    l(1) = 0.5*(e1-Sqrt(e2))

    Vd(1,1)=l(1)
    Vd(2,2)=l(2)
                                        
    call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
    call inverse2(CC,CC1)
    call PotentialP(q+0.5*k2q, Vp) !The derivative of the potential is calculated in the adiabatic state
  
    Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
    d12 = Vpd(2,1)/(l(2)-l(1))
         
    d(2,1) = -d12
    d(1,2) = d12 
    rdot =(p+0.5*k2p)/m
    pdot = -Vpd(state,state)
          
    call getcoefDot(coefdot, coef + 0.5*k2c, Vd, rdot, d)
    call getADot(adot, a+0.5*k2a, Vd, rdot, d)
    
    k3p = h * pdot
    k3q= h*rdot
    k3a =h*adot
    k3c= h*coefdot

    call Potential(q+k3q,V)!The diabatic potential matrix is evaluated
    e1= V(1,1)+V(2,2)
    e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
    l(2) = 0.5*(e1+Sqrt(e2)) !The two eigen values are evaluated 
    l(1) = 0.5*(e1-Sqrt(e2))

    Vd(1,1)=l(1)
    Vd(2,2)=l(2)
                                        
    call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed
    call inverse2(CC,CC1)
    call PotentialP(q+k3q, Vp) !The derivative of the potential is calculated in the adiabatic state
  
    Vpd = matmul(CC1,matmul(Vp,CC))!The derivative of the potential in the adiabatic state is calculated
    d12 = Vpd(2,1)/(l(2)-l(1))
         
    d(2,1) = -d12
    d(1,2) = d12 
    rdot =(p+k3p)/m
    pdot = -Vpd(state,state)
          
    call getcoefDot(coefdot, coef + k3c, Vd, rdot, d)
    call getADot(adot, a+k3a, Vd, rdot, d)
    
    k4p = h * pdot
    k4q= h*rdot
    k4a =h*adot
    k4c= h*coefdot
    
     
    p = p + (k1p + 2*k2p+ 2*k3p + k4p)/6.0
    q = q + (k1q + 2*k2q+ 2*k3q + k4q)/6.0
    a = a + (k1a + 2*k2a+ 2*k3a + k4a)/6.0
    coef = coef + (k1c + 2*k2c+ 2*k3c + k4c)/6.0

end subroutine propagate

subroutine getdelRdelPdot(delR,delP,delRdot,delPdot, V, d, Vpd, coef, m, p, state)
    implicit none
    
    complex*16, intent(inout), dimension(2,2):: delRdot, delPdot
    complex*16, intent(in), dimension(2,2):: delR, delP
    complex*16, intent(in), dimension(2):: coef
    real*8, intent(in) :: p,m
    integer, intent(in) :: state
    real*8, intent(in), dimension(2,2)::V,d,Vpd
    integer j,k,l
    complex*16 II, temp
    complex*16 , dimension(2,2) :: rho, delF
    !write(*,*) '0'
    delF = -Vpd
    do j=1,2
       delF(j,j) = delF(j,j) +Vpd(state,state)
       do k  =1,2
           rho(j,k) = conjg(coef(j))*coef(k)
       end do
    end do
     
    !write(8,*) "Inside the loop"
    II=(0,1.0)
    delRdot = -II*(matmul(V,delR) - matmul(delR,V))
    delRdot = delRdot  + delP/m
    delRdot = delRdot - (p/m)*(matmul(d,delR) - matmul(delR,d))

    delPdot = -II*(matmul(V,delP) - matmul(delP,V))
    delPdot = delPdot  + 0.5*(matmul(delF,rho) + matmul(rho,delF))
    delPdot = delPdot - (p/m)*(matmul(d,delP) - matmul(delP,d))

 end subroutine getdelRdelPdot
  


subroutine getcoefDot(coefdot, coef ,V, rdot, d)
    implicit none
    real*8, parameter:: hbar=1 !1.05457173d-34
    complex*16, intent(inout), dimension(2):: coefdot
    complex*16, intent(in), dimension(2):: coef
    real*8, intent(in), dimension(2,2)::V
    real*8, intent(in):: rdot
    real*8, intent(in), dimension(2,2)::d
    integer j,k,l
    complex*16 II

    II=(0,1.0)
    coefdot(1)=coef(1)*V(1,1)/(II *hbar) - coef(2)*rdot*d(1,2) !this is test for the propagator for the coeeficients

    coefdot(2)=coef(2)*V(2,2)/(II *hbar) + coef(1)*rdot*d(1,2)

end subroutine getcoefDot

subroutine getADot(adot, a ,V, rdot, d)
    implicit none
    real*8, parameter:: hbar=1!1.05457173d-34
    complex*16, intent(inout), dimension(2,2):: adot
    complex*16, intent(in), dimension(2,2):: a
    real*8, intent(in), dimension(2,2)::V
    real*8, intent(in):: rdot
    real*8, intent(in), dimension(2,2)::d
    integer j,k,l
    complex*16 II
  
    II=(0,1.0)
    adot=(0,0)
 
    adot(1,1) = -rdot*d(1,2)*(a(2,1)+a(1,2))
    adot(1,2) = -II*a(1,2)*(V(2,2)-V(1,1)) -(a(1,1)-a(2,2))*rdot*d(1,2)

    adot(2,2) = rdot*d(1,2)*(a(2,1)+a(1,2))

    adot(2,1) = II*a(2,1)*(V(2,2)-V(1,1)) +(a(1,1)-a(2,2))*rdot*d(1,2)
end subroutine getAdot



Subroutine Potentialp(x, Vp)
implicit none
    Real*8, intent(IN):: x
    Real*8, dimension(2,2), intent(INOUT):: Vp
    call Potentialpc(x, Vp)
END subroutine Potentialp 

Subroutine Potentialpa(x, Vp)
    IMPLICIT none
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::Vp
!
    if(x>0) then 
      Vp(1,1) = 0.01* exp(-1.6*x)*1.6
    else
      Vp(1,1) = 0.01*(1.6*exp(1.6*x))
    end if
    Vp(1,2) = -0.005*(exp(-(x**2)))*2*x
    Vp(2,1)=Vp(1,2)
    Vp(2,2) = -Vp(1,1)

END subroutine Potentialpa

Subroutine Potentialpb(x, V)
    IMPLICIT none
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B, C, D, E0
    A= 0.10
    B=0.28
    C=0.015
    D=0.06
    E0 = 0.05
!
    
    V(1,1) = 0.0
    
    V(2,2) = 2*B*A*(exp(-B*x*x))*x
    
    V(1,2) = -2*D*C*(exp(-(D*x**2)))*x
    V(2,1)=V(1,2)
    

END subroutine Potentialpb 
 
Subroutine Potentialpc(x, V)
    IMPLICIT none
!
!		declarations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B,C
    A = 0.0006
    B = 0.10
    C = 0.90
!
    if(x<0) then 
      V(1,2) = C*B*(exp(C*x))
    else
      V(1,2) = C*B*(exp(-C*x))
    end if
    V(1,1) = 0
    V(2,1)=V(1,2)
    V(2,2) = -V(1,1)

END subroutine Potentialpc 

Subroutine Potential(x, V)
!    IMPLICIT none
    implicit none
    Real*8, intent(IN):: x
    Real*8, dimension(2,2), intent(INOUT):: V
    call Potentialc(x,V)
END subroutine Potential 

Subroutine Potentiala(x, V)
    IMPLICIT none
!
!		de!larations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
!
    if(x>0) then 
      V(1,1) = 0.01*(1- exp(-1.6*x))
    else
    V(1,1) = -0.01*(1- exp(1.6*x))
    end if
    V(1,2) = 0.005*(exp(-(x**2)))
    V(2,1)=V(1,2)
    V(2,2) = -V(1,1)

END subroutine Potentiala 

Subroutine Potentialb(x, V)
    IMPLICIT none
!
!		de!larations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B, C, D, E0
    A= 0.10
    B=0.28
    C=0.015
    D=0.06
    E0 = 0.05
!
    
    V(1,1) = 0.0
    
    V(2,2) = -A*(exp(-B*x*x))+E0
    
    V(1,2) = C*(exp(-(D*x**2)))
    V(2,1)=V(1,2)
    

END subroutine Potentialb 

Subroutine Potentialc(x, V)
    IMPLICIT none
!
!		declarations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B,C
    A = 0.0006
    B = 0.10
    C = 0.90
!
    if(x<0) then 
      V(1,2) = B*(exp(C*x))
    else
      V(1,2) = B*(2- exp(-C*x))
    end if
    V(1,1) = A
    V(2,1)=V(1,2)
    V(2,2) = -V(1,1)

END subroutine Potentialc 

subroutine createEigenVector(l1, l2, CC, V)
    implicit none
    Real*8 k1,k2,d1,d2
    Real*8, intent(in):: l1,l2
    Real*8, intent(inout), dimension(2,2):: CC, V
    if(abs(v(1,2))>1d-15) then
      k1 = (l1-V(1,1))/(V(1,2))
      k2 = (l2-V(1,1))/(V(1,2))
      d1 = Sqrt(1.0+k1**2)
      d2 = Sqrt(1.0+k2**2)

      CC(1,1) = 1.0/d1
      CC(1,2) = 1.0/d2
      CC(2,1) = k1/d1
      CC(2,2) = k2/d2
    else if (V(1,1)<V(2,2)) then
      CC(1,1) = 1.0
      CC(1,2) = 0.0
      CC(2,1) = 0.0
      CC(2,2) = 1.0
    else
      CC(1,1) = 0.0
      CC(1,2) = 1.0
      CC(2,1) = 1.0
      CC(2,2) = 0.0
    
    end if  
end subroutine createEigenVector

subroutine inverse2(A, B)
    implicit none

    Real*8, dimension(2,2), intent(in)::A
    Real*8, dimension(2,2), intent(inout)::B
    
    B= Transpose(A)    
end subroutine inverse2
 subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed
function ran2(idum)
!-------------------------------------------------------------------------------
! Description:
!-------------------------------------------------------------------------------
!     From Numerical Recipes, chapter 7
!-------------------------------------------------------------------------------
!
!    Long period (> 2 * 10^18) random number generator of L'Ecuyer with 
!    Bays-Durham shuffle and added safeguards. Returns a uniform random
!    deviate between 0.0 and 1.0 exclusive the endpoint values. Call with
!    idum a negative integer to initialize; thereafter, do not alter idum
!    between successive deviates in a sequence. Variable rnmx should
!    approximate the largest floating point value that is less than 1.
!-------------------------------------------------------------------------------
  implicit none

 ! include 'fortrankinds.h'

!-------------------------------------------------------------------------------
! constants
!-------------------------------------------------------------------------------
  integer, parameter :: im1=2147483563, im2=2147483399, imm1=im1-1
  integer, parameter :: ia1=40014, ia2=40692, iq1=53668, iq2=52774
  integer, parameter :: ir1=12211, ir2=3791, ntab=32, ndiv=1+imm1/ntab
  real*8, parameter :: am=1.0/im1, eps=epsilon(1.2E-7)
  real*8, parameter :: rnmx=1.0

!-------------------------------------------------------------------------------
! Input parameters
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Input/output parameters
!-------------------------------------------------------------------------------
  integer, intent(inout) :: idum

!-------------------------------------------------------------------------------
! Output parameters
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Local variables
!-------------------------------------------------------------------------------
  integer :: j, k
  integer, save :: iv(ntab)=0, iy=0, idum2=123456789
  real*8 :: ran2

!-------------------------------------------------------------------------------
  if (idum <= 0) then
     idum = max(-idum,1)
     idum2 = idum
     do j = ntab+8, 1, -1
        k = idum/iq1
        idum = ia1*(idum-k*iq1)-k*ir1
        if (idum < 0) idum = idum+im1
        if (j <= ntab) iv(j) = idum
     end do
     iy = iv(1)
  end if
  k = idum/iq1
  idum = ia1*(idum-k*iq1)-k*ir1
  if (idum < 0) idum = idum+im1
  k = idum2/iq2
  idum2 = ia2*(idum2-k*iq2)-k*ir2
  if (idum2 < 0) idum2 = idum2+im2
  j = 1+iy/ndiv
  iy=iv(j)-idum2
  iv(j)=idum
  if (iy < 1) iy = iy+imm1
  ran2 = min(am*iy,rnmx)

end function ran2
  !open(8, File = 'gibberish.dat')
  !write(7,*) '#t, p, q, Vd(state,state), Vd(3-state,3-state), probability, pop(1)'
  !write(7,*) '#t, p, q, Vd(state,state), real(delR(1,1) - delR(2,2)), gamma_c, random'



   ! temp = delRdot(state,state)
   ! delRdot(1,1) = delR(1,1) - temp
   ! delRdot(2,2) = delR(2,2) - temp
   ! temp = delPdot(state,state)
   ! delPdot(1,1) = delP(1,1) - temp
   ! delPdot(2,2) = delP(2,2) - temp




!    coefdot(1) = -II*coef(2)
!    coefdot(2) = -II*coef(1)

    !write(8,*) 'Write the rho'
    !write(8,*) rho
    !write(8,*) 'dForce'
    !write(8,*) delF


    !write(8,*) '(matmul(delF,rho) + matmul(rho,delF))*0.5'
    !write(8,*)  (matmul(delF,rho) + matmul(rho,delF))*0.5
             !  write(7,*) t, p, q, Vd(state,state), real(delR(1,1) - delR(2,2)), gamma_c, random
             !  write(8,*) 'delR'
             !  write(8,*) delR
             !  write(8,*) '***********************************'
                          !  write(8,*) 'delP'
                          !  write(8,*) delP
                          !  write(8,*) 'gamma_c'
                          !  write(8,*) gamma_c

                          !  write(8,*) '***********************************'
                          !  write(8,*) 'Vpd'
                          !  write(8,*) Vpd

                          !  write(8,*) '***********************************'
                          !  write(8,*) coef
                          !  write(8,*) '***********************************'
