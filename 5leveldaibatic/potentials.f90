subroutine Potentiala(x, V)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::V
  real*8 A, B, C, D

  A = 0.01
  B = 1.6
  C = 0.005
  D = 1.0

  if(x>0) then 
      V(1,1) = 0.01*(1- exp(-1.6*x))
  else
      V(1,1) = -0.01*(1- exp(1.6*x))
  end if
  V(1,2) = 0.005*(exp(-(x**2)))
  V(2,1)=V(1,2)
  V(2,2) = -V(1,1)

end subroutine Potentiala 

subroutine Potentialb(x, V)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::V
  real*8 A, B, C, D, E0
  A= 0.10
  B=0.28
  C=0.015
  D=0.06
  E0 = 0.05
    
  V(1,1) = 0.0
  V(2,2) = -A*(exp(-B*x*x))+E0
    
  V(1,2) = C*(exp(-(D*x**2)))
  V(2,1)=V(1,2)
    

end subroutine Potentialb 

subroutine Potentialc(x, V)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::V
  real*8 A, B,C
  A = 0.0006
  B = 0.10
  C = 0.90
  if(x<0) then 
      V(1,2) = B*(exp(C*x))
  else
      V(1,2) = B*(2- exp(-C*x))
  end if
  V(1,1) = A
  V(2,1)=V(1,2)
  V(2,2) = -V(1,1)

end subroutine Potentialc
 
subroutine Potentialpa(x, Vp)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::Vp
  real*8 A, B, C, D

  A = 0.01
  B = 1.6
  C = 0.005
  D = 1.0
  
  if(x>0) then 
      Vp(1,1) = 0.01* exp(-1.6*x)*1.6
  else
      Vp(1,1) = 0.01*(1.6*exp(1.6*x))
  end if
  Vp(1,2) = -0.005*(exp(-(x**2)))*2*x
  Vp(2,1)=Vp(1,2)
  Vp(2,2) = -Vp(1,1)

end subroutine Potentialpa

subroutine Potentialpb(x, V)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::V
  real*8 A, B, C, D, E0
  A= 0.10
  B=0.28
  C=0.015
  D=0.06
  E0 = 0.05
    
  V(1,1) = 0.0
  V(2,2) = 2*B*A*(exp(-B*x*x))*x
   
  V(1,2) = -2*D*C*(exp(-(D*x**2)))*x
  V(2,1)=V(1,2)
    
end subroutine Potentialpb 
 
subroutine Potentialpc(x, V)
  implicit none
  real*8, intent(in) :: x
  complex*16, intent(inout),dimension(2,2) ::V
  real*8 A, B,C
  A = 0.0006
  B = 0.10
  C = 0.90

  if(x<0) then 
      V(1,2) = C*B*(exp(C*x))
  else
      V(1,2) = C*B*(exp(-C*x))
  end if
  V(1,1) = 0
  V(2,1)=V(1,2)
  V(2,2) = -V(1,1)

end subroutine Potentialpc 

subroutine modelx(x,V)
  implicit none
  
  real*8, intent(in) :: x
  complex*16, intent(inout) :: V(3,3)

  real*8 A, B, C, D

  D = 7.0

  A = 0.03
  B = 1.6
  C = 0.005

  V(1,1) = A*(tanh(B*x) + tanh(B*(x+D)))
  V(2,2) = -A*(tanh(B*x) + tanh(B*(x-D)))
  V(3,3) = -A*(tanh(B*(x+D)) - tanh(B*(x-D)))
  V(1,2) = C*exp(-(x*x))
  V(1,3) = C*exp(-(x+D)*(x+D))
  V(2,3) = C*exp(-(x-D)*(x-D))
  
  V(2,1) = V(1,2)
  V(3,1) = V(1,3)
  V(3,2) = V(2,3)

end subroutine modelx

subroutine modelxp(x,V)
  implicit none
  
  real*8, intent(in) :: x
  complex*16, intent(inout) :: V(3,3)

  real*8 A, B, C, D

  D = 7.0

  A = 0.03
  B = 1.6
  C = 0.005
  
  V(1,1) = A*B*(1/((cosh(B*x))**2) + 1/((cosh(B*(x+D)))**2))
  V(2,2) = -A*B*(1/((cosh(B*x))**2) + 1/((cosh(B*(x-D)))**2))
  V(3,3) = -A*B*(1/((cosh(B*(x+D)))**2) - 1/((cosh(B*(x-D)))**2))

  V(1,2) = -2*C*x*exp(-x*x)
  V(1,3) = -2*(x+D)*C*exp(-(x+D)*(x+D))
  V(2,3) = -2*(x-D)*C*exp(-(x-D)*(x-D))
  
  V(2,1) = V(1,2)
  V(3,1) = V(1,3)
  V(3,2) = V(2,3)
end subroutine


subroutine Holstein(q,V,mass,omega)
  implicit none
  
  real*8, intent(in) :: q(5), mass, omega
  complex*16, intent(inout) :: V(5,5)
  real*8 Vc, g, classpe, au2rcm
  integer it
  au2rcm=219474.63067d0
  V = 0
  Vc = 50/au2rcm
  g = 3091.8/au2rcm
  classpe = 0
  do it = 1,5
      classpe = classpe + 0.5*mass*omega**2*q(it)**2
      V(it,it) = q(it)*g
      if ((it-1) >0) V(it,it-1) = Vc
      if ((it+1) <6) V(it,it+1) = Vc
  end do

  do it = 1,5
     V(it,it) = V(it,it) + classpe
  end do

end subroutine


subroutine HolsteinGrad(q,V,mass,omega)
  implicit none
  
  real*8, intent(in) :: q(5), mass, omega
  complex*16, intent(inout) :: V(5,5,5)
  real*8 Vc, g, classpe, au2rcm
  integer it
  V = 0
  au2rcm=219474.63067d0
  Vc = 50/au2rcm
  g = 3091.8/au2rcm
  classpe = 0
  do it = 1,5
      V(it,it,it) = g + mass*omega**2*q(it)
  end do

end subroutine









