program frog_batch
  implicit none
  real*8 mass,tmax,timstp,p_initial,p_final,p_step,kepara,xranf
  real*8 time,scr,onedq
  real*8, allocatable :: proba(:,:), gamma_collapse(:),gamma_reset(:),p(:),q(:)

  complex*16, allocatable :: V(:,:), Vd(:), Vp(:,:,:), Vpd(:,:,:)
  complex*16, allocatable :: densmat(:,:),  delR(:,:,:),nacl(:),coef(:),densmatnew(:,:)
  complex*16, allocatable :: delP(:,:,:),c_matrix(:,:),nacv(:,:,:),eigvec(:,:),densadia(:,:)
  complex*16, allocatable :: b_matrix(:,:), b_matrix_temp(:,:), c_matrix_temp(:,:)
   
  complex*16 alpha,beta 
  integer active,it,is,irun,nruns,activeold,nstates,ndim,sh_seed, counter
  integer timecounter,pitch
  integer, allocatable :: stat(:),values(:)
  real*8 au2rcm,pi,speedoflight,au2aa,omega
  real*8 ran2, kine, pote, totale
  
  logical terminate, testmode, lafssh, ldiabatic, ldebug,lcol, ldebug3,debugadiab

  testmode = .false.
  p_initial = 0.3
  p_final = 0.5
  p_step = 0.5
  tmax = 8d3
  timstp = 5d-2
  nruns = 1
  mass = 1836.0

  nstates = 5
  ndim = 5
  
  ldebug3 = .false.
  ldebug = .false.
  terminate = .false. 

  lafssh = .false.

  ldiabatic = .true.

  alpha=(1d0,0d0)
  beta=(0d0,0d0)
  pitch=1000
  debugadiab =.false.
  lcol = .false.

  allocate(V(nstates,nstates))
  allocate(Vp(nstates,nstates,ndim))
  allocate(Vd(nstates))
  allocate(coef(nstates))
  allocate(Vpd(nstates,nstates,ndim))
  allocate(eigvec(nstates,nstates))
  allocate(nacv(nstates,nstates,ndim))
  allocate(nacl((nstates*(nstates+1))/2))
  allocate(p(ndim))
  allocate(q(ndim))
  allocate(b_matrix_temp(nstates,nstates))
  allocate(values(8))

  allocate(densmat(nstates,nstates))
  allocate(densmatnew(nstates,nstates))
  allocate(delR(nstates,nstates,ndim))
  allocate(delP(nstates,nstates,ndim))
  
  allocate(proba(nstates,nstates),gamma_collapse(nstates))
  allocate(gamma_reset(nstates))
 
  allocate(stat(4))

  au2rcm=219474.63067d0
  pi = 3.14159265359d0
  speedoflight = 137.0359990740d0

  au2rcm=219474.63067d0
  au2aa = .52917721067d0
  omega = 200.0 !cm-1
  omega = ( (omega/1d8) *au2aa )*speedoflight*2*pi

  testmode = .false.
  if(testmode) then
      q(1) = -0.2
      q(2) = 0.3
      call electronic_evaluate(mass,omega,p,q,V,Vp,Vd,Vpd,nacv,ndim,&
           & nstates,active,eigvec,.false.)
      scr = -1.0
  end if
  !write(*,'(a)')  '#p,     trans_up,       refl_up,    trans_low,     refl_low'

  !! If ldiabatic then denmat is athe diabatic density matrix, that is real time
  !! propagated. densadia is used for augmented surface hoppping and the hopping
  !! probabilities.

  if(ldiabatic) then
      allocate(densadia(nstates,nstates))
  end if

  do while(p_initial < p_final) !initial momentum loop
      stat = 0
      irun = 0 
      do while(irun<nruns) !indivudual trajectory loop
         irun = irun + 1
          write(*,*) '# run  ',irun
          write(*,*) ''
          call initialize(p,q,densmat,active,ndim,p_initial,nstates,mass,omega)
          coef = 0d0
          if(ldiabatic) coef(1) =1d0
          time = 0
          call date_and_time(values=values) 
          sh_seed = (values(8)*1000 + values(7))*60 + 1
          do it = 1,30
              xranf=ran2(sh_seed)    
          end do

          call electronic_evaluate(mass,omega,p,q,V,Vp,Vd,Vpd,nacv,ndim,&
          &    nstates,active,eigvec,.false.)
          
          if(ldiabatic) then
              do it = 1,nstates
                  if(real(eigvec(it,1)) > real(eigvec(active,1))) active = it
              end do
          end if
          write(*,'(a)')  '#time, real(densmat(1,1)), real(densmat(2,2)), real(densmat(3,3)), proba,active'

      !    write(*,'(5e18.10,i5)')  time, real(densmat(1,1)), real(densmat(2,2)), real(densmat(3,3))&
       !   &,gamma_collapse(1), active
          timecounter = 0
          terminate = .false.
          activeold = active
          do while((time<tmax).and.(.not.terminate)) 
              timecounter = timecounter + 1
              call electronic_evaluate(mass,omega,p,q,V,Vp,Vd,Vpd,nacv,ndim,&
              &    nstates,active,eigvec,.false.)

              call classical_propagate(p,q,mass,Vpd,active,&
              &activeold,nacv,kepara,timstp,Vd,nstates,ndim,nacl,kine,pote,ldebug)
              

              totale = kine+pote
                  call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates,ldebug,ldiabatic,eigvec,V)
!                  call el_pro_dia(densmat,coef,(timstp/2.0),nstates,ndim,q,eigvec)
!              write(*,'(5e18.10)') (real(densmat(it,it)), it = 1,5)
              counter = 0
              if (ldiabatic) then
                  
                  call zgemm('N','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
                       densmat,nstates,beta,densmatnew,nstates)
             
                  call zgemm('N','C',nstates,nstates,nstates,alpha,densmatnew,nstates,&
                       eigvec,nstates,beta,densadia,nstates)
                  if((mod(timecounter,pitch).eq.0).or.(lcol)) then
                      
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,'(a)') 'eigvec:'
                      if(ldebug3)    write(*,'(10e18.10)') eigvec
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,'(a)') 'densmat after in diabatic:'
                      if(ldebug3)    write(*,'(10e18.10)') densmat
                      if(ldebug3)    write(*,'(a)') 'densmat after in adiabatic:'
                      if(ldebug3)    write(*,'(10e18.10)') densadia
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                      if(ldebug3)    write(*,*) ''
                  end if
              end if

              if(lafssh) then 
                  if(ldiabatic) then 
                      call aush_propagate(densadia,p,nstates,ndim,delR,delP,Vpd, &
                         &    gamma_collapse,gamma_reset,mass,nacl,Vd,timstp,active)
                  else
                      call aush_propagate(densmat,p,nstates,ndim,delR,delP,Vpd, &
                         &    gamma_collapse,gamma_reset,mass,nacl,Vd,timstp,active)
                  end if    
              end if

              do is = 1,nstates
                  do it = is,nstates
                      counter = counter + 1 
                      b_matrix_temp(is,it) = -nacl(counter)
                      b_matrix_temp(it,is) = nacl(counter)
                  end do
              end do
              do is =1,nstates
                  b_matrix_temp(is,is) = b_matrix_temp(is,is)+(0,1)*Vd(is)
              end do
!              write(*,*) 'population 11'
!              write(*,'(a,4e18.10)') 'pop', q, real(densmat(1,1)), real(densmat(2,2)), abs(densmat(1,2))
              proba = 0
              do is=1,nstates
                  if(active.eq.is) cycle
                  if(ldiabatic) then
                      proba(active,is)=real(densadia(active,is)*b_matrix_temp(is,active))
                      proba(active,is)=(2.0)*timstp*proba(active,is)/real(densadia(active,active))
                  else
                      proba(active,is)=real(densmat(active,is)*b_matrix_temp(is,active))
                      proba(active,is)=(2.0)*timstp*proba(active,is)/real(densmat(active,active))
                  end if
              end do
!              write(*,'(5e18.10,i5)')  time, q(1), p(1), kine, pote,active
!              write(*,'(7e18.10)')  time, q(1), p(1), coef(1) ,coef(2)




              !! Main printing for data

              if(((mod(timecounter-1,pitch).eq.0).or.(lcol)).and.(.not.debugadiab)) then

                  write(*,'(5e18.10,2i5)')  time, real(densmat(1,1)), real(densmat(2,2)), real(densmat(3,3))&
                  &,gamma_collapse(1), active, lcol
                  lcol = .false.
              end if

!              if(debugadiab) then
!                  write(*,'(6e16.8,i5)')  time, real(densadia(1,1)), real(densadia(2,2)), real(densadia(3,3))&
!                  &, real(densadia(4,4)), real(densadia(5,5)), active
!              end if

              call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates,.false.,ldiabatic,eigvec,V)
              xranf = ran2(sh_seed)
              call select_newstate(active,xranf,proba,nstates)
             
              if(lafssh) then 
                  if (active.ne.activeold)  then 
                      call au_hop(delR,delP,nstates,ndim,active)
                  else
                      if(ldiabatic) then  
                          call au_collapse(gamma_reset,gamma_collapse,densadia,delR,&
                           & delP,nstates,active,lcol,ldebug3)
                          call zgemm('C','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
                               densadia,nstates,beta,densmatnew,nstates)
             
                          call zgemm('N','N',nstates,nstates,nstates,alpha,densmatnew,nstates,&
                               eigvec,nstates,beta,densmat,nstates)
                          if(lcol.and.ldebug3) then
                                 write(*,'(a)') 'densmat after in diabatic:'
                                 write(*,'(10e18.10)') densmat
                                 write(*,'(a)') 'eigvec'
                                 write(*,'(10e18.10)') eigvec
                                 write(*,*) ''
                                 write(*,*) ''
                                 write(*,*) ''
                                 write(*,*) ''
                          end if     
                      else
                          call au_collapse(gamma_reset,gamma_collapse,densmat,delR,&
                           & delP,nstates,active,lcol)
                      end if
                  end if
              else
                  call zgemm('C','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
                       densmat,nstates,beta,densmatnew,nstates)
                  call zgemm('N','N',nstates,nstates,nstates,alpha,densmatnew,nstates,&
                       eigvec,nstates,beta,densmat,nstates)
              end if    
               
              time = time + timstp
          end do! time
          if ((q(1).gt.10.0).and.(active.eq.2)) stat(1) = stat(1) + 1
          if ((q(1).lt.-10.0).and.(active.eq.2)) stat(2) = stat(2) + 1
          if ((q(1).gt.10.0).and.(active.eq.1)) stat(3) = stat(3) + 1
          if ((q(1).lt.-10.0).and.(active.eq.1)) stat(4) = stat(4) + 1
      end do ! run
!      write(*,'(e15.6,4i12)') p_initial, stat
      p_initial = p_initial + p_step
  end do! momentum

end program

subroutine select_newstate(ifstate, xranf, proba, &
& nstates)
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!! Selects the active state for the next MD step based on
!! hopping probabilities and the generated random number.
!! input/output  :: integer ifstate 
!!               [IN]   initial active state
!!               [OUT]  final active state  
!! input         :: real(kdp) xranf
!!               [IN]   Random number to be checked against probability
!! input         :: integer nstate
!!               [IN]   Total number of states
!! input         :: real(kdp) proba(nstates,nstates)
!!               [IN] Matrix of probability such that
!!                    proba(i,j) is the probability of
!!                    hopping from state i to state j.
!! Note :: Indices in proba starts from 1, state numbering
!! starts from 0. (Sorry)
!!------------------------------------------------------------------------------------
!!------------------------------------------------------------------------------------
!! Date : Jan 25, 2019.
!! Author: Saswata Roy
!!------------------------------------------------------------------------------------

implicit none

integer, intent(inout):: ifstate, nstates
real*8, intent(inout):: xranf, proba(nstates,nstates)
!Local variables
real*8 Sj 
integer ii, exopt, is,it

 Sj=0
 exopt = ifstate
 do ii = 1,nstates
     Sj=Sj+proba(exopt,ii)
     if(xranf.lt.Sj) then
         ifstate = ii
         exit
     end if
 end do
end subroutine select_newstate

subroutine initialize(p,q,densmat,active,ndim,p_initial,nstates,mass,omega)
  implicit none

  integer, intent(inout) :: active,nstates,ndim
  real*8, intent(inout) :: p(ndim),q(ndim),p_initial,mass,omega
  complex*16, intent(inout) :: densmat(nstates,nstates)
  real*8 ran2, u1,u2,sigmap,sigmax,kB,temp,pi,z1,z2
  integer idimension
  integer, allocatable :: values(:)

  

  allocate(values(8))

  active = 1

  pi = 3.14159265359d0
  temp = 575.5d0 !K
  kB = 3.166811429d-6
  sigmax = sqrt(kB*temp/(mass*omega**2))
  sigmap = sqrt(mass*kB*temp)
  call date_and_time(values=values)

  do idimension = 1,ndim
      u1 = ran2(values(8))
      u2 = ran2(values(8))
      z1 = sqrt(-2.0 * log(u1)) * cos(2*pi * u2)
      z2 = sqrt(-2.0 * log(u1)) * sin(2*pi * u2)
      p(idimension) = z1 * sigmap
      q(idimension) = z1 * sigmax
  end do

  densmat = 0d0
  densmat(active,active) = (1d0,0d0)

end subroutine initialize
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
