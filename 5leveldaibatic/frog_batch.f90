program frog_batch
  implicit none
  real*8 mass,tmax,timstp,p_initial,p_final,p_step,kepara,xranf
  real*8 time,scr,onedq
  real*8, allocatable :: proba(:,:), gamma_collapse(:),gamma_reset(:),p(:),q(:)

  complex*16, allocatable :: V(:,:), Vd(:), Vp(:,:,:), Vpd(:,:,:), vl(:,:)
  complex*16, allocatable :: densmat(:,:),  delR(:,:,:),nacl(:)
  complex*16, allocatable :: delP(:,:,:), c_matrix(:,:), nacv(:,:,:)
  complex*16, allocatable :: b_matrix(:,:), b_matrix_temp(:,:), c_matrix_temp(:,:)
  
  integer active,it,is,irun,nruns,activeold,nstates,ndim,sh_seed, counter
  integer, allocatable :: stat(:),values(:)
  real*8 ran2
  
  logical terminate, testmode, lexci

  testmode = .false.
  p_initial = 0.0
  p_final = 0.1
  p_step = 0.5
  tmax = 1000
  timstp = 1.0
  nruns = 1
  mass = 1836.0

  nstates = 5
  ndim = 5

  terminate = .false. 


  allocate(V(nstates,nstates))
  allocate(vl(nstates,nstates))
  allocate(Vp(nstates,nstates,ndim))
  allocate(Vd(nstates))
  allocate(Vpd(nstates,nstates,ndim))
  allocate(nacv(nstates,nstates,ndim))
  allocate(nacl((nstates*(nstates+1))/2))
  allocate(p(ndim))
  allocate(q(ndim))
  allocate(b_matrix_temp(nstates,nstates))
  allocate(values(8))

  allocate(densmat(nstates,nstates))
  allocate(delR(nstates,nstates,ndim))
  allocate(delP(nstates,nstates,ndim))
  
  allocate(proba(nstates,nstates),gamma_collapse(nstates))
  allocate(gamma_reset(nstates))
 
  allocate(stat(4))
  q=0
  p = 1d-3
  testmode = .false.
  if(testmode) then

  q(1) = -0.2
  q(2) = 0.3
  call electronic_evaluate(mass,p,q,V,Vp,Vd,Vpd,nacv,ndim,nstates,active,vl,.false.)

      scr = -1.0

  end if
  !write(*,'(a)')  '#p,     trans_up,       refl_up,    trans_low,     refl_low'

  do while(p_initial < p_final) !initial momentum loop
      stat = 0
      irun = 0 
      do while(irun<nruns) !indivudual trajectory loop
         irun = irun + 1
          call initialize(p,q,densmat,active,ndim,p_initial,nstates)
          terminate = .false.
          active = 1
          activeold = active
          time = 0
          call date_and_time(values=values) 
          sh_seed = (values(8)*1000 + values(7))*60 + 1
          do it = 1,30
              xranf=ran2(sh_seed)    
          end do
          lexci = .false.
          do while((time<tmax).and.(.not.terminate)) 
              if(active.eq.2) lexci = .true.
              call electronic_evaluate(mass,p,q,V,Vp,Vd,Vpd,nacv,ndim,nstates,active,vl,.false.)

              call classical_propagate(p,q,mass,Vpd,active,&
              &    activeold,nacv,kepara,timstp,Vd,nstates,ndim,nacl)

              call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates,.false.)
              
!              write(*,'(5e18.10)') (real(densmat(it,it)), it = 1,5)
              call aush_propagate(densmat,nstates,ndim,delR,delP,Vpd, &
              &    gamma_collapse,gamma_reset,mass,nacl,Vd,timstp,active)
              counter = 0
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
                  proba(active,is)=real(densmat(active,is)*b_matrix_temp(is,active))
                  proba(active,is)=(2.0)*timstp*proba(active,is)/real(densmat(active,active))
              end do
              write(*,'(5e18.10,i5)')  time, q(1), p(1), densmat(1,1),active
!              write(*,'(a)') 'proba'
!              if(real(densmat(2,2))>0.00001) write(*,*) 'hello'

              call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates,.false.)
              xranf = ran2(sh_seed)
 
              call select_newstate(active,xranf,proba,nstates)
              if (active.ne.activeold)  then 
                  call au_hop(delR,delP,nstates,ndim,active)
              else
                  call au_collapse(gamma_reset,gamma_collapse,densmat,delR,&
                       & delP,nstates,active)
              end if
              time = time + timstp
!              if ((abs(q(1))) > 11.0) then
!                  terminate = .true.
!              end if
          end do! time
          if ((q(1).gt.10.0).and.(active.eq.2)) stat(1) = stat(1) + 1
          if ((q(1).lt.-10.0).and.(active.eq.2)) stat(2) = stat(2) + 1
          if ((q(1).gt.10.0).and.(active.eq.1)) stat(3) = stat(3) + 1
          if ((q(1).lt.-10.0).and.(active.eq.1)) stat(4) = stat(4) + 1
      end do ! run
      write(*,'(e15.6,4i12)') p_initial, stat
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

subroutine initialize(p,q,densmat,active,ndim,p_initial,nstates)
  implicit none

  integer, intent(inout) :: active,nstates,ndim
  real*8, intent(inout) :: p(ndim),q(ndim),p_initial
  complex*16, intent(inout) :: densmat(nstates,nstates)

  active = 1
  p = p_initial
  q = 0 
  q = -0.00
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
