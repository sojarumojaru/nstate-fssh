program frog_batch
  implicit none
  integer ndim
  real*8 mass,tmax,timstp,p_initial,p_final,p_step,kepara,xranf
  real*8 time,scr
  real*8, allocatable :: proba(:,:), gamma_collapse(:),gamma_reset(:)
  real*8, allocatable :: postn(:), momentum(:)

  complex*16, allocatable :: V(:,:), Vd(:), Vp(:,:), Vpd(:,:), vl(:,:) 
  complex*16, allocatable :: densmat(:,:),  delR(:,:),nacl(:)
  complex*16, allocatable :: delP(:,:), c_matrix(:,:), nacv(:,:)
  complex*16, allocatable :: b_matrix(:,:), b_matrix_temp(:,:), c_matrix_temp(:,:)
  
  integer active,it,is,irun,nruns,activeold,nstates
  integer, allocatable :: stat(:)
  
  logical terminate, testmode

  testmode = .false.
  p_initial = 3.0
  p_final = 3.1
  p_step = 0.5
  tmax = 100000
  timstp = 1.0
  nruns = 1
  mass = 2000.0

  nstates = 2
  ndim = 1

  terminate = .false. 


  allocate(V(nstates,nstates))
  allocate(vl(nstates,nstates))
  allocate(Vp(nstates,nstates))
  allocate(Vd(nstates))
  allocate(Vpd(nstates,nstates))
  allocate(nacv(nstates,nstates))
  allocate(nacl((nstates*(nstates+1))/2))
  allocate(postn(ndim))
  allocate(momentum(ndim))
  allocate(b_matrix_temp(nstates,nstates))

  allocate(densmat(nstates,nstates))
  allocate(delR(nstates,nstates))
  allocate(delP(nstates,nstates))
  
  allocate(proba(nstates,nstates),gamma_collapse(nstates))
  allocate(gamma_reset(nstates))
 
  allocate(stat(4))

  testmode = .true.
  if(testmode) then
      scr = -1.0

      write(*,'(a)') "#q, real(Vd(1)), real(Vd(2)), real(nacv(1,2))"
      do q = -10.0,10.0,0.2

          call electronic_evaluate(postn,momentum,V,Vp,Vd,&
               & Vpd,nacv,ndim,nstates,active,vl)
          
          write(*,'(6e18.10)') q, real(Vd(1)), real(Vd(2)), real(nacv(1,2)),&
           & real(Vpd(1,1)), real(Vpd(2,2))
!          write(*,*) ' V '
!          write(*,'(2e18.10)') V
!          write(*,*) ' Vd '
!          write(*,'(2e18.10)') Vd
!          write(*,*) 'Vp'
!          write(*,'(2e18.10)') Vp
!          write(*,*) 'Vpd'
!          write(*,'(2e18.10)') Vpd
!          write(*,*) 'Vl'
!          write(*,'(2e18.10)') Vl
!          write(*,*) ' '
!          write(*,*) ' '
!          write(*,*) ' '
    
     end do
     
  end if

!  do while(p_initial < p_final) !initial momentum loop
!      stat = 0
!      irun = 0 
!      do while(irun<nruns) !indivudual trajectory loop
!         irun = irun + 1
!          call initialize(p,q,densmat,active,p_initial,nstates)
!          terminate = .false.
!          active = 1
!          activeold = active
!          time = 0
!          do while((time<tmax).and.(.not.terminate)) 
!              
!              call electronic_evaluate(p,q,V,Vp,Vd,Vpd,nacv,nstates,active,vl)
!
!              call classical_propagate(p,q,mass,Vpd,active,&
!              &    activeold,nacv,kepara,timstp,Vd,nstates,nacl)
!
!              call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates)
!
!              call aush_propagate(densmat,nstates,delR,delP,Vpd, &
!              &    gamma_collapse,gamma_reset,mass,nacl,Vd,timstp,active)
!              do is = 1,nstates
!                  do it = 1,nstates
!                      b_matrix_temp(is,it) = nacv(is,it)
!                  end do
!              end do
!              do is =1,nstates
!                  b_matrix_temp(is,is) = b_matrix_temp(is,is)+(0,1)*Vd(is)
!              end do
!
!              do is=1,nstates
!                  if(active.eq.is) cycle
!                  proba(active,is)=real(densmat(active,is)*b_matrix_temp(is,active))
!                  proba(active,is)=(2.0)*timstp*proba(active,is)/real(densmat(active,active))
!              end do
!!              if(real(densmat(2,2))>0.00001) write(*,*) 'hello'
!
!              call electronic_propagate(densmat,Vd,nacl,(timstp/2.0),nstates)
!
!              call random_number(xranf)
!
!              call select_newstate(active,xranf,proba,nstates)
!              if (active.ne.activeold)  then 
!                  call au_hop(delR,delP,nstates,active)
!              else
!                  call au_collapse(gamma_reset,gamma_collapse,densmat,delR,&
!                       & delP,nstates,active)
!              end if
!              time = time + timstp
!              if ((abs(q)) > 11.0) then
!                  terminate = .true.
!              end if
!          end do! time
!         if ((q.gt.10.0).and.(active.eq.2)) stat(1) = stat(1) + 1
!          if ((q.lt.-10.0).and.(active.eq.2)) stat(2) = stat(2) + 1
!          if ((q.gt.10.0).and.(active.eq.1)) stat(3) = stat(3) + 1
!          if ((q.lt.-10.0).and.(active.eq.1)) stat(4) = stat(4) + 1
!      end do ! run
!      write(*,'(e15.6,4i12)') p_initial, stat
!      p_initial = p_initial + p_step
!  end do! momentum

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
     Sj=Sj+proba(exopt+1,ii)
     if(xranf.lt.Sj) then
         ifstate = ii-1
         exit
     end if
 end do
end subroutine select_newstate

subroutine initialize(postn,momentum,densmat,active,p_initial,ndim,nstates)
  implicit none

  integer, intent(inout) :: active,nstates,ndim
  real*8, intent(inout) :: momentum(ndim),postn(ndim),p_initial
  complex*16, intent(inout) :: densmat(nstates,nstates)

  active = 1
  p = p_initial
  q = -10.0
  densmat = 0d0
  densmat(active,active) = (1d0,0d0)

end subroutine initialize
