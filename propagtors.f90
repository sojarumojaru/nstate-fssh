subroutine  electronic_evaluate(p,q,V,Vp,Vd,Vpd,d,nstates,active,vl)

  call potential(q,V)
  call potentialp(q,Vp)
    
  call zgeev('v','n',nstates,V,nstates,Vd,VL,nstates,VR,nstates,&
  & WORK,(2*(nstates)),RWORK,info )
  alpha = 1d0
  beta = 0d0
  call zgemm('N','N',nstates,nstates,nstates,alpha,vl,nstates,&
       Vp,nstates,beta,c_matrix_temp,nstates)
  call zgemm('N','C',nstates,nstates,nstates,alpha,c_matrix_temp,nstates,&
       vl,nstates,beta,Vpd,nstates)

  Vd = Vd - Vd(1)  ! setting lowest energy to zero, Vd is excitation
  do i = 1, nstates
      do j = 1,nstates
          if(i.eq.j) cycle
          if (Vd(i) - Vd(j) > thresh) d(i,j) = Vpd(i,j)/(Vd(i) - Vd(j))
      end do
  end do

end subroutine electronic_evaluate

subroutine classical_propagate(p,q,mass,force,active,activeold,d,kepara,timstp,Vd)

  if(activeold.eq.active)
      vel = p/mass
      acc = -force(active,active)/mass
      vel = vel + acc*timstp
      q = q + vel*timstp
      p = mass*vel
      kepara = 0.5*vel*vel*mass
  end if 

  if(active.ne.activeold) then
      vel = p/mass
      acc = -force(activeold)/mass
      vel = vel + acc*timstp*0.5
      kepara = 0.5*vel*vel*mass
      exci = Vd(active) - Vd(activeold)
      if(kepara.gt.exci) then
          activeold = active
          acc = -force(activeold)/mass
          vel = vel + acc*timstp*0.5 -sqrt(2.0*exci/m)
          
      else
          active = activeold
          acc = -force(activeold)/mass
          vel = vel + acc*timstp*0.5 
      end if
      q = q + vel*timstp
      p = mass*vel
  end if

end subroutine classical_propagate
         
subroutine electronic_propagate(p,q,densmat,Vd,Vpd,vl,d,timstp)

  do is=1,nstates
      c_matrix(is,is)=exp((0,-1)*Vd(is)*timstp
  end do
  alpha=1d0
  beta=0d0
  call zgemm('N','N',nstates,nstates,nstates,alpha,vl,nstates,&
       c_matrix,nstates,beta,c_matrix_temp,nstates)
  call zgemm('N','C',nstates,nstates,nstates,alpha,c_matrix_temp,nstates,&
       vl,nstates,beta,c_matrix,nstates)
  ! $\rho(t+\Delta t/2) = e^{iH\Delta t/2}\rho(t)e^{-iH\Delta t/2}
  call zgemm('N','N',nstates,nstates,nstates,alpha,c_matrix,nstates,&
       dens_mat,nstates,beta,dens_temp,nstates)
  call zgemm('N','C',nstates,nstates,nstates,alpha,dens_temp,nstates,&
       c_matrix,nstates,beta,dens_mat,nstates)

end subroutine electronic_propagation

subroutine aush_propagate(p,q,densmat,nstates,delR,delP,Vpd,gamma_collapse)

  integer rkn
  integer ierr, i, rkcycle, idime, resolution
  real(kdp) minitimstp, time, dtmax
  real(kdp) rksum
  real(kdp), allocatable :: rkcoefs(:)
  real(kdp), allocatable :: taudinv(:,:)
  complex(kdp), allocatable :: delRdot(:,:,:), delRinput(:,:), delRoutput(:,:)
  complex(kdp), allocatable :: delPdot(:,:,:), delPinput(:,:), delPoutput(:,:)
  complex(kdp), allocatable :: force(:,:)
  complex(kdp) delRPnorm, zdotc


  rkn = 4                                       !< Setting the rkn integrator. Currently it is the standard rk4
  allocate(rkcoefs(rkn))
  rkcoefs(1) = 1.0
  rkcoefs(2) = 2.0
  rkcoefs(3) = 2.0
  rkcoefs(4) = 1.0
  rksum = sum(rkcoefs)
  dtmax = 0.01/exci(nstates-1)
  allocate(delRdot(nstates,nstates,rkn), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(delPdot(nstates,nstates,rkn), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(delRinput(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(delPinput(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(delRoutput(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(delPoutput(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(force(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')
  allocate(taudinv(nstates,nstates), stat=ierr)
  if(ierr.ne.0) call escape('Allocation error in A-FSSH')

  force = 0.0  
  delRdot = 0.0
  delPdot = 0.0
  time  = 0.0
  gamma_reset = 0.0
  gamma_collapse = 0.0
  resolution = max(int(timstp/dtmax)+1,resolution_au)
  minitimstp = timstp/(dble(resolution))
  
!  write(*,'(a)') 'Propagation for augmented FSSH is done using RK4'
!  write(*,'(a,f20.14)') 'Time step in propagation of delR and delP  ', minitimstp
!  write(*,'(a,i)') 'resolution  ', resolution
  call constructdelF(nstates,natoms,nacl,nacv,ngrad,ndim,exci,fstatenow,force)
  do rkcycle = 1, resolution
      do i = 1, rkn
          delRinput = delR
          delPinput = delP
          if (i > 1) then 
              delRinput = delRinput + delRdot(:,:,:,(i-1))*minitimstp/rkcoefs(i)
              delPinput = delPinput + delPdot(:,:,:,(i-1))*minitimstp/rkcoefs(i)
          end if
          delRoutput = 0
          delPoutput = 0
          
          call delRdelPdot(nstates,exci,nacl,natoms,aumass,nacv, &
                          & ngrad,force,fstatenow, dens_mat, nowpe,&
                          & delRinput,delPinput,&
                          & delRoutput,delPoutput,ldebug)
          delRdot(:,:,i) = delRoutput
          delPdot(:,:,i) = delPoutput
      end do
      do i = 1, rkn
          delR = delR + minitimstp*rkcoefs(i)*delRdot(:,:,i)/rksum    
          delP = delP + minitimstp*rkcoefs(i)*delPdot(:,:,i)/rksum    
      end do
      call probcollapse(delR,delP,force,nstates,natoms,ndim,ngrad,nacl,nacv,exci,fstatenow,taudinv)
      do i = 1,nstates
          gamma_collapse(i) = gamma_collapse(i) + minitimstp*taudinv(i,fstatenow)
      end do
      time = time + minitimstp
  end do
end aush_propagate


subroutine constructdelF(nstates,nacl,nacv,ngrad,ndim,exci,fstatenow,force)
  implicit none
  
  ! Input parameters
  integer, intent(in)           :: nstates                                  !< number of electronic states
  integer, intent(in)           :: fstatenow                                !< current active state (convention, lowest  = 1)
  real(kdp), intent(in)         :: exci(nstates)                            !< the excitation energies, careful with convention
  real(kdp), intent(in)         :: ngrad(nstates)                           !< gradient from all the nstates, separately
  real(kdp), intent(in)         :: nacv((nstates)*(nstates+1))/2)           !< the coupling vectors
  real(kdp), intent(in)         :: nacl(((nstates)*(nstates+1))/2)          !< the coupling vector.velocity
  complex(kdp), intent(inout)   :: force(nstates,nstates)       !< The delF constructed

  !Local Variables
  
  integer ii, ij, idimension, iatoms, counter, idof

  counter = 0
  
  do ij = 1,nstates
      do ii = ij,nstates
  
  
          counter = counter+1
          do iatoms = 1,natoms
              do idimension = 1,ndim
                   idof = ndim*(iatoms-1)+idimension 
                   if (ii.eq.ij) then
                       force(ii,ij,idof) = -ngrad(idimension,iatoms,ii)
                   else
                       force(ii,ij,idof) = force(ii,ij,idof)-nacv(idimension,iatoms,counter)*(exci(ii-1)-exci(ij-1))
                 
                       force(ij,ii,idof) = force(ij,ii,idof)-nacv(idimension,iatoms,counter)*(exci(ii-1)-exci(ij-1))
                   end if
              end do
          end do
      end do
  end do    
  
  do ii = 1,nstates
      do iatoms = 1,natoms
          do idimension = 1,ndim
              idof = ndim*(iatoms-1)+idimension 
              force(ii,ii,idof) = force(ii,ii,idof) + ngrad(idimension,iatoms,fstatenow)
          end do
      end do
  end do

end subroutine constructdelF


subroutine delRdelPdot(nstates,exci,nacl,mass,nacv, &
                          & ngrad,force,fstatenow,dens_mat, &
                          & delRinput,delPinput,&
                          & delRdot,delPdot,ldebug)

  implicit none
  include 'fortrankinds.h'
  
  ! Input parameters
  integer, intent(in)         :: nstates                                  !< number of electronic states
  integer, intent(in)         :: fstatenow                                !< current active state (convention, lowest  = 1)
  real(kdp), intent(in)       :: exci(nstates)                        !< the excitation energies, careful with convention
  real(kdp), intent(in)       :: nacv(((nstates)*(nstates+1))/2) !< the coupling vectors
  real(kdp), intent(in)       :: nacl(((nstates)*(nstates+1))/2)          !< the coupling vector.velocity
  real(kdp), intent(in)       :: ngrad(nstates)                  !< gradient from all the nstates, separately
  real(kdp), intent(in)       :: aumass                           !< array of atomic mass
  complex(kdp), intent(in)    :: dens_mat(nstates,nstates)                !< electronic RDM
  complex(kdp), intent(in)    :: delRinput(nstates,nstates)      !< The distance between frozen gaussian on a state and the classical position
  complex(kdp), intent(in)    :: delPinput(nstates,nstates)      !< The distance between frozen gaussian on a state and the classical momentum 
  complex(kdp), intent(inout) :: delRdot(nstates,nstates)        !< time derivative of delR
  complex(kdp), intent(inout) :: delPdot(nstates,nstates)        !< time derivative of delP
  complex(kdp), intent(in)    :: force(nstates,nstates)        !< time derivative of delP
  logical, intent(in) :: ldebug

  !> To use this subroutine:  It requires the current delR and delP as dynamic inputs
  !! and the Hamiltonian in forms of nowpe, excitation energy and derivative coupling
  !! It also needs the gradients of all the states through ngrad(:,:) and e-RDM

  !Local variables

  integer ierr, counter, ndim, idimension, i, j
  real(kdp), allocatable :: nacmat(:,:)
  complex(kdp) alpha, beta
  complex(kdp), allocatable :: comoperator(:,:),    
  ndim = 3
  allocate(comoperator(nstates,nstates), stat = ierr)
  if(ierr.ne.0) call afail(1,'comoperator','delRdelPdot')
  allocate(nacmat(nstates,nstates), stat = ierr)
  if(ierr.ne.0) call afail(1,'nacmat','delRdelPdot')

  comoperator=0
  counter = 0
  do i = 1, nstates
      comoperator(i, i) = -(0,1)*(exci(i) + nowpe)
      do j = i,nstates
          counter = counter +1
          nacmat(j,i) = -nacl(counter)
          nacmat(i,j) = nacl(counter)
      end do
  end do
  comoperator = comoperator - nacmat
  if(ldebug) write(*,'(a)') 'comoperator'
  if(ldebug) write(*,'(6e25.14)') comoperator
  if(ldebug) write(*,'(a)') 'delF'
  if(ldebug) write(*,'(6e25.14)') force

  !  delRdot = [H,delR] + i[nacmat,delR] + delP/M
  alpha = 1d0
  beta = 0d0
  call zgemm('N','N',nstates,nstates,nstates,alpha, &
       comoperator,nstates,delRinput(1,1),nstates, &
       beta,delRdot(1,1),nstates)
  alpha = -1d0
  beta = 1d0
  call zgemm('N','N',nstates,nstates,nstates,alpha, &
      delRinput(1,1),nstates,comoperator,nstates, &
      beta,delRdot(1,1),nstates)
  
  delRdot(:,:) = delRdot(:,:) + &
           delPinput(:,:)/aumass
 
  !  delPdot = [H,delP] + i[nacmat,delP] + {rho,delF}
  alpha = 1d0
  beta = 0d0
  call zgemm('N','N',nstates,nstates,nstates,alpha, &
       comoperator,nstates,delPinput(1,1,idimension),nstates, &
       beta,delPdot(1,1),nstates)
  alpha = -1d0
  beta = 1d0
  call zgemm('N','N',nstates,nstates,nstates,alpha, &
       delPinput(1,1,idimension),nstates,comoperator,nstates,&
       beta,delPdot(1,1),nstates)
  
  alpha = 0.5d0
  beta = 1d0 
  call zgemm('T','T',nstates,nstates,nstates,alpha, &
       dens_mat(1,1),nstates,force(1,1,idimension),nstates, &
       beta,delPdot(1,1,idimension),nstates)
  
  call zgemm('T','T',nstates,nstates,nstates,alpha, &
            force(1,1),nstates,dens_mat(1,1),nstates, &
            beta,delPdot(1,1),nstates)
    
end subroutine delRdelPdot


subroutine probcollapse(delR,delP,force,nstates,natoms,ndim,ngrad,nacl,nacv,exci,fstatenow,taudinv)
  implicit none


  ! Input parameters
  integer, intent(in)         :: nstates                                  !< number of electronic states
  integer, intent(in)         :: natoms                                   !< number of classical nucleus
  integer, intent(in)         :: fstatenow                                !< current active state (convention, lowest  = 1)
  real(kdp), intent(in)       :: exci(nstates)                        !< the excitation energies, careful with convention
  real(kdp), intent(in)       :: ngrad(nstates)                  !< gradient from all the nstates, separately
  real(kdp), intent(in)       :: nacv(((nstates)*(nstates+1))/2) !< the coupling vectors
  real(kdp), intent(in)       :: nacl(((nstates)*(nstates+1))/2)          !< the coupling vector.velocity
  complex(kdp), intent(in)    :: delR(nstates,nstates)           !< The distance between frozen gaussian on a state and the classical position
  complex(kdp), intent(in)    :: delP(nstates,nstates)           !< The distance between frozen gaussian on a state and the classical momentum 
  complex(kdp), intent(in)    :: force(nstates,nstates)          !< force difference between the active state and the
                                                                          !!other states
  real(kdp), intent(inout)    :: taudinv(nstates,nstates)                 !< The collapse variables to be integrated over dt
  
  !Local Variable
  integer ii, idimension
  real(kdp) temp1, temp2, temp3
  
  taudinv = 0.0
  

  do ii = 1,nstates
          temp1 =0.5*Real(force(ii,ii)&
                        &*(delR(ii,ii) - delR(fstatenow,fstatenow))
          temp2 = Real((delR(ii,ii) - delR(fstatenow,fstatenow))*(delP(ii,ii) - delP(fstatenow,fstatenow)))
  
          temp3 = sign(temp1,temp2)
          taudinv(ii,fstatenow) = taudinv(ii,fstatenow)+temp3
                                   
  end do

end subroutine probcollapse


!> This where the fate of the wavefunction is decided. According to 
!!A-FSSH, it may be collapsed to the active state, or not. If the 
!!conditions for collapsing is satisfied, this subrountine alters the
!!density matrix and delR and delP.
subroutine au_collapse(gamma_reset,gamma_collapse,natoms,nstates,ifstatenow,&
&                      dens_mat,delR,delP,colpse,tshop,time,itotalstps)
  use type_random_state
  implicit none
  
  integer, intent(in) :: nstates                                !< The number of electronic states involved
  integer, intent(in) :: itotalstps                             !< The number of total steps executed      
  integer, intent(in) :: natoms                                 !< The number of 'classical' nuclei
  integer, intent(in) :: ifstatenow                             !< The active state (convention here, the lowest state => ifstatenow = 1)
  real(kdp), intent(in) :: gamma_collapse(nstates)              !< The array containing the probability of destruction of a state
  real(kdp), intent(in) :: gamma_reset(nstates)                 !< The array containing the propability of resetting the pseudo-dynamical variables (delR and delP)
  complex(kdp), intent(inout) :: delR(nstates,nstates,3*natoms) !< The distance between frozen gaussian on a state and the classical position
  complex(kdp), intent(inout) :: delP(nstates,nstates,3*natoms) !< The distance between frozen gaussian on a state and the classical momentum
  complex(kdp), intent(inout) :: dens_mat(nstates,nstates)      !< The electronic RDM 
  logical, intent(inout)      :: colpse                         !< The flag that prevents the propagation of the half step of RDM after collapse
  type(shopt), intent(inout) :: tshop                           !< Contains the various options (including the non-standard ones)
  real(kdp), intent(in) :: time                                 !< Current time
  !*******************************************************************************************************************************
  !*******************************************************************************************************************************
  !Local variables
  real(kdp)               xranf
  integer                 i, j, k, iatoms, iseed_temp, ierr, iaunit, icolstate
  logical                 lneedcollapse,exists,saveseed
  complex(kdp)            temp, one
  character(len=80)       filnam, filnamtemp
  type(random_state) auranobj
  !*********************************************************************!!
  !>  To use this subroutine, the current electronic RDM is to be passed with an 
  !!  array, gamma_collapse which contains the probability of collapse of each state.
  !!  (Collapse here means seeting the set to zero). The corresponding delR and delP are
  !!  also set to zero.
  !!  Options: tshop%force_gs == 3 means bypass this routine
  !!           tshop%override = .true. means every seed gives the same random number. To be
  !!                            explicit, when true, for every sh_seed in mdmaster, it creates
  !!                            an unique seed and random number. the seed is written in sh_seed
  !!                            creating a repeatable sequence of random numbers.
  
  if(mod(tshop%force_gs,3).eq.1) return
  one=(1d0,0d0)
  lneedcollapse = .false.
  filnamtemp = "aush-seed-data"

  call auranobj%initialize(filnamtemp,ldiskdump=.true.,lsaveinitial=.true.,seed=tshop%au_seed)
  call auranobj%random(xranf)
  tshop%au_seed = auranobj%idum

  write(istdout,'(a30,e18.10)') 'Random number for A-FSSH',xranf
  !> This is the block that calls the various random number generator options
  do i = 1,nstates
      
      if(gamma_collapse(i)>xranf) then
      write(*,'(a,i3)') 'Condition satisfied for state ', i
          icolstate = i
          colpse = .true.
          temp = dens_mat(i,i)
          do j = 1,nstates
              do k = 1,nstates
                  dens_mat(j,k) = dens_mat(j,k)/(one-temp)
              end do
         
              dens_mat(i,j) = 0.0
              dens_mat(j,i) = 0.0
          end do
          lneedcollapse = .true.
      end if
  end do   

  if(lneedcollapse) then
      delR = 0.0
      delP = 0.0
  end if

end subroutine au_collapse

subroutine au_hop(nstates,delR,delP,ifstate)

  implicit none

  integer, intent(in) :: nstates                                  !< number of electronic states in RDM
  integer, intent(in) :: ifstate                                   !< current active state. 
  complex(kdp), intent(inout) :: delR(nstates,nstates)     !< The distance between frozen gaussian on a state and the classical position
  complex(kdp), intent(inout) :: delP(nstates,nstates)     !< The distance between frozen gaussian on a state and the classical momentum
  !********************************************************************!
  !> This subroutine is to be called only when there is a hop in FSSH

  !Local variables
  
  integer i
  
  do i = 1,nstates
      if(i.eq.ifstate) cycle
      
      delP(i,i) = delP(i,i) - delP(ifstate,ifstate)
      delR(i,i) = delR(i,i) - delR(ifstate,ifstate)
  
  end do
  
  delP(ifstate,ifstate) = 0.0
  delR(ifstate,ifstate) = 0.0

end subroutine au_hop
