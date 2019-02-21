program frog_batch
  implicit none
  real*16 p,q,mass,tmax,timstp,p_initial,p_final,p_step,kepara,xranf
  real*16 time
  real*16, allocatable :: proba(:,:), gamma_collapse(:),gamma_reset(:)

  complex*32, allocatable :: V(:,:), Vd(:), Vp(:,:), Vpd(:,:), vl(:,:) 
  complex*32, allocatable :: densmat(:,:),  delR(:,:),nacl(:)
  complex*32, allocatable :: delP(:,:), c_matrix(:,:), nacv(:,:)
  complex*32, allocatable :: b_matrix(:,:), b_matrix_temp(:,:), c_matrix_temp(:,:)
  
  integer active,it,is,irun,nruns,activeold,nstates
  integer, allocatable :: stat(:)
  
  logical terminate

  


  allocate(V(nstates,nstates))
  allocate(vl(nstates,nstates))
  allocate(Vp(nstates,nstates))
  allocate(Vd(nstates))
  allocate(Vpd(nstates,nstates))
  allocate(nacv(nstates,nstates))
  allocate(nacl((nstates*(nstates+1))/2))

  allocate(b_matrix_temp(nstates,nstates))

  allocate(densmat(nstates,nstates))
  allocate(delR(nstates,nstates))
  allocate(delP(nstates,nstates))
  
  allocate(proba(nstates,nstates),gamma_collapse(nstates))
  allocate(gamma_reset(nstates))


  p_initial = 10.0
  p_final = 10.0
  p_step = 1.0
  tmax = 100
  timstp = 1
  nruns = 10
  mass = 2000.0

  nstates = 2
  

  stat = 0

  do while(p_initial < p_final) !initial momentum loop
      irun = 0 
      do while(irun<nruns) !indivudual trajectory loop
         irun = irun + 1
          call initialize(p,q,densmat,active,p_initial)
          terminate = .false.
          active = 1
          activeold = active
          time = 0
          do while((time<tmax).or.(terminate)) 
              
              call electronic_evaluate(p,q,V,Vp,Vd,Vpd,nacv,nstates,active,vl)

              call classical_propagate(p,q,mass,Vpd,active,&
              &    activeold,nacv,kepara,timstp,Vd,nstates,nacl)

              call electronic_propagate(p,q,densmat,Vd,Vpd,nacv,(timstp/2.0))

              call aush_propagate(densmat,delR,Vpd, &
              &    gamma_collapse,gamma_reset,nacl,Vd,timstp)
              do is = 1,nstates
                  do it = 1,nstates
                      b_matrix_temp(is,it) = nacv(is,it)
                  end do
              end do
              do is =1,nstates
                  b_matrix_temp(is,is) = b_matrix_temp(is,is)+(0,1)*Vd(is)
              end do

              do is=1,nstates
                  if(active.eq.is) cycle
                  proba(active,is)=real(densmat(active,is)*b_matrix_temp(is,active))
                  proba(active,is)=(2.0)*timstp*proba(active,is)/real(densmat(active,active))
              end do

              call electronic_propagate(p,q,densmat,Vd,Vpd,nacv,(timstp/2.0))

              call random_number(xranf)

              call select_newstate(active,xranf,proba,nstates)
              if (active.ne.activeold)  then 
                  call au_hop(delR,delP,nstates,active)
              else
                  call au_collapse(gamma_collapse,nstates,active,densmat,delR,delP,timstp)
              end if
              time = time + timstp
              if ((abs(q)) > 10) terminate = .true.
          end do! time
          if ((q.gt.0.0).and.(active.eq.1)) stat(1) = stat(1) + 1
          if ((q.lt.0.0).and.(active.eq.1)) stat(2) = stat(2) + 1
          if ((q.gt.0.0).and.(active.eq.0)) stat(3) = stat(3) + 1
          if ((q.lt.0.0).and.(active.eq.0)) stat(4) = stat(4) + 1
      end do ! run
      write(*,*) p_initial, stat
  end do! momentum

end program

