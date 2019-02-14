program frog_batch
  implicit none
  
  real*8 p, q, pdot, rdot, h, k1p, k2p, k3p, k4p, k1q, k2q, k3q, k4q,random,random2
  real*8,  allocatable :: V(:,:),d(:,:)
  complex*16, dimension(2) :: coef(:), coefdot(:)
  complex*16, dimension(2,2):: a, adot, delR, delP, delRdot, delPdot
  real*8 t,m, runtime, gamma_c, temp3
  real*8, dimension(4):: plist
  real*8, parameter:: hbar=1 !1.05457173d-34
  Integer state, counter, flag, flag2, loopcount, run, idex, sh_seed, rnd_size
  real*8 e1, e2, d12, poriginal, b12, b21, probability
  real*8, dimension(2):: l
  real*8, dimension(2,2) :: CC, CC1, Vp, Vpd, Vd
  
  inputuinit = 1
  open(inputunit,file='Input',stat=ierr)
  read(inputunit,*) dummy, nstates


  allocate(V(nstates,nstates))
  allocate(Vp(nstates,nstate))
  allocate(Vd(nstates))
  allocate(Vpd(nstates,nstate))

  p_initial = 10.0
  p_final = 40.0
  p_step = 1.0
  tmax = 100
  timstp = 1
  nruns = 10
  mass = 2000.0

  do while(p_initial < p_final) !initial momentum loop
      irun = 0 
      do while(irun<nruns) !indivudual trajectory loop
         irun = irun + 1
          call initialize(p,q,densmat,active,p_initial)
          terminate = .false.
          activeold = active
          time = 0
          do while((time<tmax).or.(terminate)) 
              
              call electronic_evaluate(p,q,V,Vp,Vd,Vpd,d,nstates,active,vl)

              call classical_propagate(p,q,mass,Vpd,active,&
              &    activeold,d,kepara,timstp, Vd)
              call electronic_propagate(p,q,densmat,Vd,Vpd,d,(timstp/2.0))

              call aush_propagate(p,q,densmat,delR,delP,Vpd,gamma_collapse)
              do is = 1,nstates
                  do it = 1,nstates
                      b_matrix_temp(is,it) = d(is,it)
                  end do
              end do
              do is =1,nstates
                  b_matrix_temp(is,is) = b_matrix_temp(is,is)+(0,1)*vl(is)
              end do

              do is=1,nstates
                  if(active.eq.is) cycle
                  proba(active,is)=real(densmatnew(exopt+1,is)*b_matrix_temp(is,exopt+1))
                  proba(active,is)=(2.0)*timstp*proba(exopt+1,is)/pop(exopt+1)
              end do

              call electronic_propagate(p,q,densmat,Vd,Vpd,d,(timstp/2.0))

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

