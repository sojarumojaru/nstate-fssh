subroutine el_pro_dia(densmat,coef,timstp,nstates,ndim,q,eigvec)

  implicit none
  integer, intent(in) :: nstates,ndim
  complex*16, intent(inout) :: coef(nstates),densmat(nstates,nstates),eigvec(nstates,nstates)
  real*8, intent(in) :: q(ndim),timstp
  
  !! Local Var
  
  complex*16, allocatable :: coefdot(:),densmatdia(:,:)
  complex*16, allocatable :: densmatnew(:,:)
  real*16 Vc, g,au2aa,au2rcm
  complex*16 imag,alpha,beta
  integer it,is

  au2aa = .52917721067d0
  au2rcm=219474.63067d0
  Vc = 50.0d0/au2rcm
  g = (3091.8/au2rcm)*au2aa
  imag = (0d0,1d0)

  allocate(coefdot(nstates))
  allocate(densmatdia(nstates,nstates))
  allocate(densmatnew(nstates,nstates))
  coefdot = 0d0
  do it = 1,nstates
      coefdot(it)  = imag*g*q(it)
      if(it>1) coefdot(it) = coefdot(it) +imag*Vc*coef(it-1)
      if(it<5) coefdot(it) = coefdot(it) +imag*Vc*coef(it+1)
  end do

  coef = coef + coefdot*timstp

  densmat = 0d0

  do it = 1,nstates
      do is = 1,nstates
          densmatdia(it,is) = conjg(coef(it))*coef(is)
      end do
  end do
  alpha=(1d0,0d0)
  beta=(0d0,0d0)
  call zgemm('C','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
           densmatdia,nstates,beta,densmatnew,nstates)

  call zgemm('N','N',nstates,nstates,nstates,alpha,densmatnew,nstates,&
           eigvec,nstates,beta,densmat,nstates)
end subroutine el_pro_dia
   

subroutine  electronic_evaluate(mass,omega,p,q,V,Vp,Vd,Vpd,nacv,ndim,&
                                & nstates,active,eigvec,ldebug)

  !! Currently evaluating things only for adiabatic basis.


  implicit none
  logical, intent(in) :: ldebug
  integer, intent(in) :: nstates,active,ndim
  real*8, intent(in) :: mass, omega, p(ndim),q(ndim)
  complex*16, intent(inout) :: V(nstates,nstates), Vd(nstates),Vp(nstates, nstates,ndim)
  complex*16, intent(inout) :: eigvec(nstates,nstates)
  complex*16, intent(inout) :: nacv(nstates,nstates,ndim), Vpd(nstates,nstates,ndim)

  !! Local variables
  
  complex*16, allocatable :: VR(:,:),c_matrix_temp(:,:), temp1(:,:)
  complex*16, allocatable :: WORK(:),vpcopy(:,:,:)
  complex*16 alphac,betac
  real*8 alpha, beta, thresh!, au2rcm,pi,speedoflight,au2aa
  real*8, allocatable ::  RWORK(:),eigen(:)
  integer i,j,info,lwork,idime,it
  logical swapped
  
  allocate(c_matrix_temp(nstates,nstates))
  allocate(Vpcopy(nstates,nstates,ndim))
  allocate(eigen(nstates))
  allocate(temp1(nstates,nstates))
  allocate(RWORK(3*nstates),WORK(2*nstates),VR(nstates,nstates))


  lwork = 2*nstates

  thresh = 1d-6

  call Holstein(q,V,mass,omega)
  call HolsteinGrad(q,Vp,mass,omega)

  if(ldebug) then
      write(*,'(a)') 'V'
      write(*,'(10e18.10)') V(:,:)
      write(*,'(a)') 'Vp'
      write(*,'(10e18.10)') Vp(:,:,:)
  end if
  eigvec = V
  Vpcopy = Vp
  Vpd = 0d0
  eigen = 0d0
  call zheev('V','U',nstates , eigvec, nstates, eigen, work,lwork,rwork,info)

  alphac = (1d0,0d0)
  betac = (0d0,0d0)
!  call sorteigen(nstates,Vd,Vl,swapped)
  if(ldebug) then
      write(*,'(a)') 'eigen vectors'
      write(*,'(10e18.10)') eigvec(:,:)
  
  end if

  do it = 1,ndim
       temp1 = Vp(:,:,it)
       call zgemm('C','N',nstates,nstates,nstates,alphac,eigvec,nstates,&
          & temp1,nstates,betac,c_matrix_temp,nstates)
      if(ldebug) then
          write(*,'(a)') 'c_matrix_temp'
          write(*,'(10e18.10)') c_matrix_temp(:,:)
      
      end if

       call zgemm('N','N',nstates,nstates,nstates,alphac,c_matrix_temp,nstates,&
          & eigvec,nstates,betac,Vpd(:,:,it),nstates)
      if(ldebug) then
          write(*,'(a)') 'inside loop Vpd'
          write(*,'(10e18.10)') Vpd(:,:,:)
      
      end if
  end do
  if(ldebug) then
      write(*,'(a)') 'Vpd'
      write(*,'(10e18.10)') Vpd(:,:,:)
  end if
  nacv = 0d0
  do i = 1, nstates
      Vd(i) = eigen(i)
  end do
  do i = 1, nstates
      do j = 1,nstates
          if(i.eq.j) cycle
          do it = 1,ndim
              if(abs(Vd(j)-Vd(i))>thresh)  nacv(i,j,it) = (Vpd(i,j,it))/(Vd(j) - Vd(i))
          end do
      end do
  end do
end subroutine electronic_evaluate

subroutine sorteigen(nstates,Vd,Vl,swapped)
  implicit none
  logical, intent(inout) :: swapped
  integer, intent(in) :: nstates
  complex*16, intent(inout) :: Vd(nstates), Vl(nstates,nstates)

  !!Local
  complex*16,allocatable:: scratch(:)
  complex*16 scratchval
  integer ii, ij, ik
  swapped = .false.
  allocate(scratch(nstates))
  do ii = 1,nstates
     do ij = 1, nstates - 1
         if(real(Vd(ij)) > real(Vd(ij+1))) then
             swapped = .true.
             scratchval = Vd(ij)
             Vd(ij) = Vd(ij+1)
             Vd(ij+1) = scratchval
             do ik = 1,nstates
                 scratchval = Vl(ij,ik)
                 Vl(ij,ik) = Vl(ij+1,ik)
                 Vl(ij+1,ik) = scratchval
             end do
         end if
     end do
 end do

end subroutine sorteigen

subroutine classical_propagate(p,q,mass,force,active,activeold,nacv,kepara,&
           & timstp,Vd,nstates,ndim,nacl,kine,pote,ldebug)
  implicit none
  
  integer, intent(in) :: nstates,ndim
  integer, intent(inout) :: active,activeold
  real*8, intent(inout) :: p(ndim),q(ndim),kepara
  real*8, intent(inout) :: kine, pote
  real*8, intent(in) :: mass,timstp
  complex*16, intent(in) :: force(nstates,nstates,ndim), nacv(nstates,nstates,ndim)
  complex*16, intent(in) :: Vd(nstates)
  complex*16, intent(inout) :: nacl((nstates*(nstates+1))/2)
  logical, intent(in) :: ldebug

  !! Local variables
  real*8, allocatable :: vel(:), acc(:),rescale_vec(:),ppara(:)
  real*8 exci,sgn,sqvec,a,b,c,vt,tsqm,dnrm2
  real*8 faktorplusp,faktorplusm,faktor
  integer ii,ij,icounter,idime
  allocate(vel(ndim))
  allocate(acc(ndim))
  allocate(rescale_vec(ndim))
  allocate(ppara(ndim))
  icounter = 0
  nacl = 0
  kine = 0d0
  pote = 0d0

  pote = real(Vd(active))

  if(activeold.eq.active) then
      vel = p/mass
      acc= -real(force(active,active,:))/mass
      vel = vel + acc*timstp
      
      do ij = 1,nstates
          do ii = ij,nstates
              icounter = icounter+1
              do idime = 1,ndim
                  nacl(icounter)=nacl(icounter)+(vel(idime)-acc(idime)*timstp/2.0)*nacv(ii,ij,idime)
              end do
          end do
      end do
      q = q + vel*timstp
      p = mass*vel
      pote = real(Vd(active))
      

  end if 
!  write(*,*) 'nacv'
!  write(*,*) nacv
  if(active.ne.activeold) then
      do idime = 1,ndim
          vel(idime) = p(idime)/mass
          acc(idime) = -force(activeold,activeold,idime)/mass
          vel(idime) = vel(idime) + acc(idime)*timstp*0.5
      end do
      icounter = 0
      do ij = 1,nstates
          do ii = ij,nstates
              icounter = icounter+1
              do idime = 1,ndim
                  nacl(icounter)=nacl(icounter)+(vel(idime))*nacv(ii,ij,idime)
              end do
          end do
      end do
      kepara = 0
      rescale_vec = real(nacv(active,activeold,:))
      sqvec = dnrm2(ndim,rescale_vec,1)
      kepara = 0

      do idime = 1,ndim
          ppara(idime) = p(idime)*rescale_vec(idime)
      end do

      ppara = ppara/sqvec

      do idime = 1,ndim
          kepara = kepara+0.5*ppara(idime)**2/mass
      end do

      exci = real(Vd(active)) - real(Vd(activeold))


      
      do idime = 1,ndim
         vt = vt + mass*vel(idime)*rescale_vec(idime)
         b = b + vel(idime)*rescale_vec(idime)
         a = a + 0.5* rescale_vec(idime)*rescale_vec(idime)/mass
         tsqm = tsqm + rescale_vec(idime)*rescale_vec(idime)/mass
      end do
      if (exci > 0) sgn = -1.0
      if (exci < 0) sgn = 1.0
      c = exci
      if((kepara.gt.exci).and.(b*b>4*a*c)) then
          
          faktorplusp = (-b + sqrt(b*b - 4*a*c))/(2*a)
          faktorplusm = (-b - sqrt(b*b - 4*a*c))/(2*a)

          if (abs(faktorplusp) <= abs(faktorplusm)) then
              faktor = faktorplusp
          else
              faktor = faktorplusm
          end if
          activeold = active
          do idime = 1,ndim
              acc(idime) = -real(force(activeold,activeold,idime))/mass
              vel(idime) = vel(idime) + acc(idime)*timstp*0.5 + faktor*rescale_vec(idime)/mass
          end do
          
      else
          if(ldebug)  write(*,*) 'Frustrated Hop'
          if(ldebug) write(*,*) 'kepara, exci  ', kepara, exci
          active = activeold
          do idime = 1,ndim
              acc(idime) = -real(force(active,active,idime))/mass
              vel(idime) = vel(idime) + acc(idime)*timstp*0.5 
          end do
      end if
      q = q + vel*timstp
      p = mass*vel
  end if
  pote = real(Vd(active))
  do idime = 1,ndim
     kine = kine + 0.5d0*mass*vel(idime)**2 
  end do
!  write(*,'(a)') '#p(1),q(1), kepara,pote,totale,active'
!  write(*,'(5e15.6,i3)') p(1),q(1), kepara,pote,totale,active

  if(isnan(q(1))) then
      write(*,*) 'faktor'
      write(*,*) faktor

      write(*,*) 'exci'
      write(*,*) exci

      write(*,*) 'kepara'
      write(*,*) kepara

      write(*,*) 'rescale_vec'
      write(*,'(10e18.10)') rescale_vec
     
      write(*,*) 'force'
      write(*,'(10e18.10)') force

      write(*,*) 'nacv'
      write(*,'(10e18.10)') nacv

      write(*,*) 'Vd'
      write(*,'(10e18.10)') Vd
      stop

  end if





end subroutine classical_propagate


subroutine electronic_propagate(densmat,Vd,nacl,timstp,nstates,ldebug,ldiabatic,eigvec,V)
  implicit none
  logical, intent(in) :: ldebug,ldiabatic
  integer, intent(in) :: nstates
  real*8, intent(in) :: timstp
  complex*16, intent(in) :: nacl((nstates*(nstates+1)/2))
  complex*16, intent(in) :: Vd(nstates)
  complex*16, intent(inout) :: densmat(nstates,nstates)
  complex*16, intent(inout) :: eigvec(nstates,nstates)
  complex*16, intent(inout) :: V(nstates,nstates)

  !! Local 
  complex*16, allocatable :: densmatnew(:,:)
  complex*16, allocatable :: b_matrix(:,:)
  complex*16, allocatable :: c_matrix(:,:)
  complex*16, allocatable :: c_matrix_temp(:,:)
  complex*16 alpha, beta
  integer it,is,counter
  complex*16, allocatable :: vr(:,:), work(:), vl(:,:), w(:)
  real*8, allocatable ::  rwork(:),w_real(:)
  integer i,j,info,lwork,ij
 
  complex*16 minusi
 
  allocate(rwork(4*nstates),work(4*nstates),vr(nstates,nstates),vl(nstates,nstates))
  lwork = 4*nstates

  allocate(b_matrix(nstates,nstates))
  allocate(c_matrix(nstates,nstates))
  allocate(c_matrix_temp(nstates,nstates))
  allocate(densmatnew(nstates,nstates))
  allocate(w(nstates))
  allocate(w_real(nstates))
  
  minusi = (0d0,-1d0)
  counter = 0
  b_matrix = 0d0
  c_matrix_temp = 0d0
  densmatnew = 0d0
  c_matrix = 0d0
  w = 0d0
  
  vl = 0d0

  alpha=(1d0,0d0)
  beta=(0d0,0d0)
  if(ldebug)  write(*,'(a)') 'nacl'
  if(ldebug)  write(*,'(2e18.10)') nacl
  if(ldiabatic) then
      if(ldebug)  write(*,'(a)') 'eigvec'
      if(ldebug)  write(*,'(2e18.10)') eigvec

      if(ldebug) write(*,*) 'densmat before transform'
      if(ldebug)  write(*,'(2e18.10)') densmat
      densmatnew = 0d0
  !    do it = 1,nstates
  !        do is = 1,nstates
  !            do ij = 1,nstates
  !                densmatnew(it,is) = densmatnew(it,is) + eigvec(it,ij)*densmat(ij,is)
  !            end do
  !        end do
  !    end do
      if(ldebug) write(*,*) 'densmatnew before transform'
      if(ldebug)  write(*,'(2e18.10)') densmatnew
!      call zgemm('N','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
!           densmat,nstates,beta,densmatnew,nstates)
!      densmat= 0d0
!      do it = 1,nstates
!          do is = 1,nstates
!              do ij = 1,nstates
!                  densmat(it,is) = densmat(it,is) + densmatnew(it,ij)*conjg(eigvec(is,ij))
!              end do
!          end do
!      end do
!      call zgemm('N','C',nstates,nstates,nstates,alpha,densmatnew,nstates,&
!           eigvec,nstates,beta,densmat,nstates)
      if(ldebug) write(*,*) 'densmat after transform'
      if(ldebug)  write(*,'(2e18.10)') densmat
      if(ldebug) write(*,*) 'V-diabatic'
      if(ldebug)  write(*,'(2e18.10)') V
      b_matrix = V
  else
      do it = 1, nstates
          do is = it, nstates
              counter=counter+1
              b_matrix(it,is) = -(0.0,1.0)*nacl(counter)
              b_matrix(is,it) = (0.0,1.0)*nacl(counter)
          end do
      end do

      do it = 1, nstates
          b_matrix(it,it) = b_matrix(it,it) + Vd(it)
      end do
  end if
!  write(*,'(a)') 'b_matrix'
!  write(*,'(2e18.10)') b_matrix

  call zheev('V','U',nstates,b_matrix,nstates,w_real,work,lwork,rwork,info)
!  call zgeev('V','N',nstates,b_matrix,nstates,w,vl,nstates,vr,&
!       & nstates,work,(2*(nstates)),rwork,info )
  if(ldebug)  write(*,'(a)') 'info'
  if(ldebug)  write(*,'(i3)') info
  if(ldebug)  write(*,*) 'w_real'
  if(ldebug)  write(*,'(2e18.10)') w_real

  do is=1,nstates
      c_matrix(is,is)=exp(minusi*w_real(is)*timstp)
  end do
  
  ! $\rho(t+\Delta t/2) = e^{iH\Delta t/2}\rho(t)e^{-iH\Delta t/2}
  call zgemm('C','N',nstates,nstates,nstates,alpha,b_matrix,nstates,&
       c_matrix,nstates,beta,c_matrix_temp,nstates)
  call zgemm('N','N',nstates,nstates,nstates,alpha,c_matrix_temp,nstates,&
       b_matrix,nstates,beta,c_matrix,nstates)

  if(ldebug) write(*,*) 'c_matrix after'
  if(ldebug)  write(*,'(2e18.10)') c_matrix
  if(ldebug) write(*,*) 'densmat before'
  if(ldebug)  write(*,'(2e18.10)') densmat
  call zgemm('N','N',nstates,nstates,nstates,alpha,c_matrix,nstates,&
       densmat,nstates,beta,densmatnew,nstates)
  
  if(ldebug) write(*,*) 'densmatnew intermed'
  if(ldebug) write(*,'(2e18.10)') densmatnew
  call zgemm('N','C',nstates,nstates,nstates,alpha,densmatnew,nstates,&
       c_matrix,nstates,beta,densmat,nstates)
  if(ldebug) write(*,*) 'densmat after'
  if(ldebug) write(*,'(2e18.10)') densmat
  if(ldiabatic) then
!      densmatnew = 0d0
!      do it = 1,nstates
!          do is = 1,nstates
!              do ij = 1,nstates
!                  
!                  densmatnew(it,is) = densmatnew(it,is) + conjg(eigvec(ij,it))*densmat(ij,is)
!              end do
!          end do
!      end do
!      call zgemm('N','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
!           densmat,nstates,beta,densmatnew,nstates)
!      densmat = 0d0
!      do it = 1,nstates
!          do is = 1,nstates
!              do ij = 1,nstates
!                  densmat(it,is) = densmat(it,is) + densmatnew(it,ij)*eigvec(ij,is)
!              end do
!          end do
!      end do
!      call zgemm('N','C',nstates,nstates,nstates,alpha,densmatnew,nstates,&
!           eigvec,nstates,beta,densmat,nstates)
!      call zgemm('C','N',nstates,nstates,nstates,alpha,eigvec,nstates,&
!           densmat,nstates,beta,densmatnew,nstates)

!      call zgemm('N','N',nstates,nstates,nstates,alpha,densmatnew,nstates,&
!           eigvec,nstates,beta,densmat,nstates)
       
      if(ldebug) write(*,*) 'densmat after back to adiabatic'
      if(ldebug) write(*,'(2e18.10)') densmat
  end if 



end subroutine electronic_propagate




subroutine aush_propagate(densmat,p,nstates,ndim,delR,delP,Vpd,gamma_collapse,&
           & gamma_reset,mass,nacl,exci,timstp,ifstatenow)

  implicit none 

  integer, intent(in) :: nstates,ndim
  

  complex*16, intent(inout) :: delR(nstates,nstates,ndim)
  complex*16, intent(inout) :: delP(nstates,nstates,ndim)
  complex*16, intent(inout) :: densmat(nstates,nstates)
  complex*16, intent(inout) :: Vpd(nstates,nstates,ndim)
  complex*16, intent(in) :: nacl((nstates*(nstates+1))/2)
  complex*16, intent(in)    :: exci(nstates)
  real*8,intent(inout) :: gamma_collapse(nstates)
  real*8, intent(inout) :: p(ndim)
  real*8,intent(inout) :: gamma_reset(nstates)
  real*8, intent(in) ::  mass,timstp
  integer, intent(in) :: ifstatenow
  


  integer rkn,it,is
  integer ierr, i, rkcycle, idime, resolution,resolution_au
  real*8 minitimstp, time, dtmax
  real*8 rksum
  real*8, allocatable :: rkcoefs(:)
  real*8, allocatable :: taudinv(:,:),taurinv(:,:)
  complex*16, allocatable :: delRdot(:,:,:,:), delRinput(:,:,:), delRoutput(:,:,:)
  complex*16, allocatable :: delPdot(:,:,:,:), delPinput(:,:,:), delPoutput(:,:,:)
  complex*16, allocatable :: force(:,:,:), vel(:)
  complex*16 delRPnorm, zdotc
  

  rkn = 4                                       !< Setting the rkn integrator. Currently it is the standard rk4
  resolution_au = 10
  allocate(rkcoefs(rkn))
  rkcoefs(1) = 1.0
  rkcoefs(2) = 2.0
  rkcoefs(3) = 2.0
  rkcoefs(4) = 1.0
  rksum = sum(rkcoefs)

  dtmax = 0.01/real(exci(nstates-1))

  allocate(delRdot(nstates,nstates,ndim,rkn), stat=ierr)
  allocate(delPdot(nstates,nstates,ndim,rkn), stat=ierr)
  allocate(delRinput(nstates,nstates,ndim), stat=ierr)
  allocate(delPinput(nstates,nstates,ndim), stat=ierr)
  allocate(delRoutput(nstates,nstates,ndim), stat=ierr)
  allocate(delPoutput(nstates,nstates,ndim), stat=ierr)
  allocate(force(nstates,nstates,ndim), stat=ierr)
  allocate(vel(ndim),stat=ierr)
  allocate(taudinv(nstates,nstates), stat=ierr)
  allocate(taurinv(nstates,nstates), stat=ierr)

  vel = p/mass
  force = -Vpd  
  do it = 1,nstates
      do is = 1,ndim
          force(it,it,is) = force(it,it,is) + Vpd(ifstatenow,ifstatenow,is)
      end do
  end do

  delRdot = 0d0
  delPdot = 0d0
  time  = 0d0
  gamma_reset = 0d0
  gamma_collapse = 0d0
  resolution = max(int(timstp/dtmax)+1,resolution_au)
  minitimstp = timstp/(dble(resolution))

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
          
          call delRdelPdot(nstates,ndim,exci,mass,nacl, &
                          & force,ifstatenow, densmat,&
                          & delRinput,delPinput,&
                          & delRoutput,delPoutput)
          delRdot(:,:,:,i) = delRoutput
          delPdot(:,:,:,i) = delPoutput
      end do
      do i = 1, rkn
          delR = delR + minitimstp*rkcoefs(i)*delRdot(:,:,:,i)/rksum    
          delP = delP + minitimstp*rkcoefs(i)*delPdot(:,:,:,i)/rksum    
      end do
      call probcollapse(vel,delR,delP,force,nstates,ndim,ifstatenow,taudinv,taurinv)
      do i = 1,nstates
          gamma_collapse(i) = gamma_collapse(i) + minitimstp*taudinv(i,ifstatenow)
          gamma_reset(i) = gamma_reset(i) + minitimstp*taurinv(i,ifstatenow)
      end do
      time = time + minitimstp
  end do
end subroutine aush_propagate

          
subroutine delRdelPdot(nstates,ndim,exci,mass,nacl, &
                          & force,fstatenow, densmat,&
                          & delRinput,delPinput,&
                          & delRdot,delPdot)


  implicit none
  
  ! Input parameters
  integer, intent(in)         :: nstates                            !< number of electronic states
  integer, intent(in)         :: ndim                               !< number of classical degrees of freedom
  integer, intent(in)         :: fstatenow                                !< current active state (convention, lowest  = 1)
  complex*16, intent(in)      :: exci(nstates)                        !< the excitation energies, careful with convention
  real*8, intent(in)          :: mass                                 !< the excitation energies, careful with convention
  complex*16, intent(in)          :: nacl(((nstates)*(nstates+1))/2)          !< the coupling vectors
  complex*16, intent(in)    :: densmat(nstates,nstates)                !< electronic RDM
  complex*16, intent(in)    :: delRinput(nstates,nstates,ndim)      !< The distance between frozen gaussian on a state and the classical position
  complex*16, intent(in)    :: delPinput(nstates,nstates,ndim)      !< The distance between frozen gaussian on a state and the classical momentum 
  complex*16, intent(inout) :: delRdot(nstates,nstates,ndim)        !< time derivative of delR
  complex*16, intent(inout) :: delPdot(nstates,nstates,ndim)        !< time derivative of delP
  complex*16, intent(in)    :: force(nstates,nstates,ndim)        !< time derivative of delP

  !> To use this subroutine:  It requires the current delR and delP as dynamic inputs
  !! and the Hamiltonian in forms of nowpe, excitation energy and derivative coupling
  !! It also needs the gradients of all the states through ngrad(:,:) and e-RDM

  !Local variables

  integer ierr, counter, i, j, idimension
  real*8, allocatable :: nacmat(:,:)
  complex*16 alpha, beta
  complex*16, allocatable :: comoperator(:,:)  

  allocate(comoperator(nstates,nstates), stat = ierr)
  allocate(nacmat(nstates,nstates), stat = ierr)

  comoperator=0
  counter = 0
  do i = 1, nstates
      comoperator(i, i) = -(0,1)*exci(i) 
      do j = i,nstates
          counter = counter +1
          nacmat(j,i) = -nacl(counter)
          nacmat(i,j) = nacl(counter)
      end do
  end do
  comoperator = comoperator - nacmat

  !  delRdot = [H,delR] + i[nacmat,delR] + delP/M
  alpha = 1d0
  beta = 0d0
  do idimension = 1,ndim
      call zgemm('N','N',nstates,nstates,nstates,alpha, &
           comoperator,nstates,delRinput(1,1,idimension),nstates, &
           beta,delRdot(1,1,idimension),nstates)
      alpha = -1d0
      beta = 1d0
      call zgemm('N','N',nstates,nstates,nstates,alpha, &
          delRinput(1,1,idimension),nstates,comoperator,nstates, &
          beta,delRdot(1,1,idimension),nstates)
      
      delRdot(:,:,idimension) = delRdot(:,:,idimension) + &
           delPinput(:,:,idimension)/mass
  end do
  !  delPdot = [H,delP] + i[nacmat,delP] + {rho,delF}
  do idimension = 1,ndim
      alpha = 1d0
      beta = 0d0
      call zgemm('N','N',nstates,nstates,nstates,alpha, &
           comoperator,nstates,delPinput(1,1,idimension),nstates, &
           beta,delPdot(1,1,idimension),nstates)
      alpha = -1d0
      beta = 1d0
      call zgemm('N','N',nstates,nstates,nstates,alpha, &
           delPinput(1,1,idimension),nstates,comoperator,nstates,&
           beta,delPdot(1,1,idimension),nstates)
      
      alpha = 0.5d0
      beta = 1d0 
      call zgemm('T','T',nstates,nstates,nstates,alpha, &
           densmat(1,1),nstates,force(1,1,idimension),nstates, &
           beta,delPdot(1,1,idimension),nstates)
      
      call zgemm('T','T',nstates,nstates,nstates,alpha, &
                force(1,1,idimension),nstates,densmat(1,1),nstates, &
                beta,delPdot(1,1,idimension),nstates)
  end do
        
end subroutine delRdelPdot


subroutine probcollapse(vel,delR,delP,force,nstates,ndim,ifstatenow,taudinv,taurinv)
  implicit none


  ! Input parameters
  integer, intent(in)       :: nstates                                  !< number of electronic states
  integer, intent(in)       :: ndim                                     !< number of classical degrees of freedom
  integer, intent(in)       :: ifstatenow                               !< current active state (convention, lowest  = 1)
  complex*16, intent(in)    :: delR(nstates,nstates,ndim)               !< The distance between frozen gaussian on a state and the classical position
  complex*16, intent(in)    :: delP(nstates,nstates,ndim)               !< The distance between frozen gaussian on a state and the classical momentum
  complex*16, intent(in)    :: force(nstates,nstates,ndim)              !< force difference between the active state and the
                                                                           !!other states
  real*8, intent(inout)    :: taudinv(nstates,nstates)                  !< The collapse variables to be integrated over dt
  real*8, intent(inout)    :: taurinv(nstates,nstates)                  !< The reset variables to be integrated over dt
  real*8, intent(inout)    :: vel(ndim)
  
  !Local Variable
  integer ii,idimension 
  real*8 temp1, temp2, temp3
  
  taudinv = 0d0
  taurinv = 0d0
  do ii = 1,nstates
      do idimension = 1,ndim
          temp1 =0.5*Real(force(ii,ii,idimension)&
                        &*(delR(ii,ii,idimension) - delR(ifstatenow,ifstatenow,idimension)))


          temp2 = abs(force(ii,ifstatenow,idimension)*(delR(ii,ii,idimension) -&
          &delR(ifstatenow,ifstatenow,idimension)*vel(idimension)))


          temp3 = Real((delR(ii,ii,idimension) - delR(ifstatenow,ifstatenow,idimension))*&
          &       (delP(ii,ii,idimension) - delP(ifstatenow,ifstatenow,idimension)))

          taudinv(ii,ifstatenow) = taudinv(ii,ifstatenow)+temp1 - 2.0 * temp2/(sum(vel*vel))
          taurinv(ii,ifstatenow) = taurinv(ii,ifstatenow)-temp1

      end do
  end do

!  do ii = 1,nstates
!      temp1 =0.5*Real(force(ii,ii)&
!            &*(delR(ii,ii) - delR(ifstatenow,ifstatenow)))
!      temp2 = Real((delR(ii,ii) - delR(ifstatenow,ifstatenow))*(delP(ii,ii) - delP(ifstatenow,ifstatenow)))
!
!      temp3 = 2.0*real(force(ii,ifstatenow)*(delR(ii,ii) - delR(ifstatenow,ifstatenow)))
!      taurinv(ii,ifstatenow) = taurinv(ii,ifstatenow)-temp2
!
!      taudinv(ii,ifstatenow) = taudinv(ii,ifstatenow) + temp2 - temp3
!
!  end do

end subroutine probcollapse


!> This where the fate of the wavefunction is decided. According to 
!!A-FSSH, it may be collapsed to the active state, or not. If the 
!!conditions for collapsing is satisfied, this subrountine alters the
!!density matrix and delR and delP.
subroutine au_collapse(gamma_reset,gamma_collapse,densmat,delR,delP,nstates,ifstatenow,lcol,ldebug3)
  implicit none
  
  logical, intent(inout) :: lcol,ldebug3
  integer, intent(in) :: nstates                                !< The number of electronic states involved
  integer, intent(in) :: ifstatenow                             !< The active state (convention here, the lowest state => ifstatenow = 1)
  real*8, intent(in) :: gamma_collapse(nstates)              !< The array containing the probability of destruction of a state
  real*8, intent(in) :: gamma_reset(nstates)                 !< The array containing the propability of resetting the pseudo-dynamical variables (delR and delP)
  complex*16, intent(inout) :: delR(nstates,nstates) !< The distance between frozen gaussian on a state and the classical position
  complex*16, intent(inout) :: delP(nstates,nstates) !< The distance between frozen gaussian on a state and the classical momentum
  complex*16, intent(inout) :: densmat(nstates,nstates)      !< The electronic RDM 

  !*******************************************************************************************************************************
  !*******************************************************************************************************************************
  !Local variables
  real*8  xranf
  integer i, j, k,icolstate
  logical  lcolpse, lreset
  complex*16  temp1, temp2, one
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
  



  one=(1d0,0d0)
  call random_number(xranf)

  
  lcol = .false.
  !> This is the block that calls the various random number generator options
  lreset = .false. 

  do i = 1,nstates
      if(gamma_reset(i)>xranf) then
          if(ldebug3) write(*,'(a,i3)') 'Reset condition satisfied for state ', i
          icolstate = i
          lreset = .true.
      end if

  end do   

  if(lreset) then
      delR = 0.0
      delP = 0.0
  end if

  if(.not.lreset) then
      do i = 1,nstates
    
          lcolpse = .false. 
          if((gamma_collapse(i)>xranf).and.(.not.lcolpse)) then
              lcol = .true.
              if(ldebug3) write(*,'(a)') 'gamma_collapse'
              if(ldebug3) write(*,'(5e18.10)') gamma_collapse
              if(ldebug3) write(*,'(a,i3)') 'Collapse condition satisfied for state ', i
              if(ldebug3) write(*,'(a)') 'densmat before in adiabatic before collapse:'
              if(ldebug3) write(*,'(10e18.10)') densmat
              if(ldebug3) write(*,'(a)') ''
              icolstate = i
              lcolpse = .true.

              temp1 = densmat(i,i)
              temp2 = densmat(ifstatenow,ifstatenow)

              densmat(ifstatenow,ifstatenow) = (temp1+temp2)
              do j = 1,nstates
                  if (j.eq.ifstatenow) cycle
                  densmat(ifstatenow,j) = densmat(ifstatenow,j) *sqrt((temp1+temp2)/temp2)
                  densmat(j,ifstatenow) = densmat(j,ifstatenow) *sqrt((temp1+temp2)/temp2)
                  densmat(i,j) = 0.0
                  densmat(j,i) = 0.0
              end do



              if(ldebug3) write(*,'(a)') 'densmat in adiabatic after collapse:'
              if(ldebug3) write(*,'(10e18.10)') densmat
              if(ldebug3) write(*,'(a)') ''
          end if
    
          if(lcolpse) then
              delR = 0.0
              delP = 0.0
          end if
      end do
  end if

end subroutine au_collapse

subroutine au_hop(delR,delP,nstates,ndim,ifstate)

  implicit none

  integer, intent(in) :: nstates                                  !< number of electronic states in RDM
  integer, intent(in) :: ndim                                     !< number of classical degrees of freedom
  integer, intent(in) :: ifstate                                   !< current active state. 
  complex*16, intent(inout) :: delR(nstates,nstates,ndim)     !< The distance between frozen gaussian on a state and the classical position
  complex*16, intent(inout) :: delP(nstates,nstates,ndim)     !< The distance between frozen gaussian on a state and the classical momentum
  !********************************************************************!
  !> This subroutine is to be called only when there is a hop in FSSH

  !Local variables
  
  integer ii,ij
  
  do ii = 1,nstates
      if(ii.eq.ifstate) cycle
      do ij = 1,ndim
          delP(ii,ii,ij) = delP(ii,ii,ij) - delP(ifstate,ifstate,ij)
          delR(ii,ii,ij) = delR(ii,ii,ij) - delR(ifstate,ifstate,ij)
      end do
  end do
  
  delP(ifstate,ifstate,:) = 0d0
  delR(ifstate,ifstate,:) = 0d0

end subroutine au_hop

subroutine createEigenVectors(Vd, CC, V)
  implicit none
  complex*16, intent(inout):: Vd(2)
  complex*16, intent(inout), dimension(2,2):: CC, V
  !! local
  complex*16 k1,k2,d1,d2,e1,e2

  e1= V(1,1)+V(2,2)
  e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)
  Vd(1) = (e1-sqrt(e2))/2.0
  Vd(2) = (e2+sqrt(e2))/2.0

  if(abs(V(1,2))>1d-15) then
      k1 = (Vd(1)-V(1,1))/(V(1,2))
      k2 = (Vd(2)-V(1,1))/(V(1,2))
      d1 = Sqrt(1.0+k1**2)
      d2 = Sqrt(1.0+k2**2)

      CC(1,1) = 1.0/d1
      CC(1,2) = 1.0/d2
      CC(2,1) = k1/d1
      CC(2,2) = k2/d2
  else if(real(V(1,1))<real(V(2,2))) then
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
end subroutine createEigenVectors
