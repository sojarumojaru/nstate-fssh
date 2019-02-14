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
