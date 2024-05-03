!
! 2D isotropic visco-acoustic modeling program
!!! modified on January 12th, 2014
  subroutine iso_visco_step(vx,vz,epsilon_m,p,r,memory_dvx_dx,memory_dvz_dz,&
             memory_dp_dx,memory_dp_dz,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
             a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
             c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

     implicit none
     integer::lv,nx,nz,nrel
     real::dlh(1:lv),dlv(1:lv)
     real::vx(-lv:nx+lv,-lv:nz+lv),vz(-lv:nx+lv,-lv:nz+lv)
     real::p(-lv:nx+lv,-lv:nz+lv),r(1:nrel,-lv:nx+lv,-lv:nz+lv)
     real::epsilon_m(-lv:nx+lv,-lv:nz+lv)

     real::memory_dvx_dx(1:nx,1:nz),memory_dvz_dz(1:nx,1:nz)
     real::memory_dp_dx(1:nx,1:nz),memory_dp_dz(1:nx,1:nz)

     real::a_x(1:nx,1:nz),a_z(1:nx,1:nz),b_x(1:nx,1:nz),b_z(1:nx,1:nz),&
           k_x(1:nx),k_z(1:nz)

     real::a_x_half(1:nx,1:nz),a_z_half(1:nx,1:nz),b_x_half(1:nx,1:nz),&
           b_z_half(1:nx,1:nz),k_x_half(1:nx),k_z_half(1:nz)

     real::c11_u(1:nx,1:nz),tau11(1:nrel,1:nx,1:nz),&
           dens(1:nx,1:nz)

     real,dimension(1:nrel,1:nx,1:nz)::x1,x2

     real::sum_r,value_dvx_dx,value_dvz_dz,value_dp_dx,&
           value_dp_dz,zero

     logical::modeling

     real::deltat
     integer::i,j,l,isource,jsource

     zero=0.e0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,value_dp_dx)
!$OMP DO SCHEDULE(DYNAMIC,1)
     do j = 1,Nz
       do i = 1,NX
 
        value_dp_dx = dlh(1)*(p(i+1,j) - p(i,j)) &
                         +dlh(2)*(p(i+2,j) - p(i-1,j)) &
                         +dlh(3)*(p(i+3,j) - p(i-2,j)) &
                         +dlh(4)*(p(i+4,j) - p(i-3,j)) 

        memory_dp_dx(i,j) = b_x_half(i,j) * memory_dp_dx(i,j) &
                          + a_x_half(i,j) * value_dp_dx

        value_dp_dx = value_dp_dx / K_x_half(i) + memory_dp_dx(i,j)

        vx(i,j) = vx(i,j) - value_dp_dx * DELTAT / dens(i,j)

       enddo
     enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,value_dp_dz)
!$OMP DO SCHEDULE(DYNAMIC,1)
     do j = 1,Nz
       do i = 1,NX
 
        value_dp_dz = dlv(1)*(p(i,j+1) - p(i,j)) & 
                    + dlv(2)*(p(i,j+2) - p(i,j-1)) & 
                    + dlv(3)*(p(i,j+3) - p(i,j-2)) & 
                    + dlv(4)*(p(i,j+4) - p(i,j-3)) 

        memory_dp_dz(i,j) = b_z_half(i,j) * memory_dp_dz(i,j) &
                          + a_z_half(i,j) * value_dp_dz

        value_dp_dz = value_dp_dz / K_z_half(j) + memory_dp_dz(i,j)

        vz(i,j) = vz(i,j) - value_dp_dz * DELTAT / dens(i,j)

       enddo
     enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,l,sum_r,value_dvx_dx,value_dvz_dz)
!$OMP DO SCHEDULE(DYNAMIC,1)
     do j = 1,nz
        do i = 1,nx

          value_dvx_dx = dlh(1)*(vx(i,j) - vx(i-1,j)) &
                     +dlh(2)*(vx(i+1,j) - vx(i-2,j)) &
                     +dlh(3)*(vx(i+2,j) - vx(i-3,j)) &
                     +dlh(4)*(vx(i+3,j) - vx(i-4,j)) 

          value_dvz_dz = dlv(1)*(vz(i,j) - vz(i,j-1)) &
                     +dlv(2)*(vz(i,j+1) - vz(i,j-2)) &
                     +dlv(3)*(vz(i,j+2) - vz(i,j-3)) &
                     +dlv(4)*(vz(i,j+3) - vz(i,j-4)) 

         memory_dvx_dx(i,j) = b_x(i,j) * memory_dvx_dx(i,j) &
                            + a_x(i,j) * value_dvx_dx
         memory_dvz_dz(i,j) = b_z(i,j) * memory_dvz_dz(i,j) &
                            + a_z(i,j) * value_dvz_dz

         value_dvx_dx = value_dvx_dx / K_x(i) + memory_dvx_dx(i,j)
         value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)

!         epsilon_m(i,j)=epsilon_m(i,j)-deltat*(value_dvx_dx+value_dvz_dz)

         p(i,j) = p(i,j) + c11_u(i,j) * (value_dvx_dx + value_dvz_dz)  
 
!         p(i,j)=-c11_u(i,j)/deltat*epsilon_m(i,j)

        if(.false.) then
!------------------------------------------------------------
! compute stress sigma and update memory variables for C-PML
!------------------------------------------------------------
         sum_r = zero

        do  l = 1, nrel
!       assign half of the sums from the previous time step (half of the
!       memory variables are stored in each r)
           sum_r = sum_r + 0.5e0*r(l,i,j)
        enddo 

        do  l = 1, nrel
!        update the stiff system of memory variables by the
!        Crank-Nicholson (time-averaged) scheme.
!        Half of each memory variable is updated and used
!           if(modeling) then
!              x1 = 1.0/(1.0 + deltat*0.50*tau_sigma(l,i,j))
!              x2 = 1.0 - deltat*0.50*tau_sigma(l,i,j)
!             else
!              x1 = 1.0/(1.0 + deltat*0.50*tau_sigma(l,i,j))
!              x2 = 1.0 - deltat*0.50*tau_sigma(l,i,j)
!           endif

             r(l,i,j) = x1(l,i,j)*(x2(l,i,j)*r(l,i,j) + tau11(l,i,j)&
                      * (value_dvx_dx + value_dvz_dz) )

             sum_r = sum_r + 0.5e0*r(l,i,j)

         enddo
 
         p(i,j) = p(i,j) - deltat*sum_r

        endif

       enddo
      enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    end subroutine iso_visco_step
