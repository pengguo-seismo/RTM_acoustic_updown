     module wave_propagator
       use precision_m
       use parameter_input
       use main_module

       implicit none

      contains 

       subroutine iso_visco_step_p

        real (fp_kind) :: value_dvx_dx, &
                          value_dvz_dz, sum_r

        integer::num_threads, omp_get_num_threads

!$acc data &
!$acc  copyout(sisvx, sisvy), &
!$acc  create( vx_1,vx_2,vx_3,&
!$acc          vy_1,vy_2,vy_3,&
!$acc          vz_1,vz_2,vz_3,&
!$acc          sigmaxx_1,sigmaxx_2,sigmaxx_3,&
!$acc          sigmayy_1,sigmayy_2,sigmayy_3,&
!$acc          sigmazz_1,sigmazz_2,sigmazz_3,&
!$acc          sigmaxy_1,sigmaxy_2,&
!$acc          sigmaxz_1,sigmaxz_3,&
!$acc          sigmayz_2,sigmayz_3), &
!$acc  copyin(dz_half_over_two,dx_over_two,dy_half_over_two, &
!$acc         dx_half_over_two,dy_over_two,dz_over_two)

     do j = 1,nz
        do i = 1,nx

          value_dvx_dx = dlh(1)*(vx(i,j) - vx(i-1,j)) &
                       + dlh(2)*(vx(i+1,j) - vx(i-2,j)) &
                       + dlh(3)*(vx(i+2,j) - vx(i-3,j)) &
                       + dlh(4)*(vx(i+3,j) - vx(i-4,j))

          value_dvz_dz = dlv(1)*(vz(i,j) - vz(i,j-1)) &
                       + dlv(2)*(vz(i,j+1) - vz(i,j-2)) &
                       + dlv(3)*(vz(i,j+2) - vz(i,j-3)) &
                       + dlv(4)*(vz(i,j+3) - vz(i,j-4))

         memory_dvx_dx(i,j) = b_x(i,j) * memory_dvx_dx(i,j) + a_x(i,j) * value_dvx_dx
         memory_dvz_dz(i,j) = b_z(i,j) * memory_dvz_dz(i,j) + a_z(i,j) * value_dvz_dz

         value_dvx_dx = value_dvx_dx / K_x(i) + memory_dvx_dx(i,j)
         value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)

!         epsilon_m(i,j)=epsilon_m(i,j)-deltat*(value_dvx_dx+value_dvz_dz)

         p(i,j) = p(i,j) + c11_u(i,j) * (value_dvx_dx + value_dvz_dz)

!         p(i,j)=-c11_u(i,j)/deltat*epsilon_m(i,j)

!------------------------------------------------------------
! compute stress sigma and update memory variables for C-PML
!------------------------------------------------------------
         sum_r = zero

!! for Q modeling
      if(.false.) then

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

             r(l,i,j) = x1(l,i,j)*(x2(l,i,j)*r(l,i,j) + tau11(l,i,j)*(value_dvx_dx + &
                       value_dvz_dz) )

             sum_r = sum_r + 0.5e0*r(l,i,j)

         enddo
       endif

         p(i,j) = p(i,j) - deltat*sum_r
       enddo
      enddo

    end subroutine iso_visco_step_p

    subroutine iso_visco_step_vxz

      use omp_lib
      real (fp_kind) :: value_dp_dx, &
                        value_dp_dz    

      integer::num_threads

     do j = 1,nz
       do i = 1,nx

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


    do j = 1,nz
       do i = 1,nx

        value_dp_dz = dlv(1)*(p(i,j+1) - p(i,j)) &
                         +dlv(2)*(p(i,j+2) - p(i,j-1)) &
                         +dlv(3)*(p(i,j+3) - p(i,j-2)) &
                         +dlv(4)*(p(i,j+4) - p(i,j-3))

        memory_dp_dz(i,j) = b_z_half(i,j) * memory_dp_dz(i,j) &
                          + a_z_half(i,j) * value_dp_dz

        value_dp_dz = value_dp_dz / K_z_half(j) + memory_dp_dz(i,j)

        vz(i,j) = vz(i,j) - value_dp_dz * DELTAT / dens(i,j)
       enddo
     enddo

     end subroutine iso_visco_step_vxz

     end module wave_propagator
