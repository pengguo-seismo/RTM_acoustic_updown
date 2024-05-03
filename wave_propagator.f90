     module wave_propagator
       use precision_m
       use parameter_input
       use main_module
       implicit none

      contains 

       subroutine para_prepare

!!     assign values for pml layers
!!     1:npoints_pml,nx-npoints_pml+1:nx
        call adding_abs_val(vp_f0,dens,tau11_eps,tau_sigma,nx,nz,&
             nrel,npoints_pml)

        open(12,file='vp_test',access='direct',form='unformatted',&
                recl=ii_bi*nz)
         do i=1,nx
            write(12,rec=i) &
                     (vp_f0(i,j),j=1,nz)
         enddo
        close(12)

        call vp_f0_2_vp(nx,nz,lv,nrel,f_ref,c11_r,vp_f0,vp,dens,&
                       tau11_eps,tau_sigma)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---
        c(1)=1.23538628085693269
        c(2)=-0.109160358903484037
        c(3)=0.246384765606012697d-1
        c(4)=-0.660472512788868975d-2

!!!!!  FD coefficients, holberg version
        c(1)=1.231666
        c(2)=-0.1041182
        c(3)=0.02063707
        c(4)=-0.003570998
!!!!!!!!!!!!!!!!! Stability Condition !!!!!!!!!!!!!!
        cons1=sqrt(1.0/deltax**2+1.0/deltaz**2)
!
        sum_a=0.0
!
        do i=1,lv
           sum_a=sum_a+abs(c(i))
        enddo
!
        sum_a=1.0/sum_a
!
        do i=1,nx
           do j=1,nz
              deltat_v=deltat*vp(i,j)
              cons2=deltat_v*cons1
              if(cons2.gt.sum_a) then
                print *,'Stability condition not reached !!!!'
                stop
              endif
           enddo
        enddo

!        print *,'Maximum allowed V', sum_a/(cons1*deltat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,lv
           dlh(i)=c(i)/deltax
           dlv(i)=c(i)/deltaz
        enddo

        call assign_parameter(nx,nz,lv,nrel,deltat,&
                             tau11_eps,tau11,tau11_new,tau_sigma,x1,x2,&
                             c11_r,dens,c11_u)
! compute cpml coefficients !
        call  cpml_coef(a_x,b_x,k_x,a_z,b_z,k_z,&
              a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
              npoints_pml,nx,nz,deltax,deltaz,deltat,vp_f0,f0)

       end subroutine para_prepare

       subroutine iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)

        real (fp_kind) :: value_dvx_dx, &
                          value_dvz_dz, sum_r

        integer::num_threads, omp_get_num_threads
        real,dimension(-lv:nx+lv,-lv:nz+lv)::vx,vz,p
        real,dimension(1:nrel,-lv:nx+lv,-lv:nz+lv)::r

        real,dimension(1:nx,1:nz)::memory_dvx_dx,memory_dvz_dz

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,l,sum_r,value_dvx_dx,value_dvz_dz)
!$OMP DO

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

         p(i,j) = p(i,j) - deltat*sum_r

       enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL 
    end subroutine iso_visco_step_p

    subroutine iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

      use omp_lib
      real (fp_kind) :: value_dp_dx, &
                        value_dp_dz    

      real,dimension(-lv:nx+lv,-lv:nz+lv)::vx,vz,p

      real,dimension(1:nx,1:nz)::memory_dp_dx,memory_dp_dz

      integer::num_threads

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,value_dp_dx)
!$OMP DO
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
!$OMP END DO
!$OMP END PARALLEL 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,value_dp_dz)
!$OMP DO
 
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
!$OMP END DO
!$OMP END PARALLEL 

     end subroutine iso_visco_step_vxz

     end module wave_propagator
