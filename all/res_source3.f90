
!res_source2(res_source_x,res_source_z,&
!                   ux0,uz0,mu_f0,nx,nz,nrel,lv,deltax,deltaz)

        subroutine res_source2(res_source_x,res_source_z,&
                   ux0,uz0,mu_f0,nx,nz,nrel,lv,dx,dz)

        implicit none        
        integer::nx,nz,nrel,lv,i,j
        real,dimension(-lv:nx+lv,-lv:nz+lv)::dila,&
                      value_dux_dx,value_duz_dz,&
                      value_duz_dx,value_dux_dz,&
                      res_source_x,res_source_z,&
                      epsilon_xz,tmp_x1,tmp_z1,tmp

        real,dimension(-lv:nx+lv,-lv:nz+lv)::ux0,uz0

        real,dimension(1:nx,1:nz)::grad_dila_x,grad_dila_z,mu_f0

        real,dimension(1:nx,1:nz)::grad_depsxx_dx,grad_depsxz_dz,&
                                   grad_depszz_dz,grad_depsxz_dx

        real,dimension(1:nx,1:nz)::grad_mu_x,grad_mu_z
        real::mu_x,mu_z

        real::dlh(1:lv),dlv(1:lv)
        real::dx,dz
        integer::ii_bi

        dila=0.0
        value_dux_dx=0.0
        value_duz_dz=0.0
        value_duz_dx=0.0
        value_dux_dz=0.0
        epsilon_xz=0.0

        tmp_x1=0.
        tmp_z1=0.
        tmp=0.

        dlh(:)=0.0
        dlh(1)=1.1382
        dlh(2)=-0.046414

      dlh(1)=1.231666e0
      dlh(2)=-0.1041182
      dlh(3)=0.02063707
      dlh(4)=-0.003570998

!        dlh(2)=-0.2
!        dlh(3)=0.038095238
!        dlh(4)=-0.0035714286

        dlh(:)=dlh(:)/dx
        dlv(:)=dlh(:)

        grad_mu_x=0.0
        grad_mu_z=0.0

        do j=1,nz-1
           do i=1,nx-1
              grad_mu_x(i,j)= dlh(1)*(mu_f0(i+1,j) - mu_f0(i,j)) &
                            + dlh(2)*(mu_f0(i+2,j) - mu_f0(i-1,j)) &
                            + dlh(3)*(mu_f0(i+3,j) - mu_f0(i-2,j)) &
                            + dlh(4)*(mu_f0(i+4,j) - mu_f0(i-3,j))

              grad_mu_z(i,j)= dlh(1)*(mu_f0(i,j+1) - mu_f0(i,j)) &
                            + dlh(2)*(mu_f0(i,j+2) - mu_f0(i,j-1)) &
                            + dlh(3)*(mu_f0(i,j+3) - mu_f0(i,j-2)) &
                            + dlh(4)*(mu_f0(i,j+4) - mu_f0(i,j-3))
           enddo 
        enddo

       ii_bi=1
       open(12,file='grad_mu_testx.dat',access='direct',&
           form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
          write(12,rec=i) (grad_mu_x(i,j),j=1,nz)
       enddo
       close(12)

       open(12,file='grad_mu_testz.dat',access='direct',&
           form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
          write(12,rec=i) (grad_mu_z(i,j),j=1,nz)
       enddo
       close(12)

        do j=1,nz
          do i=1,nx
!! integer points
            value_dux_dx(i,j) = dlh(1)*(ux0(i,j) - ux0(i-1,j)) &
                              + dlh(2)*(ux0(i+1,j) - ux0(i-2,j)) &
                              + dlh(3)*(ux0(i+2,j) - ux0(i-3,j)) &
                              + dlh(4)*(ux0(i+3,j) - ux0(i-4,j))
!! half points
            value_dux_dz(i,j) = dlh(1)*(ux0(i,j+1) - ux0(i,j)) &
                              + dlh(2)*(ux0(i,j+2) - ux0(i,j-1)) &
                              + dlh(3)*(ux0(i,j+3) - ux0(i,j-2)) &
                              + dlh(4)*(ux0(i,j+4) - ux0(i,j-3))
!! integer points
            value_duz_dz(i,j) = dlv(1)*(uz0(i,j) - uz0(i,j-1)) &
                              + dlv(2)*(uz0(i,j+1) - uz0(i,j-2)) &
                              + dlv(3)*(uz0(i,j+2) - uz0(i,j-3)) &
                              + dlv(4)*(uz0(i,j+3) - uz0(i,j-4))
!! half points
            value_duz_dx(i,j) = dlv(1)*(uz0(i+1,j) - uz0(i,j)) &
                              + dlv(2)*(uz0(i+2,j) - uz0(i-1,j)) &
                              + dlv(3)*(uz0(i+3,j) - uz0(i-2,j)) &
                              + dlv(4)*(uz0(i+4,j) - uz0(i-3,j))

            dila(i,j) = value_dux_dx(i,j) + value_duz_dz(i,j)

            epsilon_xz(i,j) = 0.5*(value_dux_dz(i,j)+value_duz_dx(i,j))

!            tmp_11(i,j)= 2.e0*c55_u(i,j)*((value_ux_dx(i,j)-dila(i,j)))
!            tmp_12(i,j)= 2.e0*c55_u(i,j)*0.5e0*(value_ux_dz(i,j)+value_uz_dx(i,j))  
!            tmp_22(i,j)= 2.e0*c55_u(i,j)*((value_uz_dz(i,j)-dila(i,j)))
          enddo
        enddo

       do j=1,nz
          do i=1,nx

!! at half x points.
             grad_depsxx_dx(i,j)= dlh(1)*(value_dux_dx(i+1,j) - value_dux_dx(i,j)) &
                                + dlh(2)*(value_dux_dx(i+2,j) - value_dux_dx(i-1,j)) &
                                + dlh(3)*(value_dux_dx(i+3,j) - value_dux_dx(i-2,j)) &
                                + dlh(4)*(value_dux_dx(i+4,j) - value_dux_dx(i-3,j))
!! at half x points
             grad_depsxz_dz(i,j)= dlh(1)*(epsilon_xz(i,j) - epsilon_xz(i,j-1)) &
                                + dlh(2)*(epsilon_xz(i,j+1) - epsilon_xz(i,j-2)) &
                                + dlh(3)*(epsilon_xz(i,j+2) - epsilon_xz(i,j-3)) &
                                + dlh(4)*(epsilon_xz(i,j+3) - epsilon_xz(i,j-4))

!! at half z points.
             grad_depszz_dz(i,j)= dlv(1)*(value_duz_dz(i,j+1) - value_duz_dz(i,j)) &
                                + dlv(2)*(value_duz_dz(i,j+2) - value_duz_dz(i,j-1)) &
                                + dlv(3)*(value_duz_dz(i,j+3) - value_duz_dz(i,j-2)) &
                                + dlv(4)*(value_duz_dz(i,j+4) - value_duz_dz(i,j-3))

             grad_depsxz_dx(i,j)= dlv(1)*(epsilon_xz(i,j) - epsilon_xz(i-1,j)) &
                                + dlv(2)*(epsilon_xz(i+1,j) - epsilon_xz(i-2,j)) &
                                + dlv(3)*(epsilon_xz(i+2,j) - epsilon_xz(i-3,j)) &
                                + dlv(4)*(epsilon_xz(i+3,j) - epsilon_xz(i-4,j))

!! at half x points
             grad_dila_x(i,j)= dlh(1)*(dila(i+1,j) - dila(i,j)) &
                             + dlh(2)*(dila(i+2,j) - dila(i-1,j)) &
                             + dlh(3)*(dila(i+3,j) - dila(i-2,j)) &
                             + dlh(4)*(dila(i+4,j) - dila(i-3,j))

!! at half z points
             grad_dila_z(i,j)= dlh(1)*(dila(i,j+1) - dila(i,j)) &
                             + dlh(2)*(dila(i,j+2) - dila(i,j-1)) &
                             + dlh(3)*(dila(i,j+3) - dila(i,j-2)) &
                             + dlh(4)*(dila(i,j+4) - dila(i,j-3))
          enddo
       enddo

       do j=1,nz-1
          do i=1,nx-1
             tmp_x1(i,j)=2.0*mu_f0(i,j)*(value_dux_dx(i,j)-dila(i,j))
             tmp_z1(i,j)=2.0*mu_f0(i,j)*(value_duz_dz(i,j)-dila(i,j))
             tmp(i,j)=2.0*mu_f0(i,j)*epsilon_xz(i,j)
          enddo 
       enddo

       do j=40,nz-1
          do i=1,nx-1
             res_source_x(i,j)=dlh(1)*(tmp_x1(i+1,j) - tmp_x1(i,j)) &
                             + dlh(2)*(tmp_x1(i+2,j) - tmp_x1(i-1,j)) &
                             + dlv(1)*(tmp(i,j) - tmp(i,j-1)) &
                             + dlv(2)*(tmp(i,j+1) - tmp(i,j-2))
                                   
             res_source_z(i,j)=dlv(1)*(tmp_z1(i,j+1) - tmp_z1(i,j)) &
                             + dlv(2)*(tmp_z1(i,j+2) - tmp_z1(i,j-1)) &
                             + dlh(1)*(tmp(i,j) - tmp(i-1,j)) &
                             + dlh(2)*(tmp(i+1,j) - tmp(i-2,j))
          enddo 
       enddo

      if(.false.) then

       do j=40,nz-1
          do i=1,nx-1

!             value_dux_dx_half=0.5e0*(value_dux_dx(i+1,j)+value_dux_dx(i,j))
!             dila_half=0.5e0*(dila(i+1,j)+dila(i,j))
             mu_z=0.5e0*(mu_f0(i,j)+mu_f0(i,j+1))
             mu_x=0.5e0*(mu_f0(i,j)+mu_f0(i+1,j))

             mu_z=mu_f0(i,j)
             mu_x=mu_f0(i,j)

!             res_source_x(i,j)=(-value_duz_dz(i,j))*2.0*grad_mu_x(i,j)+&
!                              2.0*epsilon_xz(i,j)*grad_mu_z(i,j)&
!                              +2.e0*mu_x*(grad_depsxx_dx(i,j)-grad_dila_x(i,j)&
!                              +grad_depsxz_dz(i,j))

!             res_source_x(i,j)=            

!             res_source_x(i,j)=&
!                              2.e0*mu_x*(grad_depsxx_dx(i,j)-grad_dila_x(i,j)&
!                              +grad_depsxz_dz(i,j))
!             value_duz_dz_half=0.5e0*(value_duz_dz(i,j+1)+value_duz_dz(i,j))
!             dila_half=0.5e0*(dila(i,j+1)+dila(i,j))

!             res_source_z(i,j)=(-value_dux_dx(i,j))*2.0*grad_mu_z(i,j)+&
!                              2.0*epsilon_xz(i,j)*grad_mu_x(i,j)&
!                              +2.e0*mu_z*(grad_depszz_dz(i,j)-grad_dila_z(i,j)&
!                              +grad_depsxz_dx(i,j))

!             res_source_z(i,j)=&
!                              2.e0*mu_z*(grad_depszz_dz(i,j)-grad_dila_z(i,j)&
!                              +grad_depsxz_dx(i,j))
          enddo
       enddo
  
       endif

      end
