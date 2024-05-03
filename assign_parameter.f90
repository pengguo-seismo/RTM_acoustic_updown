
!!  define intermediate parameters for the visco-acoustic propagator!!

! Coded by : Peng Guo
! Date : January 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! updated, July 2016
! --------------------------------------------------------------------

    subroutine assign_parameter(nx,nz,lv,nrel,deltat,&
                                tau11_eps,tau11,tau11_new,tau_sigma,x1,x2,&
                                c11_r,dens,c11_u)

      implicit none
      integer::nx,nz,lv,nrel
      integer::i,j,l
      integer,parameter::dp=kind(0.e0)
      real(dp)::tau11_eps(1:nrel,1:nx,1:nz),tau11_eps_tmp(1:nrel,1:nx,1:nz)
      real(dp)::tau11(1:nrel,1:nx,1:nz),tau11_new(1:nrel,1:nx,1:nz)
      real(dp)::tau11_sum(1:nx,1:nz)
      real(dp)::tau_sigma(1:nrel,1:nx,1:nz),tau_sigma_tmp(1:nrel,1:nx,1:nz)
      real(dp),dimension(1:nrel,1:nx,1:nz)::x1,x2
      real(dp)::c11_u(1:nx,1:nz),c11_r(1:nx,1:nz)
      real(dp)::dens(1:nx,1:nz)
      real(dp)::deltat

       do i=1,nx
          do j=1,nz
             do l=1,nrel
                tau11_eps_tmp(l,i,j)=tau11_eps(l,i,j)/tau_sigma(l,i,j)-1.0
             enddo 
          enddo 
       enddo 

       do i=1,nx
          do j=1,nz
             do l = 1, nrel
                tau_sigma_tmp(l,i,j) = 1.0/tau_sigma(l,i,j)
             enddo 
          enddo
       enddo 
    
       do i=1,nx
          do j=1,nz
             do l=1,nrel
                x1(l,i,j)=1.0/(1.0 + deltat*0.50*tau_sigma_tmp(l,i,j))
                x2(l,i,j)=1.0 - deltat*0.50*tau_sigma_tmp(l,i,j)
             enddo 
          enddo 
       enddo 

       tau11_sum=0.0

       do i=1,nx
          do j=1,nz
             do l = 1, nrel
                tau11_sum(i,j) = tau11_sum(i,j) + tau11_eps_tmp(l,i,j)     
             enddo  
          enddo 
       enddo 

       do i=1,nx
          do j=1,nz
            c11_u(i,j)=c11_r(i,j)*(1.0+1.0/real(nrel)*tau11_sum(i,j))
          enddo 
       enddo 

       do j=1,nz
          do i=1,nx
             do l=1,nrel
                tau11(l,i,j) = - tau11_eps_tmp(l,i,j)*c11_r(i,j)/real(nrel)&
                             * tau_sigma_tmp(l,i,j)*deltat

                tau11_new(l,i,j) = tau11_eps_tmp(l,i,j)*c11_r(i,j)/real(nrel)
             enddo 
!           Store unrelaxed modulii in their variables and include dt
                c11_u(i,j) = - c11_u(i,j)*deltat
          enddo 
       enddo 

      end
