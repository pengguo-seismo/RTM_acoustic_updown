
        subroutine adding_abs_val(vp_f0,dens,tau11_eps,tau_sigma,&
                   nx,nz,nrel,npoints_pml)
          implicit none
          integer::nx,nz,nrel,npoints_pml
          integer::i,j,l
          real,dimension(1:nx,1:nz)::vp_f0,dens
          real,dimension(1:nrel,1:nx,1:nz)::tau11_eps,tau_sigma

          do i=npoints_pml+1,nx-npoints_pml
             do j=1,npoints_pml
!! top 
                vp_f0(i,j)=vp_f0(i,npoints_pml+1)
                dens(i,j)=dens(i,npoints_pml+1)
             enddo 
             do j=nz-npoints_pml+1,nz
!! bottom
                vp_f0(i,j)=vp_f0(i,nz-npoints_pml)
                dens(i,j)=dens(i,nz-npoints_pml)
             enddo 
          enddo 
       
          do j=1,nz
             do i=1,npoints_pml
!! left 
                vp_f0(i,j)=vp_f0(npoints_pml+1,j)
                dens(i,j)=dens(npoints_pml+1,j)
             enddo 
             do i=nx-npoints_pml+1,nx
!! right
                vp_f0(i,j)=vp_f0(nx-npoints_pml,j) 
                dens(i,j)=dens(nx-npoints_pml,j) 
             enddo 
           enddo 

!! for tau11_eps,tau_sigma
          do l=1,nrel
             do i=npoints_pml+1,nx-npoints_pml
                do j=1,npoints_pml
!! top 
                   tau11_eps(l,i,j)=tau11_eps(l,i,npoints_pml+1)
                   tau_sigma(l,i,j)=tau_sigma(l,i,npoints_pml+1)
                enddo 
                do j=nz-npoints_pml+1,nz
!! bottom
                   tau11_eps(l,i,j)=tau11_eps(l,i,nz-npoints_pml)
                   tau_sigma(l,i,j)=tau_sigma(l,i,nz-npoints_pml)
                enddo 
             enddo 
          enddo 

          do l=1,nrel
             do j=1,nz
               do i=1,npoints_pml
!! left 
                  tau11_eps(l,i,j)=tau11_eps(l,npoints_pml+1,j)
                  tau_sigma(l,i,j)=tau_sigma(l,npoints_pml+1,j)
               enddo
              
               do i=nx-npoints_pml+1,nx
!! right
                  tau11_eps(l,i,j)=tau11_eps(l,nx-npoints_pml,j) 
                  tau_sigma(l,i,j)=tau_sigma(l,nx-npoints_pml,j) 
               enddo 
             enddo 
          enddo 

         end subroutine adding_abs_val                           
