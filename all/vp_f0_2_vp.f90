
        subroutine vp_f0_2_vp(nx,nz,lv,nrel,f_ref,c11_r,vp_f0,vp,rho,&
                              tau11_eps,tau11_sigma)

         implicit none
         integer::nx,nz,lv,nrel,i,j,l
         real,dimension(1:nx,1:nz)::vp_f0,vp,c11_r,tmp_11,&
                                    rho,c11_f0_r
         real,dimension(1:nrel,1:nx,1:nz)::tau11_eps,tau11_sigma
         real::f_ref,omega,pi

         pi=3.1415926
         omega=2.e0*pi*f_ref

         tmp_11(:,:)=0.0

         do i=1,nx
            do j=1,nz
               c11_f0_r(i,j)=vp_f0(i,j)*vp_f0(i,j)*rho(i,j)
            enddo 
          enddo

          do i=1,nx
             do j=1,nz
                do l=1,nrel
                   tmp_11(i,j)=tmp_11(i,j)+(1.e0+omega*omega&
                               *tau11_eps(l,i,j)*tau11_sigma(l,i,j))/&
                               (1.e0+omega*omega&
                               *tau11_sigma(l,i,j)*tau11_sigma(l,i,j))
                enddo 
             enddo  
           enddo 

           do i=1,nx 
              do j=1,nz
                 c11_r(i,j)=c11_f0_r(i,j)*real(nrel)/tmp_11(i,j)
                 vp(i,j)=sqrt(c11_r(i,j)/rho(i,j))
              enddo 
           enddo 

           end subroutine vp_f0_2_vp
              
      
