
!      author: Peng Guo
       subroutine bvr_record_wave_p(it,nt,vx,vz,p,p_l,p_r,p_t,p_b,&
                  mx,mz,lv,npoints_pml)

          implicit none 
          integer::mx,mz,lv,npoints_pml,nt,nxpmls,nzpmls

          real,dimension(1:mz,1:lv,1:nt)::vx_l,vx_r,vz_l,vz_r,p_l,p_r
          real,dimension(1:mx,1:lv,1:nt)::vx_t,vx_b,vz_t,vz_b,p_t,p_b

          integer::i,j,it

          real::p(-lv:mx+lv,-lv:mz+lv),vx(-lv:mx+lv,-lv:mz+lv),&
                vz(-lv:mx+lv,-lv:mz+lv)

!          vx_vz=.false.

           nxpmls=npoints_pml
           nzpmls=npoints_pml

           do i=1,mz
              do j=1,lv
               p_l(i,j,it)=p(nxpmls+j,i)
               p_r(i,j,it)=p(mx-nxpmls-lv+j,i)
              enddo 
           enddo 

           do i=1,mx
              do j=1,lv
               p_t(i,j,it)=p(i,nzpmls+j)
               p_b(i,j,it)=p(i,mz-nzpmls-lv+j)
              enddo 
           enddo 

       end subroutine bvr_record_wave_p
