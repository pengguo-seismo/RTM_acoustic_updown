!! wavefield decomposition code 
!! apply globally
        subroutine wavefield_decomp2(p,p_hb,sou_p,sou_p_hb,&
                   image_du,taper_x,taper_z,&
                   nx,nz,lv,nxfft,nzfft,flag)

          use fftw_module_f
          implicit none
          include "fftw3.f"
         
          integer::nx,nz,nxfft,nzfft,lv,flag
!! up/down  &&  left/right
          real,dimension(-lv:nx+lv,-lv:nz+lv)::p,p_hb,sou_p,sou_p_hb
!! 
          real,dimension(-lv:nx+lv,-lv:nz+lv)::p_lu,p_ld,p_ru,p_rd           
          real,dimension(1:nx,1:nz)::image_du 
          complex,dimension(1:nx,1:nz)::sou_p_hb_t_spa,p_hb_t_spa
 
!          complex,dimension(1:nzfft)::tmp_in_ud,tmp_out_ud,&
!                                   tmp_in_dw,tmp_out_dw,&         
!                                   tmp_in_up,tmp_out_up      
! 
!          complex,dimension(1:nxfft)::tmp_in_lr,tmp_out_lr,&
!                                   tmp_in_lf,tmp_out_lf,&         
!                                   tmp_in_rt,tmp_out_rt      
!
!          complex,dimension(1:nxfft,1:nzfft)::tmp_in_bin,tmp_out_bin,&
!                       tmp_in_lu,tmp_out_lu,tmp_in_ru,tmp_out_ru,&
!                       tmp_in_ld,tmp_out_ld,tmp_in_rd,tmp_out_rd

          real,dimension(1:nxfft)::taper_x
          real,dimension(1:nzfft)::taper_z

          integer::i,j
          real::pi,tp

       if(flag.eq.1) then 
!!  begin do i=1,nx
         do i=1,nx
            tmp_in_ud(:)=cmplx(0.0,0.0)
            tmp_out_ud(:)=cmplx(0.0,0.0)
            do j=1,nz
               tmp_in_ud(j)=cmplx(p(i,j),p_hb(i,j))
            enddo

            call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

            do j=2,nzfft/2
               tmp_out_ud(j)=2.0*cmplx(real(tmp_out_ud(j)),aimag(tmp_out_ud(j)))
            enddo
            do j=nzfft/2+2,nzfft
               tmp_out_ud(j)=cmplx(0.0,0.0)
            enddo
          if(.false.) then 
            do j=nzfft/2+2,nzfft
               tmp_out_ud(j)=2.0*cmplx(real(tmp_out_ud(j)),aimag(tmp_out_ud(j)))
            enddo
            do j=2,nzfft/2
               tmp_out_ud(j)=cmplx(0.0,0.0)
            enddo
!             do j=1,nzfft/2
          endif
!                tmp_out_ud(j)=cmplx(aimag(tmp_out_ud(j)),-real(tmp_out_ud(j)))
!             enddo
!             do j=nzfft/2+1,nzfft
!                tmp_out_ud(j)=cmplx(-aimag(tmp_out_ud(j)),real(tmp_out_ud(j)))
!             enddo

            tmp_in_ud=cmplx(0.0,0.0)
            call sfftw_execute(plan_ud_inv,tmp_out_ud,tmp_in_ud)

            do j=1,nz
               p_hb_t_spa(i,j)=tmp_in_ud(j)/real(nzfft)
            enddo 
!!  enddo i=1,nx 
          enddo 

!!  begin do i=1,nx
         do i=1,nx
            tmp_in_ud(:)=cmplx(0.0,0.0)
            tmp_out_ud(:)=cmplx(0.0,0.0)
            do j=1,nz
               tmp_in_ud(j)=cmplx(sou_p(i,j),sou_p_hb(i,j))
            enddo

            call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

!             do j=1,nzfft/2
!                tmp_out_ud(j)=cmplx(aimag(tmp_out_ud(j)),-real(tmp_out_ud(j)))
!             enddo
!             do j=nzfft/2+1,nzfft
!                tmp_out_ud(j)=cmplx(-aimag(tmp_out_ud(j)),real(tmp_out_ud(j)))
!             enddo

          if(.false.) then
            do j=2,nzfft/2
               tmp_out_ud(j)=2.0*cmplx(real(tmp_out_ud(j)),aimag(tmp_out_ud(j)))
            enddo
            do j=nzfft/2+2,nzfft
               tmp_out_ud(j)=cmplx(0.0,0.0)
            enddo
          endif
            do j=2,nzfft/2
               tmp_out_ud(j)=cmplx(0.0,0.0)
            enddo
            do j=nzfft/2+2,nzfft
               tmp_out_ud(j)=2.0*cmplx(real(tmp_out_ud(j)),aimag(tmp_out_ud(j)))
            enddo

            tmp_in_ud=cmplx(0.0,0.0)
            call sfftw_execute(plan_ud_inv,tmp_out_ud,tmp_in_ud)

            do j=1,nz
               sou_p_hb_t_spa(i,j)=tmp_in_ud(j)/real(nzfft)
            enddo
         enddo 
!!  enddo i=1,nx 
        endif

           do j=1,nz
              do i=1,nx
                 image_du(i,j)=image_du(i,j)+real(sou_p_hb_t_spa(i,j)*p_hb_t_spa(i,j))
              enddo 
           enddo  

        end subroutine wavefield_decomp2
