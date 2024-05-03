
!! apply cross-correlation imaging condition
!! Tong's Geophysics paper
! --------------------------------------------------------------------

       subroutine imaging_condition_udlr_eff(sou_p,&
        rcv_p,rcv_p_hb_t,taper_x,taper_z,image,image_du,image_ud,&
        image_dd,image_uu,image_duud,image_dduu,nx,nz,lv,&
        nxfft,nzfft,flag)

          use fftw_module_f

          implicit none
          integer::nx,nz,lv,flag
          integer::nxfft,nzfft
          real,dimension(-lv:nx+lv,-lv:nz+lv)::sou_p,rcv_p,rcv_p_hb_t
          real,dimension(1:nx,1:nz)::image,image_ud,image_du,image_duud
          real,dimension(1:nx,1:nz)::image_uu,image_dd,image_dduu
          integer::i,j

          real,dimension(1:nx,1:nz)::sou_p_hb_spa,rcv_p_hb_spa,&
                                     rcv_p_hb_t_spa

          real,dimension(1:nzfft)::taper_z
          real,dimension(1:nxfft)::taper_x

          real::cons

!! hilbert transform
!          complex,dimension(1:nzfft)::tmp_in_ud,tmp_out_ud
!          complex,dimension(1:nxfft)::tmp_in_lr,tmp_out_lr

         if(flag.eq.1) then 
!! source side
          do i=1,nx

             tmp_in_ud(:)=cmplx(0.0,0.0)
             tmp_out_ud(:)=cmplx(0.0,0.0)
             tmp_in_ud(1:nz)=cmplx(sou_p(i,1:nz),0.0)
              
             call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

             do j=1,nzfft/2
                tmp_out_ud(j)=cmplx(aimag(tmp_out_ud(j)),-real(tmp_out_ud(j)))*taper_z(j)
             enddo
             do j=nzfft/2+1,nzfft
                tmp_out_ud(j)=cmplx(-aimag(tmp_out_ud(j)),real(tmp_out_ud(j)))*taper_z(j)
             enddo

             tmp_in_ud=cmplx(0.0,0.0)

             call sfftw_execute(plan_ud_inv,tmp_out_ud,tmp_in_ud)

             sou_p_hb_spa(i,1:nz)=real(tmp_in_ud(1:nz)/nzfft)

          enddo
!! receiver side
          do i=1,nx

             tmp_in_ud(:)=cmplx(0.0,0.0)
             tmp_out_ud(:)=cmplx(0.0,0.0)
             tmp_in_ud(1:nz)=cmplx(rcv_p(i,1:nz),0.0)
              
             call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

             do j=1,nzfft/2
                tmp_out_ud(j)=cmplx(aimag(tmp_out_ud(j)),-real(tmp_out_ud(j)))*taper_z(j)
             enddo
             do j=nzfft/2+1,nzfft
                tmp_out_ud(j)=cmplx(-aimag(tmp_out_ud(j)),real(tmp_out_ud(j)))*taper_z(j)
             enddo

             tmp_in_ud=cmplx(0.0,0.0)

             call sfftw_execute(plan_ud_inv,tmp_out_ud,tmp_in_ud)

             rcv_p_hb_spa(i,1:nz)=real(tmp_in_ud(1:nz)/nzfft)

          enddo

!! hilbert-transformed receiver side 
          do i=1,nx

             tmp_in_ud(:)=cmplx(0.0,0.0)
             tmp_out_ud(:)=cmplx(0.0,0.0)
             tmp_in_ud(1:nz)=cmplx(rcv_p_hb_t(i,1:nz),0.0)
              
             call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

             do j=1,nzfft/2
                tmp_out_ud(j)=cmplx(aimag(tmp_out_ud(j)),-real(tmp_out_ud(j)))*taper_z(j)
             enddo
             do j=nzfft/2+1,nzfft
                tmp_out_ud(j)=cmplx(-aimag(tmp_out_ud(j)),real(tmp_out_ud(j)))*taper_z(j)
             enddo

             tmp_in_ud=cmplx(0.0,0.0)

             call sfftw_execute(plan_ud_inv,tmp_out_ud,tmp_in_ud)

             rcv_p_hb_t_spa(i,1:nz)=real(tmp_in_ud(1:nz)/nzfft)

          enddo

         else if (flag.eq.2) then

          do j=1,nz

             tmp_in_lr(:)=cmplx(0.0,0.0)
             tmp_in_lr(1:nx)=cmplx(sou_p(1:nx,j),0.0)
              
             call sfftw_execute(plan_lr,tmp_in_lr,tmp_out_lr)

             do i=1,nxfft/2
                tmp_out_lr(i)=cmplx(aimag(tmp_out_lr(i)),-real(tmp_out_lr(i)))
             enddo
             do i=nxfft/2+1,nxfft
                tmp_out_lr(i)=cmplx(-aimag(tmp_out_lr(i)),real(tmp_out_lr(i)))
             enddo

             tmp_in_lr=cmplx(0.0,0.0)

             call sfftw_execute(plan_lr_inv,tmp_out_lr,tmp_in_lr)

             sou_p_hb_spa(1:nx,j)=real(tmp_in_lr(1:nx)/nxfft)

          enddo

          do j=1,nz

             tmp_in_lr(:)=cmplx(0.0,0.0)
             tmp_in_lr(1:nx)=cmplx(rcv_p(1:nx,j),0.0)
              
             call sfftw_execute(plan_lr,tmp_in_lr,tmp_out_lr)

             do i=1,nxfft/2
                tmp_out_lr(i)=cmplx(aimag(tmp_out_lr(i)),-real(tmp_out_lr(i)))
             enddo
             do i=nxfft/2+1,nxfft
                tmp_out_lr(i)=cmplx(-aimag(tmp_out_lr(i)),real(tmp_out_lr(i)))
             enddo

             tmp_in_lr=cmplx(0.0,0.0)

             call sfftw_execute(plan_lr_inv,tmp_out_lr,tmp_in_lr)

             rcv_p_hb_spa(1:nx,j)=real(tmp_in_lr(1:nx)/nxfft)

          enddo

          do j=1,nz

             tmp_in_lr(:)=cmplx(0.0,0.0)
             tmp_in_lr(1:nx)=cmplx(rcv_p_hb_t(1:nx,j),0.0)
              
             call sfftw_execute(plan_lr,tmp_in_lr,tmp_out_lr)

             do i=1,nxfft/2
                tmp_out_lr(i)=cmplx(aimag(tmp_out_lr(i)),-real(tmp_out_lr(i)))
             enddo
             do i=nxfft/2+1,nxfft
                tmp_out_lr(i)=cmplx(-aimag(tmp_out_lr(i)),real(tmp_out_lr(i)))
             enddo

             tmp_in_lr=cmplx(0.0,0.0)

             call sfftw_execute(plan_lr_inv,tmp_out_lr,tmp_in_lr)

             rcv_p_hb_t_spa(1:nx,j)=real(tmp_in_lr(1:nx)/nxfft)

          enddo

         endif

         cons=0.25

          do j=1,nz
             do i=1,nx
                image(i,j)=image(i,j)+(sou_p(i,j)*rcv_p(i,j))

                image_du(i,j)=image_du(i,j)+cons*(sou_p(i,j)*rcv_p(i,j)-&
                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)-&
                sou_p(i,j)*rcv_p_hb_t_spa(i,j)-&
                sou_p_hb_spa(i,j)*rcv_p_hb_t(i,j))

                image_ud(i,j)=image_ud(i,j)+cons*(sou_p(i,j)*rcv_p(i,j)-&
                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)+&
                sou_p(i,j)*rcv_p_hb_t_spa(i,j)+&
                sou_p_hb_spa(i,j)*rcv_p_hb_t(i,j))

                image_duud(i,j)=image_duud(i,j)+cons*(2.0*(sou_p(i,j)*rcv_p(i,j)-&
                                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)))

                image_dd(i,j)=image_dd(i,j)+cons*(sou_p(i,j)*rcv_p(i,j)+&
                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)+&
                sou_p(i,j)*rcv_p_hb_t_spa(i,j)-&
                sou_p_hb_spa(i,j)*rcv_p_hb_t(i,j))
 
                image_uu(i,j)=image_uu(i,j)+cons*(sou_p(i,j)*rcv_p(i,j)+&
                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)-&
                sou_p(i,j)*rcv_p_hb_t_spa(i,j)+&
                sou_p_hb_spa(i,j)*rcv_p_hb_t(i,j))

                image_dduu(i,j)=image_dduu(i,j)+cons*(2.0*(sou_p(i,j)*rcv_p(i,j)+&
                sou_p_hb_spa(i,j)*rcv_p_hb_spa(i,j)))
             enddo 
          enddo 

        end subroutine imaging_condition_udlr_eff
          
