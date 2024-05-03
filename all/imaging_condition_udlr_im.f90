
!! apply cross-correlation imaging condition
!! Kai Gao's expanded abstract
! --------------------------------------------------------------------

       subroutine imaging_condition_udlr_im(sou_p,&
        rcv_p,taper_x,taper_z,image_duud,image_dduu,nx,nz,lv, &
        nxfft,nzfft,flag,it)

          use fftw_module_f

          implicit none
          integer::nx,nz,lv
          integer::nxfft,nzfft
          real,dimension(-lv:nx+lv,-lv:nz+lv)::sou_p,rcv_p
          real,dimension(1:nx,1:nz)::image_ud,image_du,image_duud
          real,dimension(1:nx,1:nz)::image_uu,image_dd,image_dduu
          real,dimension(-lv:nx+lv,-lv:nz+lv)::sou_p_hb,rcv_p_hb
          integer::i,j,flag
!! hilbert transform
          real,dimension(1:nzfft)::taper_z
          real,dimension(1:nxfft)::taper_x
          integer::it
          character *256::snap_part
!          complex,dimension(1:nzfft)::tmp_in_ud,tmp_out_ud
!          complex,dimension(1:nxfft)::tmp_in_lr,tmp_out_lr

         if(flag.eq.1) then 
!! source side
          do i=1,nx

             tmp_in_ud(:)=cmplx(0.0,0.0)
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

             sou_p_hb(i,1:nz)=real(tmp_in_ud(1:nz)/nzfft)

          enddo

          do i=1,nx

             tmp_in_ud(:)=cmplx(0.0,0.0)
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

             rcv_p_hb(i,1:nz)=real(tmp_in_ud(1:nz)/nzfft)

          enddo

         else if(flag.eq.2) then 
!! left and right
!! source side
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

             sou_p_hb(1:nx,j)=real(tmp_in_lr(1:nx)/nxfft)

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

             rcv_p_hb(1:nx,j)=real(tmp_in_lr(1:nx)/nxfft)

          enddo

         endif

          if(mod(it,50000).eq.0) then
              snap_part='snap_sou_p'
              call snap_shot(sou_p,it,nx,nz,lv,500,snap_part,12)
              snap_part='snap_rcv_p'
              call snap_shot(rcv_p,it,nx,nz,lv,500,snap_part,12)
              snap_part='snap_sou_ph'
              call snap_shot(sou_p_hb,it,nx,nz,lv,500,snap_part,12)
              snap_part='snap_rcv_ph'
              call snap_shot(rcv_p_hb,it,nx,nz,lv,500,snap_part,12)
           endif


          do j=1,nz
             do i=1,nx
                image_duud(i,j)=image_duud(i,j)+sou_p(i,j)*rcv_p(i,j)-&
                                sou_p_hb(i,j)*rcv_p_hb(i,j) 
!                image_dduu(i,j)=image_dduu(i,j)+sou_p_up(i,j)*rcv_p_up(i,j)&
!                               +sou_p_dw(i,j)*rcv_p_dw(i,j)
             enddo 
          enddo 

!          image_duud(:,:)=2.0*image_duud(:,:)

        end subroutine imaging_condition_udlr_im
          
