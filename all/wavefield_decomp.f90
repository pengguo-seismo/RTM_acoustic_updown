!! wavefield decomposition code 
!! apply globally
        subroutine wavefield_decomp(p,p_hb,p_up,p_dw,p_lf,&
                   p_rt,p_lu,p_ld,p_ru,p_rd,taper_x,taper_z,&
                   nx,nz,lv,nxfft,nzfft,flag)

          use fftw_module_f
          implicit none
          include "fftw3.f"
         

          integer::nx,nz,nxfft,nzfft,lv,flag
!! up/down  &&  left/right
          real,dimension(-lv:nx+lv,-lv:nz+lv)::p,p_hb,p_up,p_dw,p_lf,p_rt
!! 
          real,dimension(-lv:nx+lv,-lv:nz+lv)::p_lu,p_ld,p_ru,p_rd           
  
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
!!  up/down decomposition
!!  begin do i=1,nx
         do i=1,nx
            tmp_in_ud(:)=cmplx(0.0,0.0)
            tmp_out_ud(:)=cmplx(0.0,0.0)
            do j=1,nz
               tmp_in_ud(j)=cmplx(p(i,j),p_hb(i,j))
            enddo

            call sfftw_execute(plan_ud,tmp_in_ud,tmp_out_ud)

!! up going wave
            do j=1,nzfft/2
               tmp_in_up(j)=tmp_out_ud(j)*taper_z(j)
            enddo
            do j=nzfft/2+1,nzfft
               tmp_in_up(j)=cmplx(0,0)
            enddo
!! down going wave
            do j=1,nzfft/2
               tmp_in_dw(j)=cmplx(0,0)
            enddo
            do j=nzfft/2+1,nzfft
               tmp_in_dw(j)=tmp_out_ud(j)*taper_z(j)
            enddo

            call sfftw_execute(plan_ud1,tmp_in_up,tmp_out_up)
            call sfftw_execute(plan_ud2,tmp_in_dw,tmp_out_dw)  

            do j=1,nz
               p_up(i,j)=tmp_out_up(j)/real(nzfft)
               p_dw(i,j)=tmp_out_dw(j)/real(nzfft)
            enddo 
!!  enddo i=1,nx 
          enddo 

        endif

        if(flag.eq.2) then

!! left/right decomposition 
          do j=1,nz
 
             tmp_in_lr=cmplx(0.0,0.0) 
             tmp_out_lr=cmplx(0.0,0.0) 
             do i=1,nx
                tmp_in_lr(i)=cmplx(p(i,j),p_hb(i,j))
             enddo 
         
             call sfftw_execute(plan_lr,tmp_in_lr,tmp_out_lr)

             do i=1,nxfft/2
                tmp_in_rt(i)=tmp_out_lr(i)*taper_x(i)
             enddo

             do i=nxfft/2+1,nxfft
                tmp_in_rt(i)=cmplx(0,0)
             enddo

             do i=1,nxfft/2
                tmp_in_lf(i)=cmplx(0,0)
             enddo
             do i=nxfft/2+1,nxfft
                tmp_in_lf(i)=tmp_out_lr(i)*taper_x(i)
             enddo

             call sfftw_execute(plan_lr1,tmp_in_lf,tmp_out_lf)
             call sfftw_execute(plan_lr2,tmp_in_rt,tmp_out_rt)  

             do i=1,nx
                p_lf(i,j)=tmp_out_lf(i)/real(nxfft)
                p_rt(i,j)=tmp_out_rt(i)/real(nxfft)
             enddo 
!!  enddo j=1,nz 
          enddo 

         endif

!!  decompose the wavefield into several dominant angles           
         if (flag.eq.3) then 

             do j=1,nz
                do i=1,nx
                   tmp_in_bin(i,j)=cmplx(p(i,j),p_hb(i,j))
                enddo 
             enddo 

!! 2D fft this time
             call sfftw_execute(plan_bin,tmp_in_bin,tmp_out_bin)

!! quadrant 1 
             do j=1,nzfft/2
                do i=1,nxfft/2
                   tmp_in_rd(i,j)=tmp_out_bin(i,j)*taper_x(i)*taper_z(j)
!! other directions
                   tmp_in_ru(i,j)=0.e0
                   tmp_in_ld(i,j)=0.e0
                   tmp_in_lu(i,j)=0.e0
                 enddo 
             enddo 
!! quadrant 2 
             do j=1,nzfft/2
                do i=nxfft/2+1,nxfft
                   tmp_in_ld(i,j)=tmp_out_bin(i,j)*taper_x(i)*taper_z(j)
!! other directions
                   tmp_in_ru(i,j)=0.e0
                   tmp_in_lu(i,j)=0.e0
                   tmp_in_rd(i,j)=0.e0
                 enddo 
              enddo 
!! quadrant 3 
               do j=nzfft/2+1,nzfft
                  do i=nxfft/2+1,nxfft
                     tmp_in_lu(i,j)=tmp_out_bin(i,j)*taper_x(i)*taper_z(j)
!! other directions
                     tmp_in_ru(i,j)=0.e0
                     tmp_in_ld(i,j)=0.e0
                     tmp_in_rd(i,j)=0.e0
                   enddo 
                enddo 
!! quadrant 4 
               do j=nzfft/2+1,nzfft
                  do i=1,nxfft/2
                     tmp_in_ru(i,j)=tmp_out_bin(i,j)*taper_x(i)*taper_z(j)
!! other directions
                     tmp_in_lu(i,j)=0.e0
                     tmp_in_ld(i,j)=0.e0
                     tmp_in_rd(i,j)=0.e0
                   enddo 
                enddo 

!! inverse fft
                call sfftw_execute(plan_bin1,tmp_in_rd,tmp_out_rd)
                call sfftw_execute(plan_bin2,tmp_in_ld,tmp_out_ld)
                call sfftw_execute(plan_bin3,tmp_in_lu,tmp_out_lu)
                call sfftw_execute(plan_bin4,tmp_in_ru,tmp_out_ru)

                do j=1,nz
                   do i=1,nx
                      p_rd(i,j)=tmp_out_rd(i,j)/real(nxfft*nzfft) 
                      p_ld(i,j)=tmp_out_ld(i,j)/real(nxfft*nzfft) 
                      p_lu(i,j)=tmp_out_lu(i,j)/real(nxfft*nzfft) 
                      p_ru(i,j)=tmp_out_ru(i,j)/real(nxfft*nzfft)
                   enddo 
                enddo 

              endif 

             end subroutine wavefield_decomp 
