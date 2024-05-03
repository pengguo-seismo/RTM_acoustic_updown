
      subroutine rickerfunc_w_hb(dt,fp,length,nt,ricker,ricker_hb)

       include 'fftw3.f'

       integer length,nt,nnt,nt_fft
       integer,parameter::dp=kind(0.e0)
       real(dp)::dt,fp,tshift
!       parameter(PI=3.1415926)
       real::PI
       real(dp)::ricker(1:nt),ricker_hb(1:nt)
       integer*8::plan1,plan2
       complex,allocatable::ricker_in(:),ricker_out(:)

       PI=4.0*atan(1.0)

       nt_fft=nt*2

       allocate(ricker_in(1:nt_fft),ricker_out(1:nt_fft))

       ricker_in(:)=cmplx(0.0,0.0)
       ricker_out(:)=cmplx(0.0,0.0)

       nnt=nt

        call sfftw_plan_dft_1d &
             (plan1,nt_fft,ricker_in,ricker_out,FFTW_FORWARD,FFTW_ESTIMATE)

        call sfftw_plan_dft_1d &
             (plan2,nt_fft,ricker_out,ricker_in,FFTW_BACKWARD,FFTW_ESTIMATE)

! Hilbert transform
        ricker_in=cmplx(0.0,0.0)

        do i=1,nt
           ricker_in(i)=cmplx(ricker(i),0.0)
        enddo

        call sfftw_execute(plan1,ricker_in,ricker_out)

         do it=1,nt_fft/2
                ricker_out(it)=cmplx(aimag(ricker_out(it)),-real(ricker_out(it)))
         enddo
         do it=nt_fft/2+1,nt_fft
                ricker_out(it)=cmplx(-aimag(ricker_out(it)),real(ricker_out(it)))
         enddo

         ricker_in(:)=cmplx(0.0,0.0)

!            do it=2,nt/2
!               ricker_out(it)=2.e0*ricker_out(it)
!            enddo
!            do it=nt/2+2,nt
!               ricker_out(it)=0.e0
!            enddo
!
!            ricker_in(:)=cmplx(0.0,0.0)

        call sfftw_execute(plan2,ricker_out,ricker_in)

        call sfftw_destroy_plan(plan1)
        call sfftw_destroy_plan(plan2)

         do it=1,nt
            ricker_hb(it)=real(ricker_in(it))/nt_fft
         enddo

!       if(.false.) then 
        open(12,file='ricker_wavelet_hb_fre',access='direct',&
                form='unformatted',recl=nt)
         write(12,rec=1) (ricker_hb(i),i=1,nt)
        close(12)

!        open(12,file='ricker',access='direct',form='unformatted',&
!                recl=1000)
!         write(12,rec=1) (ricker(i),i=1,1000)
!        close(12)
!       endif

        open(12,file='ricker_wavelet_hb_fre.ascii')
         do i=1,nt
            write(12,*) i,ricker(i)
         enddo 
        close(12)
!          ricker_hb(:)=aimag(ricker_in(:))/nt_fft


      end subroutine

