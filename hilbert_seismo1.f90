
!! correct

        subroutine hilbert_seismo1(input,output,nx,nt)
         implicit none

         include 'fftw3.f'

         integer::nx,nt
         real::input(1:nx,1:nt)
         real::output(1:nx,1:nt)
         integer*8::plan1,plan2
         integer::i,j,it,nt_fft,n_of_fft
         complex,allocatable::input_fft(:),output_fft(:)

         nt_fft=n_of_fft(nt)

         allocate(input_fft(1:nt_fft),output_fft(1:nt_fft))

         print *,'nt',nt,'nt_fft',nt_fft

         call sfftw_plan_dft_1d &
              (plan1,nt_fft,input_fft,output_fft,FFTW_FORWARD,FFTW_ESTIMATE)

         call sfftw_plan_dft_1d &
              (plan2,nt_fft,output_fft,input_fft,FFTW_BACKWARD,FFTW_ESTIMATE)

! Hilbert transform
         do i=1,nx
            input_fft(:)=cmplx(0.0,0.0)
            output_fft(:)=cmplx(0.0,0.0)

            do j=1,nt
               input_fft(j)=cmplx(input(i,j),0.0)
            enddo

            call sfftw_execute(plan1,input_fft,output_fft)

            do it=2,nt_fft/2
                output_fft(it)=2.0*cmplx(real(output_fft(it)),aimag(output_fft(it)))
            enddo
            do it=nt_fft/2+2,nt_fft
                output_fft(it)=cmplx(0.0,0.0)
            enddo

            input_fft(:)=cmplx(0.0,0.0)

            call sfftw_execute(plan2,output_fft,input_fft)

            output(i,1:nt)=aimag(input_fft(1:nt)/nt_fft)
!            output(i,1:nt)=cmplx(real(input_fft(1:nt)/nt_fft),aimag(input_fft(1:nt)/nt_fft))
         enddo

         call sfftw_destroy_plan(plan1)
         call sfftw_destroy_plan(plan2)

         deallocate(input_fft,output_fft)

        end subroutine hilbert_seismo1
