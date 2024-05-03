
     real,allocatable::tap(:),tap2(:)
     
     n=9
     nn=2
!     call npad(n,nn,npaded)
     npaded=n_of_fft(n)

     print *,npaded

     allocate(tap(1:npaded))
     allocate(tap2(1:2*npaded))

     tap=0.0
     tap2=0.0

     call staper(tap,3,npaded)

     do i=1,npaded
        tap2(i)=tap(i)
     enddo
     do i=npaded+1,2*npaded
        tap2(i)=tap(i-npaded)
     enddo

     nn=14
     ntfft=n_of_fft(nn)
     print *,'ns',ntfft

     call time_fftlen(nn,ntfft)
     print *,'ns2',ntfft

     do i=1,2*npaded
        print *,i,tap2(i)
     enddo 

     end
