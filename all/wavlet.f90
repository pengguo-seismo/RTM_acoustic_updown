!     Modified by Feng Deng
!
!
!     Modification by Xinfa Zhu, Nov 15 2011
!     nprd is not easily understandable for users
!     change nprd to f0, dominant frequency of the output wavelet
!
!     nprd = 1.0/(1.4*f0*dt)
!
!     f0 in Hz, dt in second, nprd is number of time points per half period.    
!     the output is a 90-degree-phase Ricker wavelet, anti-symmetric,
!       starts from time zero (it = 1);
!       negative at the first half period (it = 1, nprd);
!       positive at the second half period (it = nprd+1, nprd*2);
!       zero for other times (it = nprd*2+1, nt)
!     Looks like negative sine function from 0 to 360 degrees.
!     The above relation is found and tested by Xinfa Zhu
!     
!
! **********************************************************************
!
!  Variable :
!
!   nstep  : # of time steps.
!
! **********************************************************************
!
!  Array :
!
!   seri    : source as function of time.
!
! **********************************************************************
!
        SUBROUTINE sourcefunc(amp,nt,nprd,dt)
        implicit none   !No undeclared variables are allowed
        INTEGER nt, nprd, it, j
        REAL*4  dt, amp(1:nt), pi, amax, fdom
!
!       !NOTE: changing nprd is equivalent to modifying the dominant 
!       !frequency of the source function
!        nprd = nint(1./(1.4*fdom*dt))
        pi   = 4.*atan(1.)
        amax = 0.       !Assign 0 to maximum value
!
!       !Choose source as a derivative (w.r.t. time) of a gaussian:
        do it=1,nt
           j = it-nprd+2
           amp(it) = -2.*j*exp(-(pi*dt*j/(nprd*dt))**2)
           if(abs(amp(it)) .gt. amax) amax = abs(amp(it))
        enddo
!
        do it=1,nt
           amp(it) = amp(it) / amax     !Normalize the source function
        enddo
!
!       !write out source funtion
!        open(7, file='source.bin', access='direct',
        open(7, form='unformatted', access='direct', recl=4*nt)
        write(7,rec=1) (amp(it), it = 1,nt)
        close(7)
!
        RETURN
        END

!! I will use this wavelet
      subroutine rickerfunc_new_seri(dt,fp,length,nt,ricker)

       include 'fftw3.f'

       integer length,nt,nnt,nt_fft
       integer,parameter::dp=kind(0.e0)
       real(dp)::dt,fp,tshift
       real::PI
       real(dp)::ricker(1:nt),ricker_hb(1:nt)
       integer*8::plan1,plan2
       integer::ii_bi

       ii_bi=4

       PI=4.0*atan(1.0)

       tshift=0.0

       ricker=0.

       nnt=nt

       amax=0.

       a=PI*PI*fp*fp
       t0 = 1.50 /fp

       do it=1,nt
          t = real(it-1)*dt
          ricker(it) = (1.e0 - 2.e0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
          if(abs(ricker(it)).gt.amax)amax=abs(ricker(it))
       enddo

       do it=1,nt
          ricker(it)=ricker(it)/amax
       enddo

       print *,'test1',nt     
 
        open(12,file='ricker_wavelet',status='replace',access='direct',&
                form='unformatted',recl=ii_bi*nt)
         write(12,rec=1) (ricker(i),i=1,nt)
        close(12)

       print *,'test2'  
   
        open(12,file='ricker_wavelet.ascii')
         do i=1,nt
            write(12,*) i,ricker(i)
         enddo
        close(12)

      end

      subroutine rickerfunc_new(dt,fp,length,nt,ricker)
       integer length,nt
       integer,parameter::dp=kind(0.e0)
       real(dp)::dt,fp
       real::PI
!       parameter(PI=3.1415926)
       real(dp)::ricker(1:nt)

       ricker=0.

       nnt=nt

       tshift=0.0

       PI=4.0*atan(1.0)

       do i=1,nnt
          t=(i-1)*dt
          ts=1.0/fp
          ag = PI*PI*fp*fp
          tau=PI*(t-1.5*ts-tshift)/(1.5*ts)
          amp=(((1.0-4.0*tau*tau)*exp(-2.0*tau*tau)))
          ricker(i)=amp
       enddo

        open(12,file='ricker_fwd',access='direct',form='unformatted',&
                recl=1000)
         write(12,rec=1) (ricker(i),i=1,1000)
        close(12)

        open(12,file='ricker_fwd.dat')
         do i=1,1000
            write(12,*) i,ricker(i)
         enddo
        close(12)

      end

      subroutine rickerfunc(dt,fp,length,nt,ricker)  
       integer length,nt
       integer,parameter::dp=kind(0.e0)
       real(dp)::dt,fp
       parameter(PI=3.1415926)
       real(dp)::ricker(1:nt)
       real(dp),allocatable::b(:)
       real::ricker2(1:nt)
       ricker=0.

       do i=-1000,0
               ddt=dt*i
               ee=PI*ddt*fp
               ax=(1.0-2.e0*ee*ee)*EXP(-ee*ee)
              if(abs(ax).lt.1.0e-5)then 
                ax=0.0
                else
                k1=i
                 exit
              endif
       enddo
       length=abs(k1)*2+1

       print *,length,'length_length_length'
       allocate(b(1:length))
       j=0
       do i=k1,0
               ddt=dt*i
               ee=PI*ddt*fp
               ax=(1.e0-2.e0*ee*ee)*EXP(-ee*ee)
           j=j+1
           b(j)=ax
       enddo
        do i=1,abs(k1)
           b(j+i)=b(j-i)
        end do

        ricker=0.0
        do i=1,length
           ricker(i)=b(i)
        enddo 

       open(222,file='ricker_0.dat')
       do i=1,1000
         write(222,*) i, ricker(i)
       enddo 
       close(222)
        do i=1,50000
            ricker2(i)=ricker(i)
        enddo 

        deallocate(b)

        open(12,file='ricker_fwd_0',access='direct',form='unformatted',&
                recl=1000)
         write(12,rec=1) (ricker(i),i=1,1000)
        close(12)      
 
      end subroutine

    
      subroutine wavelet_new(fmax,afmax,dt,nt,ricker)

!
!     starting time of the wavelet
      real::ricker(10000)
      real::der(10000)
     
      integer wavtype
      wavtype=1
      trise=0.04
      pi     = 4.0*atan(1.0)
      alpha  = pi*fmax/(sqrt(-alog(afmax)))
      alpha2 = alpha*alpha
      pisqrv = 1.0/(sqrt(pi))
      print *,'trise',trise,fmax,afmax

      do i=1,nt
        t=(i-1)*dt
        arg=alpha2*(t-trise*0.5)**2
        ricker(i)=alpha*pisqrv*exp(-arg)
        if(wavtype.eq.2)then
             der(i)=-2.*alpha2*(t-trise*0.5)*ricker(i)
             ricker(i)=der(i)
        endif
      enddo
!
       open(222,file='wavelet.dat')
       do i=1,nt
         write(222,*) ricker(i)
       enddo 
       close(222)
!      write(luwav)(ricker(i),i = 1,nt)

      if(istype.eq.3)then
         call gaussian(ricker,nt,dt)
      endif
!
      return
      end


!---------------------------------------------------------------------
!
!     Time Domain Integration
!
!---------------------------------------------------------------------
!
      subroutine gaussian(temp,npoint,dt)
      dimension temp(10000)
      totint=0.0
!
      do  i=1,npoint-1
      prtint=0.5*dt*(temp(i)+temp(i+1))
      totint=totint+prtint
      temp(i)=totint
      enddo 
      do  i=1,npoint-1
!     temp(i)=temp(i+1) - temp(1)
      temp(i)=temp(i) - temp(1)
      enddo 
      return
      end
!
!=======================================================================
