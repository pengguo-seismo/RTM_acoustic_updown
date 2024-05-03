
!**********************************************************
!* Hanning taper subroutine
!* Parameters:
!*        float *taper        1D array, the signal to be tapered
!*        int   ntaper        the length of the taper
!*
! -------------------------------------------------
!*    Code from Total E & P
! -------------------------------------------
!**********************************************************
  subroutine staper(taper,ntaper,nx)

  implicit none

  integer, intent(in)            ::nx
  integer, intent(in)            ::ntaper
  real, dimension(nx),intent(out)::taper

  integer :: ix,ntp
  real    :: pi

  pi=4*atan(1.0)

   ntp=ntaper-1
   taper=1

   if (nx.lt.ntaper) then
      write(*,*)'Warn: ntaper too big! No taper in the wave field!'
      return
   endif
   if (ntaper.lt.3) then
      write(*,*)'Warn: ntaper<3! No taper in the wave field!'
      return
   endif
   do ix=1,nx
      if (ix.le.ntaper) then
         taper(ix)=0.5-0.5*cos(pi*(ix-1)/ntp)
      else if (ix.gt.nx-ntaper) then
         taper(ix)=0.5-0.5*cos(pi*(nx-ix)/ntp)
      endif
   enddo

  return
  end


! ====================================================================
  subroutine maxvalue(input,n,lgs,ik_max)
  implicit none
  integer    :: ik_max,ik,n
  complex    :: lgs,input(n)

  ik_max=1;
  lgs=cmplx(0.0,0.0)
  do ik=1,n
     if (abs(input(ik)) .gt. abs(lgs)) then
        lgs= input(ik)
        ik_max=ik
     endif
  enddo

  end subroutine
! ===============================================================

  subroutine maxvalue_2d(local_p_fft,nxfft,nzfft,value_max,i_max,j_max)
  implicit none
  integer    :: i_max,j_max
  integer    :: nxfft,nzfft,i,j
  complex    :: value_max,local_p_fft(1:nxfft,1:nzfft)

  value_max=cmplx(0.0,0.0)
  do j=1,nzfft
     do i=1,nxfft
        if (abs(local_p_fft(i,j)) .gt. abs(value_max)) then
          value_max= local_p_fft(i,j)
          i_max=i
          j_max=j
        endif
     enddo 
  enddo

  end subroutine




