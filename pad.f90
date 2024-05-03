! =================================================================
! 
!    Subroutines are from 
!                    Total E & P
!                     June, 2012
!
! ================================================================


!!********************************************************
!!* prim - Looks for p such that n=a*m**p
!!********************************************************

logical function prim(n,m,p)

   implicit none
   integer :: n, m, p
   integer :: i

   prim=.false.
   doi: do i=1,n
      if ((m**i)>n) exit doi
      if (mod(n,m**i)==0) then
         p=i
         prim=.true.
      end if 
   end do doi

end function prim


!!***********************************************************
!!* n_of_fft - integer function to return a valid FFT integer
!!*            number, which is the combination of the power
!!*            of prime numbers 2, 3 and 5, such as:
!!*                 2^k1 * 3^k2 * 5^3
!!***********************************************************


!!$#ifdef FFT_ESSL

integer function n_of_fft(nin)

implicit none

integer, parameter :: nlist=5
integer            :: maxpower(nlist)

integer, intent(in) :: nin
integer             :: primlist(nlist), power(nlist), n, ilist
logical             :: divi(nlist), prim

maxpower(1)=25
maxpower(2)=2
maxpower(3)=1
maxpower(4)=1
maxpower(5)=1

primlist(1)=2
primlist(2)=3
primlist(3)=5
primlist(4)=7
primlist(5)=11

n_of_fft=nin
if (n_of_fft==1) return
doo: do
   do ilist=1,nlist
      divi(ilist)=prim(n_of_fft,primlist(ilist),power(ilist))
   end do
   n=1
   do ilist=1,nlist
      if (divi(ilist)) n=n*primlist(ilist)**min(power(ilist),maxpower(ilist))
   end do
   if ( n/2*2 /= n) then
      n_of_fft=n_of_fft+1
   else if (n==n_of_fft) then
      exit doo
   else
      n_of_fft=n_of_fft+1
   end if
end do doo


end function n_of_fft


!!**********************************************************
!!* npad - get a valid FFT number for a given number with
!!*        integer padding number
!!**********************************************************

subroutine npad(n,pad,npaded)

implicit none

integer :: n, pad, npaded
integer :: n_of_fft, pad_old

if (n>1) then
   npaded  = n_of_fft(n+pad)
   pad_old = pad
   if(pad/=0) pad = npaded-n
else
   npaded=n
end if

end subroutine npad


!!******************************************************************
!! MAGIC NUMBER FOR FFTW:
!!
!!          8,  16,   24,   32,   36,  48,   64,   72,   96,  108, &
!!        128, 144,  192,  216,  256, 288,  384,  432,  512,  576, &
!!        768, 864, 1024, 1152, 1536,1728, 2048, 2304, 3072, 3456, &
!!       4096, 4608  /)
!!******************************************************************

subroutine time_fftlen(nt,ntfft)

implicit none

integer, parameter  :: N = 32
integer             :: nt,ntfft,i
integer             :: magic_int(N) = (/                         &
          8,  16,   24,   32,   36,  48,   64,   72,   96,  108, &
        128, 144,  192,  216,  256, 288,  384,  432,  512,  576, &
        768, 864, 1024, 1152, 1536,1728, 2048, 2304, 3072, 3456, &
       4096, 4608  /) 

if(nt<magic_int(1) .or. nt>magic_int(N)) then
   print*,'TIME_FFTLEN: input nt wrong'
   stop
endif

do i=1,N
   if(nt<=magic_int(i)) exit
enddo

ntfft = magic_int(i)

!!write(*,*)'TIME_FFTLEN: ntfft=',ntfft

return
end subroutine time_fftlen


!!***************************************************************
!! space_fftlen - get the space Fourier transform length.
!!                this subroutine is specially designed for
!!                shifted window implementation of LCB migration.
!!                because of shifted window implementation,
!!                we must be sure nxfft+lxwin is a proper
!!                number for FFT as set in array magic_int.
!!
!!   input:
!!      nx    - input space g rid length
!!      lxwin - window length
!!      nxfft - Fourier transform length.
!!              nxfft > nx is guarranted.
!!
!!
!! MAGIC NUMBER FOR FFTW:
!!
!!          8,  16,   24,   32,   36,  48,   64,   72,   96,  108, &
!!        128, 144,  192,  216,  256, 288,  384,  432,  512,  576, &
!!        768, 864, 1024, 1152, 1536,1728, 2048, 2304, 3072, 3456, &
!!       4096, 4608  /)
!!******************************************************************

subroutine space_fftlen(nx,lxwin,nxfft)

implicit none

integer, parameter  :: N = 32
integer             :: nx,nxfft,lxwin,i

integer             :: magic_int(N) = (/                         &
          8,  16,   24,   32,   36,  48,   64,   72,   96,  108, &
        128, 144,  192,  216,  256, 288,  384,  432,  512,  576, &
        768, 864, 1024, 1152, 1536,1728, 2048, 2304, 3072, 3456, &
       4096, 4608  /)

if(nx<magic_int(1) .or. nx>magic_int(N)) then
    print*,'TIME_FFTLEN: input nx wrong'
    stop
endif

do i=1,N
    if(nx+lxwin<=magic_int(i)) exit
enddo

nxfft = magic_int(i)-lxwin

!!write(*,*)'SPACE_FFTLEN: nxfft=',nxfft

return
end subroutine space_fftlen
