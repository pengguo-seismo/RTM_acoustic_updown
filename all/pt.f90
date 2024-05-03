
      Subroutine pt &
                 (mode,ialg,nchains,nsteps,iburn,&
                 &modtemp,thigh,tlow,nbins,swaprate,&
                 &iseed0,basedir,nproc,iproc)

      use ptglobal

      Double precision, allocatable                 :: logPPD(:)
      Integer, allocatable                          :: finished(:)
      Double precision, dimension(4)                :: dmsg
      Double precision, dimension(3)                :: rmsg
      Double precision                              :: E1,E2,E
      Double precision                              :: t1,t2,T
      Double precision                              :: modtemp(*)
      Double precision                              :: thigh,tlow
      Integer, dimension(2)                         :: pair
      Integer                                       :: ipA,ipB,ipX,from
      Integer                                       :: tag,code,term
      Integer                                       :: totswapreq
      Integer                                       :: totswapacc
      Integer                                       :: totswapreqw
      Integer                                       :: totswapaccw
      Integer                                       :: nbins
      Logical, allocatable                          :: swapped(:)
      Logical                                       :: yes
      Character(len=100)                            :: filename
      Character(len=80)                             :: basedir
      Character*4                                   :: cproc

      if(mode.eq.0) then 
         iprocc = iproc
         nprocc = nproc
`


