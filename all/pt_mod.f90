!! for parallel tempering
      module ptglobal

      integer                            :: iseed,nprocc,iprocc
      Double Precision, allocatable      :: Tbins(:)
      Real, allocatable                  :: sjd(:,:)
      Integer, allocatable               :: temptrans(:,:,:)
      Integer, allocatable               :: mcmctrans(:,:)
      Integer                            :: ntemps
      Integer                            :: iout
      Logical                            :: record_temp
      Logical                            :: record_mcmc
      Logical                            :: record_temp_now
      Logical                            :: record_mcmc_now
      Logical                            :: RestrictTemp
      Logical                            :: verbose
      Logical                            :: silent
      Real                               :: TimeStart,TimeEnd,time1,time2

      end module ptglobal

      Subroutine Setuptempladder &
                 (nchains,tlow,thigh,chaintemp,seed,n_count_0)

!      use PTmod

      Double precision                   :: chaintemp(nchains)
      Double precision                   :: t,dt,tlow,thigh
      Double precision                   :: aval,bval
      Character(len=90)                     filename
      integer(kind=8)                    :: seed
      integer                            :: iseed
      Double precision                   ::perc_0
      integer                            ::n_count_0

      ! Selected temperatures randomly using log-uniform distribution  
      iseed=seed+10

      aval = log(tlow)
      bval = log(thigh)
      dx = (bval-aval)/(nchains-1)

      do it=1,nchains
         chaintemp(it) = exp(aval + ran3(iseed)*(bval-aval))
      end do

      !if(rank==1)Chaintemp(1) = tlow       ! Force some chains to be at tlow
      chaintemp(1:n_count_0) = tlow                   ! Force first chain to be at tlow

      filename = 'tlevels'
      open(15,file=filename,status='unknown')

            do it = 1,nchains
              write(15,*)it,chaintemp(it)
            end do
            close(15)
      return
      end

!
! ----------------------------------------------------------------------------
!
!           Numerical Recipes random number generator 
!
! ----------------------------------------------------------------------------
            FUNCTION ran3(idum)
            INTEGER idum
            INTEGER MBIG,MSEED,MZ
!           REAL MBIG,MSEED,MZ
            REAL ran3,FAC
            PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!           PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
            INTEGER i,iff,ii,inext,inextp,k
            INTEGER mj,mk,ma(55)
!           REAL mj,mk,ma(55)
            SAVE iff,inext,inextp,ma
            DATA iff /0/
!           write(*,*)' idum ',idum
            if(idum.lt.0.or.iff.eq.0)then
               iff=1
               mj=MSEED-iabs(idum)
               mj=mod(mj,MBIG)
               ma(55)=mj
               mk=1
               do 11 i=1,54
                  ii=mod(21*i,55)
                  ma(ii)=mk
                  mk=mj-mk
                  if(mk.lt.MZ)mk=mk+MBIG
                  mj=ma(ii)
11             continue
               do 13 k=1,4
                  do 12 i=1,55
                     ma(i)=ma(i)-ma(1+mod(i+30,55))
                     if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12                continue
13             continue
               inext=0
               inextp=31
               idum=1
            endif
            inext=inext+1
            if(inext.eq.56)inext=1
            inextp=inextp+1
            if(inextp.eq.56)inextp=1
            mj=ma(inext)-ma(inextp)
            if(mj.lt.MZ)mj=mj+MBIG
            ma(inext)=mj
            ran3=mj*FAC
            return
            END


!
! -----------------------------------------------------------
!
!       tswap_accept -> decides whether to accept a proposed temperature swap 
!                       using the Parallel Tempering extension of 
!                       Metropolis-Hastings acceptance criterion.
!
! 	tf(4) 	= T1 -> temperature 1
!		= E1 -> Energy 1 
!		= T2 -> temperature 2
!		= E2 -> Energy 2
!
!	yn	= .true. or .false.
!
!       Notes:
!             In an optimization problem Energy, e.g. E1,or E2 should
!             be the property of the model to be minimized.
!
!             In a posterior PDF sampling problem E should be
!             the negative log of the posterior PDF plus the 
!             log of the proposal distribution from state 1 to state 2.
!
!             i.e. E2 = -log[p(m2|d)] + log[q(m2|m1)]
!
!             The proposal distribution must correspond to that specified
!             in used to perturb model 1 to model 2.
!
!             E need only be evaluated up to an additive constant
!	      because the decision to accept or reject this swap
!             depends only on differences in E.
!           
!             If the prior is a constant and the proposal 
!             distribution is symmetric then E may be set to 
!             to the negative log-Likelihood, or data misfit. 
!
!             This routine optionally records the history of 
!             successful temperature swap for diagnostic purposes.
!
! ---------------------------------------------------------------
       Subroutine tswap_accept(T1,T2,E1,E2,yn)

       use ptglobal

       Double precision              :: E1,E2,delE
       Double precision              :: T1,T2,delT,delS
       Logical                       :: yn
 
       it = 0
       jt = 0
       yn = .false.
       delE = E1-E2
       delT = 1.0/T1 - 1.0/T2
       delS = delE*delT
       a = ran3(iseed)

       if(log(a).le.delS)yn = .true.      ! swap successful

       print *,iseed,'iseed at tswap',real(a),'ynynyn',yn


                                          ! find temperature bins

!       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T1,k1)
!       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T2,k2)

!       if(RestrictTemp.and.abs(k1-k2).ne.1)yn = .false. ! restrict temperature swaps to neighbours
       end

! ---------------------------------------------------------------
!
!      PT_findTbin - A utility routine to locate the bin containing a 
!                    given temperature, T.
!
! ---------------------------------------------------------------
!
! -----------------------------------------------------------
!
!       PT_McMC_accept -> decides whether to accept a proposed 
!                      Metropolis-Hastings step using the M-H 
!                      acceptance criterion. Optionally records
!                      numbers of proposed and accepted transitions
!                      as a function of chain temperature.
!
!	yn	= .true. or .false.
!
!       Notes:
!             Input variable logPPD should be
!             the negative log of the posterior PDF i.e. E2 = -log[p(m2|d)]
!
!             logQ21 is the log of the proposal distribution 
!             used to perturb from state 1 to state 2, i.e. log[q(m2|m1)]
!
!             logPPD and logQ12 etc need only be evaluated up to 
!             an additive constant because the decision to accept 
!             or reject this swap depends only on differences in logPPD.
!           
!             If the prior is a constant and the proposal 
!             distribution is symmetric then logPPD may be set to 
!             to the negative log-Likelihood, or data misfit. 
!
!             Input variable T is the temperature of the chain.
!
!             This routine optionally records the history of 
!             successful moves at each temperature for diagnostic purposes.
!
!             Since this routine is only for `within chain' transitions,
!             the temperature variable is passive and used only to identify 
!             the appropriate bin to accumulate transition data.
!      
!
! ---------------------------------------------------------------
       Subroutine PT_McMC_accept(T,logPPD1,logQ12,logPPD2,logQ21,yn)

       use ptglobal

       Double precision              :: logPPD1,logPPD2
       Double precision              :: logQ21,logQ12
       Double precision              :: delS
       Double precision              :: T
       Logical                       :: yn
 
       it = 0
       jt = 0
       yn = .false.
       delS = (logPPD1-logPPD2)/T
       delS = delS + logQ12 - logQ21
       a = ran3(iseed)

       if(log(a).le.delS)yn = .true.      ! swap successful

       !write(iout,*)T,logPPD1,logPPD2,delS,a,yn

       end












