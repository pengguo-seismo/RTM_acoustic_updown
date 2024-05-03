     module precision_m_mc
       integer, parameter, public :: singlePrecision = kind(0.0) ! Single precision
       integer, parameter, public :: doublePrecision = kind(0.0d0) ! Double precision

! Comment out one of the lines below
!       integer, parameter, public :: fp_kind_mc = singlePrecision
        integer, parameter, public :: fp_kind_mc = doublePrecision
     end module precision_m_mc

    module interface_c_struct
        use, intrinsic::iso_c_binding
        implicit none

        type, bind(c)::Layr
          real(c_double)::wid
          real(c_double)::slw
        end type Layr           

        type, bind(c)::Sample
          integer(c_int)::layrC
!          type(Layr),allocatable::layr(:)
          type(Layr)::layr(100)
        end type Sample

        type(Sample),bind(c)::theta1,theta2
    end module interface_c_struct

    module interface_c_func
     use interface_c_struct
     implicit none
!     type(Sample)::theta1,theta2
       Interface
!setLayerFitParameters(lambda,mu,sigma,minS,maxS,maxlayers,depth,seed)
         function setLayerFitParameters(lambda,mu,&
           sigma,minS,maxS,maxlayers,depth,seed) BIND(c,name="setLayerFitParameters")
           use, intrinsic::iso_c_binding,only:c_int,c_double,c_long 
           implicit none
           integer(c_int)::setLayerFitParameters
           integer(c_int)::maxlayers
           integer(c_long)::seed
           real(c_double)::lambda,mu,sigma,depth,minS,maxS
!          integer(c_float)::lambda,mu,sigma,L
         end function setLayerFitParameters

         function createInitialSample() BIND(c,name="createInitialSample")
           use, intrinsic::iso_c_binding,only:c_int
           integer(c_int)::createInitialSample
         end function createInitialSample

         function getCurrentGrid(divz,vels) BIND(c,name="getCurrentGrid")
          use, intrinsic::iso_c_binding,only:c_int,c_double
          implicit none
           integer(c_int)::divz
           real(c_double)::vels(divz)
           integer(c_int)::getCurrentGrid
         end function getCurrentGrid

         function setCurrentLikelihood(like) BIND(c,name="setCurrentLikelihood")
          use, intrinsic::iso_c_binding,only:c_int,c_double
          implicit none
          integer(c_int)::setCurrentLikelihood
          real(c_double)::like
         end function setCurrentLikelihood

         function createSuggestion() BIND(c,name="createSuggestion")
           use, intrinsic::iso_c_binding,only:c_int
           integer(c_int)::createSuggestion
         end function createSuggestion

         function getSuggestedGrid(divz,vels) BIND(c,name="getSuggestedGrid")
          use, intrinsic::iso_c_binding,only:c_int,c_double
          implicit none
           integer(c_int)::divz
           real(c_double)::vels(divz)
           integer(c_int)::getSuggestedGrid
         end function getSuggestedGrid

         function returnLikelihood(like,T) BIND(c,name="returnLikelihood")
          use, intrinsic::iso_c_binding,only:c_int,c_double
          implicit none
          integer(c_int)::returnLikelihood
          real(c_double)::like,T
         end function returnLikelihood

         function getCurrentLikelihood(lk) BIND(c,name="getCurrentLikelihood")
          use, intrinsic::iso_c_binding,only:c_int,c_double
          implicit none
          integer(c_int)::getCurrentLikelihood
          real(c_double)::lk
         end function getCurrentLikelihood
          
       end interface    
     end module interface_c_func

     module mcmc_par    
      use precision_m_mc
      implicit none

      real(fp_kind_mc)::lambda,mu,sigma,maxS,minS,depth 
      integer::rns

      contains

       subroutine par_in_mcmc     
          open(12,file='parameter_mcmc.txt')
           read(12,*) rns
           read(12,*) lambda
           read(12,*) mu
           read(12,*) sigma
           read(12,*) maxS
           read(12,*) minS
           read(12,*) depth
          close(12)
       end subroutine par_in_mcmc

     end module mcmc_par               

     program test

     use, intrinsic::iso_c_binding,only:c_int
     use interface_c_func
!, only:setLayerFitParameters,createInitialSample,getCurrentGrid
!     use interface_c_struct 
     use precision_m
     use precision_m_mc

     use main_module
     use parameter_input
     use wave_propagator

     use mcmc_par
!     use mpi

     implicit none

     include 'mpif.h'

! first define global parameters
! names are the same as in the 1D-v1 design document

!     double precision::lambda,mu,sigma,L
!     real::lambda,mu,sigma,L 
     integer::maxlayers
     integer(kind=8)::seed
     integer(kind=4)::seed_tmp
     integer(c_int)::int_setparameters,int_init,int_getgrid,&
                     int_setlikelihood
     integer(c_int)::int_creatsuggestion,int_getSuggestedGrid,&
                     int_returnLikelihood,int_getCurrent
     integer(c_int)::divz
     real(fp_kind_mc),allocatable::vels(:)
     real(fp_kind_mc)::like,T
     integer::i_rns_int
     real(fp_kind_mc)::like_current

     double precision::t1,t2,t12
     double precision::t1_a,t2_a,t12_a     

     double precision::norm_term,norm_term_sum,norm_res_data
     real,allocatable::syn_p(:,:),res_p(:,:),obs_p(:,:)
     integer::num_points

     real(fp_kind),allocatable::v_post(:,:),v_post_avg(:)

     character*90::vp_rns
!     type(Sample)::theta1,theta2

!     common /theta12/ theta1,theta2     

     integer::i_rns
      integer::clock 

     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

     call getarg(1,arg)
     arg=trim(arg)
     read(arg,'(i4)') iflag_par

     ii_bi=1

     call par_in(iflag_par)

     allocate(syn_p(1:nx_rec,1:nt),res_p(1:nx_rec,1:nt),&
              obs_p(1:nx_rec,1:nt))
     syn_p=0.0
     res_p=0.0
     obs_p=0.0

! a priori expected nt of layers
!     lambda=8.0
! a priori mean of slowness
!     mu=0.0005
!     mu=0.4
! a priori stdev of slowness
!     sigma=mu*0.3
! depth of region of interest
!     depth=2880.

     if(myid.eq.0) then 
     call system_clock(count=clock)
     seed_tmp=clock
     endif

     call mpi_bcast(seed_tmp,1,mpi_integer,0,mpi_comm_world,ierr)

     seed=seed_tmp+myid
     print *,'seed',seed,'myid',myid   

     call par_in_mcmc

     maxlayers=10
!     seed = 49132411

!(&lambda,&mu,&sigma,&minS,&maxS,&maxLayers,&L,&seed)
     int_setparameters=setLayerFitParameters(lambda,mu,sigma,minS,maxS,maxlayers,depth,seed)

     print *,'int_setparameters',int_setparameters
     print *,'seed',seed

     divz=nz_p

     allocate(vels(1:divz))

     int_init=createInitialSample() 

     int_getgrid=getCurrentGrid(divz,vels)

     print *,'divz',divz
 
!     int_test=resetSamples()
      call allocate_main_module(nx,nz,nrel,nx_rec,lv)

!! initialize velocity & density
      do j=1,divz
         do i=nx_p1,nx_p2
            vp_f0(i,j+npoints_pml)=vels(j)
         enddo 
      enddo 

      open(12,file='vp_test',access='direct',form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
         write(12,rec=i) &
                     (vp_f0(i,j),j=1,nz)
       enddo
      close(12)

      like=0
      T=1

      int_setlikelihood=setCurrentLikelihood(like)

      do i_rns=1,10000

      int_creatsuggestion=createSuggestion()

      print *,'test3',' myid ',myid,' i_rns ',i_rns

      int_returnLikelihood=returnLikelihood(like,T)

      print *,'return',' myid ',myid,'int ',int_returnLikelihood

      enddo 
!! initialize velocity & density
! get the velocity grid values of the suggested sample
      int_getSuggestedGrid=getCurrentGrid(divz,vels)
      do j=1,divz
         do i=nx_p1,nx_p2
            vp_f0(i,j+npoints_pml)=vels(j)
         enddo 
      enddo 
      print *,'test4',' myid ',myid,' i_rns ',i_rns

      vp_rns='vp'//char(myid/100+48)//&
             char(mod(myid,100)/10+48)//char(mod(myid,10)+48)

      open(12,file=vp_rns,access='direct',form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
         write(12,rec=i) &
                     (vp_f0(i,j),j=1,nz)
       enddo
      close(12)
      print *,'test5',' myid ',myid,' i_rns ',i_rns

      call MPI_FINALIZE(ierr)

      end
