
!!  head files, declaration of variables 

! Coded by : Peng Guo
! Date : May 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------

  integer,parameter::dp=kind(0.e0)

  real (dp), parameter :: factor = 1.0e7

  real (dp)::to

  integer::ii_bi

  real,allocatable::geom_x(:),geom_z(:)
  integer,allocatable::rcv_x(:),rcv_z(:)

! source

  integer::isource,jsource
  integer::isource_p,jsource_p

! value of PI
  real (dp), parameter :: PI = 3.141592653589793238462643_dp

! conversion from degrees to radians
  real (dp), parameter :: DEGREES_TO_RADIANS = PI / 180.0_dp

! zero
  real (dp), parameter :: ZERO = 0.0_dp

! large value for maximum
  real (dp), parameter :: HUGEVAL = 1.0e+30

  real (dp), parameter :: STABILITY_THRESHOLD = 1.0e+25

! main arrays

  real (dp),allocatable::vx(:,:),vz(:,:),p(:,:)

  real (dp),allocatable::p_tmp(:,:)

  integer,allocatable::recv_count(:),disp(:)

! power to compute d0 profile
  real (dp), parameter :: NPOWER = 2.0_dp

  real (dp), parameter :: K_MAX_PML = 1.0_dp ! from Gedney page 8.11

   real (dp)::alpha_max_pml

! arrays for the memory variables
  real (dp), allocatable :: &
      memory_dvx_dx(:,:), &
      memory_dvz_dz(:,:), &
      memory_dp_dx(:,:), &
      memory_dp_dz(:,:)

  real (dp) :: &
      value_dvx_dx, &
      value_dvz_dz, &
      value_dp_dx, &
      value_dp_dz

!!!!!!!

  integer::int_shot_model,int_shot_move,mx
  logical::lossless

! 1D arrays for the damping profiles
! nx
  real (dp), allocatable :: K_x(:),a_x(:,:),b_x(:,:)
  real (dp), allocatable :: K_x_half(:),a_x_half(:,:),b_x_half(:,:)

  real (dp), allocatable :: K_z(:),a_z(:,:),b_z(:,:)
  real (dp), allocatable :: K_z_half(:),a_z_half(:,:),b_z_half(:,:)

! for the source
  real (dp) :: sum_a,cons1,cons2,deltat_v

  integer :: i,j,it,nl,l

  real (dp) :: Courant_number,velocnorm

  integer::nx,nz
  integer::nx_p,nz_p
  integer::nx_rec
  integer::rcv_f,rcv_l,i_zr
  integer::rcv_f_p,rcv_l_p,i_zr_p
  integer::nx_p1,nx_p2,nz_p1,nz_p2

  integer::nx_st,nx_ed,nz_st,nz_ed
  integer::nx_st_p,nx_ed_p,nz_st_p,nz_ed_p
  integer::num_ele_d

  integer::nz_s,nz_mod

! size of a grid cell
  real (dp)::deltax,deltaz

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_ZMIN = .true.
  logical, parameter :: USE_PML_ZMAX = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 20

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gird based information
  real (dp), allocatable::c11_u(:,:),c11_r(:,:)

  real (dp), allocatable::r(:,:,:)
  real (dp), allocatable::tau11(:,:,:),tau11_new(:,:,:)
  real (dp), allocatable::tau_sigma(:,:,:)

  real (dp), allocatable::x1(:,:,:),x2(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real (dp), allocatable::tau11_eps(:,:,:)

  real (dp), allocatable::tau11_sum(:,:)
  real (dp), allocatable::dens(:,:)
  real (dp), allocatable::q_zero(:),omega_zero(:)
  real (dp), allocatable::vp(:,:),vp_f0(:,:)

  real (dp), allocatable::q11(:)

  real (dp), allocatable::c(:),dlh(:),dlv(:)

  real (dp), allocatable::depth(:)

  real (dp)::ricker(1:10000)

  real (dp), allocatable::rho(:,:)
  real (dp)::f0,f_ref

  real,allocatable::tmp_x(:,:),tmp_z(:,:)

! total number of time steps
  integer::nt,nrel,length,iflag,int_num

  integer,parameter::lv=4

  integer::nt_snap

  real(dp),allocatable::tau_sigma_tmp(:)

! time step in seconds

  real (dp)::deltat

  real (dp)::sum_r,amp

!  arrays for jacobian inversion

   real (dp),allocatable::jacobian(:,:)

   integer::iflag_par
   character *90 vpname, rhoname, tau_epsname, tau_sigmaname
   character *90 file_p,file_p_part,recons_part,snap_part
   character *90 file_snap_all_part,file_snap_all
   character *90 file_image_part,file_image_normal_part
   character *90 file_image_illu_part
   character *90::file_image,file_image_normal,file_image_illu,arg
   logical::modeling
   integer::ierr,myid,num_procs
   integer::i_shot,num_shot
   integer::num_procs_shot,num_shot_once_iter

   integer::myid_s,mpi_comm_world_s
   integer::my_source

   integer::iv,jv
   integer::recv_left,recv_right,send_left,send_right
   integer::s_size

   real,allocatable::sou_p(:,:),image(:,:),image_illu(:,:)

   integer::file_p_no,file_snap_no,file_snap_all_no,file_image_no

   real,allocatable::epsilon_m(:,:)

   real,allocatable::vx_l(:,:,:),vx_r(:,:,:),vx_t(:,:,:),vx_b(:,:,:)
   real,allocatable::vz_l(:,:,:),vz_r(:,:,:),vz_t(:,:,:),vz_b(:,:,:)
   real,allocatable::p_l(:,:,:),p_r(:,:,:),p_t(:,:,:),p_b(:,:,:)
   real,allocatable::r_l(:,:,:,:),r_r(:,:,:,:),r_t(:,:,:,:),r_b(:,:,:,:)

