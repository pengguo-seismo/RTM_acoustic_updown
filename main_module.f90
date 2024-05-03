! declare and allocate host variables
! -----------------------------------

  module main_module

  use precision_m
!  use wave_propagator,only:lv
  implicit none

  double precision::t1_a,t2_a,t12_a
  double precision::t1_a_mp,t2_a_mp,t12_a_mp

  real (fp_kind), parameter :: factor = 1.0e7

  integer::ii_bi

!  integer::nxfft,nzfft
!  integer::nxpad,nzpad

!  real,allocatable::taper_x(:),taper_z(:)
!  real,allocatable::tapx_12(:),tapz_12(:)

  real,allocatable::geom_x(:),geom_z(:)
  integer,allocatable::rcv_x(:),rcv_z(:)

! source

  integer::isource,jsource
  integer::isource_p,jsource_p

! value of PI
  real (fp_kind), parameter :: PI = 3.141592653589793238462643

! conversion from degrees to radians
  real (fp_kind), parameter :: DEGREES_TO_RADIANS = PI / 180.0

! zero
  real (fp_kind), parameter :: ZERO = 0.0

  real (fp_kind), parameter :: STABILITY_THRESHOLD = 1.0e+25

!! fd coef
  real (fp_kind), allocatable::c(:),dlh(:),dlv(:)

! 1D arrays for the damping profiles
! nx
  real (fp_kind), allocatable :: k_x(:),a_x(:,:),b_x(:,:)
  real (fp_kind), allocatable :: k_x_half(:),a_x_half(:,:),b_x_half(:,:)

  real (fp_kind), allocatable :: k_z(:),a_z(:,:),b_z(:,:)
  real (fp_kind), allocatable :: k_z_half(:),a_z_half(:,:),b_z_half(:,:)

! for the source
  real (fp_kind) :: sum_a,cons1,cons2,deltat_v

  integer :: i,j,it,nl,l

  real (fp_kind) :: tmp1,tmp2

  real (fp_kind) :: Courant_number,velocnorm

  integer::int_shot_model
  logical::modeling
! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_ZMIN = .true.
  logical, parameter :: USE_PML_ZMAX = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 20

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gird based information
  real (fp_kind), allocatable::c11_u(:,:),c11_r(:,:)

  real (fp_kind), allocatable::tau11(:,:,:),tau11_new(:,:,:)
  real (fp_kind), allocatable::tau_sigma(:,:,:)

  real (fp_kind), allocatable::x1(:,:,:),x2(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real (fp_kind), allocatable::tau11_eps(:,:,:)

  real (fp_kind), allocatable::dens(:,:)
  real (fp_kind), allocatable::vp(:,:),vp_f0(:,:)

  real (fp_kind), allocatable::ricker(:),ricker_hb(:)

  real (fp_kind), allocatable::rho(:,:)

  integer::length,iflag
  character *256::arg

  real (fp_kind)::amp,amp_hb

   integer::iflag_par
   integer::ierr,myid,num_procs
   integer::i_shot

   integer::file_p_no,file_snap_no,file_snap_all_no,file_image_no
   integer::file_p_sou_no,file_p_rcv_no
 
   character *256 file_p,file_image,file_image_normal,file_image_illu
   character *256 file_p_rcv

   contains

     subroutine allocate_main_module(nx,nz,nrel,nx_rec,lv)
      integer::nx,nz,nrel,nx_rec,lv

      allocate(c(1:lv),dlh(1:lv),dlv(1:lv))
      c=0.0
      dlh=0.0
      dlv=0.0

      allocate(geom_x(1:nx_rec+1),geom_z(1:nx_rec+1))
      geom_x=0
      geom_z=0

      allocate(rcv_x(1:nx_rec),rcv_z(1:nx_rec))
      rcv_x=0
      rcv_z=0

      allocate(k_x(1:nx),a_x(1:nx,1:nz),b_x(1:nx,1:nz))
      allocate(k_x_half(1:nx),a_x_half(1:nx,1:nz),b_x_half(1:nx,1:nz))

      k_x=zero
      a_x=zero
      b_x=zero
      k_x_half=zero
      a_x_half=zero
      b_x_half=zero

      allocate(k_z(1:nz),a_z(1:nx,1:nz),b_z(1:nx,1:nz))
      allocate(k_z_half(1:nz),a_z_half(1:nx,1:nz),b_z_half(1:nx,1:nz))

      k_z=zero
      a_z=zero
      b_z=zero
      k_z_half=zero
      a_z_half=zero
      b_z_half=zero

      allocate(c11_u(1:nx,1:nz))
      allocate(tau11(1:nrel,1:nx,1:nz),tau11_new(1:nrel,1:nx,1:nz))
      allocate(tau11_eps(1:nrel,1:nx,1:nz))
      allocate(tau_sigma(1:nrel,1:nx,1:nz))

      allocate(x1(1:nrel,1:nx,1:nz))
      allocate(x2(1:nrel,1:nx,1:nz))
      allocate(dens(1:nx,1:nz))
      allocate(vp(1:nx,1:nz))

      c11_u=zero
      tau11=zero
      tau11_new=zero
      tau11_eps=zero
      tau_sigma=zero

      x1=zero
      x2=zero
      dens=zero
      vp=zero

      allocate(vp_f0(1:nx,1:nz),c11_r(1:nx,1:nz))

      vp_f0(:,:)=zero
      c11_r(:,:)=zero

    end subroutine allocate_main_module

    subroutine deallocate_main_module()

      deallocate(c,dlh,dlv)
      deallocate(geom_x,geom_z)

      deallocate(rcv_x,rcv_z)

      deallocate(k_x,a_x,b_x)
      deallocate(k_x_half,a_x_half,b_x_half)

      deallocate(k_z,a_z,b_z)
      deallocate(k_z_half,a_z_half,b_z_half)

      deallocate(c11_u,tau11,tau11_new)
      deallocate(tau11_eps,tau_sigma)

      deallocate(x1,x2,dens,vp)

      deallocate(vp_f0,c11_r)

    end subroutine deallocate_main_module


   end module main_module






