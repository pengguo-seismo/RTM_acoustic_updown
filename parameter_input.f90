
     module parameter_input
      use precision_m
      use main_module,only:npoints_pml
      implicit none

      integer::num_shot,nrel,nx,nz,nt,nx_rec,nt_snap
      integer::ud_flag
      integer::nx_p,nz_p
      integer::nx_p1,nx_p2,nz_p1,nz_p2
      real(fp_kind)::deltax,deltaz,deltat
      real(fp_kind)::f0,f_ref
      character*256::vpname,rhoname,tau_epsname,tau_sigmaname,&
                    file_p_part,file_p_rcv_part,snap_part,recons_part
      character*256::file_image_part,file_image_normal_part,&
                    file_image_illu_part
      character*256::data_direc
      logical::lossless,multiple
      integer,parameter::lv=4

     contains

      subroutine par_in(iflag_par)
       integer::iflag_par
       character*90::par_name

       if(iflag_par.eq.0) then
        par_name='parameter.txt'
       else if (iflag_par.eq.1) then
        par_name='parameter_dir.txt'
       else if (iflag_par.eq.2) then
        par_name='parameter_rtm.txt'
       endif

      open(12,file=par_name)
       read(12,*) num_shot
       print *,'num_shot',num_shot
!!    number of relaxation mechanisms
       read(12,*) nrel
       print *,'nrel',nrel
!!    grid number -x
       read(12,*) nx_p
!!    grid number -z
       read(12,*) nz_p
!!    time step 
       read(12,*) nt
       read(12,*) nx_rec
       read(12,*) lossless
!!    snapshot output interval (time step)
       read(12,*) nt_snap
!!    grid increment - x 
       read(12,*) deltax
       print *,deltax
!!    grid increment - z
       read(12,*) deltaz
       print *,deltaz
!!    time increment
       read(12,*) deltat
       print *,deltat
!!    dominant frequency
       read(12,*) f0
       print *,f0
!!    reference frequency
!!    for velocity
       read(12,*) f_ref
!!     velocity file name
       read(12,*) vpname
       print *,vpname,'vpname'
!!     density file name
       read(12,*) rhoname
       print *,rhoname,'rhoname'
!!     strain relaxation time file name (Q)
       read(12,*) tau_epsname
       print *,tau_epsname,'tau_epsname'
!!     stress relaxation times file name (Q)
       read(12,*) tau_sigmaname
       print *,tau_sigmaname,'tau_epsname'
!!     seismogram file name
       read(12,*) file_p_part
       print *,file_p_part,'seis_p_part'
!!     cross-correlation image file name (only applicable for rtm)
       read(12,*) file_image_part
!!     normalized cross-correlation image file name (only applicable for rtm)
       read(12,*) file_image_normal_part
!!     illumination file name (only applicable for rtm)
       read(12,*) file_image_illu_part
!!     source wavefield snapshot file name
       read(12,*) snap_part
       print *,snap_part,'snap_part'
!!     receiver wavefield snapshot file name
       read(12,*) recons_part
       print *,recons_part,'recons_part'
       read(12,*) data_direc
       print *,data_direc,'data_direc'
     close(12)

      nx=nx_p+2*npoints_pml
      nz=nz_p+2*npoints_pml

      nx_p1=npoints_pml+1
      nx_p2=nx_p+npoints_pml
      nz_p1=npoints_pml+1
      nz_p2=nz_p+npoints_pml

     end subroutine par_in

      subroutine par_in_ud(iflag_par)
       integer::iflag_par
       character*90::par_name

      par_name='parameter_rtm_ud.txt'

      open(12,file=par_name)
       read(12,*) ud_flag
       read(12,*) num_shot
       print *,'num_shot',num_shot
!!    number of relaxation mechanisms
       read(12,*) nrel
       print *,'nrel',nrel
!!    grid number -x
       read(12,*) nx_p
!!    grid number -z
       read(12,*) nz_p
!!    time step 
       read(12,*) nt
       read(12,*) nx_rec
       read(12,*) lossless
!!    snapshot output interval (time step)
       read(12,*) nt_snap
!!    grid increment - x 
       read(12,*) deltax
       print *,deltax
!!    grid increment - z
       read(12,*) deltaz
       print *,deltaz
!!    time increment
       read(12,*) deltat
       print *,deltat
!!    dominant frequency
       read(12,*) f0
       print *,f0
!!    reference frequency
!!    for velocity
       read(12,*) f_ref
!!     velocity file name
       read(12,*) vpname
       print *,vpname,'vpname'
!!     density file name
       read(12,*) rhoname
       print *,rhoname,'rhoname'
!!     strain relaxation time file name (Q)
       read(12,*) tau_epsname
       print *,tau_epsname,'tau_epsname'
!!     stress relaxation times file name (Q)
       read(12,*) tau_sigmaname
       print *,tau_sigmaname,'tau_epsname'
!!     seismogram file name
       read(12,*) file_p_part
       print *,file_p_part,'seis_p_part'
!!     cross-correlation image file name (only applicable for rtm)
       read(12,*) file_image_part
!!     normalized cross-correlation image file name (only applicable for rtm)
       read(12,*) file_image_normal_part
!!     illumination file name (only applicable for rtm)
       read(12,*) file_image_illu_part
!!     source wavefield snapshot file name
       read(12,*) snap_part
       print *,snap_part,'snap_part'
!!     receiver wavefield snapshot file name
       read(12,*) recons_part
       print *,recons_part,'recons_part'
       read(12,*) data_direc
       print *,data_direc,'data_direc'
     close(12)

      nx=nx_p+2*npoints_pml
      nz=nz_p+2*npoints_pml

      nx_p1=npoints_pml+1
      nx_p2=nx_p+npoints_pml
      nz_p1=npoints_pml+1
      nz_p2=nz_p+npoints_pml

     end subroutine par_in_ud


     end module parameter_input



