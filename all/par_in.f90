
!! read in parameters, for forward modeling and RTM !!

! Coded by : Peng Guo
! Date : May 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------
    subroutine par_in(num_shot,nrel,nx,nz,nt,nx_rec,&
                      lossless,deltax,deltaz,deltat,f0,f_ref,&
                      vpname,rhoname,tau_epsname,tau_sigmaname,&
                      file_p_part,&
                      file_image_part,file_image_normal_part,&
                      file_image_illu_part,&
                      snap_part,recons_part,nt_snap,iflag_par)
    implicit none
    include 'variable.h'

    character *90::par_name

     if(iflag_par.eq.0) then
        par_name='parameter.txt'
      else if (iflag_par.eq.1) then 
        par_name='parameter_dir.txt'
      else if (iflag_par.eq.2) then 
        par_name='parameter_rtm.txt'
     endif
 
     open(12,file=par_name)
       read(12,*) num_shot
!!    number of relaxation mechanisms
       read(12,*) nrel
       print *,nrel
!!    grid number -x
       read(12,*) nx
!!    grid number -z
       read(12,*) nz
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
       print *,file_p_part,'seis_part1'
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
   close(12)

   end
