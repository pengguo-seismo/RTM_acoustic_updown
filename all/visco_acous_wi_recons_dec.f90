!
! 2D isotropic visco-acoustic RTM program
! multiple rtm 
! Coded by : Peng Guo
! Date : March 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------
! TOTAL E&P USA, 2016, 05/20/16

  program seismic_CPML_2D_iso_rtm

!            ^ z
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        sigma_xz   |
!       v_z  +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!         sigma_xx    v_x
!         sigma_zz
!

   implicit none
   include 'variable.h'
   include 'mpif.h'

!---
!--- program starts here
!---
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

     ii_bi=4

!!   read in the total number of shots
     open(12,file='num_shot.txt')
         read(12,*) num_shot
     close(12)

!      call getarg(1,arg)
!      arg=trim(arg)
!      read(arg,'(i4)') iflag_par

      file_p_sou_no=21
      file_p_rcv_no=25
      file_snap_all_no=22
      file_image_no=23
      file_snap_no=24

     do i_shot=myid+1,num_shot,num_procs 

!!!    read in parameters for rtm
        call par_in_rtm_multiple_dec_local(nrel,nx_p,nz_p,nt,nt1,mx,int_shot_model,int_shot_move,&
                    int_rcv_move,int_rcv,lossless,multiple,rcv_f_p,rcv_l_p,nx_rec,deltax,deltaz,deltat,&
                    f0,f_ref,isource_p,jsource_p,i_zr_p,&
                    propg_bin_scale,num_local_x,num_local_z,nthet,nitermax,max_decomp,&
                    vpname,rhoname,tau_epsname,tau_sigmaname, &
                    file_p_sou_part,file_p_rcv_part,file_snap_all_part, &
                    file_image_part,file_image_normal_part,&
                    file_image_illu_part,&
                    snap_part,recons_part,nt_snap,dt_snap,iflag_par)

      print *,'main','nrel',nrel
      print *,'main','nx',nx_p
      print *,'main','nz',nz_p
      print *,'main','nt',nt
      print *,'main','deltax',deltax
      print *,'main','deltaz',deltaz
      print *,'main','deltat',deltat
      print *,'main','isource',isource_p
      print *,'main','jsource',jsource_p
      print *,'main','propg_bin_scale',propg_bin_scale
      print *,'main','num_local_x',num_local_x
      print *,'main','num_local_z',num_local_z
      print *,'main','nthet',nthet
      print *,'main','max_decomp',max_decomp
      print *,'main','nt_snap',nt_snap
      print *,'main','dt_snap',dt_snap

      modeling=.true.

      snap_part1='op_u'
      snap_part2='op_v'
      snap_part3='op_u_1st'
      snap_part4='op_v_1st'

      nx=nx_p+2*npoints_pml
      nz=nz_p+2*npoints_pml

      rcv_f=rcv_f_p+npoints_pml
      rcv_l=rcv_l_p+npoints_pml
      i_zr=i_zr_p+npoints_pml

      nx_p1=npoints_pml+1
      nx_p2=nx_p+npoints_pml
      nz_p1=npoints_pml+1
      nz_p2=nz_p+npoints_pml

      isource=isource_p+npoints_pml
      jsource=jsource_p+npoints_pml

      rcv_f=rcv_f+(i_shot-1)*int_rcv_move
!!!    change shot position with a increasing pattern
      isource=isource+(i_shot-1)*int_shot_move

!!!    print files names
       print *,vpname,rhoname,tau_epsname,&
              tau_sigmaname,file_p_sou_part,snap_part,nt_snap

       alpha_max_pml=2.0*PI*(f0/2.0)

           allocate(sou_p(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_vx(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_vz(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_r(1:nrel,-lv:nx+lv,-lv:nz+lv))
           allocate(sou_memory_dvx_dx(1:nx,1:nz))
           allocate(sou_memory_dvz_dz(1:nx,1:nz))
           allocate(sou_memory_dp_dx(1:nx,1:nz))
           allocate(sou_memory_dp_dz(1:nx,1:nz))

           sou_p(:,:)=0.0
           sou_vx(:,:)=0.e0
           sou_vz(:,:)=0.e0
           sou_r(:,:,:)=0.e0
           sou_memory_dvx_dx(:,:)=0.e0
           sou_memory_dvz_dz(:,:)=0.e0
           sou_memory_dp_dx(:,:)=0.e0
           sou_memory_dp_dz(:,:)=0.e0

           allocate(sou_p_hb(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_vx_hb(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_vz_hb(-lv:nx+lv,-lv:nz+lv))
           allocate(sou_r_hb(1:nrel,-lv:nx+lv,-lv:nz+lv))
           allocate(sou_memory_dvx_dx_hb(1:nx,1:nz))
           allocate(sou_memory_dvz_dz_hb(1:nx,1:nz))
           allocate(sou_memory_dp_dx_hb(1:nx,1:nz))
           allocate(sou_memory_dp_dz_hb(1:nx,1:nz))

           sou_p_hb(:,:)=0.0
           sou_vx_hb(:,:)=0.e0
           sou_vz_hb(:,:)=0.e0
           sou_r_hb(:,:,:)=0.e0
           sou_memory_dvx_dx_hb(:,:)=0.e0
           sou_memory_dvz_dz_hb(:,:)=0.e0
           sou_memory_dp_dx_hb(:,:)=0.e0
           sou_memory_dp_dz_hb(:,:)=0.e0

           allocate(vx_l(1:nz,1:lv,1:nt))
           allocate(vx_r(1:nz,1:lv,1:nt))
           allocate(vx_t(1:nx,1:lv,1:nt))
           allocate(vx_b(1:nx,1:lv,1:nt))

           allocate(vz_l(1:nz,1:lv,1:nt))
           allocate(vz_r(1:nz,1:lv,1:nt))
           allocate(vz_t(1:nx,1:lv,1:nt))
           allocate(vz_b(1:nx,1:lv,1:nt))

           allocate(p_l(1:nz,1:lv,1:nt))
           allocate(p_r(1:nz,1:lv,1:nt))
           allocate(p_t(1:nx,1:lv,1:nt))
           allocate(p_b(1:nx,1:lv,1:nt))

           vx_l=0.e0
           vx_r=0.e0
           vx_t=0.e0
           vx_b=0.e0

           vz_l=0.e0
           vz_r=0.e0
           vz_t=0.e0
           vz_b=0.e0

           p_l=0.e0
           p_r=0.e0
           p_t=0.e0
           p_b=0.e0

           allocate(vx_l_hb(1:nz,1:lv,1:nt))
           allocate(vx_r_hb(1:nz,1:lv,1:nt))
           allocate(vx_t_hb(1:nx,1:lv,1:nt))
           allocate(vx_b_hb(1:nx,1:lv,1:nt))

           allocate(vz_l_hb(1:nz,1:lv,1:nt))
           allocate(vz_r_hb(1:nz,1:lv,1:nt))
           allocate(vz_t_hb(1:nx,1:lv,1:nt))
           allocate(vz_b_hb(1:nx,1:lv,1:nt))

           allocate(p_l_hb(1:nz,1:lv,1:nt))
           allocate(p_r_hb(1:nz,1:lv,1:nt))
           allocate(p_t_hb(1:nx,1:lv,1:nt))
           allocate(p_b_hb(1:nx,1:lv,1:nt))

           vx_l_hb=0.e0
           vx_r_hb=0.e0
           vx_t_hb=0.e0
           vx_b_hb=0.e0

           vz_l_hb=0.e0
           vz_r_hb=0.e0
           vz_t_hb=0.e0
           vz_b_hb=0.e0

           p_l_hb=0.e0
           p_r_hb=0.e0
           p_t_hb=0.e0
           p_b_hb=0.e0

           allocate(last_p(-lv:nx+lv,-lv:nz+lv))
           allocate(last_vx(-lv:nx+lv,-lv:nz+lv))
           allocate(last_vz(-lv:nx+lv,-lv:nz+lv))

           last_p(:,:)=0.0
           last_vx(:,:)=0.0
           last_vz(:,:)=0.0

           allocate(last_p_hb(-lv:nx+lv,-lv:nz+lv))
           allocate(last_vx_hb(-lv:nx+lv,-lv:nz+lv))
           allocate(last_vz_hb(-lv:nx+lv,-lv:nz+lv))

           last_p_hb(:,:)=0.0
           last_vx_hb(:,:)=0.0
           last_vz_hb(:,:)=0.0

!! snapshot of the previous step
       allocate(old_p(-lv:nx+lv,-lv:nz+lv))
       old_p(:,:)=0.0

!! optical flow
       allocate(op_u(-lv:nx+lv,-lv:nz+lv),&
                op_v(-lv:nx+lv,-lv:nz+lv))  

       op_u(:,:)=0.e0
       op_v(:,:)=0.e0

       allocate(op_u_1st(-lv:nx+lv,-lv:nz+lv),&
                op_v_1st(-lv:nx+lv,-lv:nz+lv))  

       op_u_1st(:,:)=0.e0
       op_v_1st(:,:)=0.e0
!!     
       allocate(seis_sou(1:nx_rec,1:nt),seis_rcv(1:nx_rec,1:nt))
       
       seis_sou(:,:)=0.e0
       seis_rcv(:,:)=0.e0

       allocate(seis_hb_sou(1:nx_rec,1:nt),&
                seis_hb_rcv(1:nx_rec,1:nt))
       
       seis_hb_sou(:,:)=0.e0
       seis_hb_rcv(:,:)=0.e0
!! angles 
       allocate(op_angle(1:nx,1:nz),op_angle_1st(1:nx,1:nz))
       op_angle(:,:)=0.e0
       op_angle_1st(:,:)=0.e0

       allocate(op_angle_save(1:nx,1:nz),&
                op_angle_1st_save(1:nx,1:nz))
       op_angle_save(:,:)=0.e0
       op_angle_1st_save(:,:)=0.e0

       allocate(sou_op_angle_m(1:max_decomp,1:nx,1:nz),&
                rcv_op_angle_m(1:max_decomp,1:nx,1:nz))
       sou_op_angle_m(:,:,:)=0.e0
       rcv_op_angle_m(:,:,:)=0.e0

       allocate(sou_op_angle_m_tmp(1:nx,1:nz),&
                rcv_op_angle_m_tmp(1:nx,1:nz))
       sou_op_angle_m_tmp(:,:)=0.e0
       rcv_op_angle_m_tmp(:,:)=0.e0

       allocate(sou_p_bin(1:max_decomp,-lv:nx+lv,-lv:nz+lv),&
                rcv_p_bin(1:max_decomp,-lv:nx+lv,-lv:nz+lv))
       sou_p_bin(:,:,:)=0.e0
       rcv_p_bin(:,:,:)=0.e0
 
       allocate(sou_p_bin_tmp(-lv:nx+lv,-lv:nz+lv),&
                rcv_p_bin_tmp(-lv:nx+lv,-lv:nz+lv))
       sou_p_bin_tmp(:,:)=0.e0
       rcv_p_bin_tmp(:,:)=0.e0

       allocate(max_p(1:nx,1:nz))
       max_p(:,:)=0.e0

!! up/down/left/right
       allocate(sou_p_up(-lv:nx+lv,-lv:nz+lv),sou_p_dw(-lv:nx+lv,-lv:nz+lv),&
                sou_p_lf(-lv:nx+lv,-lv:nz+lv),sou_p_rt(-lv:nx+lv,-lv:nz+lv))

       allocate(sou_p_up_old(-lv:nx+lv,-lv:nz+lv),&
                sou_p_dw_old(-lv:nx+lv,-lv:nz+lv))

       allocate(rcv_p_up(-lv:nx+lv,-lv:nz+lv),rcv_p_dw(-lv:nx+lv,-lv:nz+lv),&
                rcv_p_lf(-lv:nx+lv,-lv:nz+lv),rcv_p_rt(-lv:nx+lv,-lv:nz+lv))

       allocate(rcv_p_up_old(-lv:nx+lv,-lv:nz+lv),&
                rcv_p_dw_old(-lv:nx+lv,-lv:nz+lv))

       sou_p_up(:,:)=0.e0
       sou_p_dw(:,:)=0.e0
       sou_p_lf(:,:)=0.e0
       sou_p_rt(:,:)=0.e0
       
       sou_p_up_old(:,:)=0.e0
       sou_p_dw_old(:,:)=0.e0

       rcv_p_up(:,:)=0.e0
       rcv_p_dw(:,:)=0.e0
       rcv_p_lf(:,:)=0.e0
       rcv_p_rt(:,:)=0.e0

       rcv_p_up_old(:,:)=0.e0
       rcv_p_dw_old(:,:)=0.e0

       allocate(sou_op_angle_up(1:nx,1:nz),sou_op_angle_dw(1:nx,1:nz),&
                rcv_op_angle_up(1:nx,1:nz),rcv_op_angle_dw(1:nx,1:nz))

       sou_op_angle_up(:,:)=0.e0
       sou_op_angle_dw(:,:)=0.e0
       rcv_op_angle_up(:,:)=0.e0
       rcv_op_angle_dw(:,:)=0.e0

!!   lu,p_ld,p_ru,p_rd
       allocate(sou_p_lu(-lv:nx+lv,-lv:nz+lv),sou_p_ld(-lv:nx+lv,-lv:nz+lv),&
                sou_p_ru(-lv:nx+lv,-lv:nz+lv),sou_p_rd(-lv:nx+lv,-lv:nz+lv))

       allocate(rcv_p_lu(-lv:nx+lv,-lv:nz+lv),rcv_p_ld(-lv:nx+lv,-lv:nz+lv),&
                rcv_p_ru(-lv:nx+lv,-lv:nz+lv),rcv_p_rd(-lv:nx+lv,-lv:nz+lv))

       sou_p_lu(:,:)=0.e0
       sou_p_ld(:,:)=0.e0
       sou_p_ru(:,:)=0.e0
       sou_p_rd(:,:)=0.e0
       
       rcv_p_lu(:,:)=0.e0
       rcv_p_ld(:,:)=0.e0
       rcv_p_ru(:,:)=0.e0
       rcv_p_rd(:,:)=0.e0

!! cross-correlated image
       allocate(image(1:nx,1:nz),image_illu(1:nx,1:nz))
       allocate(image_cross(1:nx,1:nz))
       allocate(image_adcig(1:nthet,1:nx,1:nz))

       image(:,:)=0.0
       image_illu(:,:)=0.0
       image_cross(:,:)=0.0
       image_adcig(:,:,:)=0.0

       allocate(k_x(1:nx),a_x(1:nx,1:nz),b_x(1:nx,1:nz))
       allocate(k_x_half(1:nx),a_x_half(1:nx,1:nz),b_x_half(1:nx,1:nz))
 
       k_x=1.0
       a_x=zero
       b_x=zero

       k_x_half=1.0
       a_x_half=zero
       b_x_half=zero

       allocate(k_z(1:nz),a_z(1:nx,1:nz),b_z(1:nx,1:nz))
       allocate(k_z_half(1:nz),a_z_half(1:nx,1:nz),b_z_half(1:nx,1:nz))

       k_z=1.0
       a_z=zero
       b_z=zero

       k_z_half=1.0
       a_z_half=zero
       b_z_half=zero

       allocate(c11_u(1:nx,1:nz))  
       c11_u=zero

       allocate(tau11(1:nrel,1:nx,1:nz))
       allocate(tau11_eps(1:nrel,1:nx,1:nz))
       allocate(tau_sigma(1:nrel,1:nx,1:nz))

       tau11=zero
       tau11_eps=zero
       tau_sigma=zero

       allocate(x1(1:nrel,1:nx,1:nz),x2(1:nrel,1:nx,1:nz))
       allocate(dens(1:nx,1:nz))
       allocate(vp(1:nx,1:nz))

       x1=zero
       x2=zero
       dens=zero
       vp=zero

       allocate(vp_f0(1:nx,1:nz),c11_r(1:nx,1:nz))

       vp_f0(:,:)=zero
       c11_r(:,:)=zero

       allocate(vx(-lv:nx+lv,-lv:nz+lv),vz(-lv:nx+lv,-lv:nz+lv))
       allocate(p(-lv:nx+lv,-lv:nz+lv))
       allocate(r(1:nrel,-lv:nx+lv,-lv:nz+lv))
 
       vx=zero
       vz=zero
       p=zero
       r=zero

       allocate(memory_dvx_dx(1:nx,1:nz))
       allocate(memory_dvz_dz(1:nx,1:nz))
       allocate(memory_dp_dx(1:nx,1:nz))
       allocate(memory_dp_dz(1:nx,1:nz))

       memory_dvx_dx=zero
       memory_dvz_dz=zero
       memory_dp_dx=zero
       memory_dp_dz=zero

       allocate(vx_hb(-lv:nx+lv,-lv:nz+lv),vz_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(p_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(r_hb(1:nrel,-lv:nx+lv,-lv:nz+lv))
 
       vx_hb=zero
       vz_hb=zero
       p_hb=zero
       r_hb=zero

       allocate(memory_dvx_dx_hb(1:nx,1:nz))
       allocate(memory_dvz_dz_hb(1:nx,1:nz))
       allocate(memory_dp_dx_hb(1:nx,1:nz))
       allocate(memory_dp_dz_hb(1:nx,1:nz))

       memory_dvx_dx_hb=zero
       memory_dvz_dz_hb=zero
       memory_dp_dx_hb=zero
       memory_dp_dz_hb=zero

       allocate(c(1:lv),dlh(1:lv),dlv(1:lv))
       allocate(epsilon_m(-lv:nx+lv,-lv:nz+lv))
 
       c=zero
       dlh=zero
       dlv=zero 
       epsilon_m=zero

!       allocate(dp_dt(-lv:nx+lv,-lv:nz+lv))
!       dp_dt(:,:)=zero

       open(12,file=vpname,access='direct',form='unformatted',&
            recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           print *,'i','vp',i
           read(12,rec=(i_shot-1)*int_shot_model+i-nx_p1+1) &
                       (vp_f0(i,j),j=nz_p1,nz_p2)
           print *,'i','vp',i
        enddo 
       close(12)        

       open(12,file=rhoname,access='direct',form='unformatted',&
            recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           read(12,rec=(i_shot-1)*int_shot_model+i-nx_p1+1) &
                       (dens(i,j),j=nz_p1,nz_p2)
        enddo 
       close(12)        

!!!    if attenuation 
       if(.not.lossless) then 

          open(12,file=tau_epsname,access='direct',form='unformatted',&
               recl=ii_bi*nz_p)
           do l=1,nrel
              do i=nx_p1,nx_p2
                 read(12,rec=(l-1)*mx+(i_shot-1)*int_shot_model+i-nx_p1+1) &
                 (tau11_eps(l,i,j),j=nz_p1,nz_p2)
              enddo 
           enddo 
          close(12)   
    
          open(12,file=tau_sigmaname,access='direct',form='unformatted',&
               recl=ii_bi*nz_p)
           do l=1,nrel
              do i=nx_p1,nx_p2
                 read(12,rec=(l-1)*mx+(i_shot-1)*int_shot_model+i-nx_p1+1) &
                 (tau_sigma(l,i,j),j=nz_p1,nz_p2)
              enddo 
           enddo
          close(12) 
       endif

!!!    if no attenuation 
       if(lossless) then 

         tau_sigma(:,:,:)=1.e0

         do l=1,nrel
            do j=nz_p1,nz_p2
               do i=nx_p1,nx_p2
                  tau11_eps(l,i,j)=tau_sigma(l,i,j)
               enddo
            enddo
         enddo
       endif

!!     assign values for pml layers
!!     1:npoints_pml,nx-npoints_pml+1:nx
       call adding_abs_val(vp_f0,dens,tau11_eps,tau_sigma,nx,nz,nrel,npoints_pml)

       call vp_f0_2_vp(nx,nz,lv,nrel,f_ref,c11_r,vp_f0,vp,dens,&
                       tau11_eps,tau_sigma)

!! calculate intermediate variables
       call assign_parameter(nx,nz,lv,nrel,deltat,& 
                            tau11_eps,tau11,tau_sigma,x1,x2,&
                            c11_r,dens,c11_u)
!! call ricker wavelet
!       call rickerfunc(deltat,f0,length,nt,ricker)
       call rickerfunc_w_hb(deltat,f0,length,nt,ricker,ricker_hb)
 
!! display size of the model
       print *
       print *,'NX = ',NX
       print *,'Nz = ',Nz
       print *
       print *,'size of the model along X = ',(NX - 1) * DELTAX
       print *,'size of the model along Y = ',(Nz - 1) * DELTAz
       print *
       print *,'Total number of grid points = ',NX * Nz
       print *
       print *

       call cpml_coef(a_x,b_x,k_x,a_z,b_z,k_z,&
             a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
             npoints_pml,nx,nz,deltax,deltaz,deltat,vp,f0)

       c(1)=1.23538628085693269
       c(2)=-0.109160358903484037
       c(3)=0.246384765606012697d-1
       c(4)=-0.660472512788868975d-2

!!!    Holberg FD coefs
       c(1)=1.231666
       c(2)=-0.1041182
       c(3)=0.02063707
       c(4)=-0.003570998

!!!!!!!!!!!!!!!!! Stability Condition !!!!!!!!!!!!!!
      cons1=sqrt(1.0/deltax**2+1.0/deltaz**2)
!
      sum_a=0.0
!
      do i=1,lv
         sum_a=sum_a+abs(c(i))
      enddo 
!
      sum_a=1.0/sum_a
!
      do i=1,nx
         do j=1,nz
          deltat_v=deltat*vp(i,j)
          cons2=deltat_v*cons1
          if(cons2.gt.sum_a) then 
              print *,'Stability condition not reached !!!!'
              stop
          endif
         enddo
      enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do i=1,lv
          dlh(i)=c(i)/deltax
          dlv(i)=c(i)/deltaz
          print *,dlh(i),dlv(i),'dlh','dlv'
       enddo 

!!!    file names for pressure seismograms 
       file_p_sou=trim(file_p_sou_part)//char(i_shot/100+48)//&
               char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
!       file_p_sou=trim(file_p_sou_part)

!!!    file names for snapshots of the whole time history
       file_snap_all=trim(file_snap_all_part)//char((myid+1)/100+48)//&
               char(mod((myid+1),100)/10+48)//char(mod((myid+1),10)+48)

!!!    record # 
!!     nx_rec=rcv_l-rcv_f+1
!!     contains both primary and multiples
       open(file_p_sou_no,file=file_p_sou,access='direct',form='unformatted'&
            ,recl=ii_bi*nx_rec)
        do j=1,nt
          read(file_p_sou_no,rec=j) (seis_sou(i,j),i=1,nx_rec)
        enddo 
       close(file_p_sou_no)

!       open(file_snap_all_no,file=file_snap_all,access='direct',&
!            form='unformatted',recl=ii_bi*nx_p*nz_p)

       do it = nt1,nt

         if(multiple) then 
!!     Multiple rtm 
!!     read in seismograms
             do i=1,nx_rec
                p(rcv_f+int_rcv*(i-1),i_zr)=p(rcv_f+int_rcv*(i-1),i_zr)+seis_sou(i,it)
                
                p_hb(rcv_f+int_rcv*(i-1),i_zr)=p_hb(rcv_f+int_rcv*(i-1),i_zr)+seis_hb_sou(i,it)
             enddo 
           else
!!       if else
!!       traditional primary rtm
             if(it.lt.50000) then 
                amp=ricker(it)
                p(isource,jsource) = p(isource,jsource) + amp
                amp_hb=ricker_hb(it)
                p_hb(isource,jsource) = p_hb(isource,jsource) + amp_hb
              endif
           endif
!!         endif multiple
         
!!  for monitoring 
          print *,p(isource,jsource),vx(isource,jsource),&
                  vz(isource,jsource),'dt_snap',dt_snap,&
                 'haha','it',it,'myid',myid,'i_shot',i_shot

!!   wavefield reconstruction, removed since it is proved that w r is 
!!   not applicable for viscous media

!!   save boundary values for wavefield reconstruction
          call bvr_record_wave_vxz(it,nt,vx,vz,p,vx_l,vx_r,vx_t,vx_b,vz_l,vz_r,vz_t,vz_b,&
               nx,nz,lv,npoints_pml)

          call bvr_record_wave_vxz(it,nt,vx_hb,vz_hb,p_hb,vx_l_hb,vx_r_hb,vx_t_hb,vx_b_hb,&
               vz_l_hb,vz_r_hb,vz_t_hb,vz_b_hb,&
               nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
!             if(.not.snap_all) then
                 last_vx_hb(:,:)=vx_hb(:,:)
                 last_vz_hb(:,:)=vz_hb(:,:)
                 last_vx(:,:)=vx(:,:)
                 last_vz(:,:)=vz(:,:)
!             endif
          endif

!!  kernel subroutine, update wavefield variables 
!!------------------------------------------------------------
!! compute stress sigma and update memory variables for C-PML
!!------------------------------------------------------------
          call iso_visco_step(vx,vz,epsilon_m,p,r,memory_dvx_dx,memory_dvz_dz,&
               memory_dp_dx,memory_dp_dz,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

          call iso_visco_step(vx_hb,vz_hb,epsilon_m,p_hb,r_hb,memory_dvx_dx_hb,memory_dvz_dz_hb,&
               memory_dp_dx_hb,memory_dp_dz_hb,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)
          print *,p(isource,jsource),p_hb(isource,jsource),'isource'

!! save pressure
          call bvr_record_wave_p(it,nt,vx,vz,p,p_l,p_r,p_t,p_b,&
               nx,nz,lv,npoints_pml)

          call bvr_record_wave_p(it,nt,vx_hb,vz_hb,p_hb,p_l_hb,p_r_hb,p_t_hb,p_b_hb,&
               nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
!             if(.not.snap_all) then
                 last_p(:,:)=p(:,:)
                 last_p_hb(:,:)=p_hb(:,:)
!             endif
          endif

!!   write snapshots for a time interval
          if(mod(it,nt_snap).eq.0) then
             call snap_shot(p,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
!             call snap_shot(op_u,it,nx,nz,lv,nt_snap,snap_part1,file_snap_no)
!             call snap_shot(op_v,it,nx,nz,lv,nt_snap,snap_part2,file_snap_no)
!             call snap_shot(op_u_1st,it,nx,nz,lv,nt_snap,snap_part3,file_snap_no)
!             call snap_shot(op_v_1st,it,nx,nz,lv,nt_snap,snap_part4,file_snap_no)
          endif

       enddo   ! end of time loop

!      close(file_seis_no)
!       close(file_snap_all_no)

       if(.false.) then
       open(12,file='op_angle',access='direct',form='unformatted',&
               recl=ii_bi*nz)
         do i=1,nx
            write(12,rec=i) (op_angle_save(i,j),j=1,nz)
         enddo 
       close(12)

       open(12,file='op_angle_1st',access='direct',form='unformatted',&
               recl=ii_bi*nz)
         do i=1,nx
            write(12,rec=i) (op_angle_1st_save(i,j),j=1,nz)
         enddo 
       close(12)
       endif

!!     now it's rtm's turn

       modeling=.false.

       deltat=-deltat
       c11_u=-c11_u

       vx(:,:)=zero
       vz(:,:)=zero
       p(:,:)=zero
       r(:,:,:)=zero
       memory_dvx_dx=zero
       memory_dvz_dz=zero
       memory_dp_dx=zero
       memory_dp_dz=zero
   
       vx_hb(:,:)=zero
       vz_hb(:,:)=zero
       p_hb(:,:)=zero
       r_hb(:,:,:)=zero
       memory_dvx_dx_hb=zero
       memory_dvz_dz_hb=zero
       memory_dp_dx_hb=zero
       memory_dp_dz_hb=zero
   
       sou_vx(:,:)=zero
       sou_vz(:,:)=zero
       sou_p(:,:)=zero
       sou_r(:,:,:)=zero
       sou_memory_dvx_dx=zero
       sou_memory_dvz_dz=zero
       sou_memory_dp_dx=zero
       sou_memory_dp_dz=zero

       sou_vx_hb(:,:)=zero
       sou_vz_hb(:,:)=zero
       sou_p_hb(:,:)=zero
       sou_r_hb(:,:,:)=zero
       sou_memory_dvx_dx_hb=zero
       sou_memory_dvz_dz_hb=zero
       sou_memory_dp_dx_hb=zero
       sou_memory_dp_dz_hb=zero

!!     compute array dimensions
       lwx=num_local_x
       lwz=num_local_z

       padpercx=2.0
       padpercz=2.0
       tappercx=0.0
       tappercz=0.0

       num_angle=360.0

!       lwx2=2*lwx
!       lwz2=2*lwz
       lwx2=lwx
       lwz2=lwz

       nwx=max(ceiling(nx_p/real(lwx))-1,1)
       nwz=max(ceiling(nz_p/real(lwz))-1,1)
!       nwx=max(ceiling(nx_p/real(lwx2)),1)
!       nwz=max(ceiling(nz_p/real(lwz2)),1)

!! padding size
       call npad(lwx2,int(lwx2*padpercx),nxfft)
       call npad(lwz2,int(lwz2*padpercz),nzfft)

!       allocate(p_opti(1:num_bin,-lv:nx+lv,-lv:nz+lv))
!       p_opti(:,:,:)=zero

       allocate(kx_vec(1:nxfft),kz_vec(1:nzfft),k_wav(1:nxfft,1:nzfft))
       kx_vec(:)=zero
       kz_vec(:)=zero
       k_wav(:,:)=zero

!!     set wavenumber 
!!     
       call set_wavenumber(k_wav,kx_vec,kz_vec,nxfft,nzfft,deltax,deltaz)

!       propg_bin_scale=int(360/num_bin)

       allocate(filter_2d(1:nxfft,1:nzfft))
       filter_2d(:,:)=1.0

       allocate(angle_filter(1:propg_bin_scale))
       angle_filter(:)=1.0

!! build wavenumber filter
       allocate(kx_filter(1:nxfft),kz_filter(1:nzfft))
       kx_filter(:)=1.0
       kz_filter(:)=1.0

       value_eps=1.0E-6
       napp=int(nxfft/5)
       call gaussian_2d(filter_2d,kx_filter,kz_filter,nxfft,nzfft,napp,napp,value_eps)

       napp=int(propg_bin_scale/4)
       call gaussian_1d(angle_filter,propg_bin_scale,napp,value_eps)

!       open(file_snap_all_no,file=file_snap_all,access='direct',&
!            form='unformatted',recl=ii_bi*nx_p*nz_p)
!
!       nx_rec=rcv_l-rcv_f+1

       file_p_rcv=trim(file_p_rcv_part)//char(i_shot/100+48)//&
               char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
!       file_p_rcv=trim(file_p_rcv_part)

       open(file_p_rcv_no,file=file_p_rcv,access='direct',form='unformatted'&
            ,recl=ii_bi*nx_rec)
        do j=1,nt
          read(file_p_rcv_no,rec=j) (seis_rcv(i,j),i=1,nx_rec)
        enddo
       close(file_p_rcv_no)

!!  call hilbert transform
       call hilbert_seismo(seis_rcv,seis_hb_rcv,nx_rec,nt)

       do it=nt,nt1,-1
   
!!     read in seismograms 
          do i=1,nx_rec
             p(rcv_f+int_rcv*(i-1),i_zr)=p(rcv_f+int_rcv*(i-1),i_zr)+seis_rcv(i,it)

             p_hb(rcv_f+int_rcv*(i-1),i_zr)=p_hb(rcv_f+int_rcv*(i-1),i_zr)+seis_hb_rcv(i,it)
          enddo

!!     receiver adjoint operation
          call iso_visco_step(vx,vz,epsilon_m,p,r,memory_dvx_dx,memory_dvz_dz,&
               memory_dp_dx,memory_dp_dz,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

!!     hilbert receiver adjoint operation 
          call iso_visco_step(vx_hb,vz_hb,epsilon_m,p_hb,r_hb,memory_dvx_dx_hb,memory_dvz_dz_hb,&
               memory_dp_dx_hb,memory_dp_dz_hb,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

!!     source wavefield reconstruction
          call iso_visco_step(sou_vx,sou_vz,epsilon_m,sou_p,sou_r,sou_memory_dvx_dx,sou_memory_dvz_dz,&
               sou_memory_dp_dx,sou_memory_dp_dz,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

!!      hilbert
          call iso_visco_step(sou_vx_hb,sou_vz_hb,epsilon_m,sou_p_hb,sou_r_hb,&
               sou_memory_dvx_dx_hb,sou_memory_dvz_dz_hb,&
               sou_memory_dp_dx_hb,sou_memory_dp_dz_hb,dlh,dlv,a_x,b_x,k_x,a_z,b_z,k_z,&
               a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
               c11_u,x1,x2,tau11,dens,nx,nz,nrel,lv,deltat,modeling)

!!      boundary value wavefield reconstruction
            do i=1,nz
               do j=1,lv
!!      left
               sou_vx(npoints_pml+j,i)=vx_l(i,j,it)
               sou_vz(npoints_pml+j,i)=vz_l(i,j,it)
               sou_p(npoints_pml+j,i)=p_l(i,j,it)

               sou_vx_hb(npoints_pml+j,i)=vx_l_hb(i,j,it)
               sou_vz_hb(npoints_pml+j,i)=vz_l_hb(i,j,it)
               sou_p_hb(npoints_pml+j,i)=p_l_hb(i,j,it)
!!      right
               sou_vx(nx-npoints_pml-lv+j-1,i)=vx_r(i,j,it)
               sou_vz(nx-npoints_pml-lv+j-1,i)=vz_r(i,j,it)
               sou_p(nx-npoints_pml-lv+j-1,i)=p_r(i,j,it)

               sou_vx_hb(nx-npoints_pml-lv+j-1,i)=vx_r_hb(i,j,it)
               sou_vz_hb(nx-npoints_pml-lv+j-1,i)=vz_r_hb(i,j,it)
               sou_p_hb(nx-npoints_pml-lv+j-1,i)=p_r_hb(i,j,it)
               enddo
             enddo

             do i=1,nx
                do j=1,lv
!!       top
                sou_vx(i,npoints_pml+j)=vx_t(i,j,it)
                sou_vz(i,npoints_pml+j)=vz_t(i,j,it)
                sou_p(i,npoints_pml+j)=p_t(i,j,it)

                sou_vx_hb(i,npoints_pml+j)=vx_t_hb(i,j,it)
                sou_vz_hb(i,npoints_pml+j)=vz_t_hb(i,j,it)
                sou_p_hb(i,npoints_pml+j)=p_t_hb(i,j,it)
!!       bottom
                sou_vx(i,nz-npoints_pml-lv+j-1)=vx_b(i,j,it)
                sou_vz(i,nz-npoints_pml-lv+j-1)=vz_b(i,j,it)
                sou_p(i,nz-npoints_pml-lv+j-1)=p_b(i,j,it)

                sou_vx_hb(i,nz-npoints_pml-lv+j-1)=vx_b_hb(i,j,it)
                sou_vz_hb(i,nz-npoints_pml-lv+j-1)=vz_b_hb(i,j,it)
                sou_p_hb(i,nz-npoints_pml-lv+j-1)=p_b_hb(i,j,it)
                enddo
             enddo
 
            if(it.eq.nt) then 
              sou_p(:,:)=last_p(:,:)
              sou_vx(:,:)=last_vx(:,:)
              sou_vz(:,:)=last_vz(:,:)

              sou_p_hb(:,:)=last_p_hb(:,:)
              sou_vx_hb(:,:)=last_vx_hb(:,:)
              sou_vz_hb(:,:)=last_vz_hb(:,:)
            endif

!!         
            if(mod(it,nt_snap).eq.0) then 
               recons_part='snap_recons'
               call snap_shot(sou_p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_rcv'
               call snap_shot(p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
            endif

!! local wavefield directional decomposition
!! with optical flow and adcig gathers

!          if(local_decomp) then           
!            nitermax=6

!! source wavefield decomposition
           call wavefield_decomp(sou_p,sou_p_hb,sou_p_up,sou_p_dw,sou_p_lf,sou_p_rt,&
                sou_p_lu,sou_p_ld,sou_p_ru,sou_p_rd,nx,nz,lv)

!! receiver wavefield decomposition 
           call wavefield_decomp(p,p_hb,rcv_p_up,rcv_p_dw,rcv_p_lf,rcv_p_rt,&
                rcv_p_lu,rcv_p_ld,rcv_p_ru,rcv_p_rd,nx,nz,lv)

!! call optical flow for source wavefield

           call optical_flow(sou_p_up,sou_p_up_old,op_u,op_v,op_u_1st,op_v_1st,&
                alpha,nx,nz,lv,nitermax,it)

           sou_p_up_old(-lv:nx+lv,-lv:nz+lv)=sou_p_up(-lv:nx+lv,-lv:nz+lv)

           i_type=1
           call calc_angle(op_u,op_v,sou_op_angle_up,nx,nz,lv,i_type)
 
           call optical_flow(sou_p_dw,sou_p_dw_old,op_u,op_v,op_u_1st,op_v_1st,&
                alpha,nx,nz,lv,nitermax,it)

           sou_p_dw_old(-lv:nx+lv,-lv:nz+lv)=sou_p_dw(-lv:nx+lv,-lv:nz+lv)

           call calc_angle(op_u,op_v,sou_op_angle_dw,nx,nz,lv,i_type)

!! call optical flow for receiver wavefield 
           call optical_flow(rcv_p_up,rcv_p_up_old,op_u,op_v,op_u_1st,op_v_1st,&
                alpha,nx,nz,lv,nitermax,it)

           rcv_p_up_old(-lv:nx+lv,-lv:nz+lv)=rcv_p_up(-lv:nx+lv,-lv:nz+lv)

           i_type=2
           call calc_angle(op_u,op_v,rcv_op_angle_up,nx,nz,lv,i_type)
 
           call optical_flow(rcv_p_dw,rcv_p_dw_old,op_u,op_v,op_u_1st,op_v_1st,&
                alpha,nx,nz,lv,nitermax,it)

           rcv_p_dw_old(-lv:nx+lv,-lv:nz+lv)=rcv_p_dw(-lv:nx+lv,-lv:nz+lv)

           call calc_angle(op_u,op_v,rcv_op_angle_dw,nx,nz,lv,i_type)

!! organize into bin
           sou_p_bin(1,-lv:nx+lv,-lv:nz+lv)=sou_p_up(-lv:nx+lv,-lv:nz+lv)
           sou_p_bin(2,-lv:nx+lv,-lv:nz+lv)=sou_p_dw(-lv:nx+lv,-lv:nz+lv)
           rcv_p_bin(1,-lv:nx+lv,-lv:nz+lv)=rcv_p_up(-lv:nx+lv,-lv:nz+lv)
           rcv_p_bin(2,-lv:nx+lv,-lv:nz+lv)=rcv_p_dw(-lv:nx+lv,-lv:nz+lv)

           sou_op_angle_m(1,1:nx,1:nz)=sou_op_angle_up(1:nx,1:nz)
           sou_op_angle_m(2,1:nx,1:nz)=sou_op_angle_dw(1:nx,1:nz)
           rcv_op_angle_m(1,1:nx,1:nz)=rcv_op_angle_up(1:nx,1:nz)
           rcv_op_angle_m(2,1:nx,1:nz)=rcv_op_angle_dw(1:nx,1:nz)

           n_dec=2
           call imaging_adcig_udlr(sou_p_bin,rcv_p_bin,sou_op_angle_m,rcv_op_angle_m,&
                image_adcig,image_cross,nx,nz,lv,max_decomp,n_dec,nthet)           
 
           do j=nz_p1,nz_p2
              do i=nx_p1,nx_p2
                 image_illu(i,j)=image_illu(i,j)+sou_p(i,j)*sou_p(i,j)
              enddo
           enddo
 
            if(mod(it,nt_snap).eq.0) then 
               recons_part='snap_sou_p_bin1'
               sou_p_bin_tmp(:,:)=sou_p_bin(1,:,:)
               call snap_shot(sou_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_sou_p_bin2'
               sou_p_bin_tmp(:,:)=sou_p_bin(2,:,:)
               call snap_shot(sou_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_sou_p_bin3'
               sou_p_bin_tmp(:,:)=sou_p_bin(3,:,:)
               call snap_shot(sou_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_sou_p_bin4'
               sou_p_bin_tmp(:,:)=sou_p_bin(4,:,:)
               call snap_shot(sou_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     

               recons_part='snap_rcv_p_bin1'
               rcv_p_bin_tmp(:,:)=rcv_p_bin(1,:,:)
               call snap_shot(rcv_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_rcv_p_bin2'
               rcv_p_bin_tmp(:,:)=rcv_p_bin(2,:,:)
               call snap_shot(rcv_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_rcv_p_bin3'
               rcv_p_bin_tmp(:,:)=rcv_p_bin(3,:,:)
               call snap_shot(rcv_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
               recons_part='snap_rcv_p_bin4'
               rcv_p_bin_tmp(:,:)=rcv_p_bin(4,:,:)
               call snap_shot(rcv_p_bin_tmp,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
            endif

           call imaging_condition(sou_p,p,image,nx,nz,lv,lv)

       enddo 

!       close(file_snap_all_no)
       close(file_p_no)

       file_image='seis_p_image_cross'//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           write(file_image_no,rec=i-nx_p1+1) (image_cross(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       file_image='seis_p_adcig'//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           do i_angle=1,nthet
             write(file_image_no,rec=(i-nx_p1)*nthet+i_angle) &
                     (image_adcig(i_angle,i,j),j=nz_p1,nz_p2)
           enddo
        enddo
       close(file_image_no)

!!     cross-correlation 
       file_image=trim(file_image_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
!!     source-normalized cross-correlation 
       file_image_normal=trim(file_image_normal_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
!!     illumination term
       file_image_illu=trim(file_image_illu_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)

       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           write(file_image_no,rec=i-nx_p1+1) (image(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       image(nx_p1:nx_p2,nz_p1:nz_p2)=image(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       open(file_image_no,file=file_image_normal,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           write(file_image_no,rec=i-nx_p1+1) (image(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       open(file_image_no,file=file_image_illu,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           write(file_image_no,rec=i-nx_p1+1) (image_illu(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       print *
       print *,'End of the simulation'
       print *

       print *,'i_shot',i_shot

       deallocate(image,image_illu,image_cross,image_adcig)
       deallocate(k_x,a_x,b_x,k_x_half,a_x_half,b_x_half)
       deallocate(k_z,a_z,b_z,k_z_half,a_z_half,b_z_half)
       deallocate(c11_u,dens,vp)
       deallocate(tau11,tau11_eps,tau_sigma)
       deallocate(x1,x2)
       deallocate(vp_f0,c11_r)
       deallocate(vx,vz,p,r)
       deallocate(memory_dvx_dx,memory_dvz_dz)
       deallocate(memory_dp_dx,memory_dp_dz)
       deallocate(vx_hb,vz_hb,p_hb,r_hb)
       deallocate(memory_dvx_dx_hb,memory_dvz_dz_hb)
       deallocate(memory_dp_dx_hb,memory_dp_dz_hb)
       deallocate(c,dlh,dlv)
       deallocate(epsilon_m)
       deallocate(op_u,op_v,op_u_1st,op_v_1st)
       deallocate(op_angle,op_angle_1st)
       deallocate(old_p,max_p)
       deallocate(op_angle_save,op_angle_1st_save)
       deallocate(sou_op_angle_m,rcv_op_angle_m)
       deallocate(sou_op_angle_m_tmp,rcv_op_angle_m_tmp)
       deallocate(seis_sou,seis_rcv)

       deallocate(sou_p,sou_vx,sou_vz,sou_r)
       deallocate(sou_memory_dvx_dx,sou_memory_dvz_dz)
       deallocate(sou_memory_dp_dx,sou_memory_dp_dz)
       deallocate(vx_l,vx_r,vx_t,vx_b)
       deallocate(vz_l,vz_r,vz_t,vz_b)
       deallocate(p_l,p_r,p_t,p_b)

       print *,'test2'
       deallocate(sou_p_hb,sou_vx_hb,sou_vz_hb,sou_r_hb)
       deallocate(sou_memory_dvx_dx_hb,sou_memory_dvz_dz_hb)
       deallocate(sou_memory_dp_dx_hb,sou_memory_dp_dz_hb)
       deallocate(vx_l_hb,vx_r_hb,vx_t_hb,vx_b_hb)
       deallocate(vz_l_hb,vz_r_hb,vz_t_hb,vz_b_hb)
       deallocate(p_l_hb,p_r_hb,p_t_hb,p_b_hb)

       print *,'test3'
       deallocate(last_p,last_vx,last_vz)
       deallocate(filter_2d,kx_filter,kz_filter,angle_filter)
!       deallocate(dp_dt)
       deallocate(sou_p_bin,rcv_p_bin)
       deallocate(sou_p_bin_tmp,rcv_p_bin_tmp)

       deallocate(sou_p_up,sou_p_dw,sou_p_lf,sou_p_rt,&
                  rcv_p_up,rcv_p_dw,rcv_p_lf,rcv_p_rt) 
       deallocate(sou_p_up_old,sou_p_dw_old,&
                  rcv_p_up_old,rcv_p_dw_old)
       deallocate(sou_op_angle_up,sou_op_angle_dw,&
                  rcv_op_angle_up,rcv_op_angle_dw)
       deallocate(sou_p_lu,sou_p_ld,sou_p_ru,sou_p_rd)
       deallocate(rcv_p_lu,rcv_p_ld,rcv_p_ru,rcv_p_rd)

       print *,'test4'
    enddo   !! shot domain

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

  end program seismic_CPML_2D_iso_rtm

