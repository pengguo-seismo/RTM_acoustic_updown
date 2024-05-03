!

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
    use omp_lib
    use precision_m
    use main_module
    use parameter_input
    use wave_propagator    
    use wave_module
    use wave_drec_sep_module
    use wave_hilbert_module
    use fftw_module_f
!---
    implicit none
    include 'mpif.h'
    include 'fftw3.f'

    integer::ii,n_count

    integer::nxfft,nzfft
    integer::nxpad,nzpad
    integer::flag
    integer::n_of_fft
    real::padpercx,padpercz

    real,allocatable::taper_x(:),taper_z(:)
    real,allocatable::tapx_12(:),tapz_12(:)

!    integer*8::plan_ud,plan_ud1,plan_ud2
!    integer*8::plan_lr,plan_lr1,plan_lr2

!    integer*8::plan_bin,plan_bin1,plan_bin2,plan_bin3,plan_bin4

!--- program starts here
!---
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

     ii_bi=1

     flag=1

!!   read in parameter, to decide the parmeter file name in subroutine par_in
!!   0: forward modeling, parameter file: parameter.txt
!!   1: forward modeling for direct waves, parameter file: parameter_dir.txt
!!   2: rtm, parameter file: parameter_rtm.txt
      call getarg(1,arg)
      arg=trim(arg)
      read(arg,'(i4)') iflag_par

      call par_in(iflag_par)

      file_p_sou_no=21
      file_p_rcv_no=25
      file_snap_all_no=22
      file_image_no=23
      file_snap_no=24

      call allocate_main_module(nx,nz,nrel,nx_rec,lv)

      allocate(ricker(1:nt))
      ricker=0.0

      allocate(ricker_hb(1:nt))
      ricker_hb=0.0

!! begin allocate wavefield-related variables
      allocate(p(-lv:nx+lv,-lv:nz+lv))
      allocate(vx(-lv:nx+lv,-lv:nz+lv))
      allocate(vz(-lv:nx+lv,-lv:nz+lv))
      allocate(r(1:nrel,-lv:nx+lv,-lv:nz+lv))

      allocate(memory_dvx_dx(1:nx,1:nz))
      allocate(memory_dvz_dz(1:nx,1:nz))
      allocate(memory_dp_dx(1:nx,1:nz))
      allocate(memory_dp_dz(1:nx,1:nz))

      p(:,:)=0.0
      vx(:,:)=0.e0
      vz(:,:)=0.e0
      r(:,:,:)=0.e0
      memory_dvx_dx(:,:)=0.e0
      memory_dvz_dz(:,:)=0.e0
      memory_dp_dx(:,:)=0.e0
      memory_dp_dz(:,:)=0.e0

       allocate(vx_hb(-lv:nx+lv,-lv:nz+lv),vz_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(p_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(r_hb(1:nrel,-lv:nx+lv,-lv:nz+lv))

       vx_hb=0.0
       vz_hb=0.0
       p_hb=0.0
       r_hb=0.0

       allocate(memory_dvx_dx_hb(1:nx,1:nz))
       allocate(memory_dvz_dz_hb(1:nx,1:nz))
       allocate(memory_dp_dx_hb(1:nx,1:nz))
       allocate(memory_dp_dz_hb(1:nx,1:nz))

       memory_dvx_dx_hb=0.0
       memory_dvz_dz_hb=0.0
       memory_dp_dx_hb=0.0
       memory_dp_dz_hb=0.0

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

       allocate(last_p(-lv:nx+lv,-lv:nz+lv))
       allocate(last_vx(-lv:nx+lv,-lv:nz+lv))
       allocate(last_vz(-lv:nx+lv,-lv:nz+lv))

       last_p(:,:)=0.0
       last_vx(:,:)=0.0
       last_vz(:,:)=0.0

!! for hilbert wavefield
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

       allocate(last_p_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(last_vx_hb(-lv:nx+lv,-lv:nz+lv))
       allocate(last_vz_hb(-lv:nx+lv,-lv:nz+lv))

       last_p_hb(:,:)=0.0
       last_vx_hb(:,:)=0.0
       last_vz_hb(:,:)=0.0

       allocate(image(1:nx,1:nz),image_illu(1:nx,1:nz))
       image(:,:)=0.0
       image_illu(:,:)=0.0

       allocate(image_0(1:nx,1:nz),image_01(1:nx,1:nz),&
                image_illu_0(1:nx,1:nz))
 
       image_0=0.0
       image_01=0.0
       image_illu_0=0.0 

       allocate(image_0_tmp(1:nx,1:nz),image_01_tmp(1:nx,1:nz),&
                image_illu_0_tmp(1:nx,1:nz))
 
       image_0_tmp=0.0
       image_01_tmp=0.0
       image_illu_0_tmp=0.0 

       allocate(image_ud(1:nx,1:nz),image_du(1:nx,1:nz),&
                image_uu(1:nx,1:nz),image_dd(1:nx,1:nz))

       allocate(image_duud(1:nx,1:nz),image_dduu(1:nx,1:nz))

       image_ud=0.0
       image_du=0.0
       image_uu=0.0
       image_dd=0.0 
       image_duud=0.0
       image_dduu=0.0

!! summed results
       allocate(image_ud_0(1:nx,1:nz),image_du_0(1:nx,1:nz),&
                image_uu_0(1:nx,1:nz),image_dd_0(1:nx,1:nz))

       allocate(image_duud_0(1:nx,1:nz),image_dduu_0(1:nx,1:nz))

       image_ud_0=0.0
       image_du_0=0.0
       image_uu_0=0.0
       image_dd_0=0.0 
       image_duud_0=0.0
       image_dduu_0=0.0

       allocate(image_ud_0_tmp(1:nx,1:nz),image_du_0_tmp(1:nx,1:nz),&
                image_uu_0_tmp(1:nx,1:nz),image_dd_0_tmp(1:nx,1:nz))

       allocate(image_duud_0_tmp(1:nx,1:nz),image_dduu_0_tmp(1:nx,1:nz))

!! summed results (with approximate hessian conditioned)
       allocate(image_ud_01(1:nx,1:nz),image_du_01(1:nx,1:nz),&
                image_uu_01(1:nx,1:nz),image_dd_01(1:nx,1:nz))

       allocate(image_duud_01(1:nx,1:nz),image_dduu_01(1:nx,1:nz))

       image_ud_01=0.0
       image_du_01=0.0
       image_uu_01=0.0
       image_dd_01=0.0 
       image_duud_01=0.0
       image_dduu_01=0.0

       allocate(image_ud_01_tmp(1:nx,1:nz),image_du_01_tmp(1:nx,1:nz),&
                image_uu_01_tmp(1:nx,1:nz),image_dd_01_tmp(1:nx,1:nz))

       allocate(image_duud_01_tmp(1:nx,1:nz),image_dduu_01_tmp(1:nx,1:nz))

!!!!!!
       allocate(seis_sou(1:nx_rec,1:nt),seis_rcv(1:nx_rec,1:nt))
       
       seis_sou(:,:)=0.e0
       seis_rcv(:,:)=0.e0

       allocate(seis_sou_hb(1:nx_rec,1:nt),&
                seis_rcv_hb(1:nx_rec,1:nt))

       seis_sou_hb=0.0
       seis_rcv_hb=0.0

!! for wavefield decomposition
       allocate(sou_p_up(-lv:nx+lv,-lv:nz+lv),sou_p_dw(-lv:nx+lv,-lv:nz+lv),&
                sou_p_lf(-lv:nx+lv,-lv:nz+lv),sou_p_rt(-lv:nx+lv,-lv:nz+lv))

       allocate(rcv_p_up(-lv:nx+lv,-lv:nz+lv),rcv_p_dw(-lv:nx+lv,-lv:nz+lv),&
                rcv_p_lf(-lv:nx+lv,-lv:nz+lv),rcv_p_rt(-lv:nx+lv,-lv:nz+lv))

       sou_p_up(:,:)=0.e0
       sou_p_dw(:,:)=0.e0
       sou_p_lf(:,:)=0.e0
       sou_p_rt(:,:)=0.e0

       rcv_p_up(:,:)=0.e0
       rcv_p_dw(:,:)=0.e0
       rcv_p_lf(:,:)=0.e0
       rcv_p_rt(:,:)=0.e0

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

       int_shot_model=0

! call ricker wavelet
       call rickerfunc_new_seri(deltat,f0,length,nt,ricker)

       call rickerfunc_w_hb(deltat,f0,length,nt,ricker,ricker_hb)

!! padding parameters
     if(.false.) then
       padpercx=0.06
       padpercz=0.06
       call npad(nx,ceiling(nx*padpercx),nxfft)
       call npad(nz,ceiling(nx*padpercz),nzfft)
     endif

     nxfft=n_of_fft(nx)
     nzfft=n_of_fft(nz)
     if(mod(nxfft,2).eq.1) nxfft=nxfft+1
     if(mod(nzfft,2).eq.1) nzfft=nzfft+1

     allocate(tapx_12(1:nxfft/2),tapz_12(1:nzfft/2))
     allocate(taper_x(1:nxfft),taper_z(1:nzfft))

     tapx_12=0.
     tapz_12=0.
     taper_x=0.
     taper_z=0.

     call staper(tapx_12,int(0.036*nxfft),nxfft/2)
     do i=1,nxfft/2
        taper_x(i)=tapx_12(i)
     enddo 
     do i=nxfft/2+1,nxfft
        taper_x(i)=tapx_12(i-nxfft/2)
     enddo 

     call staper(tapz_12,int(0.036*nzfft),nzfft/2)
     do i=1,nzfft/2
        taper_z(i)=tapz_12(i)
     enddo 
     do i=nzfft/2+1,nzfft
        taper_z(i)=tapz_12(i-nzfft/2)
     enddo

!     taper_x=1.0
!     taper_z=1.0

     if(flag.eq.1) then
       allocate(tmp_in_ud(1:nzfft))
       allocate(tmp_out_ud(1:nzfft))
       allocate(tmp_in_dw(1:nzfft))
       allocate(tmp_out_dw(1:nzfft))
       allocate(tmp_in_up(1:nzfft))
       allocate(tmp_out_up(1:nzfft))

       call sfftw_plan_dft_1d(plan_ud,nzfft,tmp_in_ud,tmp_out_ud,FFTW_FORWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_1d(plan_ud1,nzfft,tmp_in_up,tmp_out_up,FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_1d(plan_ud2,nzfft,tmp_in_dw,tmp_out_dw,FFTW_BACKWARD,FFTW_ESTIMATE)
     else if(flag.eq.2) then 
       allocate(tmp_in_lr(1:nxfft))
       allocate(tmp_out_lr(1:nxfft))
       allocate(tmp_in_lf(1:nxfft))
       allocate(tmp_out_lf(1:nxfft))
       allocate(tmp_in_rt(1:nxfft))
       allocate(tmp_out_rt(1:nxfft))

       call sfftw_plan_dft_1d(plan_lr,nxfft,tmp_in_lr,tmp_out_lr,FFTW_FORWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_1d(plan_lr1,nxfft,tmp_in_lf,tmp_out_lf,FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_1d(plan_lr2,nxfft,tmp_in_rt,tmp_out_rt,FFTW_BACKWARD,FFTW_ESTIMATE)
     else if(flag.eq.3) then 

       allocate(tmp_in_bin(1:nxfft,1:nzfft))
       allocate(tmp_out_bin(1:nxfft,1:nzfft))
       allocate(tmp_in_lu(1:nxfft,1:nzfft))
       allocate(tmp_out_lu(1:nxfft,1:nzfft))
       allocate(tmp_in_ru(1:nxfft,1:nzfft))
       allocate(tmp_out_ru(1:nxfft,1:nzfft))
       allocate(tmp_in_ld(1:nxfft,1:nzfft))
       allocate(tmp_out_ld(1:nxfft,1:nzfft))
       allocate(tmp_in_rd(1:nxfft,1:nzfft))
       allocate(tmp_out_rd(1:nxfft,1:nzfft))

       call sfftw_plan_dft_2d(plan_bin,nxfft,nzfft,tmp_in_bin,tmp_out_bin,&
                 fftw_forward,fftw_estimate)
       call sfftw_plan_dft_2d(plan_bin1,nxfft,nzfft,tmp_in_rd,tmp_out_rd,&
                 FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_2d(plan_bin2,nxfft,nzfft,tmp_in_ru,tmp_out_ru,&
                 FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_2d(plan_bin3,nxfft,nzfft,tmp_in_ld,tmp_out_ld,&
                 FFTW_BACKWARD,FFTW_ESTIMATE)
       call sfftw_plan_dft_2d(plan_bin4,nxfft,nzfft,tmp_in_lu,tmp_out_lu,&
                 FFTW_BACKWARD,FFTW_ESTIMATE)

      endif

     do i_shot=myid+1,num_shot,num_procs 

       modeling=.true.
!!     
       open(12,file=vpname,access='direct',form='unformatted',&
            recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           read(12,rec=(i_shot-1)*int_shot_model+i-nx_p1+1) &
                       (vp_f0(i,j),j=nz_p1,nz_p2)
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
           do i=nx_p1,nx_p2
              do l=1,nrel
                 read(12,rec=((i_shot-1)*int_shot_model+i-nx_p1+1)*nrel+l) &
                 (tau11_eps(l,i,j),j=nz_p1,nz_p2)
              enddo 
           enddo 
          close(12)   
    
          open(12,file=tau_sigmaname,access='direct',form='unformatted',&
               recl=ii_bi*nz_p)
           do i=nx_p1,nx_p2
              do l=1,nrel
                 read(12,rec=((i_shot-1)*int_shot_model+i-nx_p1+1)*nrel+l) &
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

       call para_prepare

!! calculating likelihood
       open(12,file='geom_x',access='direct',form='unformatted',&
               recl=ii_bi*(nx_rec+1)) 
         read(12,rec=i_shot) (geom_x(i),i=1,nx_rec+1)
       close(12)

       open(12,file='geom_z',access='direct',form='unformatted',&
               recl=ii_bi*(nx_rec+1))
         read(12,rec=i_shot) (geom_z(i),i=1,nx_rec+1)
       close(12)

!! shot position
       isource=int(geom_x(1))+npoints_pml
       jsource=int(geom_z(1))+npoints_pml

!! receiver position
       do i=1,nx_rec
          rcv_x(i)=int(geom_x(i+1))+npoints_pml
          rcv_z(i)=int(geom_z(i+1))+npoints_pml
       enddo

       vx=0.e0
       vz=0.e0
       p=0.e0
       r=0.e0
       memory_dvx_dx=0.e0
       memory_dvz_dz=0.e0
       memory_dp_dx=0.e0
       memory_dp_dz=0.e0

       vx_hb=0.e0
       vz_hb=0.e0
       p_hb=0.e0
       r_hb=0.e0
       memory_dvx_dx_hb=0.e0
       memory_dvz_dz_hb=0.e0
       memory_dp_dx_hb=0.e0
       memory_dp_dz_hb=0.e0

       call cpu_time (t1_a)
       t1_a_mp = omp_get_wtime ()

!      taper_x=1.0 
!      taper_z=1.0 
       do it = 1,nt

          if(it.lt.50000) then
             amp=ricker(it)
             p(isource,jsource) = p(isource,jsource) + amp

             amp_hb=ricker_hb(it)
             p_hb(isource,jsource) = p_hb(isource,jsource) + amp_hb
          endif

          if(mod(it,2000).eq.0) print *,'myid ',&
                          myid,'it',it,p(isource,jsource)


          call iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)

          call iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

!! for hilbert 
          call iso_visco_step_p(vx_hb,vz_hb,p_hb,r_hb,memory_dvx_dx_hb,memory_dvz_dz_hb)

          call iso_visco_step_vxz(vx_hb,vz_hb,p_hb,memory_dp_dx_hb,memory_dp_dz_hb)

!! wavefield decomposition
!! source wavefield decomposition
           call wavefield_decomp(p,p_hb,sou_p_up,sou_p_dw,sou_p_lf,sou_p_rt,&
                sou_p_lu,sou_p_ld,sou_p_ru,sou_p_rd,taper_x,taper_z,&
                nx,nz,lv,nxfft,nzfft,flag)

!! save pressure
          call bvr_record_wave_p(it,nt,vx,vz,p,p_l,p_r,p_t,p_b,&
               nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
               last_p(:,:)=p(:,:)
          endif

!!   save boundary values for wavefield reconstruction
          call bvr_record_wave_vxz(it,nt,vx,vz,p,vx_l,vx_r,vx_t,vx_b,vz_l,vz_r,vz_t,vz_b,&
                                   nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
              last_vx(:,:)=vx(:,:)
              last_vz(:,:)=vz(:,:)
          endif

!! save pressure
          call bvr_record_wave_p(it,nt,vx_hb,vz_hb,p_hb,p_l_hb,p_r_hb,p_t_hb,p_b_hb,&
               nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
               last_p_hb(:,:)=p_hb(:,:)
          endif

!!   save boundary values for wavefield reconstruction
          call bvr_record_wave_vxz(it,nt,vx_hb,vz_hb,p_hb,vx_l_hb,vx_r_hb,&
               vx_t_hb,vx_b_hb,vz_l_hb,vz_r_hb,vz_t_hb,vz_b_hb,&
               nx,nz,lv,npoints_pml)

          if(it.eq.nt) then
              last_vx_hb(:,:)=vx_hb(:,:)
              last_vz_hb(:,:)=vz_hb(:,:)
          endif

!         if(i_shot.eq.int(num_shot/2)) then 
          if(mod(it,nt_snap).eq.0) then 
            snap_part='snap_fwd'
            call snap_shot(p,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
            snap_part='snap_fwd_up'
            call snap_shot(sou_p_up,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
            snap_part='snap_fwd_dw'
            call snap_shot(sou_p_dw,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
          endif
!         endif
      enddo !! enddo it

    if(.true.) then 

!!     now it's rtm's turn
       modeling=.false.

!   if(.false.) then 
       deltat=-deltat
       c11_u=-c11_u
!      tau11=-tau11
!   endif

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

       image(:,:)=0.0
       image_illu(:,:)=0.0
       
       image_0_tmp=0.0
       image_01_tmp=0.0
       image_illu_0_tmp=0.0       
       
       image_ud=0.0
       image_du=0.0
       image_uu=0.0
       image_dd=0.0
       image_duud=0.0
       image_dduu=0.0

       image_ud_0_tmp=0.0
       image_du_0_tmp=0.0
       image_uu_0_tmp=0.0
       image_dd_0_tmp=0.0
       image_duud_0_tmp=0.0
       image_dduu_0_tmp=0.0

       image_ud_01_tmp=0.0
       image_du_01_tmp=0.0
       image_uu_01_tmp=0.0
       image_dd_01_tmp=0.0
       image_duud_01_tmp=0.0
       image_dduu_01_tmp=0.0

!       data_direc='/flush1/guo103/wave_decomp_csiro/bp_data/'
       file_p_rcv=trim(data_direc)//trim(file_p_part)//char(i_shot/100+48)//&
               char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)

!! read in observed data
       open(file_p_rcv_no,file=file_p_rcv,access='direct',form='unformatted'&
            ,recl=ii_bi*nx_rec)
        do j=1,nt
          read(file_p_rcv_no,rec=j) (seis_rcv(i,j),i=1,nx_rec)
        enddo
       close(file_p_rcv_no)

!! calculate its hilbert transform version
       call hilbert_seismo(seis_rcv,seis_rcv_hb,nx_rec,nt)

     do it=nt,1,-1
   
!!     receiver wavefield propagation
          call iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)

          call iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

!!     hilbert receiver wavefield propagation
          call iso_visco_step_p(vx_hb,vz_hb,p_hb,r_hb,memory_dvx_dx_hb,memory_dvz_dz_hb)

          call iso_visco_step_vxz(vx_hb,vz_hb,p_hb,memory_dp_dx_hb,memory_dp_dz_hb)

!!     source wavefield reconstruction
          call iso_visco_step_vxz(sou_vx,sou_vz,sou_p,sou_memory_dp_dx,sou_memory_dp_dz)

          call iso_visco_step_p(sou_vx,sou_vz,sou_p,sou_r,&
               sou_memory_dvx_dx,sou_memory_dvz_dz)

          call iso_visco_step_vxz(sou_vx_hb,sou_vz_hb,sou_p_hb,&
               sou_memory_dp_dx_hb,sou_memory_dp_dz_hb)

          call iso_visco_step_p(sou_vx_hb,sou_vz_hb,sou_p_hb,sou_r_hb,&
               sou_memory_dvx_dx_hb,sou_memory_dvz_dz_hb)

          if(mod(it,500).eq.0)  print *,'it',it,p(isource,jsource)
!            do i=nz_p1-4,nz_p2+4
            do i=1,nz
               do j=1,lv
!!      left
               sou_vx(npoints_pml+j,i)=vx_l(i,j,it)
               sou_vz(npoints_pml+j,i)=vz_l(i,j,it)
               sou_p(npoints_pml+j,i)=p_l(i,j,it)
!!      right
               sou_vx(nx-npoints_pml-lv+j,i)=vx_r(i,j,it)
               sou_vz(nx-npoints_pml-lv+j,i)=vz_r(i,j,it)
               sou_p(nx-npoints_pml-lv+j,i)=p_r(i,j,it)
               enddo
            enddo

!             do i=nx_p1,nx_p2
              do i=1,nx
                do j=1,lv
!!       top
                sou_vx(i,npoints_pml+j)=vx_t(i,j,it)
                sou_vz(i,npoints_pml+j)=vz_t(i,j,it)
                sou_p(i,npoints_pml+j)=p_t(i,j,it)
!!       bottom
                sou_vx(i,nz-npoints_pml-lv+j)=vx_b(i,j,it)
                sou_vz(i,nz-npoints_pml-lv+j)=vz_b(i,j,it)
                sou_p(i,nz-npoints_pml-lv+j)=p_b(i,j,it)
                enddo
             enddo
 
            if(it.eq.nt) then 
              sou_p(:,:)=last_p(:,:)
              sou_vx(:,:)=last_vx(:,:)
              sou_vz(:,:)=last_vz(:,:)
            endif

!! for hilbert
            do i=1,nz
               do j=1,lv
!!      left
               sou_vx_hb(npoints_pml+j,i)=vx_l_hb(i,j,it)
               sou_vz_hb(npoints_pml+j,i)=vz_l_hb(i,j,it)
               sou_p_hb(npoints_pml+j,i)=p_l_hb(i,j,it)
!!      right
               sou_vx_hb(nx-npoints_pml-lv+j,i)=vx_r_hb(i,j,it)
               sou_vz_hb(nx-npoints_pml-lv+j,i)=vz_r_hb(i,j,it)
               sou_p_hb(nx-npoints_pml-lv+j,i)=p_r_hb(i,j,it)
               enddo
            enddo

!             do i=nx_p1,nx_p2
              do i=1,nx
                do j=1,lv
!!       top
                sou_vx_hb(i,npoints_pml+j)=vx_t_hb(i,j,it)
                sou_vz_hb(i,npoints_pml+j)=vz_t_hb(i,j,it)
                sou_p_hb(i,npoints_pml+j)=p_t_hb(i,j,it)
!!       bottom
                sou_vx_hb(i,nz-npoints_pml-lv+j)=vx_b_hb(i,j,it)
                sou_vz_hb(i,nz-npoints_pml-lv+j)=vz_b_hb(i,j,it)
                sou_p_hb(i,nz-npoints_pml-lv+j)=p_b_hb(i,j,it)
                enddo
             enddo
 
            if(it.eq.nt) then 
              sou_p_hb(:,:)=last_p_hb(:,:)
              sou_vx_hb(:,:)=last_vx_hb(:,:)
              sou_vz_hb(:,:)=last_vz_hb(:,:)
            endif

!! wavefield decomposition
!! source wavefield decomposition
           call wavefield_decomp(sou_p,sou_p_hb,sou_p_up,sou_p_dw,sou_p_lf,sou_p_rt,&
                sou_p_lu,sou_p_ld,sou_p_ru,sou_p_rd,taper_x,taper_z,&
                nx,nz,lv,nxfft,nzfft,flag)

!! receiver wavefield decomposition 
           call wavefield_decomp(p,p_hb,rcv_p_up,rcv_p_dw,rcv_p_lf,rcv_p_rt,&
                rcv_p_lu,rcv_p_ld,rcv_p_ru,rcv_p_rd,taper_x,taper_z,&
                nx,nz,lv,nxfft,nzfft,flag)

!          call snap_shot(p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)    
!           if(i_shot.eq.int(num_shot/2)) then
            if(mod(it,nt_snap).eq.0) then 
              recons_part='snap_sou'
              call snap_shot(sou_p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_sou_up'
              call snap_shot(sou_p_up,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_sou_dw'
              call snap_shot(sou_p_dw,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_rcv'
              call snap_shot(p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_rcv_up'
              call snap_shot(rcv_p_up,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_rcv_dw'
              call snap_shot(rcv_p_dw,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
            endif
!           endif

            do j=nz_p1,nz_p2
               do i=nx_p1,nx_p2
                  image_illu(i,j)=image_illu(i,j)+sou_p(i,j)*sou_p(i,j)
               enddo
            enddo

            call imaging_condition(sou_p,p,image,nx,nz,lv,lv)

            call imaging_condition_udlr(sou_p_up,sou_p_dw,&
                 rcv_p_up,rcv_p_dw,image_du,image_ud,image_dd,image_uu,&
                 image_duud,image_dduu,nx,nz,lv) 
            
!           endif

!!     read in seismograms 
            do i=1,nx_rec
               p(rcv_x(i),rcv_z(i))=p(rcv_x(i),rcv_z(i))-seis_rcv(i,it)
               p_hb(rcv_x(i),rcv_z(i))=p_hb(rcv_x(i),rcv_z(i))-seis_rcv_hb(i,it)
            enddo

!! enddo time loop
       enddo

       call cpu_time (t2_a)
       t2_a_mp = omp_get_wtime ()

       t12_a=t2_a-t1_a
       t12_a_mp=t2_a_mp-t1_a_mp
 
       if(myid.eq.0) then
         open(81,file='cmpt_time_rtm.txt')
           write(81,*) t12_a
           write(81,*) t12_a_mp
         close(81)
       endif

       close(file_p_no)

       n_count=nx*nz
       call mpi_reduce(image,image_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_ud,image_ud_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_du,image_du_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_uu,image_uu_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dd,image_dd_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_duud,image_duud_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dduu,image_dduu_0_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

       call mpi_reduce(image_illu,image_illu_0,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

!!     cross-correlation 
       file_image=trim(file_image_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)
!!     source-normalized cross-correlation 
       file_image_normal=trim(file_image_normal_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)

       ii=0
       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_duud(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dduu(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_du(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_ud(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_uu(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dd(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       image(nx_p1:nx_p2,nz_p1:nz_p2)=image(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_duud(nx_p1:nx_p2,nz_p1:nz_p2)=image_duud(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_dduu(nx_p1:nx_p2,nz_p1:nz_p2)=image_dduu(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_du(nx_p1:nx_p2,nz_p1:nz_p2)=image_du(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_ud(nx_p1:nx_p2,nz_p1:nz_p2)=image_ud(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_uu(nx_p1:nx_p2,nz_p1:nz_p2)=image_uu(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_dd(nx_p1:nx_p2,nz_p1:nz_p2)=image_dd(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       call mpi_reduce(image,image_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_ud,image_ud_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_du,image_du_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_uu,image_uu_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dd,image_dd_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_duud,image_duud_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dduu,image_dduu_01_tmp,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

       ii=0
       open(file_image_no,file=file_image_normal,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_duud(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dduu(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_du(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_ud(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_uu(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dd(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

    endif !! endif false

       print *
       print *,'End of the simulation'
       print *

     if(myid.eq.0) then 
       image_0=image_0+image_0_tmp
       image_ud_0=image_ud_0+image_ud_0_tmp
       image_du_0=image_du_0+image_du_0_tmp
       image_uu_0=image_uu_0+image_uu_0_tmp
       image_dd_0=image_dd_0+image_dd_0_tmp
       image_duud_0=image_duud_0+image_duud_0_tmp
       image_dduu_0=image_dduu_0+image_dduu_0_tmp
       image_illu_0=image_illu_0+image_illu_0_tmp
       
       image_01=image_01+image_01_tmp
       image_ud_01=image_ud_01+image_ud_01_tmp
       image_du_01=image_du_01+image_du_01_tmp
       image_uu_01=image_uu_01+image_uu_01_tmp
       image_dd_01=image_dd_01+image_dd_01_tmp
       image_duud_01=image_duud_01+image_duud_01_tmp
       image_dduu_01=image_dduu_01+image_dduu_01_tmp
     endif

!     call mpi_barrier(mpi_comm_world,ierr)

    enddo   !! shot domain

     if(myid.eq.0) then 
       file_image=trim(file_image_part)//'_sum'

       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        ii=0
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_duud_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dduu_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_du_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_ud_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_uu_0(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dd_0(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       file_image=trim(file_image_part)//'_hes'//'_sum'

       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        ii=0
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_duud_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dduu_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_du_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_ud_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_uu_01(i,j),j=nz_p1,nz_p2)
        enddo 
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_dd_01(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

       file_image=trim(file_image_part)//'_illu'//'_sum'
       open(file_image_no,file=file_image,access='direct',&
             form='unformatted',recl=ii_bi*nz_p)
        ii=0
        do i=nx_p1,nx_p2
           ii=ii+1
           write(file_image_no,rec=ii) (image_illu_0(i,j),j=nz_p1,nz_p2)
        enddo 
       close(file_image_no) 

      endif  !!endif myid.eq.0

      if(flag.eq.1) then 
          call sfftw_destroy_plan(plan_ud)
          call sfftw_destroy_plan(plan_ud1)
          call sfftw_destroy_plan(plan_ud2)
       else if (flag.eq.2) then 
          call sfftw_destroy_plan(plan_lr)
          call sfftw_destroy_plan(plan_lr1)
          call sfftw_destroy_plan(plan_lr2)       
       else if (flag.eq.3) then 
          call sfftw_destroy_plan(plan_bin)
          call sfftw_destroy_plan(plan_bin1)
          call sfftw_destroy_plan(plan_bin2)
          call sfftw_destroy_plan(plan_bin3)
          call sfftw_destroy_plan(plan_bin4) 
      endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

  end program seismic_CPML_2D_iso_rtm

