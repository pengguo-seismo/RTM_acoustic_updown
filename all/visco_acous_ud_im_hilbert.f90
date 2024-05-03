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

    use precision_m
    use main_module
    use parameter_input
    use wave_propagator    
    use wave_module
    use wave_drec_sep_module
    use wave_hilbert_module
    use fftw_module_f

!!---implementation of Kai's 2017 SEG expanded abstract

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

       allocate(image_duud(1:nx,1:nz),image_dduu(1:nx,1:nz))
       image_duud=0.0
       image_dduu=0.0

!! summed results
       allocate(image_duud_0(1:nx,1:nz),image_dduu_0(1:nx,1:nz))

       image_duud_0=0.0
       image_dduu_0=0.0

       allocate(image_duud_0_tmp(1:nx,1:nz),image_dduu_0_tmp(1:nx,1:nz))

       image_duud_0_tmp=0.0
       image_dduu_0_tmp=0.0
!! summed results (with approximate hessian conditioned)
       allocate(image_duud_01(1:nx,1:nz),image_dduu_01(1:nx,1:nz))
       image_duud_01=0.0
       image_dduu_01=0.0

       allocate(image_duud_01_tmp(1:nx,1:nz),image_dduu_01_tmp(1:nx,1:nz))
       image_duud_01_tmp=0.0
       image_dduu_01_tmp=0.0

       allocate(seis_sou(1:nx_rec,1:nt),seis_rcv(1:nx_rec,1:nt))
       
       seis_sou(:,:)=0.e0
       seis_rcv(:,:)=0.e0

       int_shot_model=0

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

!! fft
    if(flag.eq.1) then
     allocate(tmp_in_ud(1:nzfft))
     allocate(tmp_out_ud(1:nzfft))
     tmp_in_ud=0.0
     tmp_out_ud=0.0
     call sfftw_plan_dft_1d(plan_ud,nzfft,tmp_in_ud,tmp_out_ud,FFTW_FORWARD,FFTW_ESTIMATE)
     call sfftw_plan_dft_1d(plan_ud_inv,nzfft,tmp_out_ud,tmp_in_ud,FFTW_BACKWARD,FFTW_ESTIMATE)

     else if(flag.eq.2) then 
     allocate(tmp_in_lr(1:nzfft))
     allocate(tmp_out_lr(1:nzfft))
     tmp_in_lr=0.0
     tmp_out_lr=0.0
     call sfftw_plan_dft_1d(plan_lr,nxfft,tmp_in_lr,tmp_out_lr,FFTW_FORWARD,FFTW_ESTIMATE)
     call sfftw_plan_dft_1d(plan_lr_inv,nxfft,tmp_out_lr,tmp_in_lr,FFTW_BACKWARD,FFTW_ESTIMATE)
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

! call ricker wavelet
       call rickerfunc(deltat,f0,length,nt,ricker)

       call rickerfunc_w_hb(deltat,f0,length,nt,ricker,ricker_hb)

!      taper_x=1.0 
!      taper_z=1.0 
       do it = 1,nt

          if(it.lt.50000) then
             amp=ricker(it)
             p(isource,jsource) = p(isource,jsource) + amp
          endif

          if(mod(it,2000).eq.0) print *,'myid ',&
                          myid,'it',it,p(isource,jsource)

          call iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)

          call iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

!! wavefield decomposition
!! source wavefield decomposition

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

          if(mod(it,nt_snap).eq.0) then 
            snap_part='snap_fwd'
            call snap_shot(p,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
          endif

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
  
       sou_vx(:,:)=zero
       sou_vz(:,:)=zero
       sou_p(:,:)=zero
       sou_r(:,:,:)=zero
       sou_memory_dvx_dx=zero
       sou_memory_dvz_dz=zero
       sou_memory_dp_dx=zero
       sou_memory_dp_dz=zero

       image_0_tmp=0.0
       image_01_tmp=0.0
       image_illu_0_tmp=0.0

       image_duud_0_tmp=0.0
       image_dduu_0_tmp=0.0

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


     do it=nt,1,-1
   
!!     receiver wavefield propagation
          call iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)

          call iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

!!     source wavefield reconstruction
          call iso_visco_step_vxz(sou_vx,sou_vz,sou_p,sou_memory_dp_dx,sou_memory_dp_dz)

          call iso_visco_step_p(sou_vx,sou_vz,sou_p,sou_r,&
               sou_memory_dvx_dx,sou_memory_dvz_dz)

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

!          call snap_shot(p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)    
            if(mod(it,nt_snap).eq.0) then 
              recons_part='snap_recons'
              call snap_shot(sou_p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
              recons_part='snap_rtm'
              call snap_shot(p,it,nx,nz,lv,nt_snap,recons_part,file_snap_no)     
            endif

            do j=nz_p1,nz_p2
               do i=nx_p1,nx_p2
                  image_illu(i,j)=image_illu(i,j)+sou_p(i,j)*sou_p(i,j)
               enddo
            enddo

            call imaging_condition(sou_p,p,image,nx,nz,lv,lv)
           
            call imaging_condition_udlr_im(sou_p,p,taper_x,taper_z,&
                 image_duud,image_dduu,nx,nz,lv,nxfft,nzfft,flag,it) 
            
!           endif

!!     read in seismograms 
            do i=1,nx_rec
               p(rcv_x(i),rcv_z(i))=p(rcv_x(i),rcv_z(i))-seis_rcv(i,it)
            enddo

!! enddo time loop
       enddo
 
       close(file_p_no)

       n_count=nx*nz
       call mpi_reduce(image,image_0,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_duud,image_duud_0,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dduu,image_dduu_0,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

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
       close(file_image_no) 

       image(nx_p1:nx_p2,nz_p1:nz_p2)=image(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_duud(nx_p1:nx_p2,nz_p1:nz_p2)=image_duud(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       image_dduu(nx_p1:nx_p2,nz_p1:nz_p2)=image_dduu(nx_p1:nx_p2,nz_p1:nz_p2)&
                                     /image_illu(nx_p1:nx_p2,nz_p1:nz_p2)

       call mpi_reduce(image,image_01,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_duud,image_duud_01,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       call mpi_reduce(image_dduu,image_dduu_01,n_count,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

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
       close(file_image_no) 

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
       close(file_image_no) 

      endif  !!endif myid.eq.0

   endif !! endif false

       print *
       print *,'End of the simulation'
       print *

    enddo   !! shot domain

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

  end program seismic_CPML_2D_iso_rtm

