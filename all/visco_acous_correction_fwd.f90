!
! 2D isotropic visco-acoustic modeling program based on the GSLS model
! with elastic effect correction
! Peng Guo 
! CSIRO


   program seismic_CPML_2D_iso

!            ^ z
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!       v_z  +                   |
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            +---------+---------+  ---> x
!            P        v_x
!         
!
!     use mpi
     use precision_m
     use main_module
     use parameter_input
     use wave_propagator
     use wave_module
     use omp_lib
  
     implicit none
!     include 'variable.h'
     include 'mpif.h'

     double precision::t1,t2,t12 
     double precision::t_mp1,t_mp2,t_mp12 
     integer::num_threads

     integer::i_xs,i_zs

     real,allocatable::ux0(:,:),uz0(:,:)
     real,allocatable::vs_f0(:,:),mu_f0(:,:)  
     real,allocatable::vx_els(:,:),vz_els(:,:),&
                       res_source_x(:,:),res_source_z(:,:)
     character*256::vsname

!, omp_get_num_threads
!     double precision::omp_get_wtime
!---
!--- program starts here
!---
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

      ii_bi=1

!!   read in parameter, to decide the parmeter file name in subroutine par_in
!!   0: forward modeling, parameter file: parameter.txt
!!   1: forward modeling for direct waves, parameter file: parameter_dir.txt
!!   2: rtm, parameter file: parameter_rtm.txt
      call getarg(1,arg)
      arg=trim(arg)
      read(arg,'(i4)') iflag_par

      call par_in(iflag_par)

!!    allocate memory
      call allocate_main_module(nx,nz,nrel,nx_rec,lv)

      allocate(p(-lv:nx+lv,-lv:nz+lv))
      allocate(vx(-lv:nx+lv,-lv:nz+lv))
      allocate(vz(-lv:nx+lv,-lv:nz+lv))
      allocate(r(1:nrel,-lv:nx+lv,-lv:nz+lv))

      allocate(memory_dvx_dx(1:nx,1:nz))
      allocate(memory_dvz_dz(1:nx,1:nz))
      allocate(memory_dp_dx(1:nx,1:nz))
      allocate(memory_dp_dz(1:nx,1:nz))

!! correction modeling 
      allocate(sou_p(-lv:nx+lv,-lv:nz+lv))
      allocate(sou_vx(-lv:nx+lv,-lv:nz+lv))
      allocate(sou_vz(-lv:nx+lv,-lv:nz+lv))
      allocate(sou_r(1:nrel,-lv:nx+lv,-lv:nz+lv))

      allocate(sou_memory_dvx_dx(1:nx,1:nz))
      allocate(sou_memory_dvz_dz(1:nx,1:nz))
      allocate(sou_memory_dp_dx(1:nx,1:nz))
      allocate(sou_memory_dp_dz(1:nx,1:nz))

      allocate(vx_els(-lv:nx+lv,-lv:nz+lv),&
               vz_els(-lv:nx+lv,-lv:nz+lv))

      allocate(ux0(-lv:nx+lv,-lv:nz+lv),uz0(-lv:nx+lv,-lv:nz+lv))

      allocate(vs_f0(1:nx,1:nz))
      allocate(mu_f0(1:nx,1:nz))

      allocate(res_source_x(-lv:nx+lv,-lv:nz+lv),&
               res_source_z(-lv:nx+lv,-lv:nz+lv))

     do i_shot=myid+1,num_shot,num_procs

      file_p_no=21
      file_snap_no=22

      print *,i_shot,'i_shot',myid,'myid',num_procs,'num_procs'
      print *,'nx',nx,'nz',nz,'deltax',deltax,'deltaz',deltaz

!! read in geometry
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

!      isource=300
!      jsource=30

      print *,'nz_p',nz_p,'nz_p1',nz_p1,'nz_p2',nz_p2,&
                          'nx_p1',nx_p1,'nx_p2',nx_p2

      open(12,file=vpname,access='direct',form='unformatted',&
           recl=ii_bi*nz_p)
       do i=nx_p1,nx_p2
         read(12,rec=i-nx_p1+1) &
                     (vp_f0(i,j),j=nz_p1,nz_p2)
       enddo 
      close(12)        

      open(12,file=rhoname,access='direct',form='unformatted',&
           recl=ii_bi*nz_p)
       do i=nx_p1,nx_p2
         read(12,rec=i-nx_p1+1) &
                     (dens(i,j),j=nz_p1,nz_p2)
       enddo 
      close(12)        

!!!    if attenuation 
      if(.not.lossless) then
         open(12,file=tau_epsname,access='direct',form='unformatted',&
              recl=ii_bi*nz_p)
          do i=nx_p1,nx_p2
             do l=1,nrel
               read(12,rec=(i-nx_p1)*nrel+l) &
                   (tau11_eps(l,i,j),j=nz_p1,nz_p2)
             enddo 
          enddo 
         close(12)   

         open(12,file=tau_sigmaname,access='direct',form='unformatted',&
              recl=ii_bi*nz_p)
           do i=nx_p1,nx_p2
              do l=1,nrel
               read(12,rec=(i-nx_p1)*nrel+l) &
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

!!   read in S wave velocity, for the purpose of elastic correction !!
 
     vs_f0(:,:)=vp_f0(:,:)/1.7
     if(.false.) then 
      open(12,file=vsname,access='direct',form='unformatted',&
           recl=ii_bi*nz_p)
       do i=nx_p1,nx_p2
         read(12,rec=i-nx_p1+1) (vs_f0(i,j),j=nz_p1,nz_p2)
       enddo
      close(12)
     endif

       do i=nx_p1,nx_p2
!! top
           do j=1,nz_p1
              vs_f0(i,j)=vs_f0(i,nz_p1)
           enddo
!! bottom
           do j=nz_p2,nz
              vs_f0(i,j)=vs_f0(i,nz_p2)
           enddo
       enddo

       do j=1,nz
          do i=1,nx_p1
!! left 
              vs_f0(i,j)=vs_f0(nx_p1,j)
          enddo
          do i=nx_p2,nx
!! right
             vs_f0(i,j)=vs_f0(nx_p2,j)
          enddo
       enddo

       open(12,file='vs_test.dat',access='direct',form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
          write(12,rec=i) (vs_f0(i,j),j=1,nz)
       enddo 
       close(12) 


!!     assign values for pml layers
!!     1:npoints_pml,nx-npoints_pml+1:nx
       call adding_abs_val(vp_f0,dens,tau11_eps,tau_sigma,nx,nz,&
            nrel,npoints_pml)

       mu_f0(:,:)=dens(:,:)*vs_f0(:,:)*vs_f0(:,:)

       call vp_f0_2_vp(nx,nz,lv,nrel,f_ref,c11_r,vp_f0,vp,dens,&
                       tau11_eps,tau_sigma)

       call assign_parameter(nx,nz,lv,nrel,deltat,& 
                             tau11_eps,tau11,tau11_new,tau_sigma,x1,x2,&
                             c11_r,dens,c11_u)
! call ricker wavelet
       call rickerfunc(deltat,f0,length,nt,ricker)
 
       open(12,file='mu_test.dat',access='direct',form='unformatted',&
           recl=ii_bi*nz)
       do i=1,nx
          write(12,rec=i) (mu_f0(i,j),j=1,nz)
       enddo 
       close(12) 
! display size of the model
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

! compute cpml coefficients !
       call  cpml_coef(a_x,b_x,k_x,a_z,b_z,k_z,&
             a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
             npoints_pml,nx,nz,deltax,deltaz,deltat,vp_f0,f0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---
       c(1)=1.23538628085693269
       c(2)=-0.109160358903484037
       c(3)=0.246384765606012697d-1
       c(4)=-0.660472512788868975d-2

!!!!!  FD coefficients, holberg version
       c(1)=1.231666
       c(2)=-0.1041182
       c(3)=0.02063707
       c(4)=-0.003570998

! from DENISE
!       c(:)=0.0
!       c(1)=1.1382
!       c(2)=-0.046414
       
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
       enddo 

!$OMP PARALLEL 
       num_threads = omp_get_num_threads()
       print *,'num_threads',num_threads
!$OMP END PARALLEL

       call cpu_time (t1)
       t_mp1=omp_get_wtime ()

       file_p=trim(file_p_part)//char(i_shot/100+48)//&
               char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)

       open(file_p_no,file=file_p,access='direct',form='unformatted',&
            recl=ii_bi*nx_rec)

       vx=0.0
       vz=0.0
       p=0.0
       r=0.0
       memory_dvx_dx=0.e0
       memory_dvz_dz=0.e0
       memory_dp_dx=0.e0
       memory_dp_dz=0.e0     
 
       sou_vx=0.0
       sou_vz=0.0
       sou_p=0.0
       sou_r=0.0
       sou_memory_dvx_dx=0.e0
       sou_memory_dvz_dz=0.e0
       sou_memory_dp_dx=0.e0
       sou_memory_dp_dz=0.e0     
 
       ux0=0.
       uz0=0.

       i_xs=isource
       i_zs=jsource

       do it = 1,nt

         if(it.lt.50000) then 
            amp=ricker(it)
!            p(isource,jsource) = p(isource,jsource) + amp

            vz(i_xs,i_zs)=vz(i_xs,i_zs)-dlv(1)*amp
            vz(i_xs,i_zs+1)=vz(i_xs,i_zs+1)+dlv(1)*amp
            vx(i_xs,i_zs+1)=vx(i_xs,i_zs+1)+dlh(1)*amp
            vx(i_xs-1,i_zs+1)=vx(i_xs-1,i_zs+1)-dlh(1)*amp

            vz(i_xs,i_zs-1)=vz(i_xs,i_zs-1)-dlv(2)*amp
            vz(i_xs,i_zs+2)=vz(i_xs,i_zs+2)+dlv(2)*amp
            vx(i_xs+1,i_zs+1)=vx(i_xs+1,i_zs+1)+dlh(2)*amp
            vx(i_xs-2,i_zs+1)=vx(i_xs-2,i_zs+1)-dlh(2)*amp
            vz(i_xs,i_zs-2)=vz(i_xs,i_zs-2)-dlv(3)*amp
            vz(i_xs,i_zs+3)=vz(i_xs,i_zs+3)+dlv(3)*amp
            vx(i_xs+2,i_zs+1)=vx(i_xs+2,i_zs+1)+dlh(3)*amp
            vx(i_xs-3,i_zs+1)=vx(i_xs-3,i_zs+1)-dlh(3)*amp
            vz(i_xs,i_zs-3)=vz(i_xs,i_zs-3)-dlv(4)*amp
            vz(i_xs,i_zs+4)=vz(i_xs,i_zs+4)+dlv(4)*amp
            vx(i_xs+3,i_zs+1)=vx(i_xs+3,i_zs+1)+dlh(4)*amp
            vx(i_xs-4,i_zs+1)=vx(i_xs-4,i_zs+1)-dlh(4)*amp

          if(.false.) then
            vx(isource,jsource) = vx(isource,jsource) + 0.25*amp
            vx(isource-1,jsource) = vx(isource-1,jsource) - 0.25*amp
            vz(isource,jsource) = vz(isource,jsource) + 0.25*amp
            vz(isource,jsource-1) = vz(isource,jsource-1) - 0.25*amp
          endif
         endif

        if(mod(it,int(nt/10)).eq.0) then 
          print *,p(isource,jsource),vx(isource,jsource),&
                 vz(isource,jsource),amp, dens(isource,jsource),&
                 c11_u(isource,jsource),&
                'haha',amp,isource,jsource,'sigmaxx',it
         endif
!!------------------------------------------------------------
!! compute p, vx, vz and update memory variables for C-PML
!!------------------------------------------------------------
!! kernel subroutine  
           call iso_visco_step_vxz(vx,vz,p,memory_dp_dx,memory_dp_dz)

           call iso_visco_step_p(vx,vz,p,r,memory_dvx_dx,memory_dvz_dz)
           
           do j=1,nz
              do i=1,nx
                ux0(i,j)=ux0(i,j)+vx(i,j)*deltat
                uz0(i,j)=uz0(i,j)+vz(i,j)*deltat
              enddo
           enddo

!           ux0(1:nx,1:nz)=vx(1:nx,1:nz)
!           uz0(1:nx,1:nz)=vz(1:nx,1:nz)

           call res_source2(res_source_x,res_source_z,&
                   ux0,uz0,mu_f0,nx,nz,nrel,lv,deltax,deltaz)

!           call res_source(res_source_x,res_source_z,&
!                ux0,uz0,mu_f0,nx,nz,nrel,lv)

!! correction modeling
!! use it as source
          do j=40,nz 
            do i=1,nx
               sou_vx(i,j)=sou_vx(i,j)+res_source_x(i,j)/dens(i,j)*deltat
               sou_vz(i,j)=sou_vz(i,j)+res_source_z(i,j)/dens(i,j)*deltat
            enddo 
          enddo

           call iso_visco_step_vxz(sou_vx,sou_vz,sou_p,sou_memory_dp_dx,sou_memory_dp_dz) 

           call iso_visco_step_p(sou_vx,sou_vz,sou_p,sou_r,sou_memory_dvx_dx,sou_memory_dvz_dz)

           vx_els(1:nx,1:nz)=vx(1:nx,1:nz)+sou_vx(1:nx,1:nz)
           vz_els(1:nx,1:nz)=vz(1:nx,1:nz)+sou_vz(1:nx,1:nz)

           write(file_p_no,rec=it) (p(rcv_x(i),rcv_z(i)),i=1,nx_rec)

!          if(i_shot.eq.int(num_shot/2)) then
             if(mod(it,nt_snap).eq.0) then  
               snap_part='snap_ux'
               call snap_shot(ux0,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_uz'
               call snap_shot(uz0,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vx'
               call snap_shot(vx,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vz'
               call snap_shot(vz,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vx_els'
               call snap_shot(vx_els,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vz_els'
               call snap_shot(vz_els,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vx_sou'
               call snap_shot(sou_vx,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_vz_sou'
               call snap_shot(sou_vz,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_resx_sou'
               call snap_shot(res_source_x,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
               snap_part='snap_resz_sou'
               call snap_shot(res_source_z,it,nx,nz,lv,nt_snap,snap_part,file_snap_no)
             endif
!          endif

        enddo   ! end of time loop

       close(file_p_no)

       call cpu_time (t2)
       t_mp2=omp_get_wtime ()

       t12=t2-t1 

       t_mp12=t_mp2-t_mp1 

       print *,'t12, the computation time',t12
       print *,'t_mp12, the computation time',t_mp12

      enddo !! shot domain

      open(12,file='cmp_time.txt')
      write(12,*) t12
      write(12,*) t_mp12
      close(12)

     call deallocate_main_module()

!  print *
!  print *,'End of the simulation'
!  print *

      call MPI_FINALIZE(ierr)

      end
