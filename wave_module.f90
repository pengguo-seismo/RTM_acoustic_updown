
      module wave_module 
       use precision_m     

! main arrays
       real (fp_kind),allocatable::vx(:,:),vz(:,:),p(:,:)
       real (fp_kind), allocatable::r(:,:,:)

! main arrays for RTM
       real (fp_kind),allocatable::sou_vx(:,:),sou_vz(:,:),sou_p(:,:)
       real (fp_kind), allocatable::sou_r(:,:,:)

       real (fp_kind), allocatable::memory_dvx_dx(:,:),memory_dvz_dz(:,:),&
                                    memory_dp_dx(:,:),memory_dp_dz(:,:)

!
       real (fp_kind), allocatable::sou_memory_dvx_dx(:,:),sou_memory_dvz_dz(:,:),&
                                    sou_memory_dp_dx(:,:),sou_memory_dp_dz(:,:)
   
       real(fp_kind),allocatable::image(:,:),image_0(:,:),image_01(:,:),&
                                  image_illu(:,:),image_illu_0(:,:)

       real(fp_kind),allocatable::image_du(:,:),image_ud(:,:),&
                                  image_dd(:,:),image_uu(:,:),&
                                  image_duud(:,:),image_dduu(:,:)

       real(fp_kind),allocatable::image_du_0(:,:),image_ud_0(:,:),&
                                  image_dd_0(:,:),image_uu_0(:,:),&
                                  image_duud_0(:,:),image_dduu_0(:,:)

       real(fp_kind),allocatable::image_du_01(:,:),image_ud_01(:,:),&
                                  image_dd_01(:,:),image_uu_01(:,:),&
                                  image_duud_01(:,:),image_dduu_01(:,:)

        real(fp_kind),allocatable::image_0_tmp(:,:),image_01_tmp(:,:),&
                                   image_illu_0_tmp(:,:)

       real(fp_kind),allocatable::image_du_0_tmp(:,:),image_ud_0_tmp(:,:),&
                                  image_dd_0_tmp(:,:),image_uu_0_tmp(:,:),&
                                  image_duud_0_tmp(:,:),image_dduu_0_tmp(:,:)

       real(fp_kind),allocatable::image_du_01_tmp(:,:),image_ud_01_tmp(:,:),&
                                  image_dd_01_tmp(:,:),image_uu_01_tmp(:,:),&
                                  image_duud_01_tmp(:,:),image_dduu_01_tmp(:,:)

       real(fp_kind),allocatable::image_hb(:,:),image_duud_hb(:,:),image_dduu_hb(:,:)

       real(fp_kind),allocatable::seis_sou(:,:),seis_rcv(:,:) 

       real(fp_kind),allocatable::epsilon_m(:,:)

       real(fp_kind),allocatable::vx_l(:,:,:),vx_r(:,:,:),vx_t(:,:,:),vx_b(:,:,:)
       real(fp_kind),allocatable::vz_l(:,:,:),vz_r(:,:,:),vz_t(:,:,:),vz_b(:,:,:)
       real(fp_kind),allocatable::p_l(:,:,:),p_r(:,:,:),p_t(:,:,:),p_b(:,:,:)
       real(fp_kind),allocatable::r_l(:,:,:,:),r_r(:,:,:,:),r_t(:,:,:,:),r_b(:,:,:,:)
       real(fp_kind),allocatable::last_p(:,:),last_vx(:,:),last_vz(:,:)

       real(fp_kind),allocatable::tmp_x(:,:),tmp_z(:,:)

!      contains
 
!        subroutine allocate_wave_forward(nx,nz,lv)

     end module wave_module
