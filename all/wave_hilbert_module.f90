
      module wave_hilbert_module 
       use precision_m     

       real (fp_kind),allocatable::sou_vx_hb(:,:),sou_vz_hb(:,:),sou_p_hb(:,:)
       real (fp_kind), allocatable::sou_r_hb(:,:,:)

       real (fp_kind), allocatable::sou_memory_dvx_dx_hb(:,:),sou_memory_dvz_dz_hb(:,:),&
                                    sou_memory_dp_dx_hb(:,:),sou_memory_dp_dz_hb(:,:)

       real (fp_kind),allocatable::vx_hb(:,:),vz_hb(:,:),p_hb(:,:)
       real (fp_kind), allocatable::r_hb(:,:,:)

       real (fp_kind), allocatable::memory_dvx_dx_hb(:,:),memory_dvz_dz_hb(:,:),&
                                    memory_dp_dx_hb(:,:),memory_dp_dz_hb(:,:)  
 
       real(fp_kind),allocatable::seis_sou_hb(:,:),seis_rcv_hb(:,:) 

       real(fp_kind),allocatable::vx_l_hb(:,:,:),vx_r_hb(:,:,:),vx_t_hb(:,:,:),vx_b_hb(:,:,:)
       real(fp_kind),allocatable::vz_l_hb(:,:,:),vz_r_hb(:,:,:),vz_t_hb(:,:,:),vz_b_hb(:,:,:)
       real(fp_kind),allocatable::p_l_hb(:,:,:),p_r_hb(:,:,:),p_t_hb(:,:,:),p_b_hb(:,:,:)
       real(fp_kind),allocatable::last_p_hb(:,:),last_vx_hb(:,:),last_vz_hb(:,:)

     end module wave_hilbert_module
