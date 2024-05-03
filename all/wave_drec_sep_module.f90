
      module wave_drec_sep_module 
       use precision_m     

       real (fp_kind),allocatable::sou_p_up(:,:),sou_p_dw(:,:),&
                                   sou_p_lf(:,:),sou_p_rt(:,:)

       real (fp_kind),allocatable::sou_p_lu(:,:),sou_p_ld(:,:),&
                                   sou_p_ru(:,:),sou_p_rd(:,:)

       real (fp_kind),allocatable::rcv_p_up(:,:),rcv_p_dw(:,:),&
                                   rcv_p_lf(:,:),rcv_p_rt(:,:)

       real (fp_kind),allocatable::rcv_p_lu(:,:),rcv_p_ld(:,:),&
                                   rcv_p_ru(:,:),rcv_p_rd(:,:)
     end module wave_drec_sep_module
