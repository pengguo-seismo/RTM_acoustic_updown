
      module fftw_module_f 

!! 1d 
       integer*8::plan_ud,plan_ud1,plan_ud2,plan_ud_inv
       integer*8::plan_lr,plan_lr1,plan_lr2,plan_lr_inv

       complex,allocatable::tmp_in_ud(:),tmp_out_ud(:)
       complex,allocatable::tmp_in_lr(:),tmp_out_lr(:)

       complex,allocatable::tmp_in_up(:),tmp_out_up(:),&
                            tmp_in_dw(:),tmp_out_dw(:)

       complex,allocatable::tmp_in_lf(:),tmp_out_lf(:),&
                            tmp_in_rt(:),tmp_out_rt(:)

!! 2d
       integer*8::plan_bin,plan_bin1,plan_bin2,plan_bin3,plan_bin4

       complex,allocatable::tmp_in_bin(:,:),tmp_out_bin(:,:)
       complex,allocatable::tmp_in_rd(:,:),tmp_out_rd(:,:)
       complex,allocatable::tmp_in_ru(:,:),tmp_out_ru(:,:)
       complex,allocatable::tmp_in_ld(:,:),tmp_out_ld(:,:)
       complex,allocatable::tmp_in_lu(:,:),tmp_out_lu(:,:)

      end module fftw_module_f
