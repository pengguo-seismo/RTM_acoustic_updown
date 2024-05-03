
          subroutine waveform_obj(syn_p,obs_p,res_p,nx_rec,nt)
           implicit none
           integer::nx_rec,nt
           real,dimension(1:nx_rec,1:nt)::syn_p,obs_p,res_p

           res_p(:,:)=syn_p(:,:)-obs_p(:,:)

          end subroutine waveform_obj

          subroutine envelope_obj(syn_p,obs_p,res_p,nx_rec,nt)
           implicit none
           integer::nx_rec,nt
           real,dimension(1:nx_rec,1:nt)::syn_p,obs_p,res_p
           real,dimension(1:nx_rec,1:nt)::syn_p_hil,obs_p_hil
           real,dimension(1:nx_rec,1:nt)::syn_p_env,obs_p_env

           call hilbert_seismo1(syn_p,syn_p_hil,nx_rec,nt)
           call hilbert_seismo1(obs_p,obs_p_hil,nx_rec,nt)
          
           syn_p_env=sqrt(syn_p(:,:)*syn_p(:,:)+syn_p_hil(:,:)*syn_p_hil(:,:)) 
           obs_p_env=sqrt(obs_p(:,:)*obs_p(:,:)+obs_p_hil(:,:)*obs_p_hil(:,:)) 

           res_p(:,:)=syn_p_env(:,:)-obs_p_env(:,:)

          end subroutine envelope_obj

          subroutine waveform_enve_obj(syn_p,obs_p,res_p,nx_rec,nt)

          end subroutine waveform_enve_obj
