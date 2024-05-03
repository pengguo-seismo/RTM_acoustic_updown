
       real,allocatable::n_lay(:,:),v(:)
       integer::rns1,rns2,n_chain,nz
       character*256::file_v,file_dir,file_out,file_dir_part
 
       rns1=17001
       rns2=20000
       n_chain=30
       nz=181
       ii_bi=1

       allocate(n_lay(rns1:rns2,n_chain))
       n_lay=0
       allocate(v(1:nz))
       v=0.0

       file_dir_part='./'

       do i_chain=1,n_chain
 
          file_dir='modelaa_'//char(i_chain/100+48)//&
               char(mod(i_chain,100)/10+48)//char(mod(i_chain,10)+48)
          file_v=trim(file_dir)//'/vel_his'

          open(12,file=file_v,access='direct',form='unformatted',&
                  recl=ii_bi*nz)
           do i_his=rns1,rns2
              read(12,rec=i_his) (v(i),i=1,nz)
      
              num_lay=1
              do i=1,nz-1
                 if(abs(v(i+1)-v(i)).gt.1.0E-3) num_lay=num_lay+1
              enddo 
              n_lay(i_his,i_chain)=num_lay
           enddo           
          close(12)

        file_out='results/'//'his_num'//char(i_chain/100+48)//&
               char(mod(i_chain,100)/10+48)//char(mod(i_chain,10)+48)
      
         open(12,file=file_out)
          do i=rns1,rns2
             write(12,*) i,n_lay(i,i_chain)
          enddo 
         close(12)
        
        enddo !! i_chain 

        file_out='results/'//'his_num'
        open(12,file=file_out)
        do i_chain=1,n_chain
          do i=rns1,rns2
             write(12,*) i,i_chain,n_lay(i,i_chain)
          enddo
        enddo 
        close(12)

       end





