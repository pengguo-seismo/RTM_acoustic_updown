
       implicit none
       real,allocatable::image(:,:),image_shot(:,:)
       character *30::file_image,file_image_part,file_image_stack
       integer::i_shot,ii_bi
       integer::mx,nx,mz,n_shot
       integer::i,j,int_shot_model

       open(12,file='parameter_stack.txt')
        read(12,*) mx
        read(12,*) mz
        read(12,*) nx
        read(12,*) int_shot_model
        read(12,*) n_shot
        read(12,*) file_image_part
        read(12,*) file_image_stack
       close(12)

       allocate(image(1:mx,1:mz))
       allocate(image_shot(1:nx,1:mz))

       image(:,:)=0.0

       ii_bi=1
       do i_shot=1,n_shot

       file_image=trim(file_image_part)//char(int(i_shot/100)+48)//&
           char(int(mod(i_shot,100)/10)+48)//&
           char(mod(i_shot,10)+48)

        open(12,file=file_image,access='direct',form='unformatted',&
                recl=ii_bi*mz)
          do i=1,nx
             read(12,rec=i) (image_shot(i,j),j=1,mz)
          enddo 
        close(12)

        do i=1,nx
           do j=1,mz
              image((i_shot-1)*int_shot_model+i,j)=image((i_shot-1)*int_shot_model+i,j)+image_shot(i,j)
           enddo 
        enddo 

        enddo 

        open(12,file=file_image_stack,access='direct',&
                form='unformatted',&
                recl=ii_bi*mz)
          do i=1,mx
             write(12,rec=i) (image(i,j),j=1,mz)
          enddo 
        close(12)

!! high-pass filter
        call highpassfilter(mx,mz,image)

        open(12,file=trim(file_image_stack)//'lap',access='direct',&
                form='unformatted',recl=ii_bi*mz)
          do i=1,mx
             write(12,rec=i) (image(i,j),j=1,mz)
          enddo
        close(12)

       end
