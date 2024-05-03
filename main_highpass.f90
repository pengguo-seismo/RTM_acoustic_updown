
        real,allocatable::image(:,:),image_lap(:,:)

        integer::mx,mz

        mx=648
        mz=291

        ii_bi=1

        allocate(image(1:mx,1:mz),image_lap(1:mx,1:mz))

        open(12,file='seis_p_image_hes_sum',access='direct',&
                form='unformatted',recl=ii_bi*mz)
          do i=1,mx
             read(12,rec=i) (image(i,j),j=1,mz)
          enddo
        close(12)

        call highpassfilter(mx,mz,image)

        image_lap=image

        open(12,file='seis_p_image_hes_sum_lap',access='direct',&
                form='unformatted',recl=ii_bi*mz)
          do i=1,mx
             write(12,rec=i) (image_lap(i,j),j=1,mz)
          enddo
        close(12)

       end
