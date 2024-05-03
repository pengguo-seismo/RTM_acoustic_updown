      double precision,allocatable::par(:),par_c(:)
      integer::iter,ii_layer,totalpar
      totalpar=5
      num_par=6
      ii_layer=2

      allocate(par(1:num_par),par_c(1:num_par))
      par=0.d0
      totaliter=12

      open(12,file='model_par_correct.par')
      do ii=1,ii_layer
         do i=1,num_par
           read(12,*) par_c(i)
         enddo 
      enddo  

      open(12,file='parameter_hty.par')
      open(14,file='para_hty_mat.par')
      do ii=1,totaliter
         read(12,*) iter
         read(12,*) ii_layer
         do i=1,num_par
            read(12,*) par(i)
         enddo 

         write(14,'(I2,12F8.2)') iter,(par_c(i),par(i),i=1,num_par)
      enddo 
      close(12)
      close(14)
      end
