
      real,allocatable::a(:,:),b(:,:),b_i(:,:),b_r(:,:)
      complex,allocatable::a_c(:,:),b_c(:,:)
      
      nt=2600
      nx=301

      allocate(a(1:nx_rec,1:nt),b(1:nx_rec,1:nt),&
               b_i(1:nx_rec,1:nt),b_r(1:nx_rec,1:nt))
      allocate(a_c(1:nx_rec,1:nt),b_c(1:nx_rec,1:nt))

      a=0.0
      b=0.0
      b_i=0.0
      b_r=0.0

      a_c=0.0
      b_c=0.0

      pi=4.0*atan(1.0)

      open(12,file='seis_p001',access='direct',form='unformatted',&
              recl=4*nx_rec)
       do it=1,nt
          read(12,rec=it) (a(i,it),i=1,nx_rec)
       enddo 
      close(12)

      call hilbert_seismo1(a,b,nx_rec,nt)

      b_i(1:nx_rec,1:nt)=sqrt(a(1:nx_rec,1:nt)*a(1:nx_rec,1:nt)&
                        +b(1:nx_rec,1:nt)*b(1:nx_rec,1:nt))

      open(12,file='seis_p_env001',access='direct',form='unformatted',&
              recl=4*nx_rec)
       do it=1,nt
          write(12,rec=it) (b_i(i,it),i=1,nx_rec)
       enddo 
      close(12)

      end 
