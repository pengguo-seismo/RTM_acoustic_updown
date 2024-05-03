
!!  subtracting the direct waves from seismograms !!

! Coded by : Peng Guo
! Date : January 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------

      implicit none

      real,allocatable::p1(:,:),p2(:,:)

      integer::ii_bi,nx_rec,nt

      integer::it,i,n_shot,i_shot
      integer::shot_start,shot_end

      character *100 file1,file2,file3,file_p_part,&
                    file_p_dir_part,file_p_sub_part,&
                    file_vz_part,file_vz_dir_part,file_vz_sub_part

      open(12,file='parameter_sub.txt')
        read(12,*) nx_rec
        read(12,*) nt
        read(12,*) shot_start
        read(12,*) shot_end
        read(12,*) file_p_part      
        read(12,*) file_p_dir_part      
        read(12,*) file_p_sub_part    
      close(12)

      allocate(p1(1:nx_rec,1:nt),p2(1:nx_rec,1:nt))

      ii_bi=1

      do i_shot=shot_start,shot_end

      print *,i_shot,'i_shot'
      file1=trim(file_p_part)//char(int(i_shot/100)+48)//&
             char(int(mod(i_shot,100)/10)+48)//&
             char(mod(i_shot,10)+48)
      file2=trim(file_p_dir_part)//char(int(i_shot/100)+48)//&
             char(int(mod(i_shot,100)/10)+48)//&
             char(mod(i_shot,10)+48)

      open(12,file=file1,access='direct',form='unformatted',&
              recl=ii_bi*nx_rec)
         do it=1,nt
            read(12,rec=it) (p1(i,it),i=1,nx_rec)
         enddo 
      close(12)
      open(12,file=file2,access='direct',form='unformatted',&
              recl=ii_bi*nx_rec)
         do it=1,nt
            read(12,rec=it) (p2(i,it),i=1,nx_rec)
         enddo 
      close(12)

      p1(:,:)=p1(:,:)-p2(:,:)

      file3=trim(file_p_sub_part)//char(int(i_shot/100)+48)//&
             char(int(mod(i_shot,100)/10)+48)//&
             char(mod(i_shot,10)+48)

      open(12,file=file3,access='direct',form='unformatted',&
              recl=ii_bi*nx_rec)
         do it=1,nt
            write(12,rec=it) (p1(i,it),i=1,nx_rec)
         enddo 
      close(12)

      enddo 
      end
