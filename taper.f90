
!!    tapering seismograms   !!

! Coded by : Xinfa Zhu

! modified by: Peng Guo
! Date : June 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------

       implicit none 
       real,allocatable::taper(:),seis(:,:)
       integer::mx,nt,i,j,i_shot,n_shot
       integer::shot_start,shot_end
       real::ntpp,ntp1,tp,pi,rate
       character *90 file_in_part,file_out_part,file1,file2
       
       rate=0.5
!       mx=609
!       nt=4200
       open(12,file='parameter_taper.txt')
         read(12,*) mx
         read(12,*) nt
         read(12,*) shot_start
         read(12,*) shot_end
         read(12,*) file_in_part
         read(12,*) file_out_part
       close(12)

       allocate(taper(1:mx))

       taper(:)=1.0

       ntpp=mx/12.0
       ntp1=mx-ntpp
       pi=4.0*atan(1.0)
       tp=pi/ntpp

       allocate(seis(1:mx,1:nt))

         
       do i=1,mx

          if(i.le.floor(ntpp)) then 
              taper(i)=0.25*(4-3*cos(rate*float(i)*tp) )
          endif

          if(i.ge.floor(ntp1)) then 
              taper(i)=0.25*(1.0+3*cos((rate*(i-ntp1))*tp) )
          endif

       enddo      

     do i_shot=shot_start,shot_end

       file1=trim(file_in_part)//char(int(i_shot/100)+48)//&
             char(int(mod(i_shot,100)/10)+48)//&
             char(mod(i_shot,10)+48)

       open(12,file=file1,access='direct',form='unformatted',recl=4*mx)
         do i=1,nt
           read(12,rec=i) (seis(j,i),j=1,mx)
         enddo 
       close(12)

       do i=1,mx
          do j=1,nt
             seis(I,J)=seis(I,j)*taper(i)
          enddo 
       enddo 

       file2=trim(file_out_part)//char(int(i_shot/100)+48)//&
             char(int(mod(i_shot,100)/10)+48)//&
             char(mod(i_shot,10)+48)

       open(12,file=file2,access='direct',&
               form='unformatted',recl=4*mx)
         do i=1,nt
           write(12,rec=i) (seis(j,i),j=1,mx)
         enddo 
       close(12)

       enddo
       end
