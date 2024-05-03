
      subroutine highpassfilter(nz,nx,mig)

        implicit none
        integer, intent(in)    :: nz,nx
        real:: mig(1:nz,1:nx)
        integer                :: i,j,k,n1,n2
        real, allocatable      :: mig1(:,:),mig2(:,:)

        n1=2
        n2=2
        allocate(mig1(nz,nx))
        allocate(mig2(nz,nx))
        mig2=mig
        do k=1,n1
           do j=2,nz-1
              do i=1,nx
                 mig1(j,i)=0.25*mig(j-1,i)+0.5*mig(j,i)+0.25*mig(j+1,i)
              enddo
           enddo

           do i=1,nx
              mig1(1,i)=0.75*mig(1,i)+0.25*mig(2,i)
              mig1(nz,i)=0.75*mig(nz,i)+0.25*mig(nz-1,i)
           enddo
           mig=mig1
        enddo

        do k=1,n2
          do j=1,nz
             do i=2,nx-1
               mig1(j,i)=0.25*mig(j,i-1)+0.5*mig(j,i)+0.25*mig(j,i+1)
             enddo
          enddo
          do j=1,nz
             mig1(j,1)=0.75*mig(j,1)+0.25*mig(j,2)
             mig1(j,nx)=0.75*mig(j,nx)+0.25*mig(j,nx-1)
          enddo
          mig=mig1
        enddo

        mig=mig2-mig1
        deallocate(mig1,mig2)

        end

