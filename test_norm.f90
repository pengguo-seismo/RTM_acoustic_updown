
        real,allocatable::a(:,:)
        integer::seed

        seed=20171123

        nx=100
        nz=2000

        allocate(a(1:nx,1:nz))

        do j=1,nz
           do i=1,nx
              a(i,j)=r4_normal_01 ( seed )
           enddo 
        enddo 

        open(12,file='test_n')
         do j=1,nz
            do i=1,nx
               write(12,*) i,j,a(i,j)
            enddo 
         enddo 
        close(12)

        end


