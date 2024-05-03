      program test_random_seed
      implicit none
      real, allocatable :: seed(:)
      integer, allocatable :: seed_i(:)
      integer::irand
       integer :: n,i

       n=1000
!       call random_seed(size = n)
       allocate(seed(n))
       allocate(seed_i(n))
       call random_number(seed)
!       call random_seed(get=seed)
       print *,'n',n
       write (*, *) seed
       open(12,file='rand.dat')
       do i=1,n
          write(12,*) seed(i)
       enddo 
       close(12)

!       call init_random_seed()
       do i=1,n
          seed_i(i)=irand()
       enddo 
          
       open(12,file='rand_int.dat')
       do i=1,n
          write(12,*) seed_i(i)
       enddo 
       close(12)

       end program test_random_seed
