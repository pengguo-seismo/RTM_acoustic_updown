
     include 'mpif.h'

     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

     num_shot=10
     do i=myid+1,num_shot,num_procs
        tt=i*1.0
        call mpi_reduce(tt,tt_sum,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)   
     enddo 

      print *,tt_sum

     call mpi_finalize(ierr)
     end
