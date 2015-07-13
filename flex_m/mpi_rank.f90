#ifdef USE_MPI
include 'mpif.h'
#endif

function mpi_rank(a, b)
#ifdef USE_MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0
#endif /* USE_MPI */
    mpi_rank = rank
end function
