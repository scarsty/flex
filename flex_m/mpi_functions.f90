module mpi_functions
#ifdef USE_MPI
    include 'mpif.h'
#endif
contains
    integer function mpi_rank()
        IMPLICIT NONE
        integer rank
        real ierr
#ifdef USE_MPI
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
        rank = 0
#endif /* USE_MPI */
        mpi_rank = rank
    end function mpi_rank

END MODULE mpi_functions
