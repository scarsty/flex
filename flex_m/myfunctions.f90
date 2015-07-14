module myfunctions
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

    integer function sub_g2chi(a,b)
        IMPLICIT NONE
        integer a, b
        sub_g2chi = a+(b-1)*5
    end function sub_g2chi

    !n计算松原频率
    integer function calfreq(omega, fb)
    implicit none
    integer omega, fb
            if (fb/=0) then
            calfreq=2*omega-1
        else
            calfreq=2*omega
        endif
    end function calfreq

    ! k和松原频率的减法, fb: fermi(1) or bose(0)
    ! 1 - 2 -> 3
    subroutine komega_minus(k1, omega1, fb1, k2, omega2, fb2, k, k_minus, zero_k, k3, omega3, fb3)
        use Constants
        implicit none
        integer k1, omega1, fb1, k2, omega2, fb2, k3, omega3, fb3, sign_omega3, zero_k
        integer f1, f2, f3
        real, dimension (nk, 2) :: k
        integer, dimension (nk, nk) :: k_minus

        f1=calfreq(omega1, fb1)
        f2=calfreq(omega2, fb2)
        f3=f1-f2

        fb3 = abs(mod(f3, 2))
        omega3 = (f3+fb3)/2
        k3 = k_minus(k1, k2)

    end subroutine komega_minus

END MODULE myfunctions
