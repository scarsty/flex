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

    integer function sub_g2e(l,m,k,omega)
        use Constants
        IMPLICIT NONE
        integer l,m,k,omega
        integer omegat
        omegat=2*omega+1
        sub_g2e = l*nb*nk*omegat+m*nk*omegat+omega+nomega+1 ! 未完成
    end function sub_g2e

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

    function inverseAbyB(A0, B0)
        use Constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A0, B0, A, B, inverseAbyB
        integer info, lda, ldb, ipiv
        call cgesv(square_nb, square_nb, A, square_nb, ipiv, B, square_nb, info)
    end function inverseAbyB

    ! 需要测试, 考虑内存模式
    function ABA(A, B)
        use Constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A, B, ABA, C
        call ctrmm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            A, square_nb, B, square_nb, complex_0, C, square_nb)
        call ctrmm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            C, square_nb, A, square_nb, complex_0, ABA, square_nb)
    end function ABA

    ! 需要测试, 考虑内存模式
    function AB(A, B)
        use Constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A, B, AB
        call ctrmm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            A, square_nb, B, square_nb, complex_0, AB, square_nb)
    end function AB

END MODULE myfunctions
