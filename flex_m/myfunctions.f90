module myfunctions
#ifdef USE_MPI
    include 'mpif.h'
#endif
contains
    integer function mpi_rank()
        implicit none
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
        use constants, only: nb
        implicit none
        integer a, b
        sub_g2chi = a+(b-1)*nb
    end function sub_g2chi

    integer function sub_g2e(l,m,k,omega)
        use constants
        implicit none
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
        use constants
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
        use constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A0, B0, A, B, inverseAbyB
        integer info, lda, ldb, ipiv
        call cgesv(square_nb, square_nb, A, square_nb, ipiv, B, square_nb, info)
    end function inverseAbyB

    ! 需要测试, 考虑内存模式
    ! 因为都是方阵, 可考虑换函数
    function ABA(A, B)
        use constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A, B, ABA, C
        call cgemm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            A, square_nb, B, square_nb, complex_0, C, square_nb)
        call cgemm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            C, square_nb, A, square_nb, complex_0, ABA, square_nb)
    end function ABA

    ! 需要测试, 考虑内存模式
    function AB(A, B)
        use constants
        implicit none
        complex, dimension (nb*nb, nb*nb) :: A, B, AB
        call cgemm('N', 'N', square_nb, square_nb, square_nb, complex_1, &
            A, square_nb, B, square_nb, complex_0, AB, square_nb)
    end function AB

    subroutine matrixProduct(src, dft_matrix, dst, M, N, K)
        use constants
        implicit none
        integer M,N,K
        complex dft_matrix(K,N), src(M,K), dst(M,N)
        call cgemm('N', 'N', M, N, K, complex_1, &
            src, M, dft_matrix, K, complex_0, dst, M)
    end subroutine matrixProduct

    ! 生成傅里叶变换辅助矩阵dft_f, dft_b, idft_f, idft_b
    subroutine buildDFTMatrix()
        use constants
        use parameters2
        implicit none

        integer itau, tau, iomega, omega_f, omega_b

        ! 频域与时域的傅里叶变换辅助矩阵
        ! 可以一次计算一组, 或者构造大矩阵计算多组
        ! for one k-point, G_tau (行) = G_omega (行) * dft
        ! G_omega (行) = G_tau (行) * idft
        do itau = -ntau,ntau
            do iomega=-nomega,nomega
                omega_f = 2*iomega-1
                omega_b = 2*iomega
                tau = itau*1d0/ntau
                dft_f(iomega, itau) = exp(-2*pi*omega_f*tau*complex_i)
                dft_b(iomega, itau) = exp(-2*pi*omega_b*tau*complex_i)
                idft_f(itau, iomega) = exp(2*pi*omega_f*tau*complex_i)
                idft_b(itau, iomega) = exp(2*pi*omega_b*tau*complex_i)
            enddo
        enddo

        ! 这里的系数可能应该*2, 即所有频率一起考虑
        dft_f = dft_f / (2*nomega+1)
        dft_b = dft_b / (2*nomega+1)

        return
    end subroutine buildDFTMatrix

    ! dft频域到时域, fb表示为费米频率(1)或者玻色频率(0)
    subroutine dft(input, output, N, fb)
        use constants, only : nk, total_omega, total_tau, nomega, ntau
        use parameters2
        implicit none

        integer N, fb, l, m
        complex, dimension(N,N,nk,-nomega:nomega) :: input
        complex, dimension(N,N,nk,-ntau:ntau) :: output

        ! fb = mod(fb,2)
        do l=1,N; do m=1,N
            dft_omega = input(l,m,:,:)
            if (fb==0) then
                call matrixProduct(dft_omega, dft_b, dft_tau, nk, total_tau, total_omega)
            else
                call matrixProduct(dft_omega, dft_f, dft_tau, nk, total_tau, total_omega)
            endif
            output(l,m,:,:) = dft_tau
        enddo; enddo

        return

    end subroutine dft

    ! idft时域到频域
    subroutine idft(input, output, N, fb)
        use constants, only : nk, total_omega, total_tau, nomega, ntau
        use parameters2
        implicit none

        integer N, fb, l, m
        complex, dimension(N,N,nk,-nomega:nomega) :: output
        complex, dimension(N,N,nk,-ntau:ntau) :: input

        ! fb = mod(fb,2)
        do l=1,N; do m=1,N
            dft_tau= input(l,m,:,:)
            if (fb==0) then
                call matrixProduct(dft_tau, idft_b, dft_omega, nk, total_omega, total_tau)
            else
                call matrixProduct(dft_tau, idft_f, dft_omega, nk, total_omega, total_tau)
            endif
            output(l,m,:,:) = dft_omega
        enddo; enddo

        return

    end subroutine idft

#ifdef _DEBUG
    ! sometimes the linker cannot find this blas function
    real function scnrm2(N, A, incx)
        implicit none
        integer N, incx, i
        complex, dimension(N):: A
        scnrm2=0d0
        i=1
        do while (i<=N)
            scnrm2=scnrm2+abs(A(i))
        enddo
    end function scnrm2
#endif
END MODULE myfunctions
