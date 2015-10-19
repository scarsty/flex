module myfunctions
#ifdef USE_MPI
    include 'mpif.h'
#endif
contains
    ! mpi函数系列
    integer function mpi_rank()
        implicit none
        integer r
        integer ierr
#ifdef USE_MPI
        call mpi_comm_rank(MPI_COMM_WORLD, r, ierr)
#else
        r = 0
#endif /* USE_MPI */
        mpi_rank = r
    end function mpi_rank

    integer function mpi_size()
        implicit none
        integer s
        integer ierr
#ifdef USE_MPI
        call mpi_comm_size(MPI_COMM_WORLD, s, ierr)
#else
        s = 1
#endif /* USE_MPI */
        mpi_size = s
    end function mpi_size

    integer function mpi_init1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_init(ierr)
#else
        ierr = 0
#endif /* USE_MPI */
        mpi_init1 = ierr
    end function mpi_init1

    integer function mpi_finalize1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_finalize(ierr)
#else
        ierr = 0
#endif /* USE_MPI */
        mpi_finalize1 = ierr
    end function mpi_finalize1

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
        real(8), dimension (nk, 2) :: k
        integer, dimension (nk, nk) :: k_minus

        f1=calfreq(omega1, fb1)
        f2=calfreq(omega2, fb2)
        f3=f1-f2

        fb3 = abs(mod(f3, 2))
        omega3 = (f3+fb3)/2
        k3 = k_minus(k1, k2)

    end subroutine komega_minus

    function inverseAbyB(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, inverseAbyB
        integer info, lda, ldb
        integer ipiv(n)
        call zgesv(n, n, A, n, ipiv, B, n, info)
        inverseAbyB=B
        !write(stderr,*) info
    end function inverseAbyB

    ! 需要测试, 考虑内存模式
    ! 因为都是方阵, 可考虑换函数
    function ABA(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, ABA, C
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, C, n)
        call zgemm('N', 'N', n, n, n, complex_1, &
            C, n, A, n, complex_0, ABA, n)
    end function ABA

    function ABC(A, B, C, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, C, AB, ABC
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, AB, n)
        call zgemm('N', 'N', n, n, n, complex_1, &
            AB, n, C, n, complex_0, ABC, n)
    end function ABC

    ! 需要测试, 考虑内存模式
    function AB(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, AB
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, AB, n)
    end function AB

    ! 需要测试, 考虑内存模式
    function AHBA(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, AHBA, C
        call zgemm('C', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, C, n)
        call zgemm('N', 'N', n, n, n, complex_1, &
            C, n, A, n, complex_0, AHBA, n)
    end function AHBA

    function ABAH(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, ABAH, C
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, C, n)
        call zgemm('N', 'C', n, n, n, complex_1, &
            C, n, A, n, complex_0, ABAH, n)
    end function ABAH

    ! 松原频率转换, 数据结构设计如此
    function transfer_freq(freq)
        use constants
        implicit none
        integer freq, transfer_freq
        if (freq>=0) then
            transfer_freq=freq
        else
            transfer_freq=totalnomega+freq
        endif
    end function transfer_freq

    subroutine matrixProduct(src, dft_matrix, dst, M, N, K)
        use constants
        implicit none
        integer M,N,K
        complex(8) dft_matrix(K,N), src(M,K), dst(M,N)
        call zgemm('N', 'N', M, N, K, complex_1, &
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
        !        do itau = -ntau,ntau
        !            do iomega=-nomega,nomega
        !                omega_f = 2*iomega-1
        !                omega_b = 2*iomega
        !                tau = itau*1d0/ntau
        !                dft_f(iomega, itau) = exp(-2*pi*omega_f*tau*complex_i)
        !                dft_b(iomega, itau) = exp(-2*pi*omega_b*tau*complex_i)
        !                idft_f(itau, iomega) = exp(2*pi*omega_f*tau*complex_i)
        !                idft_b(itau, iomega) = exp(2*pi*omega_b*tau*complex_i)
        !            enddo
        !        enddo
        !
        !        ! 这里的系数可能应该*2, 即所有频率一起考虑
        !        dft_f = dft_f / (2*nomega+1)
        !        dft_b = dft_b / (2*nomega+1)

        return
    end subroutine buildDFTMatrix

    ! dft, direction means FORWORD(+) or BACKWORD(-)
    ! 调用fftw
    ! normal为归一化和展宽部分清零
    subroutine dft(input, output, N, direction, normal)
        use constants
        use parameters2
        use, intrinsic :: iso_c_binding
        implicit none

        include 'fftw3.f03'
        type(C_PTR) :: plan

        integer N, fb, l, m, normal, i
        integer(C_INT) direction, direction2
        complex(8), dimension(N,N,nkx,nky,0:totalnomega-1) :: input
        complex(8), dimension(N,N,nkx,nky,0:totalnomega-1) :: output

        if (direction >= 0) then
            direction2 = FFTW_FORWARD
        else
            direction2 = FFTW_BACKWARD
        endif

        do l=1,N; do m=1,N
            dft_in = input(l,m,:,:,:)

            plan=fftw_plan_dft_3d(totalnomega, nkx, nky , dft_in, dft_out, direction2, FFTW_ESTIMATE)
            call fftw_execute_dft(plan, dft_in, dft_out)
            call fftw_destroy_plan(plan)

            output(l,m,:,:,:) = dft_out
        enddo; enddo
        ! 卷积结果归一化
        if (normal /= 0) then
            output = output/nkx/nky/totalnomega
            ! 清空截断频率之外
            if (mod(normal,2)==1) then
                do i=0,2*nomega-2,2
                    output(:,:,:,:,i) = complex_0
                enddo
                do i=totalnomega-2*nomega+2,totalnomega,2
                    output(:,:,:,:,i) = complex_0
                enddo
                output(:,:,:,:,2*nomega:totalnomega-2*nomega) = complex_0
            else
                do i=1,4*nomega-3,2
                    output(:,:,:,:,i) = complex_0
                enddo
                do i=totalnomega-(4*nomega-3),totalnomega,2
                    output(:,:,:,:,i) = complex_0
                enddo
                output(:,:,:,:,4*nomega-1:totalnomega-(4*nomega-1)) = complex_0
            endif
        endif

        return

    end subroutine dft

    subroutine writematrix(A,n)
        use constants
        implicit none

        integer ix, iy, n
        complex(8), dimension(n,n) :: A

        do ix=1,n
            do iy=1,n
                write(stdout, '(A,2F7.3,A,$)') '(',A(ix,iy),' )  '
            enddo
            write(stdout,*)
        enddo

    end subroutine writematrix

    subroutine cleanError(A,n)
        use constants
        implicit none

        integer i, n
        complex(8), dimension(n) :: A

        do i=1,n
            if (abs(A(i))<real_error) then
                A(i)=0d0
            endif
        enddo

    end subroutine cleanError


    !pulay mixer 相关
    subroutine mixerInit()
        use constants
        use parameters
        use parameters2
        implicit none
        integer i

        G_mixer = 0
        mixer_pointer=1
        Pulay_A = 0
        do i=1,mix_num
            Pulay_A(mix_num+1,i)=-1
            Pulay_A(i,mix_num+1)=-1
        enddo
        Pulay_b = 0
        Pulay_b(mix_num+1) = -1
    end subroutine mixerInit

	! G1是新的, G是上一步
    subroutine mixer()
        use constants
        use parameters
        use parameters2
        implicit none
		integer n，i, ipiv, info
		
		n=mix_num+1
		
		G_mixer(:,:,:,:,:,mixer_pointer)=G1
		
		do i=1,mix_num
			Pulay_A(mixer_pointer,i)=
			Pulay_A(i,mixer_pointer)=
		enddo
        call dgesv(n, n, Pulay_A, n, ipiv, Pulay_B, n, info)
		G=complex_0
		do i=1,mix_num
			G=G+G_mixer(:,:,:,:,:,mixer_pointer)*Pulay_B(i)
		enddo
		mixer_pointer=mixer_pointer + 1
		if mixer_pointer>mix_num then
			mixer_pointer=1
		endif
    end subroutine mixer

    ! 部分废弃代码
    !subroutine buildkminus()
    !    use constants
    !    use parameters2
    !    implicit none
    ! k减法矩阵
    ! seems of no use
    !        k_minus=0
    !        do i1=1,nk
    !            do i2=1,nk
    !                temp = k(i1,:) - k(i2,:)
    !                !write(stdout, *) temp
    !                do while (abs(temp(1)+real_error)>0.5)
    !                    temp(1) = temp(1) -sign(1., temp(1))
    !                enddo
    !                do while (abs(temp(2)+real_error)>0.5)
    !                    temp(2) = temp(2) -sign(1., temp(2))
    !                enddo
    !                !write(stdout, *) temp
    !                do i = 1,nk
    !                    dis = norm2(temp-k(i,:))
    !                    if (dis<real_error) then
    !                        k_minus(i1, i2) = i
    !                        exit
    !                    endif
    !                enddo
    !                if (k_minus(i1, i2)<=0 .or. k_minus(i1, i2)>nk) then
    !                    write(stdout, *) 'Wrong k_minus at', i1, i2
    !                endif
    !            enddo
    !        enddo
    !
    !        ! k加法矩阵, k_minus(k1, k_minus(zero_k, k2))
    !        do i1=1,nk
    !            do i2=1,nk
    !                k_plus(i1, i2) = k_minus(i1, k_minus(zero_k, i2))
    !                if (k_plus(i1, i2)==0 .or. k_plus(i1, i2)>nk) then
    !                    write(stdout, *) 'Wrong k_plus at', i1, i2
    !                endif
    !            enddo
    !        enddo
    !end subroutine buildkminus


    !        h0_k_=h0_k(:,:,ikx,iky)
    !        if (k(ikx,iky,1)>=k(ikx,iky,2)) then
    !            h0_tilde_k(:,:,ikx,iky)=real(AHBA(i_plus,h0_k_,nb))
    !        else
    !            h0_tilde_k(:,:,ikx,iky)=real(AHBA(i_minus,h0_k_,nb))
    !        endif
    !        h0_tilde_k_=h0_tilde_k(:,:,ikx,iky)
    !        u_tilde_k_=h0_tilde_k_
    !        ! 这里h0_tilde变成一个实对称矩阵, 特征值全为实数, u为对应的正交变换阵
    !        call dsyev('V','U',nb,u_tilde_k_,nb,diag_h0_tilde_k_,diag_h0_tilde_k_lwork,nb*nb,info)
    !        u_tilde_k(:,:,ikx,iky)=u_tilde_k_
    !        diag_h0_tilde_k(:,ikx,iky)=diag_h0_tilde_k_
    !
    !    write(stdout,*) 'unitary matrix:'
    !    write(stdout,'(5F8.3)') u_tilde_k(:,:,ikx,iky)
    !    write(stdout,'(5F8.3)') u_tilde_k(:,:,ikx,iky)
    !    write(stdout,*) 'eigenvalue:'
    !    write(stdout,'(5F8.3)') diag_h0_tilde_k(:,ikx,iky)
    !    write(stdout,*) 'h:'
    !    call writematrix(h0_k(:,:,ikx,iky),nb)
    !
    !    write(stdout,*) 'h~:'
    !    write(stdout, '(5F8.3)') h0_tilde_k(:,:,ikx,iky)
    !
    !    write(stdout,*) 'u~ back to h~'
    !    diag_test=0d0
    !    do ix=1,nb
    !        diag_test(ix,ix)=diag_h0_tilde_k(ix,ikx,iky)
    !    enddo
    !    !write(stdout, '(5F8.3)') diag_test
    !    do ix=1,nb
    !        do iy=1,nb
    !            u_tilde_k_(ix,iy)=u_tilde_k(iy,ix,ikx,iky)
    !        enddo
    !    enddo
    !    !write(stdout, '(5F8.3)') u_tilde_k_
    !    write(stdout, '(5F8.3)') matmul(matmul(u_tilde_k_,diag_test), u_tilde_k(:,:,ikx,iky))

#ifdef _DEBUG
    ! sometimes the linker cannot find this blas function
    !    real function scnrm2(N, A, incx)
    !        implicit none
    !        integer N, incx, i
    !        complex, dimension(N):: A
    !        scnrm2=0d0
    !        i=1
    !        do while (i<=N)
    !            scnrm2=scnrm2+abs(A(i))
    !        enddo
    !    end function scnrm2
#endif
END MODULE myfunctions
