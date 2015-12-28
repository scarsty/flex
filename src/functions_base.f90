module functions_base
#ifdef USE_MPI
    include 'mpif.h'
#endif
    use constants
    real(8), external :: dznrm2

contains
    ! mpi函数系列
    integer function mpi_rank1()
        implicit none
        integer r
        integer ierr
#ifdef USE_MPI
        call mpi_comm_rank(MPI_COMM_WORLD, r, ierr)
#else
        r = 0
#endif
        mpi_rank1 = r
    end function

    integer function mpi_size1()
        implicit none
        integer s
        integer ierr
#ifdef USE_MPI
        call mpi_comm_size(MPI_COMM_WORLD, s, ierr)
#else
        s = 1
#endif
        mpi_size1 = s
    end function

    integer function mpi_init1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_init(ierr)
#else
        ierr = 0
#endif
        mpi_init1 = ierr
    end function

    integer function mpi_finalize1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_finalize(ierr)
#else
        ierr = 0
#endif
        mpi_finalize1 = ierr
    end function

    integer function mpi_reduce1(A,n)
        implicit none
        integer ierr, n
        complex(8) A(n),B(n)
#ifdef USE_MPI
        call mpi_allreduce(A,B,16,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
        ierr = 0
#endif
        mpi_reduce1 = ierr
    end function

    ! 下标变换
    integer function sub_g2chi(a,b)
        use constants, only: nb
        implicit none
        integer a, b
        sub_g2chi = a+(b-1)*nb
    end function

    ! B直接保存结果
    subroutine inverseAbyB(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B
        integer info, lda, ldb
        integer ipiv(n)
        call zgesv(n, n, A, n, ipiv, B, n, info)
        !write(stderr,*) info
    end subroutine

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
    end function

    function ABC(A, B, C, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, C, AB, ABC
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, AB, n)
        call zgemm('N', 'N', n, n, n, complex_1, &
            AB, n, C, n, complex_0, ABC, n)
    end function

    ! 需要测试, 考虑内存模式
    function AB(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, AB
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, AB, n)
    end function

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
    end function

    function ABAH(A, B, n)
        use constants
        implicit none
        integer n
        complex(8), dimension (n, n) :: A, B, ABAH, C
        call zgemm('N', 'N', n, n, n, complex_1, &
            A, n, B, n, complex_0, C, n)
        call zgemm('N', 'C', n, n, n, complex_1, &
            C, n, A, n, complex_0, ABAH, n)
    end function

    ! 松原频率转换, 数据结构设计如此
    function transfer_freq(freq)
        use constants
        implicit none
        integer freq, transfer_freq
        if (freq>=0) then
            transfer_freq=freq
        else
            transfer_freq=dft_grid+freq
        endif
    end function

    subroutine matrixProduct(src, dft_matrix, dst, M, N, K)
        use constants
        implicit none
        integer M,N,K
        complex(8) dft_matrix(K,N), src(M,K), dst(M,N)
        call zgemm('N', 'N', M, N, K, complex_1, &
            src, M, dft_matrix, K, complex_0, dst, M)
    end subroutine

    ! 交换
    subroutine swap_r8(a,b)
        real(8) a,b,c
        c=a
        a=b
        b=c
    end subroutine

    ! 快速排序
    recursive subroutine quick_sort(v,n,s,e)
        real(8) v(n),key
        integer n,s,e,l,r,m

        l=s
        r=e
        m=(s+e)/2
        if (l>=r) return
        key=v(m)
        do while (l<r)
            do while(v(l)<key)
                l=l+1
            enddo
            do while(v(r)>key)
                r=r-1
            enddo
            if (l<r) then
                call swap_r8(v(l),v(r))
            endif
        enddo
        call quick_sort(v,n,s,l-1)
        call quick_sort(v,n,r+1,e)
    end subroutine

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

    end subroutine

    subroutine cleanError_c(A,n)
        use constants
        implicit none

        integer i, n
        complex(8), dimension(n) :: A

        do i=1,n
            if (abs(A(i))<real_error) then
                A(i)=0d0
            endif
        enddo

    end subroutine

    subroutine cleanError(A,n)
        use constants
        implicit none

        integer i, n
        complex(8) A(n)
        real(8) r1,r2

        do i=1,n
            r1=real(A(i))
            r2=imag(A(i))
            if (abs(r1)<real_error) then
                r1=0d0
            endif
            if (abs(r2)<real_error) then
                r2=0d0
            endif
            A(i)=complex_1*r1+complex_i*r2
        enddo

    end subroutine


    ! 求两个G格式向量的点乘 sum conjg(a(i))*b(i)
    function GProduct(a, b)
        use constants
        implicit none
        complex(8) GProduct
        complex(8), dimension (total_grid) :: a, b
        complex(8), external :: zdotc
        !GProduct=dot_product(a,b)
        GProduct = zdotc(total_grid,a,1,b,1)
    end function


    ! dft, direction means FORWORD(+) or BACKWORD(-)
    ! 调用fftw
    ! 傅里叶变换用于计算卷积, 包含以下两种情况:
    ! (1) G*conjgG的卷积, 输入都是费米频率, 输出是玻色频率.
    !     这时需要补充2*nomega的0. 反傅里叶变换的输出结果从0开始
    ! (2) V*G的卷积, 输入是玻色频率和费米频率, 输出是费米频率.
    !     这时前者只需补一个0, 后者同(1)补2*nomega的0. 反傅里叶变换的结果从2*nomega开始.
    ! 变换的网格和傅里叶变换的输出都是4*nomega, 反傅里叶变换的输出长度随模式而定.
    ! outmodel: 0 - 从0开始输出, 大部分情况; 其他 - 从n*nomega输出, 仅适合自能.
    subroutine dft(input, output, N, length_in, length_out, direction, outmodel)
        use constants
        use parameters2
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'
        complex(8), dimension(N,N,nkx,nky,length_in) :: input
        complex(8), dimension(N,N,nkx,nky,length_out) :: output
        integer N, length_in, length_out, direction, outmodel
        integer(C_INT) direction2, l, m, begin_index
        type(C_PTR) :: plan

        if (direction >= 0) then
            direction2 = FFTW_FORWARD
        else
            direction2 = FFTW_BACKWARD
        endif

        ! 输出的起始下标
        begin_index=1
        if (outmodel/=0) then
            begin_index=2*nomega
        endif
        plan=fftw_plan_dft_3d(dft_grid, nky, nkx, dft_in, dft_out, direction2, FFTW_ESTIMATE)
        !!$omp parallel do private(m,dft_in,dft_out)
        do l=1,N; do m=1,N
            ! 前处理, 补0
            dft_in(:,:,1:length_in) = input(l,m,:,:,1:length_in)
            dft_in(:,:,length_in+1:dft_grid)=complex_0

            call fftw_execute_dft(plan, dft_in, dft_out)

            ! 后处理
            output(l,m,:,:,1:length_out) = dft_out(:,:,begin_index:length_out-begin_index+1)
        enddo; enddo
        !!$omp end parallel do
        call fftw_destroy_plan(plan)
        ! 所有情况都在反傅里叶变换时归一化
        if (direction2==FFTW_BACKWARD) then
            output = output/nk/dft_grid
        endif

        return

    end subroutine

    subroutine build_h0_k()
        use constants
        use parameters2
        implicit none

        integer l1, m1, ikx, iky

        h0_k = complex_0

        do l1=1,nb; do m1=1,nb
            do ikx=1,nkx; do iky=1,nky
                h0_k(l1,m1,ikx,iky) = - cos(k(ikx,iky,1)*pi) - cos(k(ikx,iky,2)*pi)
            enddo; enddo
        enddo; enddo
        write(stderr,*) h0_k
        return

    end subroutine build_h0_k

    subroutine get_tick(t)
        implicit none
        real(8) t
        integer(8) t1, clock_rate, clock_max
        call system_clock(t1,clock_rate,clock_max)
        t=1d0*t1/clock_rate
    end subroutine

    subroutine output_date()
        implicit none
        integer d(8)
        call date_and_time(values=d)
        write (stdout,'(I5,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I3.3)') &
            d(1),'-',d(2),'-',d(3),'/',d(5),':',d(6),':',d(7),'.',d(8)
    end subroutine

end module
