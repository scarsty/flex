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
            transfer_freq=dft_grid+freq
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
        begin_index=0
        if (outmodel/=0) then
            begin_index=2*nomega
        endif

        do l=1,N; do m=1,N
            ! 前处理, 补0
            dft_in(:,:,1:length_in) = input(l,m,:,:,1:length_in)
            dft_in(:,:,length_in+1:dft_grid)=complex_0

            plan=fftw_plan_dft_3d(dft_grid, nkx, nky, dft_in, dft_out, direction2, FFTW_ESTIMATE)
            call fftw_execute_dft(plan, dft_in, dft_out)
            call fftw_destroy_plan(plan)

            ! 后处理
            output(l,m,:,:,1:length_out) = dft_out(:,:,begin_index:length_out-begin_index+1)
        enddo; enddo
        ! 所有情况都在反傅里叶变换时归一化
        if (direction2==FFTW_BACKWARD) then
            output = output/nk/dft_grid
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

    ! pulay mixer 相关
    ! 移动指针
    function mixerIncPointer(p, n)
        use parameters2
        implicit none
        integer mixerIncPointer, p, n

        mixerIncPointer = p+n
        mixerIncPointer = mod(mixerIncPointer, mix_num)
        if (mixerIncPointer==0) mixerIncPointer=mix_num
    end function mixerIncPointer

    ! 求模的平方
    function mixerErrorProduct(a, b)
        use constants
        implicit none
        complex(8) mixerErrorProduct
        complex(8), dimension (nb, nb, nkx, nky, 0:dft_grid-1) :: a, b
        integer ib1,ib2,ikx,iky,iomegak

        mixerErrorProduct=0
        do ib1=1,nb;do ib2=1,nb
            do ikx=1,nkx;do iky=1,nky
                do iomegak=minomegaf,maxomegaf
                    mixerErrorProduct=mixerErrorProduct &
                        + conjg(a(ib1,ib2,ikx,ikx,iomegak)) &
                        * b(ib1,ib2,ikx,ikx,iomegak)
                enddo
            enddo;enddo
        enddo;enddo
    end function mixerErrorProduct

    ! 初始化混合器
    subroutine mixerInit()
        use parameters2
        implicit none
        integer i

        mixer_G = 0
        !G_mixer(:,:,:,:,:,1)=G1
        !error_mixer(:,:,:,:,:,1) = G1-G0

        mixer_pointer=1
        mixer_A = 0
        do i=1,mix_num
            mixer_A(0,i)=-1
            mixer_A(i,0)=-1
        enddo
        mixer_b = 0
        mixer_b(0) = -1
    end subroutine mixerInit

    ! G1是新的, G是上一步
    ! 混合算法
    ! http://vergil.chemistry.gatech.edu/notes/diis/node2.html
    subroutine mixer(num)
        use parameters2
        implicit none
        integer num, n, i, info, next_pointer, prev_pointer
        integer ipiv(mix_num+1)
        complex(8) zdotc
        external zdotc
        complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf):: b1,b2
        complex(8), dimension (mix_num*2) :: lwork
        complex(8) e

        !n=nb*nb*nk*totalnomega

        next_pointer=mixerIncPointer(mixer_pointer,1)
        prev_pointer=mixerIncPointer(mixer_pointer,-1)
        !write(*,*) 'kkk',mixerErrorProduct(G1,G1),mixerErrorProduct(G,G)
        if (num==1) then
            mixer_error(:,:,:,:,:,mixer_pointer)=G1-G
            !write(*,*) 'hhhh'
        else
            mixer_error(:,:,:,:,:,mixer_pointer)=G1-G!_mixer(:,:,:,:,:,prev_pointer)

            !write(*,*) 'klkl',mixerErrorProduct(error_mixer(:,:,:,:,:,mixer_pointer),error_mixer(:,:,:,:,:,mixer_pointer))
        endif
        mixer_G(:,:,:,:,:,mixer_pointer)=G1
        !sigma_mixer(:,:,:,:,:,mixer_pointer)=sigma

        ! A_ij=e_i**H*e_j
        do i=1,mix_num
            b1=mixer_error(:,:,:,:,:,mixer_pointer)
            b2=mixer_error(:,:,:,:,:,i)
            e=mixerErrorProduct(b1,b2)
            mixer_A(mixer_pointer,i)=e
            mixer_A(i,mixer_pointer)=conjg(e)
            !write(*,*)e
        enddo
        !Pulay_A(0,mixer_pointer)=-1
        !Pulay_A(mixer_pointer,0)=-1

        mixer_A1=mixer_A
        mixer_x=mixer_b
        ! 系数矩阵实际上多一行
        n=min(num+1,mix_num+1)
        call zhesv('U', n, 1, mixer_A1, mix_num+1, ipiv, mixer_x, mix_num+1, lwork, 2*mix_num, info)
        !write(*,*) info, mixer_pointer
        !call zgesv(n, 1, Pulay_A1, mix_num+1, ipiv, Pulay_x, mix_num+1, info)
        G=complex_0
        do i=1,1
            G=G+mixer_G(:,:,:,:,:,i)*real(mixer_x(i))
            !write(*,*) Pulay_x(i), Pulay_A(i,i), mixerErrorProduct(G_mixer(:,:,:,:,:,i),G_mixer(:,:,:,:,:,i))
        enddo
        mixer_pointer=next_pointer
        !call writematrix(Pulay_A,11)
        !stop
    end subroutine mixer

    subroutine outputSomething()
        implicit none
        !    输出部分结果检查
        !    write(stdout,*) dot_product(u_tilde_k_(1,:),u_tilde_k_(5,:))
        !
        !    ikx=2;iky=4
        !    write(stdout,*) 'h:'
        !    call writematrix(h0_k(:,:,ikx,iky),nb)
        !    write(stdout,*) 'unitary matrix:'
        !    call writematrix(u_h0_k(:,:,ikx,iky),nb)
        !
        !    write(stdout,*) 'eigenvalue:'
        !    write(stdout,*) ev_h0_k(:,ikx,iky)
        !
        !    write(stdout,*) 'u back to h'
        !    diag_h0_G0_=0d0
        !    do ix=1,nb
        !        !diag_h0_G0_(ix,ix)=complex_1
        !        diag_h0_G0_(ix,ix)=ev_h0_k(ix,ikx,iky)
        !        !diag_h0_G0_(ix,ix)=1-complex_i*ix
        !    enddo
        !    !write(stdout, '(5F8.3)') diag_test
        !    u_h0_k_=u_h0_k(:,:,ikx,iky)
        !
        !    !u_h0_k_=u_h0_k(:,:,ikx,iky)
        !    !call writematrix(matmul(u_h0_k_,u_h0_k(:,:,ikx,iky)),nb)
        !
        !    call writematrix(ABAH(u_h0_k_,diag_h0_G0_,nb),nb)
        !
        !
        !    write(stdout, *)'next'
        !
        !    h0_k_=h0_k(:,:,ikx,iky)
        !    h0_tilde_k_=real(AHBA(i_plus,h0_k_,nb))
        !    write(stdout, *)'h~'
        !    call writematrix(AHBA(i_plus,h0_k_,nb),nb)
        !    u_tilde_k_=h0_tilde_k_
        !    call dsyev('V','U',nb,u_tilde_k_,nb,ev_h0_k_,diag_h0_tilde_k_lwork,nb*nb,info)
        !    write(stdout,*) 'eigenvalue:'
        !    write(stdout,*) ev_h0_k_
        !    u_h0_k_=u_tilde_k_
        !
        !    u_h0_k_=AB(i_plus,u_h0_k_,nb)
        !    write(stdout,*) 'unitary matrix:'
        !    call writematrix(u_h0_k_,nb)
        !    write(stdout,*) 'u~ back to h~'
        !    call writematrix(ABAH(u_h0_k_,diag_h0_G0_,nb),nb)
    end subroutine

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
    !        enddo!        if (abs(cur_density-target_density)<density_tol) then
    !            density_conv=.true.
    !            !计算结束
    !        else
    !            ! 根据占据数调整化学势
    !            ! 第一步仅记录和猜测方向
    !            ! 第三步开始逐步抛弃较远的点, 依照线性趋势逼近
    !            ! 通常来说应保存一大一小
    !            ! 靠不太容易设计
    !            if (density_iter>=2) then
    !                if (density_iter>2) then
    !                    replaced=.false.
    !                    do i=1,2
    !                        if (abs(cur_density-target_density)<abs(density_old(i)-target_density)) then
    !                            mu_old(i)=mu
    !                            density_old(i)=cur_density
    !                            replaced=.true.
    !                            exit
    !                        endif
    !                    enddo
    !                    if (.not.replaced) then
    !                        max_diff_loc=1
    !                        if (abs(density_old(1)-target_density)<abs(density_old(2)-target_density)) then
    !                            max_diff_loc=2
    !                        endif
    !                        mu_old(max_diff_loc)=mu
    !                        density_old(max_diff_loc)=cur_density
    !                    endif
    !                else
    !                    mu_old(2)=mu
    !                    density_old(2)=cur_density
    !                endif
    !                mu=(mu_old(1)-mu_old(2))/(density_old(1)-density_old(2))*(target_density-density_old(2))+mu_old(2)
    !            elseif (density_iter==1) then
    !                mu_old(1)=mu
    !                density_old(1)=cur_density
    !                !deltamu_per_density = (mu-mu0)/(cur_density-density0)
    !                mu = mu - 1.0d-1*sign(1.0d0, (cur_density-target_density)*deltamu_per_density)
    !            endif
    !            write(stdout,*) 'modified new mu = ', mu
    !        endif
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
