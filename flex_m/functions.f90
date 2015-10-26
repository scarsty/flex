module functions
#ifdef USE_MPI
    include 'mpif.h'
#endif

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
    end function mpi_rank1

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
    end function mpi_size1

    integer function mpi_init1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_init(ierr)
#else
        ierr = 0
#endif
        mpi_init1 = ierr
    end function mpi_init1

    integer function mpi_finalize1()
        implicit none
        integer ierr
#ifdef USE_MPI
        call mpi_finalize(ierr)
#else
        ierr = 0
#endif
        mpi_finalize1 = ierr
    end function mpi_finalize1

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
    end function mpi_reduce1

    ! 下标变换
    integer function sub_g2chi(a,b)
        use constants, only: nb
        implicit none
        integer a, b
        sub_g2chi = a+(b-1)*nb
    end function sub_g2chi

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

    subroutine init()
        use parameters
        use parameters2
        implicit none
        integer i

        T_beta = 1d0/kB/T
        T_eV = kB*T

        sigma_state = 0

        ! I_chi, (nb*nb) order
        I_chi=complex_0
        do i=1,nb*nb
            I_chi(i,i)=complex_1
        enddo

        ! I_G, nb order
        I_G=complex_0
        do i=1,nb
            I_G(i,i)=complex_1
        enddo

        mu_less_count=0
        mu_more_count=0

    end subroutine

    subroutine init_Kpoints()
        use parameters
        use parameters2
        implicit none
        integer ikx, iky, irx, iry, info
        real(8) rdotk, temp(2)
        complex(8) fac

        ! 计算k点的坐标
        write(stdout, *) "Building k-points grid..."
        !zero_k = 1    ! k原点
        do ikx = 1, nkx
            do iky = 1, nky
                k(ikx, iky, 1)=1d0/nkx*(ikx-1)
                k(ikx, iky, 2)=1d0/nky*(iky-1)
                if (k(ikx, iky, 1)>0.5) k(ikx, iky, 1)=k(ikx, iky, 1)-1
                if (k(ikx, iky, 2)>0.5) k(ikx, iky, 2)=k(ikx, iky, 2)-1
                write(stdout, '(2I3,2F9.4)') ikx, iky, k(ikx,iky,:)
            enddo
        enddo
        write(stdout, *)

        ! 反傅里叶变换h0到k空间
        h0_k = complex_0
        do ikx=1,nkx; do iky=1,nky
            do irx=-rx,rx; do iry=-ry,ry
                temp=[irx,iry]
                rdotk = two_pi*dot_product(k(ikx,iky,:),temp)
                fac=exp(complex_i*rdotk)
                h0_k(:,:,ikx,iky)=h0_k(:,:,ikx,iky)+fac*h0_r(:,:,irx,iry)
            enddo; enddo

            ! 计算特征值和特征向量
            h0_k_=h0_k(:,:,ikx,iky)
            u_h0_k_=h0_k_
            call zheev('V','L',nb,u_h0_k_,nb,ev_h0_k_,ev_h0_k_lwork,nb*nb,ev_h0_k_rwork,info)
            u_h0_k(:,:,ikx,iky)=u_h0_k_
            ev_h0_k(:,ikx,iky)=ev_h0_k_
            ! write(stdout,*) diag_h0_tilde_k_
        enddo; enddo

    end subroutine

    subroutine init_U()
        use parameters
        use parameters2
        implicit none
        integer ix, iy
        ! U
        ! 能带下标ab, cd -> (a+(b-1)*nb, c+(d-1)*nb)
        ! real, dimension (nb*nb, nb*nb):: U_s, U_c, U_ud, U_uu
        U_ud = 0d0
        U_uu = 0d0
        do ix=1,nb
            do iy=1,nb
                if (ix==iy) then
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=h1_U
                else
                    U_ud(sub_g2chi(ix,ix), sub_g2chi(iy,iy))=h1_Up
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=h1_J
                    U_ud(sub_g2chi(ix,iy), sub_g2chi(iy,ix))=h1_Jp
                    U_uu(sub_g2chi(ix,ix), sub_g2chi(iy,iy))=h1_Up-h1_J
                    U_uu(sub_g2chi(ix,iy), sub_g2chi(ix,iy))=-h1_Up+h1_J
                endif
            enddo
        enddo
        U_s = U_ud-U_uu
        U_c = U_ud+U_uu
    end subroutine


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
        !!$omp parallel do private(m,dft_in,dft_out,plan)
        do l=1,N; do m=1,N
            ! 前处理, 补0
            dft_in(:,:,1:length_in) = input(l,m,:,:,1:length_in)
            dft_in(:,:,length_in+1:dft_grid)=complex_0

            plan=fftw_plan_dft_3d(dft_grid, nky, nkx, dft_in, dft_out, direction2, FFTW_ESTIMATE)
            call fftw_execute_dft(plan, dft_in, dft_out)
            call fftw_destroy_plan(plan)

            ! 后处理
            output(l,m,:,:,1:length_out) = dft_out(:,:,begin_index:length_out-begin_index+1)
        enddo; enddo
        !!$omp end parallel do
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
        complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: a, b
        integer ib1,ib2,ikx,iky,iomegak

        mixerErrorProduct=0
        do ib1=1,nb;do ib2=1,nb
            do ikx=1,nkx;do iky=1,nky
                do iomegak=minomegaf,maxomegaf
                    mixerErrorProduct=mixerErrorProduct &
                        + conjg(a(ib1,ib2,ikx,iky,iomegak)) &
                        * b(ib1,ib2,ikx,iky,iomegak)
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
        !$omp parallel do private(b1,b2,e)
        do i=1,mix_num
            b1=mixer_error(:,:,:,:,:,mixer_pointer)
            b2=mixer_error(:,:,:,:,:,i)
            e=mixerErrorProduct(b1,b2)
            mixer_A(mixer_pointer,i)=e
            mixer_A(i,mixer_pointer)=conjg(e)
            !write(*,*)e
        enddo
        !$omp end parallel do
        !Pulay_A(0,mixer_pointer)=-1
        !Pulay_A(mixer_pointer,0)=-1

        mixer_A1=mixer_A
        mixer_x=mixer_b

        n=min(num,mix_num)
        ! 系数矩阵实际上多一行
        call zhesv('U', n+1, 1, mixer_A1, mix_num+1, ipiv, mixer_x, mix_num+1, lwork, 2*mix_num, info)
        !write(*,*) info, mixer_pointer
        !call zgesv(n, 1, Pulay_A1, mix_num+1, ipiv, Pulay_x, mix_num+1, info)
        G=complex_0
        do i=1,n
            G=G+mixer_G(:,:,:,:,:,i)*real(mixer_x(i))
            !write(*,*) mixer_x(i), mixer_A(i,i)
        enddo
        mixer_pointer=next_pointer
        !call writematrix(Pulay_A,11)
        !stop
    end subroutine mixer


    ! 检查收敛点数, 注意使用时机
    ! model: 1-sigma, 2-G
    function convergence_test(iter, model)
        use parameters
        use parameters2
        implicit none
        logical convergence_test
        real(8) dznrm2
        external dznrm2
        integer ib1, ib2, ikx, iky, iomegak, conv_grid, iter, model
        real(8) norm_sigma_minus, norm_sigma, norm_sigma0, cur_sigma_tol, tol

        ! 计算sigma0与sigma的符合情况, 向量库
        ! dznrm2: 欧几里得模，行向量乘以自身转置共轭
        if (model==0) then
            sigma_minus = sigma0 - sigma
        else
            return
        endif

        norm_sigma_minus = dznrm2(total_grid, sigma_minus, 1)
        norm_sigma = dznrm2(total_grid, sigma, 1)
        tol = norm_sigma*sigma_tol
        conv_grid=0
        do ib1=1,nb; do ib2=1,nb; do ikx=1,nkx; do iky=1,nky; do iomegak=minomegaf,maxomegaf
            if (abs(sigma_minus(ib1,ib2,ikx,iky,iomegak))<tol) then
                conv_grid=conv_grid+1
            endif
        enddo; enddo; enddo; enddo; enddo;

        cur_sigma_tol = norm_sigma_minus / norm_sigma
        write(stdout,'(I7,I10,ES20.5)') iter, conv_grid, cur_sigma_tol
        convergence_test  = (conv_grid==total_grid)

#ifdef _DEBUG
        !norm_sigma0 = dznrm2(nb*nb*nkx*nky*nomegaf, sigma0, 1)
        !write(stdout,*) '0:',norm_sigma0, '1:',norm_sigma
        !write(stdout,*) '0-1:',norm_sigma_minus
#endif
    end function

    ! 1~3数组, 第一个保存最接近的值, 后面两个保存最近两次计算的值
    ! warning表示最新的值并未更加靠近, 可能存在数值问题
    subroutine modify_mu_record(count_, density_group, mu_group, density_, mu_, warning)
        use parameters
        use parameters2
        implicit none
        integer count_, warning, i
        real(8) density_, mu_, density_group(3), mu_group(3)

        count_=count_+1
        warning=0
        if (count_==1) then
            density_group(1)=density_
            mu_group(1)=mu_
        else
            if (abs(density_-target_density)<abs(density_group(1)-target_density)) then
                density_group(1)=density_
                mu_group(1)=mu_
            else
                warning=1
            endif
            i=3-mod(count_,2)
            density_group(i)=density_
            mu_group(i)=mu_
        endif

    end subroutine

    subroutine modify_mu()
        use parameters
        use parameters2
        implicit none
        integer warning, i1, i2

        if (cur_density<target_density) then
            call modify_mu_record(mu_less_count, density_less, mu_less, cur_density, mu, warning)
        else
            call modify_mu_record(mu_more_count, density_more, mu_more, cur_density, mu, warning)
        endif

        if (density_iter==1) then
            mu = mu - 1.0d-1*sign(1.0d0, (cur_density-target_density)*(mu/cur_density))
            return
        else
            if (mu_less_count/=0 .and. mu_more_count/=0) then
                ! 这里好像是在瞎搞
                ! 如果产生警告, 就用两边最新的值构造一个
                if (warning==0) then
                    mu = (mu_less(1)-mu_more(1))/(density_less(1)-density_more(1)) &
                        *(target_density-density_more(1))+mu_more(1)
                else
                    i1=3-mod(mu_less_count,2)
                    i2=3-mod(mu_more_count,2)
                    mu = (mu_less(i1)-mu_more(i2))/(density_less(i1)-density_more(i2)) &
                        *(target_density-density_more(i1))+mu_more(i2)
                endif
            elseif (mu_less_count==0) then
                mu = (mu_more(2)-mu_more(3))/(density_more(2)-density_more(3)) &
                    *(target_density-density_more(3))+mu_more(3)
            elseif (mu_more_count==0) then
                mu = (mu_less(2)-mu_less(3))/(density_less(2)-density_less(3)) &
                    *(target_density-density_less(3))+mu_less(3)
            endif
        endif
        ! write(stdout, *) mu_less_count, mu_more_count

    end subroutine


    subroutine testBand()
        use constants
        ! use parameters
        use parameters2
        implicit none


        ! 点数36: 1~11~21~36
        integer count_k, i, ik, ix, iy, fileunit
        complex(8) :: fac
        real(8) :: rdotk
        real(8), dimension (2) :: temp
        real(8), dimension (36, 2) :: k_band
        complex(8), dimension (nb, nb, 36) :: h0_k_band
        complex(8), dimension (nb,nb) :: A, B
        complex(8), dimension (nb, 36) :: ev_band
        complex(8), dimension (nb) :: alpha, beta
        complex(8), dimension (nb, nb) :: vl, vr
        complex(8), dimension (nb*2) :: work
        integer info
        real(8), dimension (nb*8) :: rwork


        ! ------------------------------------------------------------------------

        ! 测试能带正确性
        ! 组合一组高对称点
        k_band = 0d0;
        count_k=0;
        do i = 0,10
            count_k = count_k+1;
            k_band(count_k, 1) = 0d0
            k_band(count_k, 2) = i*0.5d0/10;
        enddo
        do i = 1,10
            count_k = count_k+1;
            k_band(count_k, 1) = i*0.5d0/10;
            k_band(count_k, 2) = 0.5d0;
        enddo
        do i = 1,15
            count_k = count_k+1;
            k_band(count_k, 1) = 0.5d0-i*0.5d0/15;
            k_band(count_k, 2) = 0.5d0-i*0.5d0/15;
        enddo

        h0_k = complex_0

        write(stdout,*) 'build k-points of band...'

        ! 反傅里叶变换到k空间的能带
        ! 在每个k点上对角化得到能带特征值
        do ik=1,36
            h0_k_band=complex_0
            do ix = -rx, rx
                do iy = -ry, ry
                    temp=[ix,iy]
                    rdotk = two_pi*dot_product(k_band(ik,:),temp)
                    fac=exp(complex_i*rdotk)
                    h0_k_band(:,:,ik)=h0_k_band(:,:,ik)+fac*h0_r(:,:,ix, iy)
                enddo
            enddo
            A = h0_k_band(:,:,ik)
            ! if (ik==1) then
            !    write(stdout,*) A
            ! endif
            B = complex_0
            do i = 1,nb
                B(i,i)=complex_1
            enddo
            ! write(stdout,*) 'calling cggev...'
            call zggev('N', 'N', nb, A, nb, B, nb, alpha, beta, vl, 1, vr, 1, work, 2*nb, rwork, info)
            ! write(stdout,*) 'finish state ', info
            ev_band(:,ik) = alpha/beta
        enddo

        fileunit = 9999
        open(fileunit, file='testband.dat')
        do ik=1,36
            write(fileunit, *) ik, real(ev_band(:,ik))
        enddo

        close(fileunit)

        return

    end subroutine testBand


end module functions
