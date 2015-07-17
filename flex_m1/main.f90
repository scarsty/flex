! 此程序是 flex_m 的简化, 目的是测试单带情况与现有文献对照验证正确性,

program flex_m_single_band
    use Constants
    implicit none

    ! 化学势, 占据数
    real mu, n_density

    ! 目标占据数, 收敛误差
    REAL target_density, density_tol

    ! U, U', J, J' for H1
    REAL h1_U, h1_Up, h1_J, h1_Jp

    ! 温度
    real T, T_beta

    ! 保存的自能函数
    LOGICAL sigma_input, sigma_output
    CHARACTER*128 sigma_input_file, sigma_output_file

    ! 自能函数收敛判据
    REAL sigma_tol
    INTEGER max_it
    REAL alpha
    INTEGER alpha_scheme

    ! 自旋态
    integer spin_state

    ! 数组专用

    ! 格林函数, 反常格林函数
    ! 自能函数, 反常自能函数
    complex, dimension (nk, -nomega:nomega) :: G, F, sigma, delta, G0, sigma0, G1
    complex, dimension (nk, -ntau:ntau) :: G_tau, sigma_tau, delta_tau

    ! 极化率, susceptibilities, effective interactions
    complex, dimension (nk, -nomega:nomega) :: chi_0, chi_s, chi_c, V, V_s
    complex, dimension (nk, -ntau:ntau) :: chi_0_tau, V_tau, V_s_tau

    ! 交换能, 单位矩阵
    complex:: U_s, U_c, U_ud, U_uu, I_chi

    ! H0
    complex, dimension (-2:2, -2:2) :: h0_r
    complex, dimension (nk):: h0_k

    ! k空间对应
    real, dimension (nk, 2) :: k
    ! 两个k的差对应的k
    integer, dimension (nk, nk) :: k_minus, k_plus

    ! 计算g0的辅助
    complex:: I_g0
    real:: u_g0

     ! 循环控制变量较多, 主要是为方便对照文献中公式
    integer ix, iy, iz, count_k, zero_k, ib, ib1, ib2, ik, iomega, ix1, ix2, ix3, ix4, i1, i2, i, iy1, iy2
    integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq,k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
    integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk
    integer l1,m1,l2,m2,l3,m3, n1,l,m
    integer elia1, elia2, info, lda, ldb, ipiv
    real rdotk, temp(2), dis
    complex temp_complex, fac

    ! 迭代sigma次数, density次数
    integer sigma_iter, density_iter, total_iter
    real cur_sigma_tol, cur_density
    real, dimension(2,2):: a, b

    ! 计算chi的辅助
    complex :: chi_0_, chi_c_, chi_s_, Iminuschi_0_

    ! 傅里叶变换, 反傅里叶变换的辅助, f/b means fermi and bose freq
    complex, dimension (nomega*2+1, ntau*2+1) :: dft_f, dft_b
    complex, dimension (ntau*2+1, nomega*2+1) :: idft_f, idft_b

    complex, dimension (nk, nomega*2+1) :: dft_omega
    complex, dimension (nk, ntau*2+1) :: dft_tau

    integer  omega_f, omega_b
    real tau
    !
    !complex, dimension (nb*nb*nk*(nomega*2-1), nb*nb*nk*(nomega*2-1)):: Elishaberg
    ! 变量段结束-------------------------------------------------------------------------------


    ! 代码段-----------------------------------------------------------------------------------
    call readin(T, target_density, density_tol, mu, &
        h1_U, h1_Up, h1_J, h1_Jp, &
        sigma_input, sigma_input_file, sigma_output, &
        sigma_output_file, sigma_tol, max_it, &
        alpha, alpha_scheme, h0_r)

    T_beta = 1d0/kB/T

    !计算k点的坐标
    count_k = 0
    zero_k = 1    ! k原点
    do ix = 1, kx
        do iy = 1, ky
            count_k = count_k + 1
            k(count_k, 1)=-1d0/2+1d0/kx*(ix-1)
            k(count_k, 2)=-1d0/2+1d0/ky*(iy-1)
            if ((abs(k(count_k, 1))<real_error) .and. (abs(k(count_k, 2))<real_error)) then
                zero_k=count_k
            endif
            !write(stdout, *) k(count_k,:)
        enddo
    enddo

    ! k减法矩阵
    k_minus=0
    do i1=1,nk
        do i2=1,nk
            temp = k(i1,:) - k(i2,:)
            !write(stdout, *) temp
            do while (abs(temp(1)+real_error)>0.5)
                temp(1) = temp(1) -sign(1., temp(1))
            enddo
            do while (abs(temp(2)+real_error)>0.5)
                temp(2) = temp(2) -sign(1., temp(2))
            enddo
            !write(stdout, *) temp
            do i = 1,nk
                dis = norm2(temp-k(i,:))
                if (dis<real_error) then
                    k_minus(i1, i2) = i
                    exit
                endif
            enddo
            if (k_minus(i1, i2)==0) then
                write(stdout, *) 'Wrong k_minus at', i1, i2
            endif
        enddo
    enddo

    ! 加法使用k_minus(k1, k_minus(zero_k, k2))
    do i1=1,nk
        do i2=1,nk
            k_plus(i1, i2) = k_minus(i1, k_minus(zero_k, i2))
            if (k_plus(i1, i2)==0) then
                write(stdout, *) 'Wrong k_plus at', i1, i2
            endif
        enddo
    enddo

    ! 反傅里叶变换h0到k空间
    h0_k = complex_0
    do ik=1,nk
        do ix = -2, 2
            do iy = -2, 2
                !temp(1)=ix
                !temp(2)=iy
                temp=[ix,iy]
                rdotk = two_pi*dot_product(k(ik,:),temp)
                fac=exp(complex_i*rdotk)
                h0_k(:,:,ik)=h0_k(:,:,ik)+fac*h0_r(:,:,ix, iy)
            enddo
        enddo
    enddo
    ! 好像没归一化? h0_k=h0_k/?

    ! G0
    ! 费米频率 pi*(2n-1)
    G0=complex_0
    do l1=1,nb; do m1=1,nb; do n1=1,nb
        do ik=1,nk
            do iomegak=-nomega,nomega
                G0(l1,m1,ik, iomega) = 1d0/(complex_i*pi*(2*iomegak-1)/T_beta-(h0_k(n1,n1,ik)-mu)) ! 未完成
            enddo
        enddo
    enddo;enddo;enddo
    G=G0

    ! 频域与时域的傅里叶变换辅助矩阵
    ! 可以一次计算一组, 或者构造大矩阵计算多组
    ! for one k-point, G_tau (行) = G_omega (行) * dft
    ! G_omega (行) = G_tau (行) * idft
    do itau = -ntau,ntau
        do iomega=-nomega,nomega
            omega_f = 2*iomega-1
            omega_b = 2*iomega
            tau = itau*1d0/ntau
            dft_f(iomega, itau) = exp(-2*pi * omega_f*tau)
            dft_b(iomega, itau) = exp(-2*pi * omega_b*tau)
            idft_f(itau, iomega) = exp(2*pi * omega_f*tau)
            idft_b(itau, iomega) = exp(2*pi * omega_b*tau)
        enddo
    enddo
    dft_f = dft_f / (2*nomega+1)
    dft_b = dft_b / (2*nomega+1)


    ! 迭代部分-----------------------------------------------------------

    total_iter = 0
    density_iter = 0
    cur_sigma_tol = 1d0
    cur_density = 1000d0
    do while (abs(target_density-cur_density)>density_tol)

        sigma_iter=0
        do while (cur_sigma_tol>sigma_tol)



            write(stdout, *) 'calculating chi_0...'

            ! dft G to G_tau
            ! 考虑一下是分批还是干脆一批
            do l1=1,nb; do m1=1,nb
                dft_omega = G(l1,m1,:,:)
                call matrixProduct(dft_omega, dft_f, dft_tau, nk, total_tau, total_omega)
                G_tau(l1,m1,:,:) = dft_tau
            enddo; enddo

            ! chi_0, 看起来需要并行
            ! 卷积形式, 改成减法 chi(q)= -G1(q-k)G2(-k), G(-k)=conjg(G(k)) on tau
            chi_0_tau=0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                do ikq=1,nk; do ikk=1,nk
                    k_qminusk = k_minus(ikq, ikk)
                    chi_0_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), ikq, :) &
                        = chi_0_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), ikq, :) &
                        - G_tau(l1, m1, k_qminusk, :)*conjg(G_tau(m2, l2, ikk, :))
                    !write(stdout,*) temp_complex
                enddo;enddo
                !write(stdout,*) l1, l2, m1, m2
            enddo;enddo;enddo;enddo

            ! idft chi_0_tau to chi_0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                dft_tau = chi_0_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:)
                call matrixProduct(dft_tau, idft_b, dft_omega, nk, total_omega, total_tau)
                chi_0(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:) = dft_omega
            enddo;enddo;enddo;enddo

            write(stdout, *) 'calculated chi_0'



            write(stdout, *) 'calculating chi_c, chi_s, V...'

            ! chi_c, chi_s, V, 需要并行和数学库
            ! 含有矩阵乘, 需改写
            do ikq=1,nk;do iomegaq=-nomega,nomega
                ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0

                ! chi_c = chi_0 - chi_0*chi_c
                chi_0_=chi_0(:, :, ikq, iomegaq)
                Iminuschi_0_ = I_chi + AB(chi_0_, U_c)
                chi_c_ = chi_0(:, :, ikq, iomegaq)
                call cgesv(square_nb, square_nb, Iminuschi_0_, square_nb, ipiv, chi_c_, square_nb, info)
                chi_c(:, :, ikq, iomegaq) = chi_c_

                ! chi_s = chi_0 + chi_0*chi_s
                Iminuschi_0_ = I_chi - AB(chi_0_, U_c)
                chi_s_ = chi_0(:, :, ikq, iomegaq)
                call cgesv(square_nb, square_nb, Iminuschi_0_, square_nb, ipiv, chi_s_, square_nb, info)
                chi_s(:, :, ikq, iomegaq) = chi_s_

                !chi_c(:, :, ikq, iomegaq)=inverse(I_chi + chi_0(:, :, ikq, iomegaq)*U_c) * chi_0(:, :, ikq, iomegaq)
                !chi_s(:, :, ikq, iomegaq)=inverse(I_chi - chi_0(:, :, ikq, iomegaq)*U_s) * chi_0(:, :, ikq, iomegaq)
                !V(:, :, ikq, iomegaq) = U_ud - 2*U_uu - U_ud*chi_0(:, :, ikq, iomegaq)*U_ud &
                    !    +1.5*U_s*chi_s(:, :, ikq, iomegaq)*U_s + 0.5*U_c*chi_c(:, :, ikq, iomegaq)*U_c
                V(:, :, ikq, iomegaq) = U_ud - 2*U_uu -ABA(U_ud, chi_0_) &
                    +1.5*ABA(U_s, chi_s_)+0.5*ABA(U_c, chi_c_)

            enddo;enddo

            write(stdout, *) 'calculated chi_c, chi_s, V'



            write(stdout, *) 'calculating sigma...'

            ! dft V to V_tau
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                dft_omega = V(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:)
                call matrixProduct(dft_omega, dft_b, dft_tau, nk, total_tau, total_omega)
                V_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:) = dft_tau
            enddo;enddo;enddo;enddo

            ! sigma_tau, 并行
            sigma_tau=0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                do ikk=1,nk; do ikq=1,nk;
                    k_kminusq = k_minus(ikk, ikq)
                    sigma_tau(l1, m1, ikk, :) = sigma_tau(l1, m1, ikk, :) &
                        + V_tau(sub_g2chi(l1,l2), sub_g2chi(m1,m2), ikq, :) * G_tau(l2,m2,k_kminusq,:)
                    !write(stdout,*) temp_complex
                enddo;enddo
                !write(stdout,*) l1, l2, m1, m2
            enddo;enddo;enddo;enddo

            ! idft sigma_tau to sigma
            do l1=1,nb; do m1=1,nb
                dft_tau = sigma_tau(l1,m1,:,:)
                call matrixProduct(dft_tau, idft_f, dft_omega, nk, total_omega, total_tau)
                sigma(l1,m1,:,:) = dft_omega
            enddo; enddo

            write(stdout, *) 'calculated sigma'



            write(stdout, *) 'calculating New G...'

            ! 新的G, 并行
            G1=G
            G=0
            do l1=1,nb; do m1=1,nb;
                do ikk=1,nk;do iomegak=-nomega,nomega
                    do l2=1,nb; do m2=1,nb;
                        G(l1, m1, ikk, iomegak) = G(l1, m1, ikk, iomegak) &
                            + G0(l1,l2,ikk,iomegak)*sigma(l2,m2,ikk,iomegak)*G1(m2, m1, ikk, iomegak)
                    enddo;enddo;
                enddo;enddo;
            enddo;enddo
            G=G+G0

            write(stdout, *) 'calculated New G'

            write(stdout, *) 'checking convergence of sigma...'

            !第一次迭代, 直接赋值sigma0
            if(sigma_iter==0) then
                sigma0 = sigma
            else
                !计算sigma0与sigma的符合情况, 向量库
                !cur_sigma_tol=norm2()???
            endif
            sigma_iter=sigma_iter+1;
            ! 未收敛处理G?
        enddo ! sigma loop

        ! 计算density
        cur_density=0d0
        ! 这部分应可放在循环外
        do ib=1, nb; do ik=1,nk
            cur_density=cur_density+1/(exp(T_beta*(h0_k(ib,ib,ik)-mu))-1)
        enddo; enddo

        do ib=1, nb; do ik=1,nk; do iomegak=-nomega,nomega
            cur_density=cur_density+G(ib, ib, ik, iomegak)-G0(ib, ib, ik, iomegak)
        enddo; enddo; enddo

        cur_density=cur_density*2
        if (abs(cur_density-target_density)<density_tol) then
            !计算结束
        else
            mu = mu + cur_density-target_density
        endif

    enddo ! density loop

    print *, 'Hello World!'

end program

