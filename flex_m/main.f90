program flex_m2d
    use constants
    use myfunctions
    use parameters
    use parameters2
    implicit none


    ! 循环控制变量较多, 主要是为方便对照文献中公式
    integer ix, iy, iz, count_k, zero_k, ib, ib1, ib2, ik, iq, iomega, ix1, ix2, ix3, ix4, i1, i2, i, iy1, iy2
    integer ikx, iky, irx, iry
    integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq, k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
    integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk
    integer l1,m1,l2,m2,l3,m3,n1,l,m
    integer elia1, elia2, info, lda, ldb, ipiv
    real(8) rdotk, temp(2), dis, mu0, deltamu_per_density, density0
    complex(8) temp_complex, fac

    ! 迭代sigma次数, density次数
    integer sigma_iter, density_iter, total_iter
    real(8) cur_sigma_tol, cur_density, density_base
    logical sigma_conv, density_conv
    real(8), dimension(2,2):: a, b

    integer  omega_f, omega_b
    real(8) tau

    real(8) dznrm2
    external dznrm2

    ! 变量段结束-------------------------------------------------------------------------------

    call readin()

    if (test_band) then
        call testband()
    endif

    T_beta = 1d0/kB/T

    !计算k点的坐标
    write(stdout, *) "Building k-points grid..."
    zero_k = 1    ! k原点
    do ikx = 1, nkx
        do iky = 1, nky
            k(ikx, iky, 1)=-1d0/2+1d0/nkx*(ikx-1)
            k(ikx, iky, 2)=-1d0/2+1d0/nky*(iky-1)
            ! write(stdout, *) k(ikx,iky,:)
        enddo
    enddo

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
    U_c = U_ud+u_uu

    ! write (stdout, *) U_s, U_c, U_ud, U_uu

    ! 反傅里叶变换h0到k空间
    h0_k = complex_0
    do ikx=1,nkx; do iky=1,nky
        do irx=-rx,rx; do iry=-ry,ry
            temp=[irx,iry]
            rdotk = two_pi*dot_product(k(ikx,iky,:),temp)
            fac=exp(complex_i*rdotk)
            h0_k(:,:,ikx,iky)=h0_k(:,:,ikx,iky)+fac*h0_r(:,:,irx,iry)
        enddo; enddo
    enddo; enddo
    ! 好像没归一化? a: seems it is ok

    if (nb==1) then
        call build_h0_k()
        T_beta = 0.25
    endif


   ! I_chi
    I_chi=complex_0
    do i=1,nb*nb
        I_chi(i,i)=complex_1
    enddo



    ! 迭代部分-----------------------------------------------------------
    write(stdout, *) "Temperature in K = ", 1d0/kB/T_beta
    write(stdout, *)

    write(stdout, *) "Begin to calculate FLEX"
    total_iter = 0
    density_iter = 0
    cur_sigma_tol = 1d0
    cur_density = 1000d0

    deltamu_per_density=1

    density_conv = .false.
    do while (.not. density_conv)

        sigma_conv=.false.
        sigma_iter=0

        ! G0
        ! 费米频率
        G0=complex_0
        do l1=1,nb; do m1=1,nb; !do n1=1,nb
            do ikx=1,nkx; do iky=1,nky
                do iomegak=-maxomegaf,maxomegaf,2
                    !G0(l1,m1,ikx,iky,iomegak) = &
                        !1d0/(complex_i*pi*(2*iomegak-1)/T_beta-(h0_k(n1,n1,ikx,iky)-mu)) ! 未完成
                    G0(l1,m1,ikx,iky,transfer_freq(iomegak)) = &
                        T_beta / (complex_i*pi*iomegak - (h0_k(l1,m1,ikx,iky)-mu))
                enddo
            enddo; enddo
        enddo; enddo; !enddo
        G=G0
        conjgG=conjg(G)

call testConvolution3()

        ! base density
        density_base = 0d0
        do ib=1,nb; do ikx=1,nkx; do iky=1,nky
            ! write(stdout, *) T_beta*(h0_k(ib,ib,ikx,iky)-mu)
            density_base=density_base+1/(exp(T_beta*(h0_k(ib,ib,ikx,iky)-mu))+1)
        enddo; enddo; enddo
        density_base=density_base*2/nk
        write(stdout, *) 'base density is ', density_base

        sigma_iter =0
        do while (.not. sigma_conv)

            ! write(stdout, *) 'calculating chi_0...'

            ! calculate chi_0 with chi(q)= -G1(q-k)G2(-k),

            ! dft G to G_r_tau
            call dft(G, G_r_tau, nb, 1, 0)
            call dft(conjgG, conjgG_r_tau, nb, 1, 0)

            ! chi_0, 看起来需要并行
            ! 卷积形式, 改成减法  on tau
            chi_0_r_tau=0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                chi_0_r_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
                    = - G_r_tau(l1, m1, :, :, :)*conjgG_r_tau(m2, l2, :, :, :)
            enddo; enddo; enddo; enddo

            ! idft chi_0_r_tau to chi_0
            call dft(chi_0_r_tau, chi_0, nb*nb, -1, 2)


            ! write(stdout, *) 'calculating chi_c, chi_s, V...'

            ! chi_c, chi_s, V, 需要并行和数学库
            ! 含有矩阵乘, 需改写
            do ikx=1,nkx; do iky=1,nky; do iomegaq=-maxomegab,maxomegab,2
                ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0

                ! chi_c = chi_0 - chi_0*chi_c
                chi_0_=chi_0(:, :, ikx, iky, transfer_freq(iomegaq))

                Iminuschi_0_ = I_chi + AB(chi_0_, U_c)

                chi_c_ = chi_0(:, :, ikx, iky, transfer_freq(iomegaq))
                call zgesv(square_nb, square_nb, Iminuschi_0_, square_nb, ipiv, chi_c_, square_nb, info)
                chi_c(:, :, ikx, iky, transfer_freq(iomegaq)) = chi_c_

                ! chi_s = chi_0 + chi_0*chi_s
                Iminuschi_0_ = I_chi - AB(chi_0_, U_s)
                chi_s_ = chi_0(:, :, ikx, iky, transfer_freq(iomegaq))
                call zgesv(square_nb, square_nb, Iminuschi_0_, square_nb, ipiv, chi_s_, square_nb, info)
                chi_s(:, :, ikx, iky, transfer_freq(iomegaq)) = chi_s_

                V(:, :, ikx, iky, transfer_freq(iomegaq)) = U_ud - 2*U_uu - ABA(U_ud, chi_0_) &
                    + 1.5*ABA(U_s, chi_s_) + 0.5*ABA(U_c, chi_c_)

            enddo; enddo; enddo



            ! write(stdout, *) 'calculating sigma...'

            ! dft V to V_r_tau
            call dft(V, V_r_tau, nb*nb, 1, 0)

            ! sigma_r_tau, 并行
            sigma_r_tau = complex_0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                sigma_r_tau(l1, m1, :, :, :) = sigma_r_tau(l1, m1, :, :, :) &
                    + V_r_tau(sub_g2chi(l1,l2), sub_g2chi(m1,m2),:,:,:) * G_r_tau(l2,m2,:,:,:)
            enddo; enddo; enddo; enddo

            ! idft sigma_r_tau to sigma
            call dft(sigma_r_tau, sigma, nb, -1, 1)

            ! write(stdout, *) 'checking convergence of sigma...'

            ! 第一次迭代, 直接赋值sigma0, 不测试收敛
            if (sigma_iter > 0) then
                ! 计算sigma0与sigma的符合情况, 向量库
                ! scnrm2: 欧几里得模，行向量乘以自身转置共轭
                sigma_minus = sigma0 - sigma
                cur_sigma_tol = dznrm2(nb*nb*nkx*nky*totalnomega, sigma_minus, 1) &
                    / dznrm2(nb*nb*nkx*nky*totalnomega, sigma, 1)
                write(stdout,*) 'sigma tolerance is ', cur_sigma_tol, '/', sigma_tol
                if (cur_sigma_tol>1) then
                    write(stdout,*) "sigma seems bad, please reset mu."
                    !stop
                endif
                if (cur_sigma_tol < sigma_tol) then
                    sigma_conv = .true.
                endif
            endif
            sigma0 = sigma

            ! write(stdout, *) 'calculating New G...'

            ! 新的G, dyson方程
            G1=G
            G=G0
            do l1=1,nb; do m1=1,nb;
                do ikx=1,nkx; do iky=1,nky; do iomegak=-maxomegaf,maxomegaf,2
                    do l2=1,nb; do m2=1,nb;
                        G(l1, m1, ikx, iky, transfer_freq(iomegak)) &
                            = G(l1, m1, ikx, iky, transfer_freq(iomegak)) &
                            + G0(l1, l2, ikx, iky, transfer_freq(iomegak)) &
                            * sigma(l2, m2, ikx, iky, transfer_freq(iomegak)) &
                            * G1(m2, m1, ikx, iky, transfer_freq(iomegak))
                    enddo;enddo
                enddo;enddo;enddo
            enddo;enddo
            conjgG=conjg(G)


            sigma_iter=sigma_iter+1;
            total_iter = total_iter + 1

            if (total_iter>20) then
                !write(stdout,*) sigma_minus
                !exit
            endif
        enddo

        ! sigma loop end

        ! 计算density
        density0 = cur_density
        cur_density=0d0

        do ib=1,nb; do ikx=1,nkx; do iky=1,nky; do iomegak=-maxomegaf,maxomegaf,2
            cur_density = cur_density &
                + G(ib, ib, ikx, iky, transfer_freq(iomegak)) &
                - G0(ib, ib, ikx, iky, transfer_freq(iomegak))
        enddo; enddo; enddo; enddo

        cur_density=cur_density*2/nk + density_base

        if (density_iter>0) then
            deltamu_per_density = (mu-mu0)/(cur_density-density0)
        endif

        write(stdout,*) 'density and mu are ', cur_density, '/', target_density, ', ', mu

        if (abs(cur_density-target_density)<density_tol) then
            density_conv=.true.
            !计算结束
        else
            ! 计算化学势变化与占据数变化的比值来调整新的化学势
            mu0 = mu
            if (density_iter>0) then
                !mu = mu - (cur_density-target_density)
                mu = mu - (cur_density-target_density)*deltamu_per_density
            else
                mu = mu - sign(1.0d0, (cur_density-target_density)*deltamu_per_density)
            endif
            write(stdout,*) 'modified new mu = ', mu
        endif
        density_iter = density_iter + 1

        if (total_iter>20) then
            !exit
        endif

    enddo
    ! density loop end

    ! 迭代部分结束--------------------------------------------------------------

    ! output chi_s(q,0)
    write(stdout,*) 'chi_s at omega = 0'

    do l1=1,nb; do l2=1,nb
        temp_complex=complex_0
        do ikx=1,nkx; do iky=1,nky
            temp_complex=temp_complex+chi_s(sub_g2chi(l1,l1),sub_g2chi(l2,l2),ikx,iky,0)
        enddo; enddo
        write(stdout,*) k(ikx,iky,:), temp_complex
    enddo; enddo



    if (solve_eliashberg) then
        call eliashberg()
    endif


    print *, 'good night.'
    ! return

end program

