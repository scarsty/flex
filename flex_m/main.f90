program flex_m
    use constants
    use myfunctions
    use parameters
    use parameters2
    implicit none


    ! 循环控制变量较多, 主要是为方便对照文献中公式
    integer ix, iy, iz, count_k, zero_k, ib, ib1, ib2, ik, iq, iomega, ix1, ix2, ix3, ix4, i1, i2, i, iy1, iy2
    integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq, k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
    integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk
    integer l1,m1,l2,m2,l3,m3,n1,l,m
    integer elia1, elia2, info, lda, ldb, ipiv
    real rdotk, temp(2), dis
    complex temp_complex, fac

    ! 迭代sigma次数, density次数
    integer sigma_iter, density_iter, total_iter
    real cur_sigma_tol, cur_density
    logical sigma_conv, density_conv
    real, dimension(2,2):: a, b

    integer  omega_f, omega_b
    real tau


    ! 变量段结束-------------------------------------------------------------------------------

    !call readin()

    if (test_band) then
        call testband()
    endif

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
    ! write(stdout, *) zero_k

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
            if (k_minus(i1, i2)<=0 .or. k_minus(i1, i2)>nk) then
                write(stdout, *) 'Wrong k_minus at', i1, i2
            endif
        enddo
    enddo

    ! k加法矩阵, k_minus(k1, k_minus(zero_k, k2))
    do i1=1,nk
        do i2=1,nk
            k_plus(i1, i2) = k_minus(i1, k_minus(zero_k, i2))
            if (k_plus(i1, i2)==0 .or. k_plus(i1, i2)>nk) then
                write(stdout, *) 'Wrong k_plus at', i1, i2
            endif
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

    ! 反傅里叶变换h0到k空间
    h0_k = complex_0
    do ik=1,nk
        do ix = -rx, rx
            do iy = -ry, ry
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
    enddo; enddo; enddo
    G=G0

    ! I_chi
    I_chi=complex_0
    do i=1,nb*nb
        I_chi(i,i)=complex_1
    enddo

    call buildDFTMatrix()



    ! 迭代部分-----------------------------------------------------------

    total_iter = 0
    density_iter = 0
    cur_sigma_tol = 1d0
    cur_density = 1000d0


    density_conv = .false.

    do while (.not. density_conv)

        sigma_conv=.false.
        sigma_iter=0
        do while (.not. sigma_conv)



            write(stdout, *) 'calculating chi_0...'

            ! dft G to G_tau
            ! 考虑一下是分批还是干脆一批
            call dft(G, G_tau, nb, 1)
            !            do l1=1,nb; do m1=1,nb
            !                dft_omega = G(l1,m1,:,:)
            !                call matrixProduct(dft_omega, dft_f, dft_tau, nk, total_tau, total_omega)
            !                G_tau(l1,m1,:,:) = dft_tau
            !            enddo; enddo

            ! chi_0, 看起来需要并行
            ! 卷积形式, 改成减法 chi(q)= -G1(q-k)G2(-k), G(-k)=conjg(G(k)) on tau
            chi_0_tau=0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                do ikq=1,nk; do ikk=1,nk
                    k_qminusk = k_minus(ikq, ikk)
                    chi_0_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), ikq, :) &
                        = chi_0_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), ikq, :) &
                        - G_tau(l1, m1, k_qminusk, :)*conjg(G_tau(m2, l2, ikk, :))
                enddo;enddo
            enddo;enddo;enddo;enddo

            ! idft chi_0_tau to chi_0
            call idft(chi_0_tau, chi_0, nb*nb, 0)

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

                ! chi_c(:, :, ikq, iomegaq)=inverse(I_chi + chi_0(:, :, ikq, iomegaq)*U_c) * chi_0(:, :, ikq, iomegaq)
                ! chi_s(:, :, ikq, iomegaq)=inverse(I_chi - chi_0(:, :, ikq, iomegaq)*U_s) * chi_0(:, :, ikq, iomegaq)
                ! V(:, :, ikq, iomegaq) = U_ud - 2*U_uu - U_ud*chi_0(:, :, ikq, iomegaq)*U_ud &
                    ! +1.5*U_s*chi_s(:, :, ikq, iomegaq)*U_s + 0.5*U_c*chi_c(:, :, ikq, iomegaq)*U_c
                V(:, :, ikq, iomegaq) = U_ud - 2*U_uu - ABA(U_ud, chi_0_) &
                    + 1.5*ABA(U_s, chi_s_) + 0.5*ABA(U_c, chi_c_)

            enddo;enddo

            write(stdout, *) 'calculated chi_c, chi_s, V'



            write(stdout, *) 'calculating sigma...'

            ! dft V to V_tau
            call dft(V, V_tau, nb*nb, 0)

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
            call idft(sigma_tau, sigma, nb, 1)

            write(stdout, *) 'calculated sigma'



            write(stdout, *) 'calculating New G...'

            ! 新的G, dython方程
            G1=G
            G=G0
            do l1=1,nb; do m1=1,nb;
                do ikk=1,nk;do iomegak=-nomega,nomega
                    do l2=1,nb; do m2=1,nb;
                        G(l1, m1, ikk, iomegak) = G(l1, m1, ikk, iomegak) &
                            + G0(l1,l2,ikk,iomegak)*sigma(l2,m2,ikk,iomegak)*G1(m2, m1, ikk, iomegak)
                    enddo;enddo;
                enddo;enddo;
            enddo;enddo

            write(stdout, *) 'calculated New G'


            write(stdout, *) 'checking convergence of sigma...'

            ! 第一次迭代, 直接赋值sigma0
            if(sigma_iter==0) then
                sigma0 = sigma
            else
                ! 计算sigma0与sigma的符合情况, 向量库
                ! scnrm2: 欧几里得模，行向量乘以自身转置共轭
                sigma_minus = sigma0 - sigma
                cur_sigma_tol=scnrm2(nb*nb*nk*total_omega, sigma_minus, 1)
                if (cur_sigma_tol<sigma_tol) then
                    sigma_conv = .true.
                endif
            endif
            sigma_iter=sigma_iter+1;
            ! 未收敛处理G?
        enddo ! sigma loop

        ! 计算density
        cur_density=0d0

        ! 这部分应可放在循环外
        do ib=1,nb; do ik=1,nk
            cur_density=cur_density+1/(exp(T_beta*(h0_k(ib,ib,ik)-mu))-1)
        enddo; enddo

        do ib=1,nb; do ik=1,nk; do iomegak=-nomega,nomega
            cur_density=cur_density+G(ib, ib, ik, iomegak)-G0(ib, ib, ik, iomegak)
        enddo; enddo; enddo

        cur_density=cur_density*2
        write(stdout,*) 'density ', cur_density, '/', target_density

        if (abs(cur_density-target_density)<density_tol) then
            density_conv=.true.
            !计算结束
        else
            mu = mu + cur_density-target_density
        endif

    enddo ! density loop

    ! 迭代部分结束--------------------------------------------------------------

    ! output chi_s(q,0)
    write(stdout,*) 'chi_s at omega = 0'

    do l1=1,nb; do l2=1,nb
        temp_complex=complex_0
        do iq=1,nk
            temp_complex=temp_complex+chi_s(sub_g2chi(l1,l1),sub_g2chi(l2,l2),iq,0)
        enddo
        write(stdout,*) k(iq,:), temp_complex
    enddo; enddo



    if (solve_eliashberg) then
        call eliashberg()
    endif


    print *, 'good night.'
    return

end program

