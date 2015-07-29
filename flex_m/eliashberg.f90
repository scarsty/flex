subroutine eliashberg()
    use constants
    use myfunctions
    use parameters
    use parameters2
    implicit none

    ! 未完成
    ! 厄立希伯格方程, 直接应用上面得到的组合

    ! 自旋态不是3就是1
    ! 含矩阵乘, 需改写

    integer elia1, elia2, sub_g2chi1, sub_g2chi2, sub_g2chi3, maxindex
    integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq,k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
    integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk
    integer l1,m1,l2,m2,l3,m3
    real lambda, lambda0
    logical elia_conv
    ! complex, dimension (nb*nb*nk*(nomega*2-1), nb*nb*nk*(nomega*2-1)):: Elishaberg

    ! ---------------------------------------------------------------------------------------------

    write(stdout,*) 'Solving Eliashberg equation...'
    ! 文献中自旋3态有区别, 需自行推导
    do ikq=1,nk; do iomegaq=-nomega,nomega
        chi_c_ = chi_0(:, :, ikq, iomegaq)
        chi_s_ = chi_0(:, :, ikq, iomegaq)
        if (spin_state==3) then
            V_s(:, :, ikq, iomegaq) = U_ud - 0.5*ABA(U_s,chi_s_) &
                - 0.5*ABA(U_c, chi_c_)
        else
            V_s(:, :, ikq, iomegaq) = U_ud + 1.5*ABA(U_s,chi_s_) &
                -0.5*ABA(U_c, chi_c_)
        endif
    enddo; enddo

    ! dft G to G_tau (seems no use)
    ! do l1=1,nb; do m1=1,nb
    !     dft_omega = G(l1,m1,:,:)
    !     call matrixProduct(dft_omega, dft_f, dft_tau, nk, total_tau, total_omega)
    !     G_tau(l1,m1,:,:) = dft_tau
    ! enddo; enddo

    ! dft V_s to V_s_tau
    do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
        dft_omega = V_s(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:)
        call matrixProduct(dft_omega, dft_b, dft_tau, nk, total_tau, total_omega)
        V_s_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2),:,:) = dft_tau
    enddo; enddo; enddo; enddo


    ! 规格化幂法求解LEV
    ! u=v
    ! while not converge do
    !     v=Au
    !     lambda1=abs_max(v)
    !     u=v/lambda1
    ! end


    ! 初始值

    delta0 = complex_1 / (nb*nb*nk*(2*nomega+1))
    lambda0=1d0
    elia_conv = .false.


    do while (.not.elia_conv)

        ! 用上一个值计算出新的delta_tau

        ! 频域上G*G*delta, 下标l2,m2,l3,m3
        GGdelta = complex_0
        do l2=1,nb; do m2=1,nb
            sub_g2chi2 = sub_g2chi(l2,m2)
            do l3=1,nb; do m3=1,nb
                sub_g2chi3 = sub_g2chi(l3,m3)
                do ikk=1,nk; do iomegak = -nomega,nomega
                    GGdelta(sub_g2chi2,sub_g2chi3,nk,iomegak) = GGdelta(sub_g2chi2,sub_g2chi3,nk,iomegak) &
                        + G(l3,l2,ikk,iomegak)*conjg(G(m3,m2,ikk,iomegak))*delta0(l2,m2,ikk,iomegak)
                enddo; enddo
            enddo; enddo
        enddo; enddo

        ! 变换至时域上G*G*delta
        call dft(GGdelta, GGdelta_tau, nb*nb, 0)

        ! 原方程包含负号, 使用减法
        delta_tau = complex_0
        do l1=1,nb; do m1=1,nb
            do ikk1=1,nk; do ikk2=1,nk;
                k_kminusk = k_minus(ikk1, ikk2)
                do l2=1,nb; do m2=1,nb
                    do l3=1,nb; do m3=1,nb
                        do itau=-ntau,ntau

                            delta_tau(l1, m1, ikk1, itau) = delta_tau(l1, m1, ikk1, itau) &
                                - &
                                V_s_tau(sub_g2chi(l1,l3), sub_g2chi(m3,m1), k_kminusk, itau) &
                                * GGdelta_tau(sub_g2chi(l2,l3), sub_g2chi(m2,m3), ikk2, itau)

                            !Elishaberg(elia1, elia2) = Elishaberg(elia1, elia2) &
                                !   - &
                                !   V_s(sub_g2chi(l1,l3), sub_g2chi(m3,m1), k_kminusk, omega_kminusk) &
                                !   *G(l3,l2,ikk2,iomegak2)*conjg(G(m3,m2,ikk2,iomegak2))

                        enddo
                    enddo; enddo
                enddo; enddo
            enddo; enddo
        enddo; enddo

        ! 变换回频域
        call idft(delta_tau, delta, nb, 1)

        ! 规格化
        lambda = 0
        do l1=1,nb; do m1=1,nb
            do ikk=1,nk; do iomegak=-nomega,nomega
                lambda = max(lambda, abs(delta(l1,m1,ikk,iomegak)))
            enddo; enddo
        enddo; enddo
        delta = delta/lambda

        ! 检测收敛性, 计算lambda
        if (abs(lambda0 - lambda) < 1e-5) then
            elia_conv=.true.
        endif

        ! 下一步
        delta0 = delta
        lambda0 = lambda
    enddo

    write(stdout,*) 'Solving ended, the eigenvalue is ', lambda

    ! output delta_nn (gap function)


    ! 求特征值和特征向量, 调用数学库, (取消)
    !call ()

end subroutine eliashberg
