subroutine eliashberg()
    ! 变量传递的目的是复用内存
    use constants
    use functions
    use parameters
    use parameters2
    implicit none

    ! 求解厄立希伯格方程以及对应特征向量
    ! 规格化幂法求解LEV
    ! u=v
    ! while not converge do
    !     v=Au
    !     lambda1=abs_max(v)
    !     u=v/lambda1
    ! end

    complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: &
        F, GGdelta, delta, delta0 ! Eliashberg方程所需变量

    integer iomegak,iomegaq
    integer ikx, iky
    integer l1,m1,l2,m2,l3,m3
    real(8) lambda, lambda0, lambda_list(3)
    logical elia_conv
    complex(8) temp_complex

    ! ---------------------------------------------------------------------------------------------


    ! 文献中自旋3态有区别, 需自行推导
    ! 自旋态不是3就是1
    ! 注意: 为节省内存, 复用V表示V_s
    do spin_state=1,3,2
        write(stdout,*)
        write(stdout,*) 'Solving Eliashberg equation for spin ', spin_state
        V = complex_0
        do ikx=1,nkx; do iky=1,nky; do iomegaq=minomegab,maxomegab
            call cal_chi_cs(ikx,iky,iomegaq)
            if (spin_state==3) then
                V(:, :, ikx, iky, iomegaq) &
                    = U_ud - 0.5*ABA(U_s,chi_s_,nb*nb) - 0.5*ABA(U_c, chi_c_,nb*nb)
            else
                V(:, :, ikx, iky, iomegaq) &
                    = U_ud + 1.5*ABA(U_s,chi_s_,nb*nb) - 0.5*ABA(U_c, chi_c_,nb*nb)
            endif
        enddo; enddo; enddo

        ! dft V_s to V_s_r_tau
        call dft(V, r_tau_sqr, nb*nb, nomegab, dft_grid, 1, 0)

        ! 初始值

        delta0 = complex_1 / (nb*nb*nk*nomegaf)
        lambda0=1d0
        elia_conv = .false.


        do while (.not.elia_conv)

            ! 用上一个值计算出新的delta_r_tau

            ! 频域上G*G*delta, 下标l2,m2,l3,m3
            GGdelta = complex_0
            do l2=1,nb; do m2=1,nb
                do l3=1,nb; do m3=1,nb
                    do ikx=1,nkx; do iky=1,nky; do iomegak = minomegaf,maxomegaf
                        GGdelta(l3,m3,ikx,iky,iomegak) &
                            = GGdelta(l3,m3,ikx,iky,iomegak) &
                            + G(l3,l2,ikx,iky,iomegak) &
                            *conjgG(m3,m2,ikx,iky,iomegak) &
                            *delta0(l2,m2,ikx,iky,iomegak)
                    enddo; enddo; enddo
                enddo; enddo
            enddo; enddo

            ! 变换至时域上G*G*delta
            call dft(GGdelta, r_tau1, nb, nomegaf, dft_grid, 1, 0)
            !GGdelta_r_tau = GGdelta_r_tau

            ! 原方程包含负号, 使用减法
            r_tau2 = complex_0
            do l1=1,nb; do m1=1,nb
                do l3=1,nb; do m3=1,nb
                    r_tau2(l1,m1,:,:,:) = r_tau2(l1,m1,:,:,:) &
                        - r_tau_sqr(sub_g2chi(l1,l3),sub_g2chi(m3,m1),:,:,:) &
                        * r_tau1(l3,m3,:,:,:)
                enddo; enddo
            enddo; enddo

            ! 变换回频域
            call dft(r_tau2, delta, nb, dft_grid, nomegaf, -1, 1)
            delta=T_eV/nk*delta
            ! 规格化
            lambda = 0
            do l1=1,nb; do m1=1,nb
                do ikx=1,nkx; do iky=1,nky; do iomegak=minomegaf,maxomegaf
                    lambda = max(lambda, abs(delta(l1,m1,ikx,iky,iomegak)))
                enddo; enddo; enddo
            enddo; enddo
            delta = delta/lambda
            write(stdout, *) lambda
            ! 检测收敛性, 计算lambda
            if (abs(lambda0 - lambda) < 1e-5) then
                elia_conv=.true.
                lambda_list(spin_state)=lambda
            endif

            ! 下一步
            delta0 = delta
            lambda0 = lambda
        enddo

        write(stdout,*) 'Solving ended'
        write(stdout,*)

        ! output delta_nn (gap function)
        write(stdout,*) 'gap function'
        write(stdout,*) '       kx        ky        delta(real and imag)'
        write(stdout,*) '---------------------------------------------------'
        do ikx=1,nkx; do iky=1,nky
            temp_complex=complex_0
            do l1=1,nb
                temp_complex=temp_complex+delta(l1,l1,ikx,iky,1)
            enddo
            write(stdout, '(2F10.4,2F14.8)') k(ikx,iky,:), temp_complex
        enddo; enddo

        write(stdout,*)
        write(stdout,*) 'The maximum eigenvalue is ', lambda, 'for spin', spin_state
        write(stdout,*)
    enddo

    write(stdout,*)
    write(stdout,*) 'Summary:'
    write(stdout,*)
    write(stdout,*) 'Temperature is ', T
    write(stdout,*)
    write(stdout,*) '  spin      max.eigenvalue'
    write(stdout,*) '------------------------------'

    do spin_state=1,3,2
        write(stdout,'(I7, F20.8)') spin_state, lambda_list(spin_state)
    enddo

end subroutine eliashberg
