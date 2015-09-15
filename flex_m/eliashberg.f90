subroutine eliashberg()
    use constants
    use myfunctions
    use parameters
    use parameters2
    implicit none

    ! δ���
    ! ����ϣ���񷽳�, ֱ��Ӧ������õ������

    ! ����̬����3����1
    ! �������, ���д

    integer elia1, elia2, sub_g2chi1, sub_g2chi2, sub_g2chi3, maxindex
    integer ikk,ikq,iomegak,iomegaq, k_kplusq, omega_kplusq,k_kminusq, omega_kminusq, itau, k_0minusq, k_qminusk
    integer ikk1, ikk2, iomegak1, iomegak2, k_kminusk, omega_kminusk, ikx, iky
    integer l1,m1,l2,m2,l3,m3
    real(8) lambda, lambda0
    logical elia_conv
    ! complex, dimension (nb*nb*nk*(nomega*2-1), nb*nb*nk*(nomega*2-1)):: Eliashberg

    ! ---------------------------------------------------------------------------------------------

    write(stdout,*)
    write(stdout,*) 'Solving Eliashberg equation...'
    ! ����������3̬������, �������Ƶ�
    V_s = complex_0
    do ikx=1,nkx; do iky=1,nky; do iomegaq=-maxomegab,maxomegab,2
        chi_c_ = chi_c(:, :, ikx, iky, transfer_freq(iomegaq))
        chi_s_ = chi_s(:, :, ikx, iky, transfer_freq(iomegaq))
        if (spin_state==3) then
            V_s(:, :, ikx, iky, transfer_freq(iomegaq)) &
            = U_ud - 0.5*ABA(U_s,chi_s_) - 0.5*ABA(U_c, chi_c_)
        else
            V_s(:, :, ikx, iky, transfer_freq(iomegaq)) &
            = U_ud + 1.5*ABA(U_s,chi_s_) - 0.5*ABA(U_c, chi_c_)
        endif
    enddo; enddo; enddo

    ! dft G to G_r_tau (seems no use)
    ! do l1=1,nb; do m1=1,nb
    !     dft_omega = G(l1,m1,:,:)
    !     call matrixProduct(dft_omega, dft_f, dft_r_tau, nk, total_r_tau, total_omega)
    !     G_r_tau(l1,m1,:,:) = dft_r_tau
    ! enddo; enddo

    ! dft V_s to V_s_r_tau
    call dft(V_s, V_s_r_tau, nb*nb, 1,0)

    ! ����ݷ����LEV
    ! u=v
    ! while not converge do
    !     v=Au
    !     lambda1=abs_max(v)
    !     u=v/lambda1
    ! end


    ! ��ʼֵ

    delta0 = complex_1 / (nb*nb*nk*totalnomega)
    lambda0=1d0
    elia_conv = .false.


    do while (.not.elia_conv)

        ! ����һ��ֵ������µ�delta_r_tau

        ! Ƶ����G*G*delta, �±�l2,m2,l3,m3
        GGdelta = complex_0
        do l2=1,nb; do m2=1,nb
            sub_g2chi2 = sub_g2chi(l2,m2)
            do l3=1,nb; do m3=1,nb
                sub_g2chi3 = sub_g2chi(l3,m3)
                do ikx=1,nkx; do iky=1,nky; do iomegak = -maxomegaf,maxomegaf,2
                    GGdelta(sub_g2chi2,sub_g2chi3,ikx,iky,transfer_freq(iomegak)) &
                        = GGdelta(sub_g2chi2,sub_g2chi3,ikx,iky,transfer_freq(iomegak)) &
                        + &
                        G(l3,l2,ikx,iky,transfer_freq(iomegak)) &
                        *conjgG(m3,m2,ikx,iky,transfer_freq(iomegak)) &
                        *delta0(l2,m2,ikx,iky,transfer_freq(iomegak))
                enddo; enddo; enddo
            enddo; enddo
        enddo; enddo

        ! �任��ʱ����G*G*delta
        call dft(GGdelta, GGdelta_r_tau, nb*nb, 1, 0)

        ! ԭ���̰�������, ʹ�ü���
        delta_r_tau = complex_0
        do l1=1,nb; do m1=1,nb
            do l3=1,nb; do m3=1,nb
                delta_r_tau(l1,m1,:,:,:) = delta_r_tau(l1,m1,:,:,:) &
                    - &
                    V_s_r_tau(sub_g2chi(l1,l3),sub_g2chi(m3,m1),:,:,:) &
                    * GGdelta_r_tau(l3,m3,:,:,:)
            enddo; enddo
        enddo; enddo

        ! �任��Ƶ��
        call dft(delta_r_tau, delta, nb, -1, 1)

        ! ���
        lambda = 0
        do l1=1,nb; do m1=1,nb
            do ikx=1,nkx; do iky=1,nky; do iomegak=-maxomegaf,maxomegaf,2
                lambda = max(lambda, abs(delta(l1,m1,ikx,iky,transfer_freq(iomegak))))
            enddo; enddo; enddo
        enddo; enddo
        delta = delta/lambda
    write(stdout, *) lambda
        ! ���������, ����lambda
        if (abs(lambda0 - lambda) < 1e-5) then
            elia_conv=.true.
        endif

        ! ��һ��
        delta0 = delta
        lambda0 = lambda
    enddo

    write(stdout,*) 'Solving ended, the eigenvalue is ', lambda

    ! output delta_nn (gap function)


end subroutine eliashberg
