module functions_g_chi
    use functions_base
    use parameters
    use parameters2

contains

    subroutine cal_chi_cs(kx,ky,omegaq)
        implicit none
        integer kx,ky,omegaq

        chi_0_=chi_0(:, :, kx, ky, omegaq)

        ! chi_c = chi_0 - chi_0*U_c*chi_c
        Iminuschi_0_ = I_chi + AB(chi_0_,U_c,nb*nb)
        !Iminuschi_0_ = I_chi + AB(U_c,chi_0_,nb*nb)
        chi_c_ = chi_0_
        call inverseAbyB(Iminuschi_0_,chi_c_,nb*nb)

        ! chi_s = chi_0 + chi_0*U_s*chi_s
        Iminuschi_0_ = I_chi - AB(chi_0_,U_s,nb*nb)
        !Iminuschi_0_ = I_chi - AB(U_s,chi_0_,nb*nb)
        chi_s_ = chi_0_
        call inverseAbyB(Iminuschi_0_,chi_s_,nb*nb)

    end subroutine


    !从G计算G_out
    subroutine cal_G_out()
        implicit none

        integer ikx, iky
        integer iomegak,iomegaq
        integer l1,m1,l2,m2

        complex(8), dimension (nb, nb) :: G_, G0_, sigma_

        !call cleanError(G,total_grid)

        call dft(G, r_tau1, nb, nomegaf, dft_grid, 1, 0)
        conjgG=conjg(G)
        call dft(conjgG, r_tau2, nb, nomegaf, dft_grid, 1, 0)

        ! chi_0, 并行
        ! 卷积形式, 改成减法 on tau
        r_tau_sqr=0
        !omp parallel do private(l2,m1,m2)
        do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
            r_tau_sqr(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
                = - r_tau1(l1, m1, :, :, :)*r_tau2(m2, l2, :, :, :)
        enddo; enddo; enddo; enddo
        !omp end parallel do

        ! idft chi_0_r_tau to chi_0
        call dft(r_tau_sqr, chi_0, nb*nb, dft_grid, nomegab, -1, 0)
        chi_0 = T_ev/nk*chi_0

        ! chi_c, chi_s, V
        ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0
        !$omp parallel do private(ikx,iky,Iminuschi_0_,chi_0_,chi_c_,chi_s_) firstprivate(I_chi)
        do iomegaq=minomegab,maxomegab; do ikx=1,nkx; do iky=1,nky;

            !call cal_chi_cs(ikx,iky,iomegaq)
            ! 上面的过程并行会导致问题, 这里直接展开

            chi_0_=chi_0(:, :, ikx, iky, iomegaq)
            ! chi_c = chi_0 - chi_0*U_c*chi_c
            Iminuschi_0_ = I_chi + AB(chi_0_,U_c,nb*nb)
            !Iminuschi_0_ = I_chi + AB(U_c,chi_0_,nb*nb)
            chi_c_ = chi_0_
            call inverseAbyB(Iminuschi_0_,chi_c_,nb*nb)

            ! chi_s = chi_0 + chi_0*U_s*chi_s
            Iminuschi_0_ = I_chi - AB(chi_0_,U_s,nb*nb)
            !Iminuschi_0_ = I_chi - AB(U_s,chi_0_,nb*nb)
            chi_s_ = chi_0_
            call inverseAbyB(Iminuschi_0_,chi_s_,nb*nb)

            V(:, :, ikx, iky, iomegaq) = U_ud - 2*U_uu &
                - ABA(U_ud, chi_0_, nb*nb) &
                + 1.5*ABA(U_s, chi_s_, nb*nb) &
                + 0.5*ABA(U_c, chi_c_, nb*nb)

        enddo; enddo; enddo
        !$omp end parallel do

        !sigma(k) = V(k-k')*G(k')

        ! dft V to V_r_tau
        call dft(V, r_tau_sqr, nb*nb, nomegab, dft_grid, 1, 0)

        ! sigma_r_tau, 并行
        r_tau2 = complex_0
        !omp parallel do private(l2,m1,m2) reduction (+:sigma_r_tau)
        do l1=1,nb; do m1=1,nb;
            do l2=1,nb; do m2=1,nb
                r_tau2(l1, m1, :, :, :) = r_tau2(l1, m1, :, :, :) &
                    + r_tau_sqr(sub_g2chi(l1,l2), sub_g2chi(m1,m2),:,:,:) * r_tau1(l2,m2,:,:,:)
            enddo; enddo;
        enddo; enddo
        !omp end parallel do

        ! idft sigma_r_tau to sigma
        call dft(r_tau2, sigma, nb, dft_grid, nomegaf, -1, 1)
        ! write(*,*) sigma(1,1,1,1,1)

        !call testConvolution3sigma()
        sigma=T_eV/nk*sigma

        !call convergence_test(G_conv)
        !if (G_conv) then
        !    exit
        !endif

        ! 新的G, dyson方程
        ! G=G0+G0*sigma*G, then we have G=(I-G0*sigma)**(-1)*G0
        !$omp parallel do private(ikx,iky,G0_,sigma_,G_)
        do iomegak=minomegaf,maxomegaf;do ikx=1,nkx;do iky=1,nky
            if (iter_method==0)then
                G0_=G0(:,:,ikx,iky,iomegak)
                sigma_=sigma(:,:,ikx,iky,iomegak)
                G_=AB(G0_,sigma_,nb)
                G_=I_G - G_
                call inverseAbyB(G_,G0_,nb)
                G_out(:,:,ikx,iky,iomegak) = G0_
            else
                G(:,:,ikx,iky,iomegak) &
                    = &
                    G0(:,:,ikx,iky,iomegak) &
                    + ABC(G0(:,:,ikx,iky,iomegak), &
                    sigma(:,:,ikx,iky,iomegak), &
                    G_out(:,:,ikx,iky,iomegak),nb)
            endif
        enddo;enddo;enddo
        !$omp end parallel do

    end subroutine

end module
