program flex_2d

    !   Compute Green function for superconductive properties
    !   author  : Sun TY (scarsty@gmail.com)
    !   status  : constructing
    !   version : none

    use constants
    use functions
    use parameters
    use parameters2
    implicit none


    ! 循环控制变量较多, 主要是为方便对照文献中公式
    integer ib
    integer ikx, iky
    integer iomegak,iomegaq
    integer l1,m1,l2,m2
    !integer mpiinfo
    complex(8) temp_complex

    ! 变量段结束-------------------------------------------------------------------------------

    mpi_info = mpi_init1()
    mpi_rank = mpi_rank1()
    mpi_size = mpi_size1()

    call readin()

    if (test_band) then
        call testband()
    endif

    call init()

    call init_Kpoints()

    call init_U()

    if (nb==1) then
        call build_h0_k()
        !T_beta = 0.25
    endif


    ! 迭代部分---------------------------------------------------------------------------------
    write(stdout, *) "Temperature in K = ", T
    write(stdout, *) "Temperature in eV = ", T_eV
    write(stdout, *) "Temperature in beta = ", T_beta
    write(stdout, *)

    write(stdout, *) "Begin to calculate FLEX"
    total_iter = 1
    density_iter = 1
    cur_density = 1000d0

    density_conv = .false.

    do while (.not. density_conv)

        sigma_conv=.false.
        sigma_iter=0

        ! G0
        ! 费米频率
        G0=complex_0

        do ikx=1,nkx; do iky=1,nky
            do iomegak=minomegaf,maxomegaf
                diag_h0_G0_=complex_0
                do ib=1,nb
                    diag_h0_G0_(ib,ib) = &
                        1/(complex_i*(2*iomegak-1)*pi*T_eV-(ev_h0_k(ib,ikx,iky)-mu))
                    !write(stdout,*)diag_h0_G0_(ib,ib)
                enddo
                u_h0_k_=u_h0_k(:,:,ikx,iky)
                diag_h0_G0_=ABAH(u_h0_k_,diag_h0_G0_,nb)
                G0(:,:,ikx,iky,iomegak)=diag_h0_G0_
                !G0(l1,m1,ikx,iky,transfer_freq(iomegak)) = &
                    !T_beta / (complex_i*pi*iomegak - (h0_k(l1,m1,ikx,iky)-mu))
            enddo
        enddo; enddo
        G=G0

        call mixerInit()

        !call testConvolution()
        !call testConvolution3()
        !call testConvolution3G()

        ! base density
        density_base = 0d0
        do ib=1,nb; do ikx=1,nkx; do iky=1,nky
            !write(stdout, *) T_beta*(h0_k(ib,ib,ikx,iky)-mu)
            density_base=density_base+1/(exp(T_beta*(real(h0_k(ib,ib,ikx,iky))-mu))+1)
        enddo; enddo; enddo
        density_base=density_base*2/nk
        write(stdout, *) 'base density is ', density_base

        sigma_iter = 1
        write(stdout,*) '  iter  conv.pts           sigma.tol'
        write(stdout,*) '--------------------------------------'

        ! sigma迭代中使用openmp并行

        do while (.not. sigma_conv)

            ! calculate chi_0 with chi(q)= -G1(q-k)G2(-k), the same to -G1(q-k)G2(k)**H
            ! dft G to G_r_tau
            call dft(G, G_r_tau, nb, nomegaf, dft_grid, 1, 0)

            conjgG=conjg(G)
            call dft(conjgG, conjgG_r_tau, nb, nomegaf, dft_grid, 1, 0)

            ! chi_0, 并行
            ! 卷积形式, 改成减法 on tau
            chi_0_r_tau=0
            !$omp parallel do private(l2,m1,m2)
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                chi_0_r_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
                    = - G_r_tau(l1, m1, :, :, :)*conjgG_r_tau(m2, l2, :, :, :)
            enddo; enddo; enddo; enddo
            !$omp end parallel do

            ! idft chi_0_r_tau to chi_0
            call dft(chi_0_r_tau, chi_0, nb*nb, dft_grid, nomegab, -1, 0)
            chi_0 = T_ev/nk*chi_0

            ! chi_c, chi_s, V
            ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0
            !$omp parallel do private(iky,iomegaq,Iminuschi_0_,chi_0_,chi_c_,chi_s_)
            do ikx=1,nkx; do iky=1,nky; do iomegaq=minomegab,maxomegab

                ! chi_c = chi_0 - chi_0*U_c*chi_c
                Iminuschi_0_ = I_chi + AB(chi_0_,U_c,nb*nb)
                chi_0_=chi_0(:, :, ikx, iky, iomegaq)
                !Iminuschi_0_ = I_chi + AB(U_c,chi_0_,nb*nb)
                chi_c_ = chi_0(:, :, ikx, iky, iomegaq)
                chi_c(:, :, ikx, iky, iomegaq)= inverseAbyB(Iminuschi_0_,chi_c_,nb*nb)

                ! chi_s = chi_0 + chi_0*U_s*chi_s
                Iminuschi_0_ = I_chi - AB(chi_0_,U_s,nb*nb)
                !Iminuschi_0_ = I_chi - AB(U_s,chi_0_,nb*nb)
                chi_s_ = chi_0(:, :, ikx, iky, iomegaq)
                chi_s(:, :, ikx, iky, iomegaq) = inverseAbyB(Iminuschi_0_,chi_s_,nb*nb)

                V(:, :, ikx, iky, iomegaq) = U_ud - 2*U_uu &
                    - ABA(U_ud, chi_0(:, :, ikx, iky, iomegaq), nb*nb) &
                    + 1.5*ABA(U_s, chi_s(:, :, ikx, iky, iomegaq), nb*nb) &
                    + 0.5*ABA(U_c, chi_c(:, :, ikx, iky, iomegaq), nb*nb)

            enddo; enddo; enddo
            !$omp end parallel do


            !sigma(k) = V(k-k')*G(k')

            ! dft V to V_r_tau
            call dft(V, V_r_tau, nb*nb, nomegab, dft_grid, 1, 0)

            ! sigma_r_tau, 并行
            sigma_r_tau = complex_0
            !!$omp parallel do private(l2,m1,m2) reduction(+:sigma_r_tau)
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                sigma_r_tau(l1, m1, :, :, :) = sigma_r_tau(l1, m1, :, :, :) &
                    + V_r_tau(sub_g2chi(l1,l2), sub_g2chi(m1,m2),:,:,:) * G_r_tau(l2,m2,:,:,:)
            enddo; enddo; enddo; enddo
            !!$omp end parallel do

            ! idft sigma_r_tau to sigma
            call dft(sigma_r_tau, sigma, nb, dft_grid, nomegaf, -1, 1)

            !call testConvolution3sigma()

            !if (sigma_iter > 1) then
            sigma_conv = convergence_test(sigma_iter, 0)
            if (sigma_conv) then
                exit
            endif
            !endif
            sigma0 = sigma

            sigma=T_eV/nk*sigma

            ! write(stdout, *) 'calculating New G...'

            ! 新的G, dyson方程
            ! G=G0+G0*sigma*G, then we have G=(I-G0*sigma)**(-1)*G0
            !$omp parallel do private(iky,iomegak,G0_,sigma_,G_)
            do ikx=1,nkx;do iky=1,nky;do iomegak=minomegaf,maxomegaf
                if (sigma_state==0)then
                    G0_=G0(:,:,ikx,iky,iomegak)
                    sigma_=sigma(:,:,ikx,iky,iomegak)
                    G_=AB(G0_,sigma_,nb)
                    G_=I_G - G_
                    G_=inverseAbyB(G_,G0_,nb)
                    G1(:,:,ikx,iky,iomegak) = G_
                else
                    G(:,:,ikx,iky,iomegak) &
                        = &
                        G0(:,:,ikx,iky,iomegak) &
                        + ABC(G0(:,:,ikx,iky,iomegak), &
                        sigma(:,:,ikx,iky,iomegak), &
                        G1(:,:,ikx,iky,iomegak),nb)
                endif
            enddo;enddo;enddo
            !$omp end parallel do
            !sigma_conv=convergence_test(sigma_iter, 1)

            select case (mixer_method)
                case (0)
                    G=G1
                case (1)
                    G=mixer_beta*G1+(1-mixer_beta)*G
                case (2)
                    call mixer(sigma_iter)
            end select

            !G2=G

            sigma_iter=sigma_iter+1;
            total_iter = total_iter + 1

            if (total_iter>20) then
                !write(stdout,*) sigma_minus
            endif
        enddo

        ! sigma loop end

        ! 计算density
        cur_density=0d0
        !$omp parallel do private(ikx,iky,iomegak) reduction(+:cur_density)
        do ib=1,nb; do ikx=1,nkx; do iky=1,nky; do iomegak=minomegaf,maxomegaf
            cur_density = cur_density &
                + real(G(ib, ib, ikx, iky, iomegak)) &
                - real(G0(ib, ib, ikx, iky, iomegak))
        enddo; enddo; enddo; enddo
        !$omp end parallel do
        cur_density=cur_density*2*T_eV/nk + density_base

        write(stdout,*) 'density and mu: ', cur_density,'/', mu
        write(stdout,*)


        if (abs(cur_density-target_density)<density_tol) then
            density_conv=.true.
            !计算结束
        else
            call modify_mu()
            write(stdout,*) 'modified new mu = ', mu
        endif

        density_iter = density_iter + 1

        if (total_iter>20) then
            !exit
        endif

    enddo
    ! density loop end

    ! 迭代部分结束--------------------------------------------------------------------------

    ! output chi_s(q,0)
    write(stdout,*) 'chi_s at omega = 0'
    write(stdout,*) 'kx, ky, chi_s(real and imag)'
    do ikx=1,nkx; do iky=1,nky
        temp_complex=complex_0
        do l1=1,nb; do m1=1,nb
            temp_complex=temp_complex+chi_s(sub_g2chi(l1,l1),sub_g2chi(m1,m1),ikx,iky,0)
        enddo; enddo
        write(stdout, '(2F10.4,2F14.8)') k(ikx,iky,:), temp_complex
    enddo; enddo



    if (solve_eliashberg) then
        call eliashberg()
    endif

    write(stdout,*)
    write(stdout,*) 'good night.'
    write(stdout,*)

    mpi_info = mpi_finalize1()


end program

