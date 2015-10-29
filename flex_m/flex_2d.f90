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
logical G_conv
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
        write(stdout,*) '  iter   iter  conv.pts           sigma.tol'
        write(stdout,*) '-----------------------------------------------'

        ! sigma迭代中使用openmp并行

        do while (.not. sigma_conv)

            ! calculate chi_0 with chi(q)= -G1(q-k)G2(-k), the same to -G1(q-k)G2(k)**H
            ! dft G to G_r_tau
            call dft(G, r_tau1, nb, nomegaf, dft_grid, 1, 0)

            conjgG=conjg(G)
            call dft(conjgG, r_tau2, nb, nomegaf, dft_grid, 1, 0)

            ! chi_0, 并行
            ! 卷积形式, 改成减法 on tau
            r_tau_sqr=0
            !$omp parallel do private(l2,m1,m2)
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                r_tau_sqr(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
                    = - r_tau1(l1, m1, :, :, :)*r_tau2(m2, l2, :, :, :)
            enddo; enddo; enddo; enddo
            !$omp end parallel do

            ! idft chi_0_r_tau to chi_0
            call dft(r_tau_sqr, chi_0, nb*nb, dft_grid, nomegab, -1, 0)
            chi_0 = T_ev/nk*chi_0

            ! chi_c, chi_s, V
            ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0
            !$omp parallel do private(ikx,iky,Iminuschi_0_,chi_0_,chi_c_,chi_s_)
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

            call convergence_test(sigma_conv)
            if (sigma_conv) then
                exit
            endif

            ! 新的G, dyson方程
            ! G=G0+G0*sigma*G, then we have G=(I-G0*sigma)**(-1)*G0
            !$omp parallel do private(ikx,iky,G0_,sigma_,G_)
            do iomegak=minomegaf,maxomegaf;do ikx=1,nkx;do iky=1,nky
                if (sigma_state==0)then
                    G0_=G0(:,:,ikx,iky,iomegak)
                    sigma_=sigma(:,:,ikx,iky,iomegak)
                    G_=AB(G0_,sigma_,nb)
                    G_=I_G - G_
                    call inverseAbyB(G_,G0_,nb)
                    G1(:,:,ikx,iky,iomegak) = G0_
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

            call convergence_testG(G_conv)

            select case (mixer_method)
                case (0)
                    G=G1
                case (1)
                    G=mixer_beta*G1+(1-mixer_beta)*G
                case (2:3)
                    call mixer(sigma_iter,mixer_method)
                    !G=mixer_beta*G1+(1-mixer_beta)*G
            end select

            !G2=G
            sigma0 = sigma

            sigma_iter=sigma_iter+1;
            total_iter = total_iter + 1

            if (total_iter>20) then
                !write(stdout,*) sigma_minus
            endif
        enddo

        ! sigma loop end
        call convergence_testG(G_conv)
        ! 计算density
        cur_density=0d0
        !$omp parallel do private(ikx,iky,ib) reduction(+:cur_density)
        do iomegak=minomegaf,maxomegaf; do ib=1,nb; do ikx=1,nkx; do iky=1,nky;
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
        !stop
    enddo
    ! density loop end

    ! 迭代部分结束--------------------------------------------------------------------------

    ! output chi_s(q,0), 未完成, 需要计算chi_s
    write(stdout,*) 'chi_s at omega = 0'
    write(stdout,*) '       kx        ky        chi_s(real and imag)'
    write(stdout,*) '---------------------------------------------------'
    do ikx=1,nkx; do iky=1,nky
        temp_complex=complex_0
        call cal_chi_cs(ikx,iky,0)
        do l1=1,nb; do m1=1,nb
            temp_complex=temp_complex+chi_s_(l1,m1)
        enddo; enddo
        write(stdout, '(2F10.4,2F14.8)') k(ikx,iky,:), temp_complex
    enddo; enddo


    if (solve_eliashberg) then
        call eliashberg()
    endif


    write(stdout,*)
    write(stdout,*) ' good night.'
    write(stdout,*)


    mpi_info = mpi_finalize1()


end program

