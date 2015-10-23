program flex_m2d

    !   Compute Green function for superconductive properties
    !   author  : Sun TY (scarsty@gmail.com)
    !   status  : constructing
    !   version : none

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
    real(8) cur_density, density_base
    logical sigma_conv, density_conv
    real(8), dimension(2,2):: a, b

    integer  omega_f, omega_b
    real(8) tau

    ! 变量段结束-------------------------------------------------------------------------------

    call readin()

    if (test_band) then
        call testband()
    endif

    T_beta = 1d0/kB/T
    T_eV = kB*T

    sigma_state = 0

    ! 计算k点的坐标
    write(stdout, *) "Building k-points grid..."
    zero_k = 1    ! k原点
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

    ! 辅助变换
    i_minus = complex_0
    i_plus = complex_0
    do ix=1,nb
        i_plus(ix,ix)=complex_1
        i_minus(ix,ix)=complex_1
        if (ix==2 .or. ix==3) then
            i_plus(ix,ix)=complex_i
            i_minus(ix,ix)=-complex_i
        endif
    enddo

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

    if (nb==1) then
        call build_h0_k()
        !T_beta = 0.25
    endif

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




    ! 迭代部分-----------------------------------------------------------
    write(stdout, *) "Temperature in K = ", T
    write(stdout, *) "Temperature in eV = ", T_eV
    write(stdout, *) "Temperature in beta = ", T_beta
    write(stdout, *)

    write(stdout, *) "Begin to calculate FLEX"
    total_iter = 1
    density_iter = 1
    cur_density = 1000d0

    deltamu_per_density=1

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
        do while (.not. sigma_conv)

            ! write(stdout, *) 'calculating chi_0...'

            ! calculate chi_0 with chi(q)= -G1(q-k)G2(-k),

            ! dft G to G_r_tau
            call dft(G, G_r_tau, nb, nomegaf, dft_grid, 1, 0)

            conjgG=conjg(G)
            call dft(conjgG, conjgG_r_tau, nb, nomegaf, dft_grid, 1, 0)

            ! chi_0, 看起来需要并行
            ! 卷积形式, 改成减法 on tau
            chi_0_r_tau=0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                chi_0_r_tau(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
                    = - G_r_tau(l1, m1, :, :, :)*conjgG_r_tau(m2, l2, :, :, :)
            enddo; enddo; enddo; enddo

            ! idft chi_0_r_tau to chi_0
            call dft(chi_0_r_tau, chi_0, nb*nb, dft_grid, nomegab, -1, 0)
            chi_0 = T_ev/nk*chi_0
            ! call cleanError(chi_0, nb**4*nk*totalnomega)

            ! write(stdout, *) 'calculating chi_c, chi_s, V...'

            ! chi_c, chi_s, V, 需要并行和数学库
            ! 含有矩阵乘, 需改写
            do ikx=1,nkx; do iky=1,nky; do iomegaq=minomegab,maxomegab
                ! the same to solve AX=B, where A = (I +(c)/-(s) chi_0) and B = chi_0

                ! chi_c = chi_0 - chi_0*chi_c
                chi_0_=chi_0(:, :, ikx, iky, iomegaq)

                !Iminuschi_0_ = I_chi + AB(chi_0_,U_c,nb*nb)
                Iminuschi_0_ = I_chi + AB(U_c,chi_0_,nb*nb)

                chi_c_ = chi_0(:, :, ikx, iky, iomegaq)
                chi_c(:, :, ikx, iky, iomegaq)= inverseAbyB(Iminuschi_0_,chi_c_,nb*nb)

                ! chi_s = chi_0 + chi_0*chi_s
                !Iminuschi_0_ = I_chi - AB(chi_0_,U_s,nb*nb)
                Iminuschi_0_ = I_chi - AB(U_s,chi_0_,nb*nb)
                chi_s_ = chi_0(:, :, ikx, iky, iomegaq)
                chi_s(:, :, ikx, iky, iomegaq) = inverseAbyB(Iminuschi_0_,chi_s_,nb*nb)

                V(:, :, ikx, iky, iomegaq) = U_ud - 2*U_uu &
                    - ABA(U_ud, chi_0(:, :, ikx, iky, iomegaq), nb*nb) &
                    + 1.5*ABA(U_s, chi_s(:, :, ikx, iky, iomegaq), nb*nb) &
                    + 0.5*ABA(U_c, chi_c(:, :, ikx, iky, iomegaq), nb*nb)

            enddo; enddo; enddo
            !write(stdout,*) V(:, :, 1, 1, 0)
            !stop

            ! write(stdout, *) 'calculating sigma...'

            ! dft V to V_r_tau
            call dft(V, V_r_tau, nb*nb, nomegab, dft_grid, 1, 0)

            ! sigma_r_tau, 并行
            sigma_r_tau = complex_0
            do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
                sigma_r_tau(l1, m1, :, :, :) = sigma_r_tau(l1, m1, :, :, :) &
                    + V_r_tau(sub_g2chi(l1,l2), sub_g2chi(m1,m2),:,:,:) * G_r_tau(l2,m2,:,:,:)
            enddo; enddo; enddo; enddo

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

            !sigma_conv=convergence_test(sigma_iter, 1)

            select case (mixer_method)
                case (0)
                    G=G1
                case (1)
                    G=mixer_beta*G1+(1-mixer_beta)*G
                case (2)
                    call mixer(sigma_iter)
            end select

            G2=G
            sigma_iter=sigma_iter+1;
            total_iter = total_iter + 1

            if (total_iter>20) then
                !write(stdout,*) sigma_minus
            endif
        enddo

        ! sigma loop end

        ! 计算density
        density0 = cur_density
        cur_density=0d0

        do ib=1,nb; do ikx=1,nkx; do iky=1,nky; do iomegak=minomegaf,maxomegaf
            cur_density = cur_density &
                + real(G(ib, ib, ikx, iky, iomegak)) &
                - real(G0(ib, ib, ikx, iky, iomegak))
        enddo; enddo; enddo; enddo

        cur_density=cur_density*2*T_eV/nk + density_base

        write(stdout,*) 'density and mu: ', cur_density,'/', mu
        write(stdout,*)


        ! 下面是两个调整方法, 测试一下哪个稳定
        if (abs(cur_density-target_density)<density_tol) then
            density_conv=.true.
            !计算结束
        else
            if (density_iter>0) then
                deltamu_per_density = (mu-mu0)/(cur_density-density0)
            endif
            ! 计算化学势变化与占据数变化的比值来调整新的化学势
            mu0 = mu
            if (density_iter>1) then
                !mu = mu - (cur_density-target_density)
                mu = mu - (cur_density-target_density)*deltamu_per_density
            else
                mu = mu - 1.0d-1*sign(1.0d0, (cur_density-target_density)*deltamu_per_density)
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


end program

