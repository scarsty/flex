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
    integer conv_count
    !integer mpiinfo
    complex(8) temp_complex
    real(8) T_ev0

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

    mu=sum(eigen_value)/nb

    do while (.not. density_conv)

        G_conv=.false.
        G_iter=0

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
        if (density_iter==1) then
            G=G0
        endif
        call mixer_init()

        !call testConvolution()
        !call testConvolution3()
        !call testConvolution3G()
        write(stdout,*) 'mu = ', mu, t_ev/t_ev0
        ! base density
        density_base = 0d0
        do ib=1,nb; do ikx=1,nkx; do iky=1,nky
            !write(stdout, *) T_beta*(h0_k(ib,ib,ikx,iky)-mu)
            density_base=density_base+1/(exp(T_beta*(real(h0_k(ib,ib,ikx,iky))-mu))+1)
        enddo; enddo; enddo
        density_base=density_base*2/nk
        write(stdout, *) 'base density is ', density_base

        G_iter = 1
        write(stdout,*) '  iter   iter  conv.pts          norm.error'
        write(stdout,*) '-----------------------------------------------'

        ! sigma迭代中使用openmp并行
        do while (.not. G_conv)
            ! calculate chi_0 with chi(q)= -G1(q-k)G2(-k), the same to -G1(q-k)G2(k)**H
            ! dft G to G_r_tau

            call cal_G_out()
            call G_conv_test(G_conv)
            if (G_conv) then
                exit
            endif

            select case (mixer_method)
                case (0)
                    G=G_out
                case (1)
                    call mixer_linear()
                case (2:3)
                    call mixer_Pulay()
                case (4)
                    call mixer_Broyden()
            end select

            !if (mod(G_iter,100)==0) mixer_beta=mixer_beta/2
            !G2=G
            !sigma0 = sigma

            G_iter=G_iter+1;
            total_iter = total_iter + 1

            if (total_iter>max_iter) then
                !write(stdout,*) sigma_minus
            endif
        enddo

        ! sigma loop end
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
        endif

        density_iter = density_iter + 1

        if (total_iter>max_iter) then
            !exit
        endif
        stop
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

    call destroy()

end program

