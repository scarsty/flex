
subroutine testFunctions()
    implicit none

    real(8), dimension (512) :: k1, k2



    return

end subroutine testFunctions

function regionParameter(p,a,b)
    implicit none
    integer p,a,b,regionParameter,c

    c=p
    do while (c<a .or. c>b)
        if (c>b) c=c-(b-a+1)
        if (c<a) c=c+(b-a+1)
    enddo
    regionParameter=c

end function regionParameter

! 卷积测试手动傅立叶变换
subroutine testConvolution()
    use constants
    implicit none

    integer, parameter:: n=8,n1=n-1
    integer, parameter:: m=4,m1=m-1
    complex(8), dimension (-n1:n) :: f, g
    complex(8), dimension (-m1:m) :: f_ft

    integer i, x, t, i1, i2, x_minus_t


    ! calculate convolution g(x)=f(t)f(x-t) with 3 methods

    f=complex_0
    do i=-n1,n
        f(i)=cmplx(i,i)
    enddo

    g=0
    do x=-n1,n
        do t=-n1,n
            if (x-t>=-n1 .and. x-t<=n) then
                g(x)=g(x)+f(t)*f(x-t)
            endif
        enddo
    enddo

    write(0, *) 'directly calculated: '
    write(0, *) g

    g=0
    do x=-n1,n
        do t=-n1,n
            x_minus_t = x-t
            do while (x_minus_t<-n1 .or. x_minus_t>n)
                x_minus_t=x_minus_t-sign(1, x_minus_t)*(2*n)
            enddo
            g(x)=g(x)+f(t)*f(x_minus_t)
        enddo
    enddo

    write(0, *) 'directly calculated with period: '
    write(0, *) g
    f_ft=complex_0
    do i1=-m1,m;do i2=-n1,n
        f_ft(i1) = f_ft(i1) + exp(2*pi*i1*i2/(2*m)*complex_i)*f(i2)
    enddo; enddo

    f_ft = f_ft*f_ft

    g=complex_0
    do i1=-n1,n;do i2=-m1,m
        g(i1) = g(i1) + exp(-2*pi*i1*i2/(2*m)*complex_i)*f_ft(i2)
    enddo; enddo
    g=g/(m*2)

    write(0, *) 'calculated by ft: '
    write(0, *) g

    return

end subroutine testConvolution


! 卷积测试快速傅立叶变换2D
subroutine testConvolution2()
    use constants
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    integer i, x, t, i1, i2, j1, j2, i_minus_j1, i_minus_j2, num
    type(C_PTR) :: plan
    real(8) :: start, finish, summary
    integer,parameter:: n=8,n1=n-1

    complex(8), dimension (0:n-1,0:n-1) :: f, g, g1
    complex(8), dimension (0:n-1,0:n-1) :: f_ft

    ! calculate 2D convolution g(x)=f(t)f(x-t) with 2 methods

    do i1=0,n1; do i2=0,n1
        f(i1,i2)=cmplx(i1*i2,i1+i2)
    enddo; enddo

    call cpu_time(start)
    do num=1,1
        g=0
        do i1=0,n1; do i2=0,n1
            do j1=0,n1; do j2=0,n1
                i_minus_j1 = i1-j1
                i_minus_j2 = i2-j2
                do while (i_minus_j1<0 .or. i_minus_j1>n1)
                    i_minus_j1=i_minus_j1-sign(1, i_minus_j1)*n
                enddo
                do while (i_minus_j2<0 .or. i_minus_j2>n1)
                    i_minus_j2=i_minus_j2-sign(1, i_minus_j2)*n
                enddo
                g(i1,i2)=g(i1,i2)+f(j1,j2)*f(i_minus_j1,i_minus_j2)
            enddo; enddo
        enddo; enddo
    enddo
    call cpu_time(finish)
    write (0,*) finish-start

    write(0, *) 'directly calculated with period: '
    do i=0,n1
        write(0, *) g(i,:)
    enddo

    g1=g
    call cpu_time(start)
    do num=1,1
        f_ft=0
        plan=fftw_plan_dft_2d(n, n, f, f_ft, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, f, f_ft)
        call fftw_destroy_plan(plan)
        !write(0, *) f_ft

        f_ft=f_ft*f_ft

        plan=fftw_plan_dft_2d(n, n, f_ft, g, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, f_ft, g)
        call fftw_destroy_plan(plan)
        g=g/n/n
    enddo

    call cpu_time(finish)
    write (0,*) finish-start

    write(0, *) 'calculated by fft: '
    do i=0,n1
        write(0, *) g(i,:)
    enddo

    g=g-g1
    summary=0d0
    do i1=0,n1; do i2=0,n1
        summary=summary+abs(g(i1,i2))
    enddo; enddo

    write(0, *) 'error is ', summary

    return

end subroutine testConvolution2


! 卷积测试快速傅立叶变换3D
subroutine testConvolution3()
    use constants
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    integer i, x, t, i1, i2, i3, j1, j2, j3, i_minus_j1, i_minus_j2, i_minus_j3, num
    type(C_PTR) :: plan
    real(8) :: start, finish, summary
    integer,parameter:: n=2,n1=3,n2=1, n3=1

    !complex(8), dimension (n,n,n) :: f, g, g1, conjg_f
    !complex(8), dimension (n,n,n1) :: f_ft, conjg_f_ft
    complex(8), dimension (0:n1,0:n2,0:n3) :: f, g, g1, conjg_f
    complex(8), dimension (0:n1,0:n2,0:n3) :: f_ft, conjg_f_ft

    integer regionParameter
    external regionParameter

    ! calculate 2D convolution g(x)=f(t)f(x-t) with 2 methods
    f=1
    do i1=0,n1; do i2=0,n2; do i3=0,n3
        f(i1,i2,i3)=cmplx(-i1/(1+i2*i3),i1-i2+i3)
    enddo; enddo; enddo

    conjg_f=conjg(f)

    call cpu_time(start)
    do num=1,1
        g=complex_0
        !do i1=1,n; do i2=1,n; do i3=1,n
        !do j1=1,n; do j2=1,n; do j3=1,n
        do i1=0,n1; do i2=0,n2; do i3=0,n3
            do j1=0,n1; do j2=0,n2; do j3=0,n3
                i_minus_j1 = regionParameter(i1-j1,0,n1)
                i_minus_j2 = regionParameter(i2-j2,0,n2)
                i_minus_j3 = regionParameter(i3-j3,0,n3)
                g(i1,i2,i3)=g(i1,i2,i3)+f(j1,j2,j3)*f(i_minus_j1,i_minus_j2,i_minus_j3)
            enddo; enddo; enddo
        enddo; enddo; enddo
    enddo
    g1=g
    call cpu_time(finish)
    !write (0,*) finish-start

    write(0, *) 'directly calculated with period: '
    write(0, *) g

    call cpu_time(start)
    do num=1,1
        !f_ft=0
        plan=fftw_plan_dft_3d(n3+1, n2+1, n1+1, f, f_ft, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, f, f_ft)
        call fftw_destroy_plan(plan)

        plan=fftw_plan_dft_3d(n1+1, n2+1, n3+1, conjg_f, conjg_f_ft, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, conjg_f, conjg_f_ft)
        call fftw_destroy_plan(plan)
        !write(0, *) f_ft

        f_ft=f_ft*f_ft

        plan=fftw_plan_dft_3d(n3+1, n2+1, n1+1, f_ft, g, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, f_ft, g)
        call fftw_destroy_plan(plan)
        g=g/(n1+1)/(n2+1)/(n3+1)
    enddo

    call cpu_time(finish)
    !write (0,*) finish-start

    write(0, *) 'calculated by fft: '
    write(0, *) g

    g=g-g1
    summary=0d0
    do i1=0,n1; do i2=0,n2; do i3=0,n3
        summary=summary+abs(g(i1,i2,i3))
    enddo; enddo; enddo

    write(0, *) 'error is ', summary
    write(0, *)

    return

end subroutine testConvolution3


! 实际卷积测试
subroutine testConvolution3G()
    use constants
    use parameters2
    use functions
    implicit none

    integer l1, l2, m1, m2
    integer ikx1, iky1, ikx2, iky2, iomega1, iomega2, kxplus, kyplus, omegaplus
    complex(8) temp_complex

    integer regionParameter
    external regionParameter

    complex(8) a(1,1,1,1,8),b(1,1,1,1,8)

    a=0
    a(1,1,1,1,1)=1
    a(1,1,1,1,2)=2
    a(1,1,1,1,3)=3
    a(1,1,1,1,4)=4

    call dft(a,b,1,4,8,1,0)
    !write(stderr,*) b
    b=b*b
    call dft(b,a,1,8,7,-1,0)
    !write(stderr,*) a

    !stop

    ! dft G to G_r_tau
    write(stderr,*)nomegaf,nomegab
    call dft(G, r_tau1, nb, nomegaf, dft_grid, 1, 0)
    call dft(conjgG, r_tau2, nb, nomegaf, dft_grid, 1, 0)
    ! 卷积形式, 改成减法
    r_tau_sqr=0
    do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
        r_tau_sqr(sub_g2chi(l1, l2), sub_g2chi(m1, m2), :, :, :) &
            = - r_tau1(l1, m1, :, :, :)*r_tau2(m2, l2, :, :, :)
    enddo; enddo; enddo; enddo

    ! idft chi_0_r_tau to chi_0
    call dft(r_tau_sqr, chi_0, nb*nb, dft_grid, nomegab, -1, 0)
    write(stderr, *) chi_0(:,:,:,:,:)

    V = chi_0
    chi_0=complex_0
    do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
        do ikx1=1,nkx;do iky1=1,nky;do iomega1=minomegaf,maxomegaf
            do ikx2=1,nkx;do iky2=1,nky;do iomega2=minomegab,maxomegab
                kxplus=regionParameter(ikx1+ikx2-1,1,nkx)
                kyplus=regionParameter(iky1+iky2-1,1,nky)
                omegaplus=((2*iomega1-1+2*iomega2)+1)/2
                if (omegaplus>=minomegaf .and. omegaplus<=maxomegaf) then
                    chi_0(sub_g2chi(l1,l2), sub_g2chi(m1,m2), ikx2, iky2, iomega2) &
                        = chi_0(sub_g2chi(l1,l2), sub_g2chi(m1,m2), ikx2, iky2, iomega2) &
                        - G(l1, m1, kxplus, kyplus, omegaplus)*G(m2, l2, ikx1, iky1, iomega1)
                endif
            enddo;enddo;enddo
        enddo;enddo;enddo
    enddo;enddo;enddo;enddo
    write(stderr, *) chi_0(:,:,:,:,:)

    V = chi_0-V
    write(stderr, *) dznrm2(nb**4*nk*nomegab, V, 1) &
        / dznrm2(nb**4*nk*nomegab, chi_0, 1)
    stop
    !write(stderr, *) G0(1,1,1,1,1), G0(1,1,1,1,-1)

end subroutine testConvolution3G

subroutine testConvolution3sigma()
    use constants
    use parameters2
    use functions
    implicit none

    integer l1, l2, m1, m2
    integer ikx1, iky1, ikx2, iky2, iomega1, iomega2, kxplus, kyplus, omegaplus
    complex(8) temp_complex

    integer regionParameter
    external regionParameter

    sigma0=sigma
    sigma=0
    do l1=1,nb; do l2=1,nb; do m1=1,nb; do m2=1,nb
        do ikx1=1,nkx;do iky1=1,nky;do iomega1=minomegaf,maxomegaf
            do ikx2=1,nkx;do iky2=1,nky;do iomega2=minomegab,maxomegab
                kxplus=regionParameter(ikx1-ikx2-1,1,nkx)
                kyplus=regionParameter(iky1-iky2-1,1,nky)
                omegaplus=((2*iomega1-1-2*iomega2)+1)/2
                if (omegaplus>=minomegaf .and. omegaplus<=maxomegaf) then
                    sigma(l1,m1, ikx1, iky1, iomega1) &
                        = sigma(l1,m1, ikx1, iky1, (iomega1)) &
                        + V(sub_g2chi(l1,l2), sub_g2chi(m1,m2), ikx2, iky2,iomega2) &
                        * G(l2, m2, kxplus, kyplus, omegaplus)
                endif
            enddo;enddo;enddo
        enddo;enddo;enddo
    enddo;enddo;enddo;enddo


    write(stderr, *) sigma0
    write(stderr, *) sigma
    write(stderr, *) dznrm2(nb**2*nk*nomegaf, sigma-sigma0, 1)
    stop
    !write(stderr, *) G0(1,1,1,1,1), G0(1,1,1,1,-1)

end subroutine testConvolution3sigma

subroutine build_h0_k()
    use constants
    use parameters2
    implicit none

    integer l1, m1, ikx, iky

    h0_k = complex_0

    do l1=1,nb; do m1=1,nb
        do ikx=1,nkx; do iky=1,nky
            h0_k(l1,m1,ikx,iky) = - cos(k(ikx,iky,1)*pi) - cos(k(ikx,iky,2)*pi)
        enddo; enddo
    enddo; enddo
    write(stderr,*) h0_k
    return

end subroutine build_h0_k

!subroutine outputSomething()
!    implicit none
!    输出部分结果检查
!    write(stdout,*) dot_product(u_tilde_k_(1,:),u_tilde_k_(5,:))
!
!    ikx=2;iky=4
!    write(stdout,*) 'h:'
!    call writematrix(h0_k(:,:,ikx,iky),nb)
!    write(stdout,*) 'unitary matrix:'
!    call writematrix(u_h0_k(:,:,ikx,iky),nb)
!
!    write(stdout,*) 'eigenvalue:'
!    write(stdout,*) ev_h0_k(:,ikx,iky)
!
!    write(stdout,*) 'u back to h'
!    diag_h0_G0_=0d0
!    do ix=1,nb
!        !diag_h0_G0_(ix,ix)=complex_1
!        diag_h0_G0_(ix,ix)=ev_h0_k(ix,ikx,iky)
!        !diag_h0_G0_(ix,ix)=1-complex_i*ix
!    enddo
!    !write(stdout, '(5F8.3)') diag_test
!    u_h0_k_=u_h0_k(:,:,ikx,iky)
!
!    !u_h0_k_=u_h0_k(:,:,ikx,iky)
!    !call writematrix(matmul(u_h0_k_,u_h0_k(:,:,ikx,iky)),nb)
!
!    call writematrix(ABAH(u_h0_k_,diag_h0_G0_,nb),nb)
!
!
!    write(stdout, *)'next'
!
!    h0_k_=h0_k(:,:,ikx,iky)
!    h0_tilde_k_=real(AHBA(i_plus,h0_k_,nb))
!    write(stdout, *)'h~'
!    call writematrix(AHBA(i_plus,h0_k_,nb),nb)
!    u_tilde_k_=h0_tilde_k_
!    call dsyev('V','U',nb,u_tilde_k_,nb,ev_h0_k_,diag_h0_tilde_k_lwork,nb*nb,info)
!    write(stdout,*) 'eigenvalue:'
!    write(stdout,*) ev_h0_k_
!    u_h0_k_=u_tilde_k_
!
!    u_h0_k_=AB(i_plus,u_h0_k_,nb)
!    write(stdout,*) 'unitary matrix:'
!    call writematrix(u_h0_k_,nb)
!    write(stdout,*) 'u~ back to h~'
!    call writematrix(ABAH(u_h0_k_,diag_h0_G0_,nb),nb)
! end subroutine

! 部分废弃代码
!subroutine buildkminus()
!    use constants
!    use parameters2
!    implicit none
! k减法矩阵
! seems of no use
!        k_minus=0
!        do i1=1,nk
!            do i2=1,nk
!                temp = k(i1,:) - k(i2,:)
!                !write(stdout, *) temp
!                do while (abs(temp(1)+real_error)>0.5)
!                    temp(1) = temp(1) -sign(1., temp(1))
!                enddo
!                do while (abs(temp(2)+real_error)>0.5)
!                    temp(2) = temp(2) -sign(1., temp(2))
!                enddo
!                !write(stdout, *) temp
!                do i = 1,nk
!                    dis = norm2(temp-k(i,:))
!                    if (dis<real_error) then
!                        k_minus(i1, i2) = i
!                        exit
!                    endif
!                enddo
!                if (k_minus(i1, i2)<=0 .or. k_minus(i1, i2)>nk) then
!                    write(stdout, *) 'Wrong k_minus at', i1, i2
!                endif
!            enddo
!        enddo
!
!        ! k加法矩阵, k_minus(k1, k_minus(zero_k, k2))
!        do i1=1,nk
!            do i2=1,nk
!                k_plus(i1, i2) = k_minus(i1, k_minus(zero_k, i2))
!                if (k_plus(i1, i2)==0 .or. k_plus(i1, i2)>nk) then
!                    write(stdout, *) 'Wrong k_plus at', i1, i2
!                endif
!            enddo
!        enddo
!end subroutine buildkminus


!        h0_k_=h0_k(:,:,ikx,iky)
!        if (k(ikx,iky,1)>=k(ikx,iky,2)) then
!            h0_tilde_k(:,:,ikx,iky)=real(AHBA(i_plus,h0_k_,nb))
!        else
!            h0_tilde_k(:,:,ikx,iky)=real(AHBA(i_minus,h0_k_,nb))
!        endif
!        h0_tilde_k_=h0_tilde_k(:,:,ikx,iky)
!        u_tilde_k_=h0_tilde_k_
!        ! 这里h0_tilde变成一个实对称矩阵, 特征值全为实数, u为对应的正交变换阵
!        call dsyev('V','U',nb,u_tilde_k_,nb,diag_h0_tilde_k_,diag_h0_tilde_k_lwork,nb*nb,info)
!        u_tilde_k(:,:,ikx,iky)=u_tilde_k_
!        diag_h0_tilde_k(:,ikx,iky)=diag_h0_tilde_k_
!
!    write(stdout,*) 'unitary matrix:'
!    write(stdout,'(5F8.3)') u_tilde_k(:,:,ikx,iky)
!    write(stdout,'(5F8.3)') u_tilde_k(:,:,ikx,iky)
!    write(stdout,*) 'eigenvalue:'
!    write(stdout,'(5F8.3)') diag_h0_tilde_k(:,ikx,iky)
!    write(stdout,*) 'h:'
!    call writematrix(h0_k(:,:,ikx,iky),nb)
!
!    write(stdout,*) 'h~:'
!    write(stdout, '(5F8.3)') h0_tilde_k(:,:,ikx,iky)
!
!    write(stdout,*) 'u~ back to h~'
!    diag_test=0d0
!    do ix=1,nb
!        diag_test(ix,ix)=diag_h0_tilde_k(ix,ikx,iky)
!    enddo
!    !write(stdout, '(5F8.3)') diag_test
!    do ix=1,nb
!        do iy=1,nb
!            u_tilde_k_(ix,iy)=u_tilde_k(iy,ix,ikx,iky)
!        enddo!        if (abs(cur_density-target_density)<density_tol) then
!            density_conv=.true.
!            !计算结束
!        else
!            ! 根据占据数调整化学势
!            ! 第一步仅记录和猜测方向
!            ! 第三步开始逐步抛弃较远的点, 依照线性趋势逼近
!            ! 通常来说应保存一大一小
!            ! 靠不太容易设计
!            if (density_iter>=2) then
!                if (density_iter>2) then
!                    replaced=.false.
!                    do i=1,2
!                        if (abs(cur_density-target_density)<abs(density_old(i)-target_density)) then
!                            mu_old(i)=mu
!                            density_old(i)=cur_density
!                            replaced=.true.
!                            exit
!                        endif
!                    enddo
!                    if (.not.replaced) then
!                        max_diff_loc=1
!                        if (abs(density_old(1)-target_density)<abs(density_old(2)-target_density)) then
!                            max_diff_loc=2
!                        endif
!                        mu_old(max_diff_loc)=mu
!                        density_old(max_diff_loc)=cur_density
!                    endif
!                else
!                    mu_old(2)=mu
!                    density_old(2)=cur_density
!                endif
!                mu=(mu_old(1)-mu_old(2))/(density_old(1)-density_old(2))*(target_density-density_old(2))+mu_old(2)
!            elseif (density_iter==1) then
!                mu_old(1)=mu
!                density_old(1)=cur_density
!                !deltamu_per_density = (mu-mu0)/(cur_density-density0)
!                mu = mu - 1.0d-1*sign(1.0d0, (cur_density-target_density)*deltamu_per_density)
!            endif
!            write(stdout,*) 'modified new mu = ', mu
!        endif
!    enddo
!    !write(stdout, '(5F8.3)') u_tilde_k_
!    write(stdout, '(5F8.3)') matmul(matmul(u_tilde_k_,diag_test), u_tilde_k(:,:,ikx,iky))

! 辅助变换
!    i_minus = complex_0
!    i_plus = complex_0
!    do ix=1,nb
!        i_plus(ix,ix)=complex_1
!        i_minus(ix,ix)=complex_1
!        if (ix==2 .or. ix==3) then
!            i_plus(ix,ix)=complex_i
!            i_minus(ix,ix)=-complex_i
!        endif
!    enddo

! sometimes the linker cannot find this blas function
!    real function scnrm2(N, A, incx)
!        implicit none
!        integer N, incx, i
!        complex, dimension(N):: A
!        scnrm2=0d0
!        i=1
!        do while (i<=N)
!            scnrm2=scnrm2+abs(A(i))
!        enddo
!    end function scnrm2


! 生成傅里叶变换辅助矩阵dft_f, dft_b, idft_f, idft_b
!    subroutine buildDFTMatrix()
!        use constants
!        use parameters2
!        implicit none

!integer itau, tau, iomega, omega_f, omega_b

! 频域与时域的傅里叶变换辅助矩阵
! 可以一次计算一组, 或者构造大矩阵计算多组
! for one k-point, G_tau (行) = G_omega (行) * dft
! G_omega (行) = G_tau (行) * idft
!        do itau = -ntau,ntau
!            do iomega=-nomega,nomega
!                omega_f = 2*iomega-1
!                omega_b = 2*iomega
!                tau = itau*1d0/ntau
!                dft_f(iomega, itau) = exp(-2*pi*omega_f*tau*complex_i)
!                dft_b(iomega, itau) = exp(-2*pi*omega_b*tau*complex_i)
!                idft_f(itau, iomega) = exp(2*pi*omega_f*tau*complex_i)
!                idft_b(itau, iomega) = exp(2*pi*omega_b*tau*complex_i)
!            enddo
!        enddo
!
!        ! 这里的系数可能应该*2, 即所有频率一起考虑
!        dft_f = dft_f / (2*nomega+1)
!        dft_b = dft_b / (2*nomega+1)

!return
!end subroutine buildDFTMatrix

!        integer function sub_g2e(l,m,k,omega)
!        use constants
!        implicit none
!        integer l,m,k,omega
!        integer omegat
!        omegat=2*omega+1
!        sub_g2e = l*nb*nk*omegat+m*nk*omegat+omega+nomega+1 ! 未完成
!    end function sub_g2e
!
!    !n计算松原频率
!    integer function calfreq(omega, fb)
!        implicit none
!        integer omega, fb
!        if (fb/=0) then
!            calfreq=2*omega-1
!        else
!            calfreq=2*omega
!        endif
!    end function calfreq
!
!    ! k和松原频率的减法, fb: fermi(1) or bose(0)
!    ! 1 - 2 -> 3
!    subroutine komega_minus(k1, omega1, fb1, k2, omega2, fb2, k, k_minus, zero_k, k3, omega3, fb3)
!        use constants
!        implicit none
!        integer k1, omega1, fb1, k2, omega2, fb2, k3, omega3, fb3, sign_omega3, zero_k
!        integer f1, f2, f3
!        real(8), dimension (nk, 2) :: k
!        integer, dimension (nk, nk) :: k_minus
!
!        f1=calfreq(omega1, fb1)
!        f2=calfreq(omega2, fb2)
!        f3=f1-f2
!
!        fb3 = abs(mod(f3, 2))
!        omega3 = (f3+fb3)/2
!        k3 = k_minus(k1, k2)
!
!    end subroutine komega_minus
