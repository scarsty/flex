subroutine testBand()
    use constants
    use myfunctions
    ! use parameters
    use parameters2
    implicit none


    ! 自旋态不是3就是1
    ! 含矩阵乘, 需改写

    ! 点数36: 1~11~21~36
    integer count_k, i, ik, ix, iy, fileunit
    complex :: fac
    real :: rdotk
    real, dimension (2) :: temp
    real, dimension (36, 2) :: k_band
    complex, dimension (nb, nb, 36) :: h0_k_band
    complex, dimension (nb,nb) :: A, B
    complex, dimension (nb, 36) :: ev_band
    complex, dimension (nb) :: alpha, beta
    complex, dimension (nb, nb) :: vl, vr
    complex, dimension (nb*2) :: work
    integer lwork, info
    real, dimension (nb*8) :: rwork


    ! ------------------------------------------------------------------------

    ! 测试能带正确性
    ! 组合一组高对称点
    k_band = 0d0;
    count_k=0;
    do i = 0,10
        count_k = count_k+1;
        k_band(count_k, 1) = 0d0
        k_band(count_k, 2) = i*0.5d0/10;
    enddo
    do i = 1,10
        count_k = count_k+1;
        k_band(count_k, 1) = i*0.5d0/10;
        k_band(count_k, 2) = 0.5d0;
    enddo
    do i = 1,15
        count_k = count_k+1;
        k_band(count_k, 1) = 0.5d0-i*0.5d0/15;
        k_band(count_k, 2) = 0.5d0-i*0.5d0/15;
    enddo

    h0_k = complex_0

    write(stdout,*) 'build k-points of band...'

    ! 反傅里叶变换到k空间的能带
    ! 在每个k点上对角化得到能带特征值
    do ik=1,36
        h0_k_band=complex_0
        do ix = -rx, rx
            do iy = -ry, ry
                temp=[ix,iy]
                rdotk = two_pi*dot_product(k_band(ik,:),temp)
                fac=exp(complex_i*rdotk)
                h0_k_band(:,:,ik)=h0_k_band(:,:,ik)+fac*h0_r(:,:,ix, iy)
            enddo
        enddo
        A = h0_k_band(:,:,ik)
        ! if (ik==1) then
        !    write(stdout,*) A
        ! endif
        B = complex_0
        do i = 1,nb
            B(i,i)=complex_1
        enddo
        ! write(stdout,*) 'calling cggev...'
        call cggev('N', 'N', nb, A, nb, B, nb, alpha, beta, vl, 1, vr, 1, work, 2*nb, rwork, info)
        ! write(stdout,*) 'finish state ', info
        ev_band(:,ik) = alpha/beta
    enddo

    fileunit = 9999
    open(fileunit, file='testband.dat')
    do ik=1,36
        write(fileunit, *) ik, real(ev_band(:,ik))
    enddo

    close(fileunit)

    return

end subroutine testBand


subroutine testFunctions()
    implicit none

    real, dimension (512) :: k1, k2



    return

end subroutine testFunctions


subroutine testConvolution()
    use constants
    implicit none

    integer, parameter:: n=8,n1=n-1
    integer, parameter:: m=4,m1=m-1
    complex, dimension (-n1:n) :: f, g
    complex, dimension (-m1:m) :: f_ft

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


subroutine testConvolution2()
    use constants
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3.f03'

    integer i, x, t, i1, i2, j1, j2, i_minus_j1, i_minus_j2, num
    type(C_PTR) :: plan
    real :: start, finish, summary
    integer,parameter:: n=8,n1=n-1

    complex, dimension (0:n-1,0:n-1) :: f, g, g1
    complex, dimension (0:n-1,0:n-1) :: f_ft

    ! calculate 2D convolution g(x)=f(t)f(x-t) with 2 methods

    do i1=0,n1; do i2=0,n1
        f(i1,i2)=cmplx(i1*i2,i1+i2)
    enddo; enddo


    call cpu_time(start)
    do num=1,10000
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
        !write(0, *) g(i,:)
    enddo

    g1=g
    call cpu_time(start)
    do num=1,10000
        f_ft=0
        plan=fftwf_plan_dft_2d(n, n, f, f_ft, FFTW_BACKWARD, FFTW_ESTIMATE)
        call fftwf_execute_dft(plan, f, f_ft)
        call fftwf_destroy_plan(plan)
        !write(0, *) f_ft

        f_ft=f_ft*f_ft

        plan=fftwf_plan_dft_2d(n, n, f_ft, g, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftwf_execute_dft(plan, f_ft, g)
        call fftwf_destroy_plan(plan)
        g=g/n/n
    enddo

    call cpu_time(finish)
    write (0,*) finish-start

    write(0, *) 'calculated by fft: '
    do i=0,n1
        !write(0, *) g(i,:)
    enddo

    g=g-g1
    summary=0d0
    do i1=0,n1; do i2=0,n1
        summary=summary+abs(g(i1,i2))
    enddo; enddo

    write(0, *) 'error is ', summary

    return

end subroutine testConvolution2

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

    return

end subroutine build_h0_k

