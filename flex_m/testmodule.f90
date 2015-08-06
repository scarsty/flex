subroutine testBand()
    use constants
    use myfunctions
    ! use parameters
    use parameters2
    implicit none

    ! 未完成
    ! 厄立希伯格方程, 直接应用上面得到的组合

    ! 自旋态不是3就是1
    ! 含矩阵乘, 需改写

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

    complex, dimension (0:3) :: f, g
    complex, dimension (0:3) :: f_ft

    integer i, x, t, i1, i2, x_minus_t


    ! calculate convolution g(x)=f(t)f(x-t) with 3 methods

    do i=0,3
        f(i)=i
    enddo

    g=0
    do x=0,3
        do t=0,3
            if (x-t>=0 .and. x-t<=3) then
                g(x)=g(x)+f(t)*f(x-t)
            endif
        enddo
    enddo

    write(0, *) 'directly calculated: ', g

    g=0
    do x=0,3
        do t=0,3
            x_minus_t = x-t
            do while (x_minus_t<0 .or. x_minus_t>3)
                x_minus_t=x_minus_t-sign(1, x_minus_t)*4
            enddo
            g(x)=g(x)+f(t)*f(x_minus_t)
        enddo
    enddo

    write(0, *) 'directly calculated with period: ', g

    f_ft=complex_0
    do i1=0,3;do i2=0,3
        f_ft(i1) = f_ft(i1) + exp(-2*pi*i1*i2/4*complex_i)*f(i2)
    enddo; enddo

    f_ft = f_ft*f_ft

    g=complex_0
    do i1=0,3;do i2=0,3
        g(i1) = g(i1) + exp(2*pi*i1*i2/4*complex_i)*f_ft(i2)
    enddo; enddo
    g=g/4

    write(0, *) 'calculated by ft: ', g


    return

end subroutine testConvolution

subroutine build_h0_k()
    use constants
    use parameters2
    implicit none

    integer l1, m1, ik, iomega

    do l1=1,nb; do m1=1,nb
        do ik=1,nk
            h0_k(l1,m1,ik) = - cos(k(ik,1)*two_pi) - cos(k(ik,2)*two_pi)
        enddo
    enddo; enddo

    return

end subroutine build_h0_k
