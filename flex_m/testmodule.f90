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
        do ix = -2, 2
            do iy = -2, 2
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
    implicit none

    real, dimension (512) :: k1, k2



    return
end subroutine testConvolution
