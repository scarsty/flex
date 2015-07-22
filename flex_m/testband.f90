subroutine eliashberg()
    use constants
    use myfunctions
    use parameters
    use parameters2
    implicit none

    ! 未完成
    ! 厄立希伯格方程, 直接应用上面得到的组合

    ! 自旋态不是3就是1
    ! 含矩阵乘, 需改写

    integer count_k, i, ik, ix, iy, fileunit
    real, dimension (31, 2) :: k_band
    complex, dimension (nb,nb,31) :: h0_k_band
    complex, dimension (nb,nb,31) :: A, B
    real, dimension(nb, 31) :: ev_band
    complex, dimension(nb) :: alpha, beta
    complex, dimension(nb, nb) :: vl, vr
    complex, dimension(nb*2) :: work
    integer lwork

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
    do i = 1,10
        count_k = count_k+1;
        k_band(count_k, 1) = 0.5d0-i*0.5d0/10;
        k_band(count_k, 2) = 0.5d0-i*0.5d0/10;
    enddo

    h0_k = complex_0


    do ik=1,30
        do ix = -2, 2
            do iy = -2, 2
                temp=[ix,iy]
                rdotk = two_pi*dot_product(k_band(ik,:),temp)
                fac=exp(complex_i*rdotk)
                h0_k_band(:,:,ik)=h0_k_band(:,:,ik)+fac*h0_r(:,:,ix, iy)
            enddo
        enddo
        A = h0_k_band(:,:,ik)
        B = complex_0
        do ik = 1,nb
            B(ik,ik)=complex_1
        enddo
        call cggbev('N', 'N', nb, A, nb, B, nb, alpha, beta, vl, 1, vr, 1)
    enddo

    fileunit = 9999
    open(fileunit, file='testband.dat')
    do ik=1,30
        write(fileunit, *) ik, ev_band(:,ik)
    enddo

close(fileunit)


    return

end subroutine testband
