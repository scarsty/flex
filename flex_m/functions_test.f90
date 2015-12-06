module functions_test
    use functions_base
    use parameters
    use parameters2

contains
    ! 检查sigma收敛点数, 未使用, 因为用自能判断不太可靠
    ! 两个相同的输入G会导致sigma相同, 但是此时不能保证G收敛
    subroutine sigma_convergence_test(conv)

        implicit none
        !real(8), external :: dznrm2
        logical conv
        integer ib1, ib2, ikx, iky, iomegak, conv_grid
        real(8) norm_sigma_minus, norm_sigma, norm_sigma0, cur_G_tol, tol
        real(8) cur_error, total_error
        complex(8) sigma_one, sigma0_one

        norm_sigma = dznrm2(total_grid, sigma, 1)
        tol = norm_sigma*G_tol/total_grid
        conv_grid=0
        total_error=0
        !conv_grid=count(abs(sigma-sigma0)<=G_tol*abs(sigma))
        cur_G_tol = total_error/norm_sigma
        write(stdout,'(I7,I7,I10,ES20.5)') density_iter, G_iter, conv_grid, cur_G_tol
        conv = (conv_grid==total_grid)
    end subroutine

    subroutine conv_test(test1, test2, conv, output)
        implicit none

        complex(8), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: test1, test2
        !real(8), external :: dznrm2
        logical conv, output
        integer ib1, ib2, ikx, iky, iomegak, conv_grid
        real(8) norm_G, cur_G_tol, tol
        real(8) cur_error, total_error
        complex(8) G_one, G_out_one, G_error_one

        norm_G = dznrm2(total_grid, test1, 1)

        G_error=test2-test1
        !conv_grid=count(abs(G_error)<=G_tol*abs(G) .or. abs(G)<=real_error)
        conv_grid=0
        !total_error=0
        do ib1=1,nb; do ib2=1,nb; do ikx=1,nkx; do iky=1,nky; do iomegak=minomegaf,maxomegaf
            G_one=test1(ib1,ib2,ikx,iky,iomegak)
            G_error_one=G_error(ib1,ib2,ikx,iky,iomegak)
            cur_error=abs(G_error_one)
            !total_error = total_error + cur_error
            if (cur_error<=G_tol*abs(G_one) .or. abs(G_one)<=real_error) then
                !if (cur_error<=G_tol*abs(G_one)) then
                conv_grid=conv_grid+1
            endif
        enddo; enddo; enddo; enddo; enddo;

        cur_G_tol = dznrm2(total_grid, G_error, 1)!/norm_G
        conv = ((conv_grid==total_grid) .or. cur_G_tol<1d-8)
        !if (conv .or. mod(G_iter,20)==0) then
        if (output) then
            write(stdout,'(I7,I7,I10,ES20.5)') density_iter, G_iter, conv_grid, cur_G_tol
        endif

        !endif

    end subroutine


    subroutine testBand()
        implicit none

        ! 点数36: 1~11~21~36
        integer count_k, i, ik, ix, iy, fileunit
        complex(8) :: fac
        real(8) :: rdotk
        real(8), dimension (2) :: temp
        real(8), dimension (36, 2) :: k_band
        complex(8), dimension (nb, nb, 36) :: h0_k_band
        complex(8), dimension (nb,nb) :: A, B
        complex(8), dimension (nb, 36) :: ev_band
        complex(8), dimension (nb) :: alpha, beta
        complex(8), dimension (nb, nb) :: vl, vr
        complex(8), dimension (nb*2) :: work
        integer info
        real(8), dimension (nb*8) :: rwork


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
            call zggev('N', 'N', nb, A, nb, B, nb, alpha, beta, vl, 1, vr, 1, work, 2*nb, rwork, info)
            ! write(stdout,*) 'finish state ', info
            ev_band(:,ik) = alpha/beta
        enddo

        fileunit = 9999
        open(fileunit, file='testband.dat')
        do ik=1,36
            write(fileunit, *) ik, real(ev_band(:,ik))
        enddo

        close(fileunit)

    end subroutine

end module
