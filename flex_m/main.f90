#define stdout 6
#define stderr 0

program flex_m
    use Constants
    use myfunctions
    implicit none
    include "parameters.F90"
    include "parameters2.F90"

    integer ix, iy, iz, count_k, zero_k, ib1, ib2, ik, iomega, ix1, ix2, ix3, ix4
    real rdotk, fac, temp(2)

    call readin(T, target_density, density_tol, mu, &
        h1_U, h1_Up, h1_J, h1_Jp, &
        sigma_input, sigma_input_file, sigma_output, &
        sigma_output_file, sigma_tol, max_it, &
        alpha, alpha_scheme, h0_r)

    T_beta = 1d0/kb/T

    !计算k点的坐标
    count_k = 0
    ! k原点
    zero_k = 1
    do ix = 1, kx
        do iy = 1, ky
            count_k = count_k + 1
            k(count_k, 1)=-1d0/2+1d0/kx*ix
            k(count_k, 2)=-1d0/2+1d0/ky*iy
            if ((abs(k(count_k, 1))<real_error) .and. (abs(k(count_k, 2))<real_error)) then
                zero_k=count_k
            endif
            !write(stdout, *) k(count_k,:)
        enddo
    enddo

    ! U
    ! 能带下标ab, cd -> (a+(b-1)*5, c+(d-1)*5)
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
    U_c = U_ud+u_uu

    ! 反傅里叶变换h0到k空间
    h0_k = cmplx_0
    do ik=1,nk
        do ix = -2, 2
            do iy = -2, 2
                temp(1)=ix
                temp(2)=iy
                rdotk = two_pi*dot_product(k(ik,:),temp)
                fac=exp(cmplx_i*rdotk)
                h0_k(:,:,ik)=h0_k(:,:,ik)+fac*h0_r(:,:,ix, iy)
            enddo
        enddo
    enddo
    ! 好像没归一化? h0_k=h0_k/?

    ! G0
    ! 费米频率 pi*(2n-1)
    G0=cmplx_0
    do ik=1,nk
        do iomega=1,nomega
            !G0(:,:,ik, iomega) = 1d0/(cmplx_i*pi*(2*iomega-1)-(h0_k(ik,ik,k)-mu))
        enddo
    enddo

    ! chi
    ! chi_c
    ! chi_s
    ! 玻色频率 pi*(2n-2), n从1开始
    print *, 'end.'
    return

end program

