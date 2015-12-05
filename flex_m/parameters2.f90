!数组专用
module parameters2

    use constants
    implicit none

    public
    ! 格林函数, 反常格林函数
    ! 自能函数, 反常自能函数
    complex(4), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: &
        G0, G, conjgG, G_out, G_error, sigma, sigma0
        !G_in0, G_out0, G_error0, G_prev, G1_prev, deltaG

    complex(4), dimension (nb, nb, nkx, nky, dft_grid) :: &
        r_tau1, r_tau2

    ! 极化率, susceptibilities, effective interactions
    complex(4), dimension (nb*nb, nb*nb, nkx, nky, minomegab:maxomegab) :: &
        chi_0, V !,V_s, chi_s, chi_c,

    complex(4), dimension (nb*nb, nb*nb, nkx, nky, dft_grid) :: &
        r_tau_sqr


    ! 交换能, 单位矩阵
    complex(4), dimension (nb*nb, nb*nb) :: U_s, U_c, U_ud, U_uu, I_chi
    ! dyson方程的辅助计算矩阵
    complex(4), dimension (nb, nb) :: I_G, G_, G0_, sigma_

    ! H0
    ! u_h0_k是计算G0所使用的辅助幺正阵, G0=u_h0_k**H*diag(1/(i*omage-e+mu))*u_h0_k
    complex(4), dimension (nb, nb, -rx:rx, -ry:ry) :: h0_r
    complex(4), dimension (nb, nb, nkx, nky) :: h0_k
    real(4), dimension (nb, nb, nkx, nky) :: u_tilde_k, h0_tilde_k

    ! k空间对应
    real(4), dimension (nkx, nky, 2) :: k
    ! 两个k的差对应的k
    integer, dimension (nk, nk) :: k_minus, k_plus

    ! 计算g0的辅助
    complex(4), dimension (nb, nb, nkx, nky) :: u_h0_k
    complex(4), dimension (nb, nb) :: I_g0, i_plus, i_minus, h0_k_, diag_h0_G0_, u_h0_k_, ev_h0_k_lwork
    real(4), dimension (nb, nb) :: u_g0, u_tilde_k_, h0_tilde_k_, diag_h0_tilde_k_lwork, diag_test, ev_h0_k_rwork
    real(4), dimension (nb, nkx, nky) :: ev_h0_k
    real(4), dimension (nb) :: ev_h0_k_

    ! 计算chi的辅助
    complex(4), dimension (nb*nb, nb*nb) :: chi_0_, chi_c_, chi_s_, Iminuschi_0_

    ! 傅里叶变换, 反傅里叶变换的辅助
    complex(4), dimension (nkx, nky, dft_grid) :: dft_in
    complex(4), dimension (nkx, nky, dft_grid) :: dft_out

    ! Pulay mixer 相关
    integer, parameter :: mix_num  = 10, mix_keep = 5
    integer mixer_pointer, mixer_order, mixer_pointer2
    complex(4), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf, mix_num) :: mixer_G, mixer_error
    complex(4), dimension (nb, nb, nkx, nky, minomegaf:maxomegaf) :: mixer_error_, mixer_error2_
    ! 根据推导, 这部分应该都是实矩阵和向量
    real(4), dimension (0:mix_num, 0:mix_num) :: mixer_A, mixer_A1
    real(4), dimension (0:mix_num) :: mixer_x, mixer_b

    complex(4), allocatable, dimension(:,:) :: Jacobian

    ! 占据数相关
    integer, parameter :: mu_num = 100
    real(4), dimension (0:mu_num) :: mu_history, mu_error
    real(4), dimension (0:mu_num, 0:mu_num) :: mu_A, mu_A1
    real(4), dimension (0:mu_num) :: mu_x, mu_b

    ! 迭代G次数, density次数
    integer G_iter, density_iter, total_iter
    real(4) cur_density, density_base
    logical G_conv, density_conv

end module parameters2
