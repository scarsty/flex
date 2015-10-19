!数组专用
module parameters2

    use constants
    implicit none

    public
    ! 格林函数, 反常格林函数
    ! 自能函数, 反常自能函数
    complex(8), dimension (nb, nb, nkx, nky, 0:totalnomega-1) :: G, F, sigma, delta, G0, sigma0, G1, delta0, sigma_minus, conjgG, GGdelta
    complex(8), dimension (nb, nb, nkx, nky, 0:totalnomega-1) :: G_r_tau, sigma_r_tau, delta_r_tau, delta_r_tau0, conjgG_r_tau, GGdelta_r_tau

    ! 极化率, susceptibilities, effective interactions
    complex(8), dimension (nb*nb, nb*nb, nkx, nky, 0:totalnomega-1) :: chi_0, chi_s, chi_c, V, V_s
    complex(8), dimension (nb*nb, nb*nb, nkx, nky, 0:totalnomega-1) :: chi_0_r_tau, V_r_tau, V_s_r_tau

    ! 交换能, 单位矩阵
    complex(8), dimension (nb*nb, nb*nb) :: U_s, U_c, U_ud, U_uu, I_chi
    ! dyson方程的辅助计算矩阵
    complex(8), dimension (nb, nb) :: I_G, G_, G0_, sigma_

    ! H0
    ! u_h0_k是计算G0所使用的辅助幺正阵, G0=u_h0_k**H*diag(1/(i*omage-e+mu))*u_h0_k
    complex(8), dimension (nb, nb, -rx:rx, -ry:ry) :: h0_r
    complex(8), dimension (nb, nb, nkx, nky) :: h0_k
    real(8), dimension(nb, nb, nkx, nky) :: u_tilde_k, h0_tilde_k

    ! k空间对应
    real(8), dimension (nkx, nky, 2) :: k
    ! 两个k的差对应的k
    integer, dimension (nk, nk) :: k_minus, k_plus

    ! 计算g0的辅助
    complex(8), dimension (nb, nb, nkx, nky) :: u_h0_k
    complex(8), dimension (nb, nb) :: I_g0, i_plus, i_minus, h0_k_, diag_h0_G0_, u_h0_k_, ev_h0_k_lwork
    real(8), dimension (nb, nb) :: u_g0, u_tilde_k_, h0_tilde_k_, diag_h0_tilde_k_lwork, diag_test, ev_h0_k_rwork
    real(8), dimension (nb, nkx, nky) :: ev_h0_k
    real(8), dimension (nb) :: ev_h0_k_

    ! 计算chi的辅助
    complex(8), dimension (nb*nb, nb*nb) :: chi_0_, chi_c_, chi_s_, Iminuschi_0_

    ! 傅里叶变换, 反傅里叶变换的辅助, f/b means fermi and bose freq
    !complex, dimension (-nomega:nomega, -ntau:ntau) :: dft_f, dft_b
    !complex, dimension (-ntau:ntau, -nomega:nomega) :: idft_f, idft_b

    complex(8), dimension (nkx, nky, 0:totalnomega-1) :: dft_in
    complex(8), dimension (nkx, nky, 0:totalnomega-1) :: dft_out

    integer, parameter :: mix_num  = 10
    complex(8), dimension (nb, nb, nkx, nky, 0:totalnomega-1, mix_num) :: G_mixer
    real(8), dimension (mix_num+1, mix_num+1) :: Pulay_A
    real(8), dimension(mix_num+1) :: Pulay_x, Pulay_b

end module parameters2
