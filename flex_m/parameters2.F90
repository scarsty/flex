!数组专用

! 格林函数, 反常格林函数
! 自能函数, 反常自能函数
complex, dimension (nb, nb, nk, -nomega:nomega) :: G, F, sigma, delta, G0, sigma0, G1
complex, dimension (nb, nb, nk, -ntau:ntau) :: G_tau, sigma_tau, delta_tau

! 极化率, susceptibilities, effective interactions
complex, dimension (nb*nb, nb*nb, nk, -nomega:nomega) :: chi_0, chi_s, chi_c, V, V_s
complex, dimension (nb*nb, nb*nb, nk, -ntau:ntau) :: chi_0_tau, V_tau, V_s_tau

! 交换能, 单位矩阵
complex, dimension (nb*nb, nb*nb):: U_s, U_c, U_ud, U_uu, I_chi

! H0
complex, dimension (nb, nb, -2:2, -2:2) :: h0_r
complex, dimension (nb, nb, nk):: h0_k

! k空间对应
real, dimension (nk, 2) :: k
! 两个k的差对应的k
integer, dimension (nk, nk) :: k_minus, k_plus

! 计算g0的辅助
complex, dimension (nb, nb) :: I_g0
real, dimension (nb, nb) :: u_g0


