!数组专用

! 格林函数, 反常格林函数
! 自能函数, 反常自能函数
complex, dimension (nb, nb, nk, nomega) :: G, F, sigma, delta, G0

! 极化率, susceptibilities, effective interactions
complex, dimension (nb*nb, nb*nb, nk, nomega) :: chi, chi_s, chi_c, V, V_s

! 交换能
real, dimension (nb*nb, nb*nb):: U_s, U_c, U_ud, U_uu

! H0
complex, dimension (nb, nb, -2:2, -2:2) :: h0_r
complex, dimension (nb, nb, nk):: h0_k

! k空间对应
real, dimension (nk, 2) :: k
! 两个k的差对应的k
integer, dimension (nk, nk) :: k_minus
