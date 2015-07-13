! 此文件偷懒专用, include文本包含, 本身不参与编译

! MPI variables
INTEGER rank, size
! Timing variables
Real start_time, end_time
Real last_it_time, this_it_time

! 格林函数, 反常格林函数
! 自能函数, 反常自能函数
complex, dimension (nb, nb, nk, nomega) :: G, F, sigma, delta

! 极化率, susceptibilities, effective interactions
complex, dimension (nb*nb, nb*nb, nk, nomega) :: chi, chi_s, chi_c, V, V_s

! 交换能
real, dimension (nb*nb, nb*nb):: U_s, U_c, U_ud, U_uu

! 化学势, 占据数
real mu, n_density

! H0
complex, dimension (-2:2, -2:2, nb, nb) :: h0_r

! 目标占据数, 收敛误差
REAL target_density, density_tol

! U, U', J, J' for H1
REAL h1_U, h1_Up, h1_J, h1_Jp

! 温度
real T

! 保存的自能函数
LOGICAL read_input, write_output
CHARACTER*128 sigma_input_file, sigma_output_file

! 自能函数收敛判据
REAL sigma_tol
INTEGER max_it
REAL alpha
INTEGER alpha_scheme

#ifdef USE_MPI
INTEGER ierr
#endif /* USE_MPI */
