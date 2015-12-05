module parameters
    implicit none

    ! 化学势, 占据数
    real(4) mu, n_density

    ! 目标占据数, 收敛误差
    real(4) target_density, density_tol

    ! U, U', J, J' for H1
    real(4) h1_U, h1_Up, h1_J, h1_Jp

    ! 温度
    real(4) T, T_beta, T_eV

    ! 保存的自能函数
    logical sigma_input, sigma_output
    character*128 sigma_input_file, sigma_output_file

    ! 自能函数收敛判据
    real(4) G_tol
    integer max_iter
    real(4) alpha
    integer alpha_scheme

    integer iter_method ! 0: dyson方程, others: G_n+1=G0+G0*sigma*G_n

    integer mixer_method
    real(4) mixer_beta, mixer_beta0
    ! 自旋态
    integer :: spin_state

    logical :: solve_eliashberg = .true.
    logical :: test_band = .false.

    integer :: from_high_T=9

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        G_tol, max_iter, alpha, alpha_scheme, spin_state,&
        mixer_method, mixer_beta, &
        solve_eliashberg, test_band, &
        from_high_T

    real(4) eigen_value(5)
    namelist /band/ eigen_value

    ! MPI variables
    integer mpi_rank, mpi_size, mpi_info
    ! Timing variables
    real(4) start_time, end_time
    real(4) last_it_time, this_it_time

    integer alloc_error

end module parameters
