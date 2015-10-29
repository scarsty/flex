module parameters
    implicit none

    ! 化学势, 占据数
    real(8) mu, n_density

    ! 目标占据数, 收敛误差
    real(8) target_density, density_tol

    ! U, U', J, J' for H1
    real(8) h1_U, h1_Up, h1_J, h1_Jp

    ! 温度
    real(8) T, T_beta, T_eV

    ! 保存的自能函数
    logical sigma_input, sigma_output
    character*128 sigma_input_file, sigma_output_file

    ! 自能函数收敛判据
    real(8) sigma_tol
    integer max_it
    real(8) alpha
    integer alpha_scheme

    integer sigma_state ! 0: dyson方程, others: G_n+1=G0+G0*sigma*G_n

    integer mixer_method
    real(8) mixer_beta
    ! 自旋态
    integer spin_state

    logical solve_eliashberg
    logical test_band

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        sigma_tol, max_it, alpha,  alpha_scheme, spin_state,&
        mixer_method, mixer_beta, &
        solve_eliashberg, test_band

   ! MPI variables
    INTEGER mpi_rank, mpi_size, mpi_info
    ! Timing variables
    Real(8) start_time, end_time
    Real(8) last_it_time, this_it_time

end module parameters
