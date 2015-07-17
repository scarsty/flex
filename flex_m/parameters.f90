module parameters
    implicit none

    ! 化学势, 占据数
    real mu, n_density

    ! 目标占据数, 收敛误差
    REAL target_density, density_tol

    ! U, U', J, J' for H1
    REAL h1_U, h1_Up, h1_J, h1_Jp

    ! 温度
    real T, T_beta

    ! 保存的自能函数
    LOGICAL sigma_input, sigma_output
    CHARACTER*128 sigma_input_file, sigma_output_file

    ! 自能函数收敛判据
    REAL sigma_tol
    INTEGER max_it
    REAL alpha
    INTEGER alpha_scheme

    ! 自旋态
    integer spin_state

	logical, parameter :: solve_eliashberg = .false.


#ifdef USE_MPI
    INTEGER ierr
#endif /* USE_MPI */

    ! MPI variables
    INTEGER rank, size
    ! Timing variables
    Real start_time, end_time
    Real last_it_time, this_it_time

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        sigma_tol, max_it, alpha,  alpha_scheme, spin_state

end module parameters
