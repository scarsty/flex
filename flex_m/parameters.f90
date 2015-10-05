module parameters
    implicit none

    ! 化学势, 占据数
    real(8) mu, n_density

    ! 目标占据数, 收敛误差
    REAL(8) target_density, density_tol

    ! U, U', J, J' for H1
    REAL(8) h1_U, h1_Up, h1_J, h1_Jp

    ! 温度
    real(8) T, T_beta, T_eV

    ! 保存的自能函数
    LOGICAL sigma_input, sigma_output
    CHARACTER*128 sigma_input_file, sigma_output_file

    ! 自能函数收敛判据
    REAL(8) sigma_tol
    INTEGER max_it
    REAL(8) alpha
    INTEGER alpha_scheme

    ! 自旋态
    integer spin_state

    logical, parameter :: solve_eliashberg = .true.
    logical, parameter :: test_band = .true.

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        sigma_tol, max_it, alpha,  alpha_scheme, spin_state



#ifdef USE_MPI
    INTEGER ierr
#endif /* USE_MPI */

    ! MPI variables
    INTEGER rank, size
    ! Timing variables
    Real(8) start_time, end_time
    Real(8) last_it_time, this_it_time

end module parameters
