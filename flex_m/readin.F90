
subroutine readin(T, target_density, density_tol, mu, &
        h1_U, h1_Up, h1_J, h1_Jp, &
        read_input, sigma_input_file, write_output, &
        sigma_output_file, sigma_tol, max_it, &
        alpha, alpha_scheme, h0_r)

    USE CONSTANTS
    use ifport
    use mpi_functions
    IMPLICIT NONE

    include "main_defs.f90"

    rank=mpi_rank()
    return
end subroutine readin
