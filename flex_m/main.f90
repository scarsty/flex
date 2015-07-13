
program flex_m
     implicit none
    use Constants


#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

#include "main_defs.F90"

      call readin(T, target_density, density_tol, mu, &
       h1_U, h1_Up, h1_J, h1_Jp, &
       read_input, sigma_input_file, write_output, &
       sigma_output_file, sigma_tol, max_it, &
       alpha, alpha_scheme, h0_r)


    print *, 'Hello World!'

end program

