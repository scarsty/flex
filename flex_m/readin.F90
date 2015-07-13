#define stdin  5
#define stdout 6
#define stderr 0

subroutine readin(T, target_density, density_tol, mu, &
        h1_U, h1_Up, h1_J, h1_Jp, &
        sigma_input, sigma_input_file, sigma_output, &
        sigma_output_file, sigma_tol, max_it, &
        alpha, alpha_scheme, h0_r)

    USE CONSTANTS
    use ifport
    use mpi_functions
    IMPLICIT NONE

    include "parameters.F90"

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        sigma_tol, max_it, alpha,  alpha_scheme

    real eigen_value(5)
    namelist /band/ eigen_value

    integer ix, iy, iz, ib1, ib2, p1, p2
    character (len=100) :: text
    character (len=100) :: text1

    rank=mpi_rank()
    !读取

    read(stdin, nml=basic)
    read(stdin, nml=band)


    write(stdout, *) "Temperature in K = ", T

    ! 返回去重新读, 后面的格式不可出错
    rewind(stdin)
    do while (.true.)
        read(stdin, "(A)") text
        !write(stdout, *) text
        if (index(text, "HOPPING") > 0) then
            exit
        endif
    enddo
    !哈密顿量
    do ix = -2, 2
        do iy = -2, 2
            !write(text1, "(A1), (I3), (I3), (A1), (A8)") "[", ix, iy, "]", "hopping"
            !write(stdout, *) text1
            read(stdin,*)
            do ib1 = 1, nb
                do ib2 = 1, nb
                    read(stdin,*) p1, p2, h0_r(ix, iy, ib1, ib2)
                    !write(stdout,*) ix,iy,ib1,ib2,h0_r(ix,iy,ib1,ib2)
                enddo
            enddo
            !read(stdin,*)
        enddo
    enddo

    return
end subroutine readin
