#define stdin  5
#define stdout 6
#define stderr 0

subroutine readin(T, target_density, density_tol, mu, &
        h1_U, h1_Up, h1_J, h1_Jp, &
        sigma_input, sigma_input_file, sigma_output, &
        sigma_output_file, sigma_tol, max_it, &
        alpha, alpha_scheme, h0_r)

    USE CONSTANTS
    !use ifport
    use mpi_functions
    IMPLICIT NONE

    include "parameters.F90"
    complex, dimension (nb, nb, -2:2, -2:2) :: h0_r

    namelist /basic/ T, target_density, density_tol, mu,&
        h1_U, h1_Up, h1_J, h1_Jp,&
        sigma_input,  sigma_input_file,&
        sigma_output, sigma_output_file,&
        sigma_tol, max_it, alpha,  alpha_scheme

    real eigen_value(5)
    namelist /band/ eigen_value

    integer ix, iy, iz, ib1, ib2, tempb1, tempb2
    character (len=100) :: text
    character (len=100) :: text1

    integer fileunit;

    complex temp_h0

    rank=mpi_rank()
    !读取
#ifdef _DEBUG
    fileunit = 10
    open(fileunit, file='input_file')
#else
    fileunit = stdin
#endif
    read(fileunit, nml=basic)
    read(fileunit, nml=band)


    write(stdout, *) "Temperature in K = ", T

    ! 返回去重新读, 后面的格式不可出错
    rewind(fileunit)
    do while (.true.)
        read(fileunit, "(A)") text
        !write(stdout, *) text
        if (index(text, "HOPPING") > 0) then
            exit
        endif
    enddo
    !哈密顿量
    h0_r = cmplx_0
    do ib1 = 1, nb
        h0_r(ib1, ib2, 0, 0)=eigen_value(ib1)
    enddo
    do ix = -2, 2
        do iy = -2, 2
            !write(text1, "(A1), (I3), (I3), (A1), (A8)") "[", ix, iy, "]", "hopping"
            !write(stdout, *) text1
            read(fileunit,*)
            do ib1 = 1, nb
                do ib2 = 1, nb
                    read(fileunit,*) tempb1, tempb2, temp_h0
                    h0_r(tempb1, tempb2, ix, iy)=temp_h0
                    !write(stdout,*) ix,iy,ib1,ib2,h0_r(ix,iy,ib1,ib2)
                enddo
            enddo
            !read(stdin,*)
        enddo
    enddo

    return
end subroutine readin
