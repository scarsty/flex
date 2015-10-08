 subroutine readin()
    use constants
    use parameters
    use parameters2
    use myfunctions
    IMPLICIT NONE

    real(8) eigen_value(5)
    namelist /band/ eigen_value

    integer ix, iy, iz, ib1, ib2, tempb1, tempb2
    character (len=100) :: text
    character (len=100) :: text1

    integer fileunit;

    complex(8) temp_h0

    !rank=mpi_rank()
    !读取
#ifdef _DEBUG
    fileunit = 10
    open(fileunit, file='input_file')
#else
    fileunit = stdin
#endif
    read(fileunit, nml=basic)
    read(fileunit, nml=band)

    if (nb/=5) then
        return
    endif

    ! 返回去重新读, 后面的格式不可出错
    rewind(fileunit)
    do while (.true.)
        read(fileunit, "(A)") text
        !write(stdout, *) text
        if (index(text, "HOPPING") > 0) then
            exit
        endif
    enddo

    ! w-s原胞点阵上的哈密顿量
    write(stdout,*)eigen_value
    h0_r = complex_0
    do ib1 = 1, nb
        !h0_r(ib1, ib1, 0, 0)=eigen_value(ib1)
    enddo
    do ix = -rx, rx
        do iy = -ry, ry
            !write(text1, "(A1), (I3), (I3), (A1), (A8)") "[", ix, iy, "]", "hopping"
            !write(stdout, *) text1
            read(fileunit,*)
            do ib1 = 1, nb
                do ib2 = 1, nb
                    read(fileunit,*) tempb1, tempb2, temp_h0
                    !if (ix/=0 .or. iy/=0) then
                        h0_r(tempb1, tempb2, ix, iy)=temp_h0
                    !endif
                    !write(stdout,*) ix,iy,ib1,ib2,h0_r(ix,iy,ib1,ib2)
                enddo
            enddo
            !read(stdin,*)
        enddo
    enddo

    return
end subroutine readin
