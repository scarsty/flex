MODULE Constants

    implicit none

    public
    integer, parameter :: nb=5 ! number of bands

    ! Lattice and frequency grid parameters. (all must be 2^n!)
    ! 松原频率的n
    integer, parameter :: nomega=8
    ! k点
    integer, parameter :: kx = 8   ! lattice dimension >=lcx
    integer, parameter :: ky = 8
    integer, parameter :: kz = 1
    integer, parameter :: nk = kx*ky*kz

    ! MPI related constants.
    integer, parameter :: np = 16 !number of processes
    !integer, parameter :: mp = m/np  ! must be an integer
    !integer, parameter :: ncp = nc/np  ! must be an integer

    ! Convenience constants
    !integer, parameter :: mp1 = mp - 1
    !integer, parameter :: m1 = m - 1

    ! Mathematical and physical constants
    REAL, parameter :: kB = 8.61734315d-05  ! eV/K, 玻尔兹曼常数
    REAL, parameter :: mub = 5.788381755d-5 ! eV/T, 磁导率, 此处无用
    REAL, parameter :: gs = 2.002319 !不知道

    real, parameter, public    :: pi=3.141592653589793238462643383279
    real, parameter, public    :: two_pi = 2*pi

    complex, parameter, public :: cmplx_i = (0.0,1.0)
    complex, parameter, public :: cmplx_0 = (0.0,0.0)
    complex, parameter, public :: cmplx_1 = (1.0,0.0)
    real, parameter, public::real_error = 1d-5
END MODULE CONSTANTS
