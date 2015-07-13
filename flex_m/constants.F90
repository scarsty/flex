MODULE Constants

    implicit none

    integer, parameter :: nb=5 ! number of bands

    ! Lattice and frequency grid parameters. (all must be 2^n!)
    ! 松原频率的n
    integer, parameter :: nomega=512
    ! k点
    integer, parameter :: kx = 32   ! lattice dimension >=lcx
    integer, parameter :: ky = 32
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
    REAL, parameter :: kb = 8.61734315d-05  ! eV/K, 波尔兹曼常数
    REAL, parameter :: mub = 5.788381755d-5 ! eV/T, 暂时不知道
    REAL, parameter :: gs = 2.002319

END MODULE CONSTANTS
