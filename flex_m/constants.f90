MODULE constants

    implicit none

    public
    integer, parameter :: nb = 5 ! number of bands
    integer, parameter :: square_nb=nb*nb ! order of chi

    ! Lattice and frequency grid parameters. (all must be 2^n!)

    ! 松原频率
    ! 这里的数据结构比较特殊, 对于格林函数和自能, 玻色频率上都是0, 对于极化率, 费米频率上都是0
    ! nomega 为正的费米频率个数
    ! 总费米频率个数 2*nomega, 总玻色频率个数4*nomega-2
    ! 计算范围较大值以避开周期
    ! 费米和玻色频率的循环中, 直接计算频率的值, 2*n-1和2*n
    integer, parameter :: nomega=4
    integer, parameter :: dft_grid=4*nomega  ! 这里实际上是fft的网格数
    integer, parameter :: minomegaf=-(nomega-1), maxomegaf=nomega
    integer, parameter :: minomegab=-2*(nomega-1), maxomegab=2*(nomega-1)
    integer, parameter :: nomegaf = 2*nomega, nomegab = 4*nomega-2

    ! 虚时间离散的点数
    ! integer, parameter :: ntau=256

    ! k点
    integer, parameter :: nkx = 4
    integer, parameter :: nky = 4
    integer, parameter :: nkz = 1
    integer, parameter :: nk = nkx*nky*nkz

    ! 归一化常数
    integer, parameter :: NN = nk

    ! r点, 表示r坐标 -r to r
    integer, parameter :: rx = 2
    integer, parameter :: ry = 2
    integer, parameter :: rz = 0

    ! MPI related constants.
    integer, parameter :: np = 16 !number of processes
    !integer, parameter :: mp = m/np  ! must be an integer
    !integer, parameter :: ncp = nc/np  ! must be an integer

    ! Convenience constants
    !integer, parameter :: mp1 = mp - 1
    !integer, parameter :: m1 = m - 1

    ! Mathematical and physical constants
    real(8), parameter :: kB = 8.61733035d-05  ! eV/K, 玻尔兹曼常数
    real(8), parameter :: mub = 5.788381755d-5 ! eV/T, 磁导率, 此处无用
    real(8), parameter :: gs = 2.002319 !不知道
    real(8), parameter :: h_bar  =6.582119514d-16 ! eV s, 普朗克常数

    real(8), parameter, public    :: pi=3.141592653589793238462643383279d0
    real(8), parameter, public    :: two_pi = 2*pi

    complex(8), parameter, public :: complex_i = (0.0d0,1.0d0)
    complex(8), parameter, public :: complex_0 = (0.0d0,0.0d0)
    complex(8), parameter, public :: complex_1 = (1.0d0,0.0d0)
    real(8), parameter, public :: real_error = 1d-5


    integer, parameter :: stdin =  5
    integer, parameter :: stdout = 6
    integer, parameter :: stderr = 0

END MODULE constants
