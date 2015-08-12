MODULE constants

	implicit none

	public
	integer, parameter :: nb=1 ! number of bands
	integer, parameter :: square_nb=nb*nb ! number of bands

	! Lattice and frequency grid parameters. (all must be 2^n!)

	! 松原频率个数, effective is -(nomega1)~(nomega1)
	integer, parameter :: nomega=8, nomega1 = nomega-1, nomega2 = 2*nomega
	!
	! 虚时间离散的点数
	! integer, parameter :: ntau=256

	integer, parameter :: total_omega=nomega*2+1
	integer, parameter :: total_tau=ntau*2+1

	! k点
	integer, parameter :: nkx = 8   ! lattice dimension >=lcx
	integer, parameter :: nky = 8
	integer, parameter :: nkz = 1
	integer, parameter :: nk = nkx*nky*nkz

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
	REAL, parameter :: kB = 8.61734315d-05  ! eV/K, 玻尔兹曼常数
	REAL, parameter :: mub = 5.788381755d-5 ! eV/T, 磁导率, 此处无用
	REAL, parameter :: gs = 2.002319 !不知道

	real, parameter, public    :: pi=3.141592653589793238462643383279
	real, parameter, public    :: two_pi = 2*pi
    ! real, public    :: pi=3.141592653589793238462643383279
	! real, public    :: two_pi

	complex, parameter, public :: complex_i = (0.0,1.0)
	complex, parameter, public :: complex_0 = (0.0,0.0)
	complex, parameter, public :: complex_1 = (1.0,0.0)
	real, parameter, public :: real_error = 1d-5


	integer, parameter :: stdin =  5
	integer, parameter :: stdout = 6
	integer, parameter :: stderr = 0

END MODULE constants
