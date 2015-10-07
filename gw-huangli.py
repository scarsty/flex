#!/usr/bin/env python

""" gw_dmft_new_exp.py

    author  : Li HUANG (email : li.huang@unifr.ch)
    status  : unstable
    version : Ver.2015/01/07D
"""

# python standard library
import os
import sys
import time
import math
import shutil
import ctypes

# python 3rd-party library
import numpy
import scipy.weave
import scipy.interpolate

# python binding for ALPS
import pyalps.mpi
import pyalps.cthyb

def gw_dmft_config():
    """ setup the config parameters for the gw + edmft code.
        note: don't modify this function. if you want to change
        the config parameters, please edit the in.param file
    """
    # setup the parameters
    # param : a dict object, used to record all the necessary parameters
    param = {
        'U'     : [2.20  , "Onsite Coulomb interaction"             ],
        'V'     : [0.00  , "Intersite Coulomb interaction"          ],
        't'     : [0.25  , "Hopping parameter for 2D Hubbard model" ],
        'mune'  : [1.10  , "Chemical potential"                     ],
        'beta'  : [100.  , "Inverse temperature"                    ],
        'occup' : [1.00  , "Impurity occupation number"             ],
        'nkx'   : [81    , "Number of k-points along x axis"        ],
        'nky'   : [81    , "Number of k-points along y axis"        ],
        'niter' : [40    , "Maximum iteration number"               ],
        'alpha' : [0.7   , "Mixing parameter"                       ],
        'nspin' : [2     , "Number of (spin-)orbitals"              ],
        'ntime' : [1024  , "Number of imaginary time points"        ],
        'nffrq' : [512   , "Number of fermionic frequencies"        ],
        'nbfrq' : [512   , "Number of bosonic frequencies"          ],
        'start' : [True  , "Whether start from scratch"             ],
        'sy_ph' : [True  , "Whether enforce symmetrization (PH)"    ],
        'sy_pm' : [True  , "Whether enforce symmetrization (PM)"    ],
        'gw_ty' : [0     , "Whether the GW or FLEX part is on"      ],
        'gw_it' : [10000 , "After n-iter the GW or FLEX part is on" ],
        'dc_ty' : [1     , "Scheme for double counting term"        ],
        'alps'  : [False , "Whether the ct-hyb in ALPS is used"     ],
        'sskip' : [False , "Whether skip the ct-qmc impurity solver"],
    }

    # read config parameters from in.param file
    # only master node can do it
    if ( pyalps.mpi.rank == 0 ) :
        # check whether file in.param exists
        if os.path.isfile("in.param"):
            in_param = {}
            # read data, pay attention to their data type
            f_in = open("in.param", "r")
            for line in f_in:
                A = line.split()
                if   A[2] == 'float':
                    in_param[A[0]] = float(A[4])
                elif A[2] == 'int'  :
                    in_param[A[0]] = int(A[4])
                elif A[2] == 'bool' :
                    in_param[A[0]] = eval(A[4])
                else :
                    sys.exit("FATAL ERROR : wrong data type")
            # close the file handler
            f_in.close()

            # now the config parameters are stored in in_param dict,
            # we need to transfer them to param dict
            for key in in_param.iterkeys():
                param[key][0] = in_param[key]
        # if in.param does not exist, we have to stop the code
        else:
            sys.exit("FATAL ERROR : in.param does not exist")

    # broadcast the param dict to all children processes
    param = pyalps.mpi.broadcast(pyalps.mpi.world, param, 0)
    pyalps.mpi.world.barrier()

    # dump the parameters into the stdout
    if pyalps.mpi.rank == 0:
        print (">>> Listing parameters:")
        for k in sorted(param.keys()):
            print ("    " + param[k][1].ljust(42) + " : " + k.rjust(5) + " : " + str(param[k][0]).rjust(5))
        print ("")

    return param

def gw_dmft_sconfig(p):
    """ extract some key parameters which are useful for the impurity solver
        p : param dictionary
    """
    if p["alps"][0]: # just for alps ctqmc impurity solver
        s_para = {
            'U'           : p["U"    ][0],     # Coulomb interaction
            'BETA'        : p["beta" ][0],     # inversion temperature
            'MU'          : p["mune" ][0],     # chemical potential
            'N_MATSUBARA' : p["nffrq"][0],     # number of matsubara frequencies
            'N_TAU'       : p["ntime"][0] - 1, # number of imaginary time points
        }
    return s_para

def gw_dmft_header():
    """ print the header to standard terminal
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    print ("GE_2d_new :: MBPT + DMFT")
    print ("Ver.2015/01/03D (accelerated by mpi + scipy.weave)")
    print ("Developed by Li HUANG (li.huang@unifr.ch)")
    print ("Department of Physics, Fribourg University, Switzerland")
    print ("")

    print (">>> Starting in " + time.asctime())
    print (">>> Using %i processors, current rank is %i" % (pyalps.mpi.size, pyalps.mpi.rank))
    print ("")

def gw_dmft_iter(p, curr_iter):
    """ print the iteration information to standard terminal
        p : param dictionary
        curr_iter : current iteration number
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    for case in switch( p["gw_ty"][0] ):
        if case(5):
            print ("-"*60)
            print (">>> SOPT + DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

        if case(4):
            print ("-"*60)
            print (">>> TMA + DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

        if case(3):
            print ("-"*60)
            print (">>> FLEX_pp + DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

        if case(2):
            print ("-"*60)
            print (">>> FLEX_ph + DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

        if case(1):
            print ("-"*60)
            print (">>> SGW + DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

        if case(0):
            print ("-"*60)
            print (">>> DMFT ITER : current => %i  maximum => %i" % (curr_iter, p["niter"][0]))
            print ("-"*60)
            print ("")
            break

def gw_dmft_footer(gw_type):
    """ print the footer to standard terminal
        gw_type : the calculating scheme
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    for case in switch( gw_type ):
        if case(5):
            print ("-"*60)
            print (">>> SOPT + DMFT ITER : Finished")
            print ("-"*60)
            break

        if case(4):
            print ("-"*60)
            print (">>> TMA + DMFT ITER : Finished")
            print ("-"*60)
            break

        if case(3):
            print ("-"*60)
            print (">>> FLEX_pp + DMFT ITER : Finished")
            print ("-"*60)
            break

        if case(2):
            print ("-"*60)
            print (">>> FLEX_ph + DMFT ITER : Finished")
            print ("-"*60)
            break

        if case(1):
            print ("-"*60)
            print (">>> SGW + DMFT ITER : Finished")
            print ("-"*60)
            break

        if case(0):
            print ("-"*60)
            print (">>> DMFT ITER : Finished")
            print ("-"*60)
            break
    print (">>> Ending in " + time.asctime())

class switch(object):
    """ This class provides the functionality we want. You only need to look
        at this if you want to know how this works. It only needs to be defined
        once, no need to muck around with its internals.
    """
    def __init__(self, value):
        """ class constructor for the switch
        """
        self.__value = value
        self.__fall = False

    def __iter__(self):
        """ return the match method once, then stop
        """
        yield self.match
        raise StopIteration

    def match(self, *args):
        """ indicate whether or not to enter a case suite
        """
        if self.__fall or not args:
            return True
        elif self.__value in args:
            self.__fall = True
            return True
        else:
            return False

class Symm(object):
    """ class Symm is used to symmetrize the input data over imaginary time
        or spin. it can be also used to get rid of the real part or imaginary
        part of the input data.
    """
    @staticmethod
    def symm_over_time(time_data):
        """ perform symmetrization over imaginary time
            time_data : imaginary time data
        """
        # get number of elements
        n = len(time_data)
        for i in range(n):
            raux = 0.5 * ( time_data[i] + time_data[n-1-i] )
            time_data[i] = raux
            time_data[n-1-i] = raux

    @staticmethod
    def symm_over_spin(time_data):
        """ perform symmetrization over spin
            time_data : imaginary time data
        """
        raux = 0.5 * ( time_data[0] + time_data[1] )
        time_data[0] = raux.copy()
        time_data[1] = raux.copy()
        del raux

    @staticmethod
    def symm_over_real(freq_data):
        """ get rid of the imaginary part of the input data
            freq_data : matsubara frequency data
        """
        # get number of elements
        n = len(freq_data)
        for i in range(n):
            freq_data[i] = freq_data[i].real + 0j

    @staticmethod
    def symm_over_imag(freq_data):
        """ get rid of the real part of the input data
            freq_data : matsubara frequency data
        """
        # get number of elements
        n = len(freq_data)
        for i in range(n):
            freq_data[i] = 0.0 + 1j * freq_data[i].imag

    @staticmethod
    def symm_smth_data(data, naver):
        """ used to smooth 1-d float data
            data : 1-d array
            naver : used to define the smooth region
        """
        # determine the length of 1-d data
        n = len(data)

        # create a temporary array
        data_s = numpy.zeros_like(data, dtype = numpy.float)

        for i in range(n):
            # bracketing the smooth region [i1, i2]
            i1 = i - naver
            if i1 < 0: i1 = 0

            i2 = i + naver
            if i2 > n - 1: i2 = n - 1

            # accumulate the data
            count = 0
            for j in range(i1, i2):
                count = count + 1
                data_s[i] = data_s[i] + data[j]

            # calculate the average
            data_s[i] = data_s[i] / float(count)

        # copy smoothed data to original data
        for i in range(n):
            data[i] = data_s[i]

        # release memory
        del data_s

class Mesh(object):
    """ class Mesh is used to define fermionic and bosonic meshes, imaginary
        time mesh is also defined in this class.
        note : the k-point mesh is defined in the Lattice class.
    """
    def __init__(self, p):
        """ class constructor for the Mesh
            p : param dictionary
            beta  : inversion temperature
            ntime : number of imaginary time points
            nffrq : number of frequencies for fermionic mesh
            nbfrq : number of frequencies for bosonic mesh
        """
        self.__beta  = p["beta"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

    def create_mesh(self, mesh_type):
        """ construct imaginary time, fermionic or bosonic mesh according
            to the mesh_type variable.
            mesh_type : which type of mesh should be constructed
        """
        if   mesh_type == 1: # construct imaginary time mesh
            tmesh = numpy.zeros(self.__ntime, dtype = numpy.float)
            for i in range(len(tmesh)):
                tmesh[i] = float(i) * self.__beta / ( self.__ntime - 1 )
            return tmesh

        elif mesh_type == 2: # construct fermionic mesh
            fmesh = numpy.zeros(self.__nffrq, dtype = numpy.float)
            for i in range(len(fmesh)):
                fmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / self.__beta
            return fmesh
        else:                # construct bosonic mesh
            bmesh = numpy.zeros(self.__nbfrq, dtype = numpy.float)
            for i in range(len(bmesh)):
                bmesh[i] = ( 2.0 * i + 0.0 ) * math.pi / self.__beta
            return bmesh

    @staticmethod
    def create_time_mesh(beta, ntime):
        """ construct the imaginary time mesh
            beta  : inversion temperature
            ntime : number of imaginary time points
        """
        tmesh = numpy.zeros(ntime, dtype = numpy.float)
        for i in range(len(tmesh)):
            tmesh[i] = float(i) * beta / ( ntime - 1)
        return tmesh

    @staticmethod
    def create_fermion_mesh(beta, nffrq):
        """ construct the fermionic mesh
            beta  : inversion temperature
            nffrq : number of frequencies for fermionic mesh
        """
        fmesh = numpy.zeros(nffrq, dtype = numpy.float)
        for i in range(len(fmesh)):
            fmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / beta
        return fmesh

    @staticmethod
    def create_boson_mesh(beta, nbfrq):
        """ create the bosonic mesh
            beta  : inversion temperature
            nbfrq : number of frequencies for bosonic mesh
        """
        bmesh = numpy.zeros(nbfrq, dtype = numpy.float)
        for i in range(len(bmesh)):
            bmesh[i] = ( 2.0 * i + 0.0 ) * math.pi / beta
        return bmesh

class Dump(object):
    """ class Dump is used to dump the key variables to the disk file.
    """
    @staticmethod
    def dump_time_data1(file_name, time_mesh, time_data):
        """ dump the imaginary time data to the disk file
            file_name : file name
            time_mesh : imaginary time mesh
            time_data : array in imaginary time mesh
        """
        f = open(file_name, "w")
        for i in range(len(time_mesh)):
            f.write( "%6u %16.8s %16.8E\n" %
                   (i + 1, time_mesh[i], time_data[i]) )
        f.close()

    @staticmethod
    def dump_time_data2(file_name, time_mesh, time_data):
        """ dump the imaginary time data to the disk file, spin-resolved version
            file_name : file name
            time_mesh : imaginary time mesh
            time_data : array in imaginary time mesh
        """
        f = open(file_name, "w")
        for i in range(len(time_mesh)):
            f.write( "%6u %16.8s %16.8E %16.8E\n" %
                   (i + 1, time_mesh[i], time_data[0][i], time_data[1][i]) )
        f.close()

    @staticmethod
    def dump_freq_data1(file_name, freq_mesh, freq_data):
        """ dump the matsubara frequency data to the disk file
            file_name : file name
            freq_mesh : matsubara frequency mesh
            freq_data : array in matsubara frequency mesh
        """
        f = open(file_name, "w")
        for i in range(len(freq_mesh)):
            f.write( "%6u %16.8s %16.8E %16.8E\n" %
                   (i + 1, freq_mesh[i], freq_data[i].real, freq_data[i].imag) )
        f.close()

    @staticmethod
    def dump_freq_data2(file_name, freq_mesh, freq_data):
        """ dump the matsubara frequency data to the disk file, spin-resolved version
            file_name : file name
            freq_mesh : matsubara frequency mesh
            freq_data : array in matsubara frequency mesh
        """
        f = open(file_name, "w")
        for i in range(len(freq_mesh)):
            f.write( "%6u %16.8s %16.8E %16.8E %16.8E %16.8E\n" %
                   (i + 1, freq_mesh[i], freq_data[0][i].real, freq_data[0][i].imag,
                                         freq_data[1][i].real, freq_data[1][i].imag) )
        f.close()

    @staticmethod
    def dump_htau_data(file_name, time_mesh, htau_data):
        """ dump the hybridization function \Delta(\tau) to the disk file,
            spin-resolved version, this function is designed for the ctqmc
            impurity solver contained in ALPS package
            file_name : file name
            time_mesh : imaginary time mesh
            htau_data : array in imaginary time mesh
        """
        f = open(file_name, "w")
        for i in range(len(time_mesh)):
            f.write( "%6u %16.8E %16.8E\n" % (i, htau_data[0][i], htau_data[1][i]) )
        f.close()

    @staticmethod
    def dump_kdep_data(file_name, nfreq, nkx, nky, kdep_data):
        """ dump k-dependent data in matsubara frequency axis to the disk file
            file_name : file name
            nfreq     : number of frequency points
            nkx       : number of k-points
            nky       : number of k-points
            kdep_data : array in matsubara frequency axis
        """
        is_complex = isinstance(kdep_data[0,0,0], complex)
        fk = open(file_name, "w")
        for f in range(nfreq):
            if is_complex:
                # write real part
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky].real ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")

                # write imaginary part
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky].imag ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")
            else:
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky] ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")
        fk.close()

class Fourier(object):
    """ class Fourier is used to perform forward and backward fourier
        transformations between imaginary time and matsubara frequency
        spaces.
    """
    def __init__(self, p, ftype):
        """ class constructor for the Fourier
            p : param dictionary
            ftype : fourier type
            beta  : inversion temperature
            ntime : number of imaginary time points
            nfreq : number of frequencies for fermionic/bosonic mesh
        """
        self.__beta  = p["beta"][0]
        self.__ntime = p["ntime"][0]

        if ftype == 1: # for fermionic system
            # set nfreq to number of fermionic frequency
            self.__nfreq = p["nffrq"][0]
        else:          # for bosonic system
            # set nfreq to number of bosonic frequency
            self.__nfreq = p["nbfrq"][0]

    def __tails(self, freq_mesh, freq_data):
        """ calculate high frequency tails using K. Haule's trick
            freq_mesh : matsubara frequency grid
            freq_data : function on matsubara frequency space
        """
        Sn = 0.0
        Sx = 0.0
        Sy = 0.0

        Sxx = 0.0
        Sxy = 0.0

        for j in range(self.__nfreq - 128, self.__nfreq):
            Sn = Sn + 1.0
            Sx = Sx + 1.0 / freq_mesh[j]**2
            Sy = Sy + freq_data[j].imag * freq_mesh[j]
            Sxx = Sxx + 1.0 / freq_mesh[j]**4
            Sxy = Sxy + freq_data[j].imag * freq_mesh[j] / freq_mesh[j]**2

        return (Sx * Sxy - Sxx * Sy) / (Sn * Sxx - Sx * Sx)

    def freq_to_time_F(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ backward fourier transformation from matsubara frequency to imaginary time,
            this function is only suitable for the fermionic system
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # calculate high frequency tails need to be subtracted
        tail = self.__tails(freq_mesh, freq_data)

        # perform backward fourier transformation
        if fast: # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq
            beta  = self.__beta
            ftail = float(tail) # attention: tail is in numpy.float64 type, while ftail is in float type

            # define c++ code
            code = """
            #include <complex>

            double raux;
            for ( int i = 0; i < ntime; i++ ) {
                raux = 0.0;
                for ( int j = 0; j < nfreq; j++ ) {
                    raux = raux + cos( freq_mesh(j) * time_mesh(i) ) *   real(freq_data(j));
                    raux = raux + sin( freq_mesh(j) * time_mesh(i) ) * ( imag(freq_data(j)) + ftail / freq_mesh(j));
                }
                time_data(i) = 2.0 * raux / beta - 0.5 * ftail;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'beta', 'ftail', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else: # pure python + numpy version
            for i in range(self.__ntime):
                raux = 0.0
                for j in range(self.__nfreq):
                    raux = raux + math.cos( freq_mesh[j] * time_mesh[i] ) *   freq_data[j].real
                    raux = raux + math.sin( freq_mesh[j] * time_mesh[i] ) * ( freq_data[j].imag + tail / freq_mesh[j])
                time_data[i] = 2.0 * raux / self.__beta - 0.5 * tail

        # corrections for the boundary point
        raux = freq_data[self.__nfreq-1].real * freq_mesh[self.__nfreq-1] / math.pi
        time_data[0] = time_data[0] + raux
        time_data[self.__ntime-1] = time_data[self.__ntime-1] - raux

    def time_to_freq_F(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the fermionic system
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # perform fourier transformation
        if fast : # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq

            # define c++ code
            code = """
            #include <complex>

            std::complex<double> ci(0.0, 1.0);
            for ( int i = 0; i < nfreq; i++ ) {
                double sre = 0.0;
                double sim = 0.0;
                for ( int j = 0; j < ntime - 1; j++ ) {
                    double c0 = cos( time_mesh(j)   * freq_mesh(i) );
                    double c1 = cos( time_mesh(j+1) * freq_mesh(i) );
                    double s0 = sin( time_mesh(j)   * freq_mesh(i) );
                    double s1 = sin( time_mesh(j+1) * freq_mesh(i) );
                    double g0 = time_data(j);
                    double g1 = time_data(j+1);
                    double dg = ( g1 - g0 ) / ( time_mesh(j+1) - time_mesh(j) );
                    sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / freq_mesh(i) ) / freq_mesh(i);
                    sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / freq_mesh(i) ) / freq_mesh(i);
                }
                freq_data(i) = sre + ci*sim;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else : # pure python + numpy version
            for i in range(self.__nfreq):
                sre = 0.0
                sim = 0.0
                for j in range(self.__ntime-1):
                    c0 = math.cos( time_mesh[j]   * freq_mesh[i] )
                    c1 = math.cos( time_mesh[j+1] * freq_mesh[i] )
                    s0 = math.sin( time_mesh[j]   * freq_mesh[i] )
                    s1 = math.sin( time_mesh[j+1] * freq_mesh[i] )
                    g0 = time_data[j]
                    g1 = time_data[j+1]
                    dg = ( g1 - g0 ) / ( time_mesh[j+1] - time_mesh[j] )
                    sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / freq_mesh[i] ) / freq_mesh[i]
                    sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / freq_mesh[i] ) / freq_mesh[i]
                freq_data[i] = sre + 1j*sim

    def freq_to_time_B(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ backward fourier transformation from matsubara frequency to imaginary time,
            this function is only suitable for the bosonic function
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # perform backward fourier transformation
        if fast: # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq
            beta  = self.__beta

            # define c++ code
            code = """
            #include <complex>

            std::complex<double> ci(0.0, 1.0);
            std::complex<double> caux;
            for ( int t = 0; t < ntime; t++ ) {
                caux = freq_data(0);
                for ( int f = 1; f < nfreq; f++ ) {
                    caux = caux + 2.0 * freq_data(f) * exp(-ci * freq_mesh(f) * time_mesh(t) );
                }
                time_data(t) = real(caux) / beta;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'beta', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else: # pure python + numpy version
            for t in range(self.__ntime):
                caux = freq_data[0]
                for f in range(1,self.__nfreq):
                    caux = caux + 2.0 * freq_data[f] * numpy.exp(-1j * freq_mesh[f] * time_mesh[t])
                time_data[t] = caux.real / self.__beta

    def time_to_freq_B(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the bosonic function
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        fit = scipy.interpolate.InterpolatedUnivariateSpline(time_mesh, time_data)
        ntime_dense = 4 * self.__ntime # used to build a dense imaginary time mesh
        # denser time mesh
        time_mesh_dense = Mesh.create_time_mesh(self.__beta, ntime_dense)
        # denser time data
        time_data_dense = fit(time_mesh_dense)
        for im in range(self.__nfreq):
            faux = time_data_dense * numpy.exp(1j * freq_mesh[im] * time_mesh_dense)
            # calculate \int^{\beta}_{0} f(\tau) e^{i \omega \tau} d\tau
            # now faux = f(\tau) e^{i \omega \tau}
            freq_data[im] = numpy.trapz(faux, time_mesh_dense)

    def time_to_freq_C(self, time_mesh, time_data, freq_mesh, freq_data):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the bosonic function \chi
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
        """
	# spline interpolation to evaluate the high-frequency tail
        fit = scipy.interpolate.InterpolatedUnivariateSpline(time_mesh, time_data)
        deriv_0 = fit.derivatives(0.0)
        deriv_beta = fit.derivatives(self.__beta)
        c1 = deriv_beta[0] - deriv_0[0]
        c2 = -(deriv_beta[1] - deriv_0[1])
        c3 = deriv_beta[2] - deriv_0[2]
        for im in range(1,self.__nfreq): # im = 0 will cause invalid arithmetic opeartion
            freq_data[im] = c1 / (1j*freq_mesh[im])    \
                          + c2 / (1j*freq_mesh[im])**2 \
                          + c3 / (1j*freq_mesh[im])**3

	# contribution from the rest part
        ntime_dense = 2 * self.__ntime # used to build a dense imaginary time mesh
        time_mesh_dense = Mesh.create_time_mesh(self.__beta, ntime_dense)
        time_data_dense = fit(time_mesh_dense)
        for i in range(ntime_dense):
            time_data_dense[i] = time_data_dense[i] \
                               - c1 * ( time_mesh_dense[i] - self.__beta/2 ) / self.__beta \
                               + c2 * ( time_mesh_dense[i] - self.__beta/2 )**2 / (2.0*self.__beta) - c2 * self.__beta / 24.0 \
                               - c3 * ( time_mesh_dense[i] - self.__beta/2 )**3 / (6.0*self.__beta)

        cutoff = min(self.__ntime, self.__nfreq)
        for im in range(cutoff):
            faux = time_data_dense * numpy.exp(1j * freq_mesh[im] * time_mesh_dense)
            # calculate \int^{\beta}_{0} f(\tau) e^{i \omega \tau} d\tau
            # now faux = f(\tau) e^{i \omega \tau}
            freq_data[im] += numpy.trapz(faux, time_mesh_dense)

class Lattice(object):
    """ class Lattice is used to define a square lattice.
        please pay attention to the definition of k-mesh
    """
    def __init__(self, p):
        """ class constructor for the Lattice
            p : param dictionary
            U   : onsite interaction
            V   : intersite interaction
            t   : hopping parameter
            nkx : k-mesh densities in x-axis
            nky : k-mesh densities in y-axis
        """
        self.__U = p["U"][0]
        self.__V = p["V"][0]
        self.__t = p["t"][0]

        self.__nkx = p["nkx"][0]
        self.__nky = p["nky"][0]

    def create_kxky(self, bz_ty = 1):
        """ build k-mesh along x and y-axis
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        kx = numpy.zeros(self.__nkx, dtype = numpy.float)
        for i in range(self.__nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx[i] = i * math.pi / ( self.__nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx[i] = i * ( math.pi * 2 ) / ( self.__nkx - 1 )

        ky = numpy.zeros(self.__nky, dtype = numpy.float)
        for j in range(self.__nky):
            if bz_ty == 0: # ky \in [0,1\pi]
                ky[j] = j * math.pi / ( self.__nky - 1 )
            else:          # ky \in [0,2\pi]
                ky[j] = j * ( math.pi * 2 ) / ( self.__nky - 1 )

        # return them as a tuple
        return (kx, ky)

    @staticmethod
    def create_kxky(nkx, nky, bz_ty = 1):
        """ build k-mesh along x and y-axis
            nkx : k-mesh densities in x-axis
            nky : k-mesh densities in y-axis
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        kx = numpy.zeros(nkx, dtype = numpy.float)
        for i in range(nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx[i] = i * math.pi / ( nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx[i] = i * ( math.pi * 2 ) / ( nkx - 1 )

        ky = numpy.zeros(nky, dtype = numpy.float)
        for j in range(nky):
            if bz_ty == 0: # ky \in [0,1\pi]
                ky[j] = j * math.pi / ( nky - 1 )
            else:          # ky \in [0,2\pi]
                ky[j] = j * ( math.pi * 2 ) / ( nky - 1 )

        # return them as a tuple
        return (kx, ky)

    def create_ek(self, bz_ty = 1):
        """ build the band dispersion for the square lattice
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        ek = numpy.zeros((self.__nkx, self.__nky), dtype = numpy.float)
        for i in range(self.__nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx = i * math.pi / ( self.__nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx = i * ( math.pi * 2 ) / ( self.__nkx - 1 )
            for j in range(self.__nky):
                if bz_ty == 0: # ky \in [0,1\pi]
                    ky = j * math.pi / ( self.__nky - 1 )
                else:          # ky \in [0,2\pi]
                    ky = j * ( math.pi * 2 ) / ( self.__nky - 1 )
                ek[i,j] = -2 * self.__t * ( math.cos(kx) + math.cos(ky) )
        return ek

    @staticmethod
    def create_ek(t, nkx, nky, bz_ty = 1):
        """ build the band dispersion for the square lattice
            t   : hopping parameter
            nkx : k-mesh densities in x-axis
            nky : k-mesh densities in y-axis
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        ek = numpy.zeros((nkx, nky), dtype = numpy.float)
        for i in range(nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx = i * math.pi / ( nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx = i * ( math.pi * 2 ) / ( nkx - 1 )
            for j in range(nky):
                if bz_ty == 0: # ky \in [0,1\pi]
                    ky = j * math.pi / ( nky - 1 )
                else:          # ky \in [0,2\pi]
                    ky = j * ( math.pi * 2 ) / ( nky - 1 )
                ek[i,j] = -2 * t * ( math.cos(kx) + math.cos(ky) )
        return ek

    def create_vk(self, bz_ty = 1):
        """ build the general interaction
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        vk = numpy.zeros((self.__nkx, self.__nky), dtype = numpy.float)
        for i in range(self.__nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx = i * math.pi / ( self.__nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx = i * ( math.pi * 2 ) / ( self.__nkx - 1 )
            for j in range(self.__nky):
                if bz_ty == 0: # ky \in [0,1\pi]
                    ky = j * math.pi / ( self.__nky - 1 )
                else:          # ky \in [0,2\pi]
                    ky = j * ( math.pi * 2 ) / ( self.__nky - 1 )
                vk[i,j] = self.__U + 2 * self.__V * ( math.cos(kx) + math.cos(ky) )
        return vk

    @staticmethod
    def create_vk(U, V, nkx, nky, bz_ty = 1):
        """ build the general interaction
            U   : onsite interaction
            V   : intersite interaction
            nkx : k-mesh densities in x-axis
            nky : k-mesh densities in y-axis
            bz_ty : type of k-mesh distribution in brillouin zone
        """
        vk = numpy.zeros((nkx, nky), dtype = numpy.float)
        for i in range(nkx):
            if bz_ty == 0: # kx \in [0,1\pi]
                kx = i * math.pi / ( nkx - 1 )
            else:          # kx \in [0,2\pi]
                kx = i * ( math.pi * 2 ) / ( nkx - 1 )
            for j in range(nky):
                if bz_ty == 0: # ky \in [0,1\pi]
                    ky = j * math.pi / ( nky - 1 )
                else:          # ky \in [0,2\pi]
                    ky = j * ( math.pi * 2 ) / ( nky - 1 )
                vk[i,j] = U + 2 * V * ( math.cos(kx) + math.cos(ky) )
        return vk

class Mixer(object):
    """ class Mixer is used to mix the old and new vectors to produce a
        newer one.
    """
    def __init__(self, p):
        """ class constructor for the Mixer
            p : param dictionary
            alpha : mixing parameter
        """
        self.__alpha = p["alpha"][0]

    def linear_mixing(self, vec1, vec2):
        """ perform the simplest linear mixing
            vec1 : old vector on input
            vec2 : new vector on input
        """
        return vec1 * ( 1 - alpha ) + vec2 * alpha

    @staticmethod
    def linear_mixing(alpha, vec1, vec2):
        """ perform the simplest linear mixing
            alpha : mixing parameter
            vec1  : old vector on input
            vec2  : new vector on input
        """
        return vec1 * ( 1 - alpha ) + vec2 * alpha

class AlpsSolver(object):
    """ class AlpsSolver is used to drive the hybridization expansion impurity
        solver contained in the ALPS package.
    """
    def __init__(self, params):
        """ class constructor for the AlpsSolver
            params : user-supply input parameters
        """
        # setup default parameters for the quantum impurity solver
        # note : MAX_TIME IS TOO SMALL. IT SHOULD BE MODIFIED BEFORE CALCULATION.
        self.__params = {
            'N_MEAS'                    : 100        ,
            'SWEEPS'                    : 200000000  ,
            'THERMALIZATION'            : 1000       ,
            'MAX_TIME'                  : 1200       ,
            'N_ORBITALS'                : 2          ,
            'N_TAU'                     : 1023       ,
            'N_nn'                      : 256        ,
            'N_MATSUBARA'               : 512        ,
            'N_LEGENDRE'                : 48         ,
            'N_HISTOGRAM_ORDERS'        : 128        ,
            'MU'                        : 4.0        ,
            'U'                         : 8.0        ,
            'BETA'                      : 100.0      ,
            'MEASURE_time'              : 1          ,
            'MEASURE_freq'              : 1          ,
            'MEASURE_nn'                : 1          ,
            'MEASURE_nnt'               : 16         ,
            'MEASURE_legendre'          : 1          ,
            'TEXT_OUTPUT'               : 1          ,
            'DELTA'                     : "D.dat"    ,
            'SEED'                      : int( time.time() ) + ( pyalps.mpi.rank + 1.0 ) * 1981
        }

        # setup additional keys and values
        for key in params.iterkeys():
            self.__params[key] = params[key]

    def setup(self):
        """ prepare necessary inputs and perform check to ensure everything is OK
        """
        # only master node can do it
        if pyalps.mpi.rank != 0:
            return

        # check whether the delta.dat is available
        if not self.__params.has_key("DELTA"):
            sys.exit("FATAL ERROR : please provide the hybridization function")

        if not os.path.isfile(self.__params["DELTA"]):
            try:
                raise IOError(self.__params["DELTA"])
            except IOError as e:
                print ('IOError : ' + str(e) + " does not exist!")
                sys.exit("FATAL ERROR : please check it and then run this code again")

    def start(self, curr_iter, params = {}):
        """ invoke the quantum impurity solver
            curr_iter : current iteration number
            params    : user-supply input parameters
        """
        # update the parameters dynamically if necessary
        for key in params.iterkeys():
            self.__params[key] = params[key]

        # dump the parameters to external file
        if pyalps.mpi.rank == 0:
            f = open("ctqmc.param." + str(curr_iter), "w")
            for key in self.__params.iterkeys():
                print >> f, key, ":", self.__params[key]
            f.close()

        # mpi barrier
        pyalps.mpi.world.barrier()
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver BEGIN ...")

        # just do it
        pyalps.cthyb.solve(self.__params)

        # mpi barrier
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver END ...")
        pyalps.mpi.world.barrier()

    def final(self, curr_iter):
        """ finialize the quantum impurity solver
            curr_iter : current iteration number
        """
        if pyalps.mpi.rank != 0:
            return

        # rename the output data of the quantum impurity solver
        # only the master node can do it
        try:
            shutil.move("Gt.dat", "t_Gt.dat." + str(curr_iter))
            shutil.move("Gw.dat", "t_Gw.dat." + str(curr_iter))
            shutil.move("Ft.dat", "t_Ft.dat." + str(curr_iter))
            shutil.move("Fw.dat", "t_Fw.dat." + str(curr_iter))
            shutil.move("Sw.dat", "t_Sw.dat." + str(curr_iter))

            shutil.move("nnt.dat", "t_nnt.dat." + str(curr_iter))

            if self.__params["MEASURE_legendre"] == 1: shutil.move("Gtl.dat", "t_Gtl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Gwl.dat", "t_Gwl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Ftl.dat", "t_Ftl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Fwl.dat", "t_Fwl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Swl.dat", "t_Swl.dat." + str(curr_iter))
        except Exception:
            print ("NON-FATAL ERROR : shutil move error")
            pass

        # deal with the other files
        try:
            shutil.move("observables.dat", "observables.dat." + str(curr_iter))
            shutil.move("orders.dat", "orders.dat." + str(curr_iter))
            shutil.move("simulation.dat", "simulation.dat." + str(curr_iter))
        except Exception:
            print ("NON-FATAL ERROR : shutil move error")
            pass

    def extract(self, obs_type, p = 0):
        """ get the simulated results when the quantum impurity solver is finished
            obs_type : used to specify the object we want to extract
            p : param dictionary
        """
        if pyalps.mpi.rank != 0:
            return

        # open results.out.h5, which is in HDF5 format
        # only the master node can do it
        ar = pyalps.ngs.h5ar('results.out.h5','r')

        for case in switch(obs_type):
            if case(1): # extract the green's function in matsubara frequency space
                gw_0 = ar['G_omega/0/mean/value'] # spin up part
                gw_1 = ar['G_omega/1/mean/value'] # spin dn part
                obs_data = numpy.zeros((2, len(gw_0)), numpy.complex)
                obs_data[0] = gw_0
                obs_data[1] = gw_1
                return obs_data
                break

            if case(2): # extract the self-energy function in matsubara frequency space
                sw_0 = ar['S_omega/0/mean/value'] # spin up part
                sw_1 = ar['S_omega/1/mean/value'] # spin dn part
                obs_data = numpy.zeros((2, len(sw_0)), numpy.complex)
                obs_data[0] = sw_0
                obs_data[1] = sw_1
                return obs_data
                break

            if case(3): # extract the green's function in matsubara frequency space
                gw_l_0 = ar['G_l_omega/0/mean/value'] # spin up part
                gw_l_1 = ar['G_l_omega/1/mean/value'] # spin dn part
                obs_data = numpy.zeros((2, len(gw_l_0)), numpy.complex)
                obs_data[0] = gw_l_0
                obs_data[1] = gw_l_1
                return obs_data
                break

            if case(4): # extract the self-energy function in matsubara frequency space
                sw_l_0 = ar['S_l_omega/0/mean/value'] # spin up part
                sw_l_1 = ar['S_l_omega/1/mean/value'] # spin dn part
                obs_data = numpy.zeros((2, len(sw_l_0)), numpy.complex)
                obs_data[0] = sw_l_0
                obs_data[1] = sw_l_1
                return obs_data
                break

            if case(100):
                # get orbital occupation
                occ_0 = ar['simulation/results/density_0/mean/value']
                occ_1 = ar['simulation/results/density_1/mean/value']
                return occ_0 + occ_1
                break

            if case():  # default, could also just omit condition or 'if True'
                sys.exit("FATAL ERROR : this obs_type parameter is not supported yet.")

        # erase the HDF5 object
        del ar

class LocVar(object):
    """ class LocVar is used to store and manipulate the local objects.
    """
    def __init__(self, p):
        """ class constructor for the LocVar
            p     : param dictionary
            nspin : number of (spin-)orbitals
            nffrq : number of fermionic frequencies
        """
        self.__nspin = p["nspin"][0]
        self.__nffrq = p["nffrq"][0]

        # allocate memory for G(i\omega), \Sigma(i\omega), \mathcal{G}(i\omega)
        self.__gloc = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)
        self.__sloc = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)
        self.__bloc = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

        # allocate memory for self-energy function in the Hubbard-I approximation
        # it is spin-degenerated
        self.__shub = numpy.zeros(self.__nffrq, dtype = numpy.complex)

        # allocate memory for G(i\omega) and \Sigma(i\omega) used as backup
        self.__gloc_save = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)
        self.__sloc_save = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

    @property
    def gloc(self):
        """ getter for __gloc property
        """
        return self.__gloc

    @property
    def sloc(self):
        """ getter for __sloc property
        """
        return self.__sloc

    @property
    def bloc(self):
        """ getter for __bloc property
        """
        return self.__bloc

    @gloc.setter
    def gloc(self, gloc):
        """ setter for __gloc property
        """
        self.__gloc = gloc

    @sloc.setter
    def sloc(self, sloc):
        """ setter for __sloc property
        """
        self.__sloc = sloc

    @bloc.setter
    def bloc(self, bloc):
        """ setter for __bloc property
        """
        self.__bloc = bloc

    def reload_sloc(self):
        """ reload local self-energy function from in.sloc file
        """
        # read data from in.sloc file, only master node can do it
        if pyalps.mpi.rank == 0 :
            print ("    reloading sloc data from in.sloc ...")
            if os.path.isfile("in.sloc"):
                f_sig = open("in.sloc", 'r')
                for f in range(self.__nffrq):
                    line = f_sig.readline()
                    spl = map(float, line.split())
                    self.__sloc[0][f] = spl[2] + 1j * spl[3] # build complex numbers
                    self.__sloc[1][f] = spl[4] + 1j * spl[5]
                f_sig.close()
            else:
                sys.exit("FATAL ERROR : in.sloc does not exist")

        # broadcast the self.__sloc data to children processes
        self.__sloc = pyalps.mpi.broadcast(pyalps.mpi.world, self.__sloc, 0)
        pyalps.mpi.world.barrier()

        # copy sloc to sloc_save
        self.__sloc_save = self.__sloc.copy()

    def mix_gloc(self, curr_iter, alpha):
        """ check the convergence of gloc, and then mix it with gloc_save
            curr_iter : current iteration number
            alpha : mixing parameter
        """
        error = (numpy.absolute(self.__gloc - self.__gloc_save)).max()
        # for the first iteration, we do not mix them
        if curr_iter > 1 :
            self.__gloc  = Mixer.linear_mixing(alpha, self.__gloc_save, self.__gloc)
        self.__gloc_save = self.__gloc.copy()
        return error

    def mix_sloc(self, curr_iter, alpha):
        """ check the convergence of sloc, and then mix it with sloc_save
            curr_iter : current iteration number
            alpha : mixing parameter
        """
        error = (numpy.absolute(self.__sloc - self.__sloc_save)).max()
        # for the first iteration, we do not mix them
        if curr_iter > 1 :
            self.__sloc  = Mixer.linear_mixing(alpha, self.__sloc_save, self.__sloc)
        self.__sloc_save = self.__sloc.copy()
        return error

    def cal_ksum_2d(self, nkx, nky, v_lat):
        """ core subroutine used to calculate the k-summing
            note : this function is used internally
            nkx  : number of k-points in x axis
            nky  : number of k-points in y axis
            v_lat : k-dependent function
        """
        # define the boundary
        kx_start = 1
        kx_end = nkx - 1
        ky_start = 1
        ky_end = nky - 1

        # local variable
        v_loc = 0.0

        # for core part
        v_loc = numpy.sum( v_lat[kx_start:kx_end,ky_start:ky_end] )

        # for boundary part
        for kx in range(1,kx_end): # x-axis
            v_loc = v_loc + ( v_lat[kx,0] + v_lat[kx,ky_end] ) * 0.5

        for ky in range(1,ky_end): # y-axis
            v_loc = v_loc + ( v_lat[0,ky] + v_lat[kx_end,ky] ) * 0.5

        # for the four corners
        v_loc = v_loc + ( v_lat[0,0] + v_lat[0,ky_end] + v_lat[kx_end,0] + v_lat[kx_end,ky_end] ) * 0.25

        # return the normalized value
        return v_loc / float(kx_end * ky_end)

    @staticmethod
    def cal_ksum_2d(nkx, nky, v_lat):
        """ core subroutine used to calculate the k-summing
            note: this function is called from outside
            nkx  : number of k-points in x axis
            nky  : number of k-points in y axis
            v_lat : k-dependent function
        """
        # define the boundary
        kx_start = 1
        kx_end = nkx - 1
        ky_start = 1
        ky_end = nky - 1

        # local variable
        v_loc = 0.0

        # for core part
        v_loc = numpy.sum( v_lat[kx_start:kx_end,ky_start:ky_end] )

        # for boundary part
        for kx in range(1,kx_end): # x-axis
            v_loc = v_loc + ( v_lat[kx,0] + v_lat[kx,ky_end] ) * 0.5

        for ky in range(1,ky_end): # y-axis
            v_loc = v_loc + ( v_lat[0,ky] + v_lat[kx_end,ky] ) * 0.5

        # for the four corners
        v_loc = v_loc + ( v_lat[0,0] + v_lat[0,ky_end] + v_lat[kx_end,0] + v_lat[kx_end,ky_end] ) * 0.25

        # return the normalized value
        return v_loc / float(kx_end * ky_end)

    def cal_gloc_by_ksum(self, nkx, nky, glat):
        """ calculate G_{loc} by G_{lat}
            nkx  : number of k-points in x axis
            nky  : number of k-points in y axis
            glat : lattice green's function
        """
        for s in range(self.__nspin):
            for f in range(self.__nffrq):
                self.__gloc[s,f] = self.cal_ksum_2d( nkx, nky, glat[s,f,:,:] )

    def cal_sloc_by_ksum(self, nkx, nky, slat):
        """ calculate \Sigma_{loc} by \Sigma_{lat}
            nkx  : number of k-points in x axis
            nky  : number of k-points in y axis
            slat : lattice self-energy function
        """
        for s in range(self.__nspin):
            for f in range(self.__nffrq):
                self.__sloc[s,f] = self.cal_ksum_2d( nkx, nky, slat[s,f,:,:] )

    def cal_bloc_by_dyson(self):
        """ calculate \mathcal{G}_{loc} by using the dyson equation
        """
        self.__bloc = 1.0 / ( 1.0 / self.__gloc + self.__sloc )

    def cal_sloc_by_dyson(self):
        """ calculate \Sigma_{loc} by using the dyson equation
        """
        self.__sloc = 1.0 / self.__bloc - 1.0 / self.__gloc

    def cal_sloc_by_dmft(self, fmesh, mune, hybf):
        """ calculate \Sigma_{loc} by using the self-consistent equation
            fmesh : fermionic mesh
            mune : chemical potential
            hybf : hybridization function
        """
        for s in range(self.__nspin):
            self.__sloc[s,:] = 1j*fmesh + mune - hybf[s,:] - 1.0 / self.__gloc[s,:]

    def cal_shub_by_hub1(self, p, fmesh):
        """ calculate self-energy function in atomic limit
            p : param dictionary
            fmesh : fermionic mesh
        """
        beta = p["beta"][0]
        mu = p["mune"][0]
        U = p["U"][0]
        N1 = ( math.exp(beta * mu) + math.exp(beta * (2.0*mu - U)) )
        N2 = ( 1.0 + 2.0 * math.exp(beta * mu) + math.exp(beta * (2.0*mu - U)) )
        N3 = N1/N2
        for i in range(self.__nffrq):
            # G is the atomic green's function
            G = ( 1.0 - N3 ) / ( 1j*fmesh[i] + mu ) + N3 / ( 1j*fmesh[i] + mu - U )
            self.__shub[i] = 1j*fmesh[i] + mu - 1.0 / G

    def cal_sloc_by_hub1(self):
        """ combine low-frequency self-energy function with its high-frequency
            part within the hubbard-1 approximation
        """
        # determine boundary
        left = 128 - 20
        right = 128 + 20
        # get the imaginary part of original self-energy function
        im_up = self.__sloc[0].imag
        im_dn = self.__sloc[1].imag
        # smooth them
        for i in range(3):
            Symm.symm_smth_data(im_up[left:right], 3)
            Symm.symm_smth_data(im_dn[left:right], 3)
        # calculation the deviations
        de_up = self.__shub[128].imag - im_up[128]
        de_dn = self.__shub[128].imag - im_dn[128]
        # combine hubbard-1 self-energy function with the low frequency part
        # we only focus on the imaginary part
        for i in range(128, self.__nffrq):
            self.__sloc[0,i] = self.__sloc[0,i].real + (self.__shub[i].imag - de_up) * 1j
            self.__sloc[1,i] = self.__sloc[1,i].real + (self.__shub[i].imag - de_dn) * 1j

    def sym_gloc(self, p):
        """ symmetrize the gloc over spin, we also remove its real part
            p : param dictionary
        """
        if p["sy_pm"][0]: Symm.symm_over_spin(self.__gloc)
        if p["sy_ph"][0]: Symm.symm_over_imag(self.__gloc)

    def sym_sloc(self, p):
        """ symmetrize the sloc over spin
            p : param dictionary
        """
        if p["sy_pm"][0]: Symm.symm_over_spin(self.__sloc)

    def sym_bloc(self, p):
        """ symmetrize the bloc over spin
            p : param dictionary
        """
        if p["sy_pm"][0]: Symm.symm_over_spin(self.__bloc)

class LatVar(object):
    """ class LatVar is used to store and manipulate the lattice K-dependent objects.
    """
    def __init__(self, p):
        """ class constructor for the LocVar
            p     : param dictionary
            nspin : number of (spin-)orbitals
            nffrq : number of fermionic frequencies
            nkx   : number of k-points in x axis
            nky   : number of k-points in y axis
        """
        self.__nspin = p["nspin"][0]
        self.__nffrq = p["nffrq"][0]

        self.__nkx   = p["nkx"][0]
        self.__nky   = p["nky"][0]

        # allocate memory for G_{lat}(i\omega), \Sigma_{lat}(i\omega)
        self.__glat = numpy.zeros((self.__nspin, self.__nffrq, self.__nkx, self.__nky), dtype = numpy.complex)
        self.__slat = numpy.zeros((self.__nspin, self.__nffrq, self.__nkx, self.__nky), dtype = numpy.complex)

    @property
    def glat(self):
        """ getter for __glat property
        """
        return self.__glat

    @property
    def slat(self):
        """ getter for __slat property
        """
        return self.__slat

    @glat.setter
    def glat(self, glat):
        """ setter for __glat property
        """
        self.__glat = glat

    @slat.setter
    def slat(self, slat):
        """ setter for __slat property
        """
        self.__slat = slat

    def sloc_to_slat(self, sloc):
        """ update slat with sloc, ignore the k-dependence
            sloc : k-independent self-energy function
        """
        for ikx in range(self.__nkx):
            for iky in range(self.__nky):
                self.__slat[:,:,ikx,iky] = sloc.copy()

    def cal_glat_by_dyson(self, fmesh, mune, ek):
        """ calculate G_{lat}(i\omega) by the dyson equation
            fmesh : matsubara frequency
            mune  : chemical potential
            ek    : band dispersion
        """
        for s in range(self.__nspin):
            for f in range(self.__nffrq):
                self.__glat[s,f,:,:] = 1.0 / ( 1j * fmesh[f] + mune - ek - self.__slat[s,f,:,:] )

class DMFT(object):
    """ class DMFT is used to defined the self-consistent equation.
    """
    def __init__(self, p):
        """ class constructor for the DMFT
            p     : param dictionary
            beta  : inversion temperature
            nspin : number of (spin-)orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
        """
        self.__beta  = p["beta"][0]

        self.__nspin = p["nspin"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]

        # allocate memory for \Delta(i\omega) and \Delta(\tau)
        self.__hybf = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)
        self.__htau = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)

    @property
    def hybf(self):
        """ getter for __hybf property
        """
        return self.__hybf

    @property
    def htau(self):
        """ getter for __htau property
        """
        return self.__htau

    @hybf.setter
    def hybf(self, hybf):
        """ setter for __hybf property
        """
        self.__hybf = hybf

    @htau.setter
    def htau(self, htau):
        """ setter for __htau property
        """
        self.__htau = htau

    def cal_hybf_by_dmft(self, mune, fmesh, sloc, gloc):
        """ calculate \Delta(i\omega) by using the self-consistent equation
            mune  : chemical potential
            fmesh : fermionic mesh
            sloc  : \Sigma(i\omega)
            gloc  : G(i\omega)
        """
        for s in range(self.__nspin):
            self.__hybf[s,:] = 1j*fmesh + mune - sloc[s,:] - 1.0 / gloc[s,:]
            # suppress the real part of hybridization function
            for f in range(self.__nffrq):
                self.__hybf[s,f] = self.__hybf[s,f] - self.__hybf[s, self.__nffrq - 1].real

    def cal_htau_by_fourier(self, p, tmesh, fmesh):
        """ fourier \Delta(i\omega) to \Delta(\tau)
            p : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
        """
        fourier = Fourier(p, 1) # just for fermionic system
        for s in range(self.__nspin):
            fourier.freq_to_time_F(tmesh, self.__htau[s,:], fmesh, self.__hybf[s,:])
        del fourier

        # check its casuality
        for s in range(self.__nspin):
            for t in range(self.__ntime):
                self.__htau[s,t] = min(self.__htau[s,t], -0.00001)

        # symmetrize over imaginary time
        if p["sy_ph"][0]:
            for s in range(self.__nspin):
                Symm.symm_over_time(self.__htau[s])

        # symmetrize over (spin-)orbital index
        if p["sy_pm"][0]:
            Symm.symm_over_spin(self.__htau)

class Diagram(object):
    @staticmethod
    def cal_frq_to_tau_f(p, tmesh, tau, fmesh, frq):
        """ fourier transformation for fermionic system
        """
        nspin = p["nspin"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]

        # create fourier instance
        fourier = Fourier(p, 1) # for the fermionic system

        # loop over k-mesh
        nkpt = 0 # k-point counter
        for ikx in range(nkx):
            for iky in range(nky):
                nkpt = nkpt + 1
                # distribute the job among children processes
                if pyalps.mpi.rank == nkpt % pyalps.mpi.size:
                    # loop over spin
                    for s in range(nspin):
                        # backward fourier grnf to gtau
                        fourier.freq_to_time_F(tmesh, tau[s,:,ikx,iky], fmesh, frq[s,:,ikx,iky])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_frq_to_tau_f_l(p, tmesh, tau, fmesh, frq):
        """ fourier transformation for fermionic system
        """
        nspin = p["nspin"][0]

        # create fourier instance
        fourier = Fourier(p, 1) # for the fermionic system

        # loop over spin
        for s in range(nspin):
            # backward fourier grnf to gtau
            fourier.freq_to_time_F(tmesh, tau[s,:], fmesh, frq[s,:])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_tau_to_frq_f(p, tmesh, tau, fmesh, frq):
        """ fourier transformation for fermionic system
        """
        nspin = p["nspin"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]

        # create fourier instance
        fourier = Fourier(p, 1) # for the fermionic system

        # loop over k-mesh
        nkpt = 0 # k-point counter
        for ikx in range(nkx):
            for iky in range(nky):
                nkpt = nkpt + 1
                # distribute the job among children processes
                if pyalps.mpi.rank == nkpt % pyalps.mpi.size:
                    # loop over spin
                    for s in range(nspin):
                        fourier.time_to_freq_F(tmesh, tau[s,:,ikx,iky], fmesh, frq[s,:,ikx,iky])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_tau_to_frq_f_l(p, tmesh, tau, fmesh, frq):
        """ fourier transformation for fermionic system
        """
        nspin = p["nspin"][0]

        # create fourier instance
        fourier = Fourier(p, 1) # for the fermionic system

        # loop over spin
        for s in range(nspin):
            fourier.time_to_freq_F(tmesh, tau[s,:], fmesh, frq[s,:])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_frq_to_tau_b(p, tmesh, tau, bmesh, frq):
        """ fourier transformation for bosonic system
        """
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]

        # create fourier instance
        fourier = Fourier(p, 2) # for the bosonic system

        # loop over k-mesh
        nkpt = 0 # k-point counter
        for ikx in range(nkx):
            for iky in range(nky):
                nkpt = nkpt + 1
                # distribute the job among children processes
                if pyalps.mpi.rank == nkpt % pyalps.mpi.size:
                    fourier.freq_to_time_B(tmesh, tau[:,ikx,iky], bmesh, frq[:,ikx,iky])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_frq_to_tau_b_l(p, tmesh, tau, bmesh, frq):
        """ fourier transformation for bosonic system
        """
        # create fourier instance
        fourier = Fourier(p, 2) # for the bosonic system

        fourier.freq_to_time_B(tmesh, tau[:], bmesh, frq[:])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_tau_to_frq_b(p, tmesh, tau, bmesh, frq):
        """ fourier transformation for bosonic system
        """
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]

        # create fourier instance
        fourier = Fourier(p, 2) # for the bosonic system

        # loop over k-mesh
        nkpt = 0 # k-point counter
        for ikx in range(nkx):
            for iky in range(nky):
                nkpt = nkpt + 1
                # distribute the job among children processes
                if pyalps.mpi.rank == nkpt % pyalps.mpi.size:
                    fourier.time_to_freq_B(tmesh, tau[:,ikx,iky], bmesh, frq[:,ikx,iky])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_tau_to_frq_b_l(p, tmesh, tau, bmesh, frq):
        """ fourier transformation for bosonic system
        """
        # create fourier instance
        fourier = Fourier(p, 2) # for the bosonic system

        fourier.time_to_freq_B(tmesh, tau[:], bmesh, frq[:])

        # delete fourier instance
        del fourier

    @staticmethod
    def cal_particle_hole_bubble(p, gtau, bubble):
        """ \Pi(k, \tau) = \sum_q G(q, \tau) G(q - k, -\tau)
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over k-mesh
        for ikx in range(nkx):
            for iky in range(nky):
                # loop over time
                for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                    vt[0] = 0.0

                    # define c++ code
                    code = """
                    // loop over k-mesh
                    for ( int jkx = 0; jkx < nkx; jkx++ ) {
                        double fx = 1.0;
                        // meet the boundary
                        if ( jkx == 0 || jkx == nkx - 1) fx = 0.5;
                        // loop over k-mesh
                        for ( int jky = 0; jky < nky; jky++ ) {
                            double fy = 1.0;
                            // meet the boundary
                            if ( jky == 0 || jky == nky - 1 ) fy = 0.5;

                            int mkx = jkx - ikx;
                            int mky = jky - iky;
                            if ( mkx <= 0 ) mkx = mkx + (nkx - 1);
                            if ( mkx >= nkx - 1 ) mkx = mkx - (nkx - 1);
                            if ( mky <= 0 ) mky = mky + (nky - 1);
                            if ( mky >= nky - 1 ) mky = mky - (nky - 1);
                            // loop over spin
                            for ( int s = 0; s < nspin; s++ ) {
                                if ( t == 0 ) {
                                    vt(0) = vt(0) - fx * fy * gtau(s,0,jkx,jky) * (gtau(s,0,mkx,mky) + 1.0);
                                } else {
                                    vt(0) = vt(0) + fx * fy * gtau(s,t,jkx,jky) * gtau(s,ntime-1-t,mkx,mky);
                                }
                            }
                        }
                    }
                    """

                    # inline c++ code using weave
                    scipy.weave.inline(code, ['ikx', 'iky', 'nkx', 'nky', 'nspin', 'ntime', 'gtau', 't', 'vt'],
                                       type_converters=scipy.weave.converters.blitz)

                    bubble[t,ikx,iky] = -vt / 2.0
        del vt

    @staticmethod
    def cal_particle_hole_bubble_l(p, gtau, bubble):
        """ \Pi(\tau) = G(\tau) G(-\tau)
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over time
        for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
            vt[0] = 0.0

            # define c++ code
            code = """
                // loop over spin
                for ( int s = 0; s < nspin; s++ ) {
                    if ( t == 0 ) {
                        vt(0) = vt(0) - gtau(s,0) * (gtau(s,0) + 1.0);
                    } else {
                        vt(0) = vt(0) + gtau(s,t) * gtau(s,ntime-1-t);
                    }
                }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['nspin', 'ntime', 'gtau', 't', 'vt'],
                               type_converters=scipy.weave.converters.blitz)

            bubble[t] = -vt / 2.0
        del vt

    @staticmethod
    def cal_particle_hole_function(p, afun, bfun, cfun):
        """ C(k, \tau) = - \sum_q A(q, \tau) B(q - k, -\tau), B is spin-dependent
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over k-mesh
        for ikx in range(nkx):
            for iky in range(nky):
                # loop over spin
                for s in range(nspin):
                    # loop over time
                    for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                        vt[0] = 0.0

                        # define c++ code
                        code = """
                        // loop over k-mesh
                        for ( int jkx = 0; jkx < nkx; jkx++ ) {
                            double fx = 1.0;
                            // meet the boundary
                            if ( jkx == 0 || jkx == nkx - 1) fx = 0.5;
                            // loop over k-mesh
                            for ( int jky = 0; jky < nky; jky++ ) {
                                double fy = 1.0;
                                // meet the boundary
                                if ( jky == 0 || jky == nky - 1 ) fy = 0.5;

                                int mkx = jkx - ikx;
                                int mky = jky - iky;
                                if ( mkx <= 0 ) mkx = mkx + (nkx - 1);
                                if ( mkx >= nkx - 1 ) mkx = mkx - (nkx - 1);
                                if ( mky <= 0 ) mky = mky + (nky - 1);
                                if ( mky >= nky - 1 ) mky = mky - (nky - 1);
                                if ( t == 0 ) {
                                    vt(0) = vt(0) + fx * fy * afun(0,jkx,jky) * (bfun(s,0,mkx,mky) + 1.0);
                                } else {
                                    vt(0) = vt(0) - fx * fy * afun(t,jkx,jky) * bfun(s,ntime-1-t,mkx,mky);
                                }
                            }
                        }
                        """

                        # inline c++ code using weave
                        scipy.weave.inline(code, ['ikx', 'iky', 'nkx', 'nky', 's', 't', 'ntime', 'bfun', 'afun', 'vt'],
                                           type_converters=scipy.weave.converters.blitz)

                        cfun[s,t,ikx,iky] = -vt[0]
        del vt

    @staticmethod
    def cal_particle_hole_function_l(p, afun, bfun, cfun):
        """ C(\tau) = - A(\tau) B(-\tau), B is spin-dependent
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over spin
        for s in range(nspin):
            # loop over time
            for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                vt[0] = 0.0

                # define c++ code
                code = """
                    double fx = 1.0;
                    double fy = 1.0;
                    if ( t == 0 ) {
                        vt(0) = vt(0) + afun(0) * (bfun(s,0) + 1.0);
                    } else {
                        vt(0) = vt(0) - afun(t) * bfun(s,ntime-1-t);
                    }
                """

                # inline c++ code using weave
                scipy.weave.inline(code, ['s', 't', 'ntime', 'bfun', 'afun', 'vt'],
                                   type_converters=scipy.weave.converters.blitz)

                cfun[s,t] = -vt[0]
        del vt

    @staticmethod
    def cal_particle_particle_bubble(p, gtau, bubble):
        """ \Pi(k, \tau) = \sum_q G(q, \tau) G(k - q, \tau)
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over k-mesh
        for ikx in range(nkx):
            for iky in range(nky):
                # loop over time
                for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                    vt[0] = 0.0

                    # define c++ code
                    code = """
                    // loop over k-mesh
                    for ( int jkx = 0; jkx < nkx; jkx++ ) {
                        double fx = 1.0;
                        // meet the boundary
                        if ( jkx == 0 || jkx == nkx - 1) fx = 0.5;
                        // loop over k-mesh
                        for ( int jky = 0; jky < nky; jky++ ) {
                            double fy = 1.0;
                            // meet the boundary
                            if ( jky == 0 || jky == nky - 1 ) fy = 0.5;

                            int mkx = ikx - jkx;
                            int mky = iky - jky;
                            if ( mkx <= 0 ) mkx = mkx + (nkx - 1);
                            if ( mkx >= nkx - 1 ) mkx = mkx - (nkx - 1);
                            if ( mky <= 0 ) mky = mky + (nky - 1);
                            if ( mky >= nky - 1 ) mky = mky - (nky - 1);
                            // loop over spin
                            for ( int s = 0; s < nspin; s++ ) {
                                vt(0) = vt(0) + fx * fy * gtau(0,t,jkx,jky) * gtau(0,t,mkx,mky);
                            }
                        }
                    }
                    """

                    # inline c++ code using weave
                    scipy.weave.inline(code, ['ikx', 'iky', 'nkx', 'nky', 'nspin', 'ntime', 'gtau', 't', 'vt'],
                                       type_converters=scipy.weave.converters.blitz)

                    bubble[t,ikx,iky] = vt[0] / 2.0
        del vt

    @staticmethod
    def cal_particle_particle_bubble_l(p, gtau, bubble):
        """ \Pi(\tau) = G(\tau) G(\tau)
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over time
        for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
            vt[0] = 0.0

            # define c++ code
            code = """
                // loop over spin
                for ( int s = 0; s < nspin; s++ ) {
                    vt(0) = vt(0) + gtau(0,t) * gtau(0,t);
                }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['nspin', 'ntime', 'gtau', 't', 'vt'],
                               type_converters=scipy.weave.converters.blitz)

            bubble[t] = vt[0] / 2.0
        del vt

    @staticmethod
    def cal_particle_particle_function(p, afun, bfun, cfun):
        """ C(k, \tau) = - \sum_q A(q, \tau) B(k - q, \tau), A is spin-dependent
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over k-mesh
        for ikx in range(nkx):
            for iky in range(nky):
                # loop over spin
                for s in range(nspin):
                    # loop over time
                    for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                        vt[0] = 0.0

                        # define c++ code
                        code = """
                        // loop over k-mesh
                        for ( int jkx = 0; jkx < nkx; jkx++ ) {
                            double fx = 1.0;
                            // meet the boundary
                            if ( jkx == 0 || jkx == nkx - 1) fx = 0.5;
                            // loop over k-mesh
                            for ( int jky = 0; jky < nky; jky++ ) {
                                double fy = 1.0;
                                // meet the boundary
                                if ( jky == 0 || jky == nky - 1 ) fy = 0.5;

                                int mkx = ikx - jkx;
                                int mky = iky - jky;
                                if ( mkx <= 0 ) mkx = mkx + (nkx - 1);
                                if ( mkx >= nkx - 1 ) mkx = mkx - (nkx - 1);
                                if ( mky <= 0 ) mky = mky + (nky - 1);
                                if ( mky >= nky - 1 ) mky = mky - (nky - 1);
                                vt(0) = vt(0) + fx * fy * afun(s,t,jkx,jky) * bfun(t,mkx,mky);
                            }
                        }
                        """

                        # inline c++ code using weave
                        scipy.weave.inline(code, ['ikx', 'iky', 'nkx', 'nky', 's', 't', 'afun', 'bfun', 'vt'],
                                           type_converters=scipy.weave.converters.blitz)

                        cfun[s,t,ikx,iky] = -vt[0]
        del vt

    @staticmethod
    def cal_particle_particle_function_l(p, afun, bfun, cfun):
        """ C(\tau) = - A(\tau) B(\tau), A is spin-dependent
        """
        nspin = p["nspin"][0]
        ntime = p["ntime"][0]
        vt = numpy.zeros(1, dtype=numpy.float)

        # loop over spin
        for s in range(nspin):
            # loop over time
            for t in range(pyalps.mpi.rank, ntime, pyalps.mpi.size):
                vt[0] = 0.0

                # define c++ code
                code = """
                    vt(0) = vt(0) + afun(s,t) * bfun(t);
                """

                # inline c++ code using weave
                scipy.weave.inline(code, ['s', 't', 'afun', 'bfun', 'vt'],
                                   type_converters=scipy.weave.converters.blitz)

                cfun[s,t] = -vt[0]
        del vt

    @staticmethod
    def cor_particle_hole_bubble(p, cph):
        nbfrq = p["nbfrq"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        U     = p["U"][0]

        eta   = 0.001
        cor   = False  
        for ikx in range(nkx):
            for iky in range(nky):
                for frq in range(nbfrq):
                    re = cph[frq,ikx,iky].real
                    im = cph[frq,ikx,iky].imag
                    if  re > ( 1.0 - eta ) / U:
                        cor = True
        if cor:
            cph[:,:,:] = cph[:,:,:] * 0.5
            if pyalps.mpi.rank == 0:
                print ("    correct the susceptibility ...") 

    @staticmethod
    def cor_particle_hole_bubble_l(p, cph):
        nbfrq = p["nbfrq"][0]
        U     = p["U"][0]

        eta   = 0.001
        cor   = False  
        for frq in range(nbfrq):
            re = cph[frq].real
            im = cph[frq].imag
            if  re > ( 1.0 - eta ) / U:
                cor = True
        if cor:
            cph[:] = cph[:] * 0.5
            if pyalps.mpi.rank == 0:
                print ("    correct the susceptibility ...")

    @staticmethod
    def cal_final_slat(p, s_mbpt, sloc, slat, sdc):
        nspin = p["nspin"][0]
        nffrq = p["nffrq"][0]
        nkx   = p["nkx"][0]
        nky   = p["nky"][0]
        sskip = p["sskip"][0]
        alpha = p["alpha"][0]

        # calculate local contribution to self-energy function
        s_mbpt_loc = numpy.zeros_like(sloc, dtype = numpy.complex)
        for s in range(nspin):
            for f in range(nffrq):
                s_mbpt_loc[s,f] = LocVar.cal_ksum_2d( nkx, nky, s_mbpt[s,f,:,:] )

        # calculate nonlocal contribution to self-energy function
        for s in range(nspin):
            for f in range(nffrq):
                if p["dc_ty"][0] == 1:
                    s_mbpt[s,f,:,:] = s_mbpt[s,f,:,:] - s_mbpt_loc[s,f]
                else:
                    s_mbpt[s,f,:,:] = s_mbpt[s,f,:,:] - sdc[s,f]

        # calculate new self-energy function
        for s in range(nspin):
            for f in range(nffrq):
                if sskip:
                    slat[s,f,:,:] = alpha * (s_mbpt[s,f,:,:] + s_mbpt_loc[s,f]) + (1 - alpha) * slat[s,f,:,:]
                else:
                    slat[s,f,:,:] = alpha * (s_mbpt[s,f,:,:] + sloc[s,f]) + (1 - alpha) * slat[s,f,:,:]

        # dump the local part and double counting term
        if pyalps.mpi.rank == 0:
            if p["dc_ty"][0] > 1:
                print ("    double counting scheme is using ...")
            Dump.dump_freq_data2("MBPT_S_loc.dat", fmesh, s_mbpt_loc)
            Dump.dump_freq_data2("MBPT_S_dc.dat", fmesh, sdc)

        # release memory
        del s_mbpt_loc

class SOPT(object):
    def __init__(self, p):
        """ class constructor for the SOPT
            p     : param dictionary
            U     : Coulomb interaction
            nspin : number of (spin-)orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            nkx   : number of k-points in x axis
            nky   : number of k-points in y axis
        """
        self.__U     = p["U"][0]

        self.__nspin = p["nspin"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

        self.__nkx   = p["nkx"][0]
        self.__nky   = p["nky"][0]

        # allocate memory for G_{lat}(k, \tau), \Sigma_{SOPT}(k, \tau)
        self.__glat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__slat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # plat2 means \chi_{ph}(k, \tau)
        self.__plat2 = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # allocate memory for G_{loc}(\tau), \Sigma_{SOPT}(\tau)
        self.__gloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)
        self.__sloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)

        # ploc2 means \chi_{ph}(\tau)
        self.__ploc2 = numpy.zeros((self.__ntime), dtype = numpy.float)

        # double counting for self-energy function
        self.__sdc   = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

    def cal_glat2_by_fourier(self, p, tmesh, fmesh, glat):
        """ calculate G_{lat}(k, \tau) from G_{lat}(k, i\omega) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            glat  : G_{lat}(k, i\omega)
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : glat // iw -> tau, please waiting ...")

        self.__glat2 = self.__glat2 * 0.0
        Diagram.cal_frq_to_tau_f(p, tmesh, self.__glat2, fmesh, glat)
        pyalps.mpi.world.barrier()
        self.__glat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__glat2, lambda x,y: x + y)

    def cal_plat2_by_sopt(self, p):
        """ calculate the k-dependent bosonic self-energy function by SOPT approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    sopt approximation : plat // time mesh, please waiting ...")

        self.__plat2 = self.__plat2 * 0.0
        Diagram.cal_particle_hole_bubble(p, self.__glat2, self.__plat2)
        pyalps.mpi.world.barrier()
        self.__plat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__plat2, lambda x,y: x + y)
        self.__plat2 = self.__plat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_slat2_by_sopt(self, p):
        """ calculate the k-dependent self-energy function by SOPT approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    sopt approximation : slat // time mesh, please waiting ...")

        self.__slat2 = self.__slat2 * 0.0
        Diagram.cal_particle_particle_function(p, self.__glat2, self.__plat2, self.__slat2)
        pyalps.mpi.world.barrier()
        self.__slat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__slat2, lambda x,y: x + y)
        self.__slat2 = self.__slat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_sdc_by_sopt(self, p, tmesh, fmesh, gloc):
        """ calculate double counting term for self-energy function
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            gloc  : G_{loc}(i\omega)
        """
        self.__gloc2 = self.__gloc2 * 0.0
        Diagram.cal_frq_to_tau_f_l(p, tmesh, self.__gloc2, fmesh, gloc)
        
        self.__ploc2 = self.__ploc2 * 0.0
        Diagram.cal_particle_hole_bubble_l(p, self.__gloc2, self.__ploc2)
        pyalps.mpi.world.barrier()
        self.__ploc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__ploc2, lambda x,y: x + y)

        self.__sloc2 = self.__sloc2 * 0.0
        Diagram.cal_particle_particle_function_l(p, self.__gloc2, self.__ploc2, self.__sloc2)
        pyalps.mpi.world.barrier()
        self.__sloc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__sloc2, lambda x,y: x + y)

        self.__sdc = self.__sdc * 0.0
        Diagram.cal_tau_to_frq_f_l(p, tmesh, self.__sloc2, fmesh, self.__sdc)
        for s in range(self.__nspin):
            self.__sdc[s] = self.__sdc[s] * self.__U * self.__U + self.__U * (self.__gloc2[s,0] + 1.0)

    def cal_new_slat_by_sopt(self, _iter_, p, tmesh, fmesh, gloc, sloc, slat):
        """ calculate nonlocal contribution to \Sigma, and then combine nonlocal
            \Sigma and local \Sigma to produce new lattice \Sigma
            _iter_: current iteration number
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            sloc  : local self-energy function, \Sigma_{loc}
            slat  : lattice self-energy function, \Sigma_{lat}
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : s_sopt // tau -> iw, please waiting ...")

        # build s_sopt: \Sigma_{SOPT}(k, i\omega)
        s_sopt = numpy.zeros_like(slat, dtype = numpy.complex)
        Diagram.cal_tau_to_frq_f(p, tmesh, self.__slat2, fmesh, s_sopt)
        pyalps.mpi.world.barrier()
        s_sopt = pyalps.mpi.all_reduce(pyalps.mpi.world, s_sopt, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print ("    sopt approximation : slat // non-local, please waiting ...")

        # calculate the hartree term
        s_hartree = numpy.zeros(self.__nspin, dtype = numpy.float)
        for s in range(self.__nspin):
            s_hartree[s] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, self.__glat2[s,0,:,:] + 1.0 )
            s_sopt[s] = s_sopt[s] * self.__U * self.__U + self.__U * s_hartree[s]

        # calculate double counting term
        if p["dc_ty"][0] == 3:
            for s in range(self.__nspin):
                for f in range(self.__nffrq):
                    gk = 1.0 / ( 1j * fmesh[f] + p["mune"][0] - ek - s_sopt[s,f,:,:] )
                    gloc[s,f] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, gk )
        self.cal_sdc_by_sopt(p, tmesh, fmesh, gloc)

        # determine new self-energy function, s_sopt and slat will be updated
        Diagram.cal_final_slat(p, s_sopt, sloc, slat, self.__sdc)

        # dump the nonlocal self-energy function to external file
        if pyalps.mpi.rank == 0:
            Dump.dump_kdep_data("SOPT_S_up.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_sopt[0])
            Dump.dump_kdep_data("SOPT_S_dn.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_sopt[1])

        # dump the full self-energy function
        if pyalps.mpi.rank == 0:
            slat.tofile("SOPT_S_lat.dat." + str(_iter_))

        # release memory
        del s_sopt
        del s_hartree

class TMA(object):
    def __init__(self, p):
        """ class constructor for the TMA
            p     : param dictionary
            U     : Coulomb interaction
            nspin : number of (spin-)orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            nkx   : number of k-points in x axis
            nky   : number of k-points in y axis
        """
        self.__U     = p["U"][0]

        self.__nspin = p["nspin"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

        self.__nkx   = p["nkx"][0]
        self.__nky   = p["nky"][0]

        # allocate memory for G_{lat}(k, \tau), \Sigma_{TMA}(k, \tau)
        self.__glat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__slat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # plat2 means \chi_{pp}(k, \tau)
        self.__plat2 = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # vpp_t means T(k, \tau)
        self.__vpp_t = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # vpp_v means T(k, i\nu)
        self.__vpp_v = numpy.zeros((self.__nbfrq, self.__nkx, self.__nky), dtype = numpy.complex)

        # allocate memory for G_{loc}(\tau), \Sigma_{TMA}(\tau)
        self.__gloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)
        self.__sloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)

        # ploc2 means \chi_{pp}(\tau)
        self.__ploc2 = numpy.zeros((self.__ntime), dtype = numpy.float)

        # vqq_t means T(\tau)
        self.__vqq_t = numpy.zeros((self.__ntime), dtype = numpy.float)

        # vqq_v means T(i\nu)
        self.__vqq_v = numpy.zeros((self.__nbfrq), dtype = numpy.complex)

        # double counting for self-energy function
        self.__sdc   = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

    def cal_glat2_by_fourier(self, p, tmesh, fmesh, glat):
        """ calculate G_{lat}(k, \tau) from G_{lat}(k, i\omega) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            glat  : G_{lat}(k, i\omega)
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : glat // iw -> tau, please waiting ...")

        self.__glat2 = self.__glat2 * 0.0
        Diagram.cal_frq_to_tau_f(p, tmesh, self.__glat2, fmesh, glat)
        pyalps.mpi.world.barrier()
        self.__glat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__glat2, lambda x,y: x + y)

    def cal_plat2_by_tma(self, p):
        """ calculate the k-dependent particle-particle bubble function by TMA approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    tma approximation : plat // time mesh, please waiting ...")

        self.__plat2 = self.__plat2 * 0.0
        Diagram.cal_particle_particle_bubble(p, self.__glat2, self.__plat2)
        pyalps.mpi.world.barrier()
        self.__plat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__plat2, lambda x,y: x + y)
        self.__plat2 = self.__plat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_vpp_v_by_tma(self, p, tmesh, bmesh):
        """ calculate T-matrix T(k, i\nu)
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : vppv // tau -> iw, please waiting ...")

        self.__vpp_v = self.__vpp_v * 0.0
        Diagram.cal_tau_to_frq_b(p, tmesh, self.__plat2, bmesh, self.__vpp_v)
        pyalps.mpi.world.barrier()
        self.__vpp_v = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vpp_v, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print "    U chi_pp((pi, pi), 0) :", self.__U * self.__vpp_v[0, (self.__nkx - 1)/2,(self.__nky - 1)/2]
       
        # build vpp_v 
        self.__vpp_v = -self.__vpp_v  * self.__vpp_v * self.__U / (1.0 + self.__U * self.__vpp_v)

    def cal_vpp_t_by_tma(self, p, tmesh, bmesh):
        """ calculate T(k, \tau) from T(k, i\nu) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : vppt // iw -> tau, please waiting ...")

        self.__vpp_t = self.__vpp_t * 0.0
        Diagram.cal_frq_to_tau_b(p, tmesh, self.__vpp_t, bmesh, self.__vpp_v)
        pyalps.mpi.world.barrier()
        self.__vpp_t = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vpp_t, lambda x,y: x + y)

    def cal_slat2_by_tma(self, p):
        """ calculate the k-dependent self-energy function by TMA approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    tma approximation : slat // time mesh, please waiting ...")

        self.__slat2 = self.__slat2 * 0.0
        Diagram.cal_particle_hole_function(p, self.__vpp_t + self.__plat2, self.__glat2, self.__slat2)
        pyalps.mpi.world.barrier()
        self.__slat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__slat2, lambda x,y: x + y)
        self.__slat2 = self.__slat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_sdc_by_tma(self, p, tmesh, fmesh, bmesh, gloc):
        """ calculate double counting term for self-energy function
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            bmesh : bosonic mesh
            gloc  : G_{loc}(i\omega)
        """
        self.__gloc2 = self.__gloc2 * 0.0
        Diagram.cal_frq_to_tau_f_l(p, tmesh, self.__gloc2, fmesh, gloc)

        self.__ploc2 = self.__ploc2 * 0.0
        Diagram.cal_particle_particle_bubble_l(p, self.__gloc2, self.__ploc2)
        pyalps.mpi.world.barrier()
        self.__ploc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__ploc2, lambda x,y: x + y)

        self.__vqq_v = self.__vqq_v * 0.0
        Diagram.cal_tau_to_frq_b_l(p, tmesh, self.__ploc2, bmesh, self.__vqq_v)

        # build vqq_v 
        self.__vqq_v = -self.__vqq_v  * self.__vqq_v * self.__U / (1.0 + self.__U * self.__vqq_v)

        self.__vqq_t = self.__vqq_t * 0.0
        Diagram.cal_frq_to_tau_b_l(p, tmesh, self.__vqq_t, bmesh, self.__vqq_v)

        self.__sloc2 = self.__sloc2 * 0.0
        Diagram.cal_particle_hole_function_l(p, self.__vqq_t + self.__ploc2, self.__gloc2, self.__sloc2)
        pyalps.mpi.world.barrier()
        self.__sloc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__sloc2, lambda x,y: x + y)

        self.__sdc = self.__sdc * 0.0
        Diagram.cal_tau_to_frq_f_l(p, tmesh, self.__sloc2, fmesh, self.__sdc)
        for s in range(self.__nspin):
            self.__sdc[s] = self.__sdc[s] * self.__U * self.__U + self.__U * (self.__gloc2[s,0] + 1.0)

    def cal_new_slat_by_tma(self, _iter_, p, tmesh, fmesh, bmesh, gloc, sloc, slat):
        """ calculate nonlocal contribution to \Sigma, and then combine nonlocal
            \Sigma and local \Sigma to produce new lattice \Sigma
            _iter_: current iteration number
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            sloc  : local self-energy function, \Sigma_{loc}
            slat  : lattice self-energy function, \Sigma_{lat}
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : s_tma // tau -> iw, please waiting ...")

        # build s_tma: \Sigma_{TMA}(k, i\omega)
        s_tma = numpy.zeros_like(slat, dtype = numpy.complex)
        Diagram.cal_tau_to_frq_f(p, tmesh, self.__slat2, fmesh, s_tma)
        pyalps.mpi.world.barrier()
        s_tma = pyalps.mpi.all_reduce(pyalps.mpi.world, s_tma, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print ("    tma approximation : slat // non-local, please waiting ...")

        # calculate the hartree term
        s_hartree = numpy.zeros(self.__nspin, dtype = numpy.float)
        for s in range(self.__nspin):
            s_hartree[s] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, self.__glat2[s,0,:,:] + 1.0 )
            s_tma[s] = s_tma[s] * self.__U * self.__U + self.__U * s_hartree[s]

        # calculate double counting term
        if p["dc_ty"][0] == 3:
            for s in range(self.__nspin):
                for f in range(self.__nffrq):
                    gk = 1.0 / ( 1j * fmesh[f] + p["mune"][0] - ek - s_tma[s,f,:,:] )
                    gloc[s,f] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, gk )
        self.cal_sdc_by_tma(p, tmesh, fmesh, bmesh, gloc)

        # determine new self-energy function, s_tma and slat will be updated
        Diagram.cal_final_slat(p, s_tma, sloc, slat, self.__sdc)

        # dump the nonlocal self-energy function to external file
        if pyalps.mpi.rank == 0:
            Dump.dump_kdep_data("TMA_S_up.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_tma[0])
            Dump.dump_kdep_data("TMA_S_dn.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_tma[1])

        # dump the full self-energy function
        if pyalps.mpi.rank == 0:
            slat.tofile("TMA_S_lat.dat." + str(_iter_))

        # release memory
        del s_tma
        del s_hartree

class FLEX(object):
    """ class FLEX is used to calculate the (bosonic) self-energy function by
        using the FLEX approximation
    """
    def __init__(self, p):
        """ class constructor for the FLEX
            p     : param dictionary
            U     : Coulomb interaction
            nspin : number of (spin-)orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            nkx   : number of k-points in x axis
            nky   : number of k-points in y axis
        """
        self.__U     = p["U"][0]

        self.__nspin = p["nspin"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

        self.__nkx   = p["nkx"][0]
        self.__nky   = p["nky"][0]

        # allocate memory for G_{lat}(k, \tau), \Sigma_{FLEX}(k, \tau)
        self.__glat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__slat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # chiph means \chi_{ph}(k, \tau), and chipp means \chi_{pp}(k, \tau)
        self.__chiph = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__chipp = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # vph_t means V_{ph}(k, \tau), and vpp_t means V_{pp}(k, \tau)
        self.__vph_t = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__vpp_t = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # vph_v means V_{ph}(k, i\nu), and vpp_v means V_{pp}(k, i\nu)
        self.__vph_v = numpy.zeros((self.__nbfrq, self.__nkx, self.__nky), dtype = numpy.complex)
        self.__vpp_v = numpy.zeros((self.__nbfrq, self.__nkx, self.__nky), dtype = numpy.complex)

        # allocate memory for G_{loc}(\tau), \Sigma_{FLEX}(\tau)
        self.__gloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)
        self.__sloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)

        # chiqh means \chi_{ph}(\tau), and chiqq means \chi_{pp}(\tau)
        self.__chiqh = numpy.zeros((self.__ntime), dtype = numpy.float)
        self.__chiqq = numpy.zeros((self.__ntime), dtype = numpy.float)

        # vqh_t means V_{ph}(\tau), and vqq_t means V_{pp}(\tau)
        self.__vqh_t = numpy.zeros((self.__ntime), dtype = numpy.float)
        self.__vqq_t = numpy.zeros((self.__ntime), dtype = numpy.float)

        # vqh_v means V_{ph}(i\nu), and vqq_v means V_{pp}(i\nu)
        self.__vqh_v = numpy.zeros((self.__nbfrq), dtype = numpy.complex)
        self.__vqq_v = numpy.zeros((self.__nbfrq), dtype = numpy.complex)

        # double counting for self-energy function
        self.__sdc   = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

    def cal_glat2_by_fourier(self, p, tmesh, fmesh, glat):
        """ calculate G_{lat}(k, \tau) from G_{lat}(k, i\omega) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            glat  : G_{lat}(k, i\omega)
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : glat // iw -> tau, please waiting ...")

        self.__glat2 = self.__glat2 * 0.0
        Diagram.cal_frq_to_tau_f(p, tmesh, self.__glat2, fmesh, glat)
        pyalps.mpi.world.barrier()
        self.__glat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__glat2, lambda x,y: x + y)

    def cal_chiph_by_flex(self, p):
        """ calculate the k-dependent particle-hole bubble function by FLEX approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    flex approximation : c_ph // time mesh, please waiting ...")

        self.__chiph = self.__chiph * 0.0
        Diagram.cal_particle_hole_bubble(p, self.__glat2, self.__chiph)
        pyalps.mpi.world.barrier()
        self.__chiph = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__chiph, lambda x,y: x + y)
        self.__chiph = -self.__chiph / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_chipp_by_flex(self, p):
        """ calculate the k-dependent particle-particle bubble function by FLEX approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    flex approximation : c_pp // time mesh, please waiting ...")

        self.__chipp = self.__chipp * 0.0
        Diagram.cal_particle_particle_bubble(p, self.__glat2, self.__chipp)
        pyalps.mpi.world.barrier()
        self.__chipp = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__chipp, lambda x,y: x + y)
        self.__chipp = self.__chipp / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_vph_v_by_flex(self, p, tmesh, bmesh):
        """ calculate V_{ph}(k, i\nu)
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : vphv // tau -> iw, please waiting ...")

        self.__vph_v = self.__vph_v * 0.0
        Diagram.cal_tau_to_frq_b(p, tmesh, self.__chiph, bmesh, self.__vph_v)
        pyalps.mpi.world.barrier()
        self.__vph_v = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vph_v, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            Dump.dump_kdep_data("FLEX_chi_ph.dat", self.__nbfrq, self.__nkx, self.__nky, self.__U * self.__vph_v)
            print "    U chi_ph((pi, pi), 0) :", self.__U * self.__vph_v[0, (self.__nkx - 1)/2,(self.__nky - 1)/2]

        # fix vph_v to avoid instability
        Diagram.cor_particle_hole_bubble(p, self.__vph_v)

        if pyalps.mpi.rank == 0:
            print "    U chi_ph((pi, pi), 0) :", self.__U * self.__vph_v[0, (self.__nkx - 1)/2,(self.__nky - 1)/2]

        # build vph_v
        chi_sp = self.__vph_v / (1.0 - self.__U * self.__vph_v) # spin channel
        chi_ch = self.__vph_v / (1.0 + self.__U * self.__vph_v) # charge channel
        chi_00 = self.__vph_v * 1.0                             # GW channel
        self.__vph_v = 1.5 * chi_sp + 0.5 * chi_ch - chi_00
        del chi_sp, chi_ch, chi_00

    def cal_vpp_v_by_flex(self, p, tmesh, bmesh):
        """ calculate V_{pp}(k, i\nu)
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : vppv // tau -> iw, please waiting ...")

        self.__vpp_v = self.__vpp_v * 0.0
        Diagram.cal_tau_to_frq_b(p, tmesh, self.__chipp, bmesh, self.__vpp_v)
        pyalps.mpi.world.barrier()
        self.__vpp_v = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vpp_v, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print "    U chi_pp((pi, pi), 0) :", self.__U * self.__vpp_v[0, (self.__nkx - 1)/2,(self.__nky - 1)/2]

        # build vpp_v
        chi_00 = self.__vpp_v * 1.0
        self.__vpp_v = self.__U * chi_00 * chi_00 / (1.0 + self.__U * chi_00)
        del chi_00

    def cal_vph_t_by_flex(self, p, tmesh, bmesh):
        """ calculate V_{ph}(k, \tau) from V_{ph}(k, i\nu) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : vpht // iw -> tau, please waiting ...")

        self.__vph_t = self.__vph_t * 0.0
        Diagram.cal_frq_to_tau_b(p, tmesh, self.__vph_t, bmesh, self.__vph_v)
        pyalps.mpi.world.barrier()
        self.__vph_t = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vph_t, lambda x,y: x + y)

    def cal_vpp_t_by_flex(self, p, tmesh, bmesh):
        """ calculate V_{pp}(k, \tau) from V_{pp}(k, i\nu) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : vppt // iw -> tau, please waiting ...")

        self.__vpp_t = self.__vpp_t * 0.0
        Diagram.cal_frq_to_tau_b(p, tmesh, self.__vpp_t, bmesh, self.__vpp_v)
        pyalps.mpi.world.barrier()
        self.__vpp_t = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__vpp_t, lambda x,y: x + y)

    def cal_slat2_by_flex(self, p):
        """ calculate the k-dependent self-energy function by FLEX approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    flex approximation : slat // time mesh, please waiting ...")

        # for particle-hole channel
        slat_ph = self.__slat2 * 0.0
        Diagram.cal_particle_particle_function(p, self.__glat2, self.__vph_t, slat_ph)
        slat_ph = -1.0 * pyalps.mpi.all_reduce(pyalps.mpi.world, slat_ph, lambda x,y: x + y)
        slat_ph = slat_ph / ((self.__nkx - 1) * (self.__nky - 1))

        # for particle-particle channel
        slat_pp = self.__slat2 * 0.0
        if p["gw_ty"][0] == 3:
            Diagram.cal_particle_hole_function(p, self.__vpp_t, self.__glat2, slat_pp)
        slat_pp = -1.0 * pyalps.mpi.all_reduce(pyalps.mpi.world, slat_pp, lambda x,y: x + y)
        slat_pp = slat_pp / ((self.__nkx - 1) * (self.__nky - 1))

        # add the contributions from particle-hole and particle-particle channels
        self.__slat2 = slat_pp + slat_ph

    def cal_sdc_by_flex(self, p, tmesh, fmesh, bmesh, gloc):
        """ calculate double counting term for self-energy function
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            bmesh : bosonic mesh
            gloc  : G_{loc}(i\omega)
        """
        self.__gloc2 = self.__gloc2 * 0.0
        Diagram.cal_frq_to_tau_f_l(p, tmesh, self.__gloc2, fmesh, gloc)

        self.__chiqh = self.__chiqh * 0.0
        Diagram.cal_particle_hole_bubble_l(p, self.__gloc2, self.__chiqh)
        pyalps.mpi.world.barrier()
        self.__chiqh = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__chiqh, lambda x,y: x + y)
        self.__chiqh = -self.__chiqh

        self.__chiqq = self.__chiqq * 0.0
        Diagram.cal_particle_particle_bubble_l(p, self.__gloc2, self.__chiqq)
        pyalps.mpi.world.barrier()
        self.__chiqq = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__chiqq, lambda x,y: x + y)

        self.__vqh_v = self.__vqh_v * 0.0
        Diagram.cal_tau_to_frq_b_l(p, tmesh, self.__chiqh, bmesh, self.__vqh_v)

        # fix vqh_v to avoid instability
        Diagram.cor_particle_hole_bubble_l(p, self.__vqh_v)

        # build vqh_v
        chi_sp = self.__vqh_v / (1.0 - self.__U * self.__vqh_v) # spin channel
        chi_ch = self.__vqh_v / (1.0 + self.__U * self.__vqh_v) # charge channel
        chi_00 = self.__vqh_v * 1.0                             # GW channel
        self.__vqh_v = 1.5 * chi_sp + 0.5 * chi_ch - chi_00
        del chi_sp, chi_ch, chi_00

        self.__vqq_v = self.__vqq_v * 0.0
        Diagram.cal_tau_to_frq_b_l(p, tmesh, self.__chiqq, bmesh, self.__vqq_v)

        # build vqq_v
        chi_00 = self.__vqq_v * 1.0
        self.__vqq_v = self.__U * chi_00 * chi_00 / (1.0 + self.__U * chi_00)
        del chi_00

        self.__vqh_t = self.__vqh_t * 0.0
        Diagram.cal_frq_to_tau_b_l(p, tmesh, self.__vqh_t, bmesh, self.__vqh_v)

        self.__vqq_t = self.__vqq_t * 0.0
        Diagram.cal_frq_to_tau_b_l(p, tmesh, self.__vqq_t, bmesh, self.__vqq_v)

        # for particle-hole channel
        slat_qh = self.__sloc2 * 0.0
        Diagram.cal_particle_particle_function_l(p, self.__gloc2, self.__vqh_t, slat_qh)
        slat_qh = -1.0 * pyalps.mpi.all_reduce(pyalps.mpi.world, slat_qh, lambda x,y: x + y)

        # for particle-particle channel
        slat_qq = self.__sloc2 * 0.0
        if p["gw_ty"][0] == 3:
            Diagram.cal_particle_hole_function_l(p, self.__vqq_t, self.__gloc2, slat_qq)
        slat_qq = -1.0 * pyalps.mpi.all_reduce(pyalps.mpi.world, slat_qq, lambda x,y: x + y)

        # add the contributions from particle-hole and particle-particle channels
        self.__sloc2 = slat_qq + slat_qh

        self.__sdc = self.__sdc * 0.0
        Diagram.cal_tau_to_frq_f_l(p, tmesh, self.__sloc2, fmesh, self.__sdc)
        for s in range(self.__nspin):
            self.__sdc[s] = self.__sdc[s] * self.__U * self.__U + self.__U * (self.__gloc2[s,0] + 1.0)

    def cal_new_slat_by_flex(self, _iter_, p, tmesh, fmesh, bmesh, gloc, sloc, slat):
        """ calculate nonlocal contribution to \Sigma, and then combine nonlocal
            \Sigma and local \Sigma to produce new lattice \Sigma
            _iter_: current iteration number
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            sloc  : local self-energy function, \Sigma_{loc}
            slat  : lattice self-energy function, \Sigma_{lat}
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : s_flex // tau -> iw, please waiting ...")

        # build s_flex: \Sigma_{FLEX}(k, i\omega)
        s_flex = numpy.zeros_like(slat, dtype = numpy.complex)
        Diagram.cal_tau_to_frq_f(p, tmesh, self.__slat2, fmesh, s_flex)
        pyalps.mpi.world.barrier()
        s_flex = pyalps.mpi.all_reduce(pyalps.mpi.world, s_flex, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print ("    flex approximation : slat // non-local, please waiting ...")

        # calculate the hartree term
        s_hartree = numpy.zeros(self.__nspin, dtype = numpy.float)
        for s in range(self.__nspin):
            s_hartree[s] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, self.__glat2[s,0,:,:] + 1.0 )
            s_flex[s] = s_flex[s] * self.__U * self.__U + self.__U * s_hartree[s]

        # calculate double counting term
        if p["dc_ty"][0] == 3:
            for s in range(self.__nspin):
                for f in range(self.__nffrq):
                    gk = 1.0 / ( 1j * fmesh[f] + p["mune"][0] - ek - s_flex[s,f,:,:] )
                    gloc[s,f] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, gk )
        self.cal_sdc_by_flex(p, tmesh, fmesh, bmesh, gloc)

        # determine new self-energy function, s_flex and slat will be updated
        Diagram.cal_final_slat(p, s_flex, sloc, slat, self.__sdc)

        # dump the nonlocal self-energy function to external file
        if pyalps.mpi.rank == 0:
            Dump.dump_kdep_data("FLEX_S_up.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_flex[0])
            Dump.dump_kdep_data("FLEX_S_dn.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_flex[1])

        # dump the full self-energy function
        if pyalps.mpi.rank == 0:
            slat.tofile("FLEX_S_lat.dat." + str(_iter_))

        # release memory
        del s_flex
        del s_hartree

class GW(object):
    """ class GW is used to calculate the (bosonic) self-energy function by
        using the GW approximation
    """
    def __init__(self, p):
        """ class constructor for the GW
            p     : param dictionary
            U     : Coulomb interaction
            nspin : number of (spin-)orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            nkx   : number of k-points in x axis
            nky   : number of k-points in y axis
        """
        self.__U     = p["U"][0]

        self.__nspin = p["nspin"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

        self.__nkx   = p["nkx"][0]
        self.__nky   = p["nky"][0]

        # allocate memory for G_{lat}(k, \tau), \Sigma_{GW}(k, \tau)
        self.__glat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__slat2 = numpy.zeros((self.__nspin, self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # allocate memory for W_{lat}(k, \tau), \Pi_{GW}(k, \tau)
        self.__wlat2 = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)
        self.__plat2 = numpy.zeros((self.__ntime, self.__nkx, self.__nky), dtype = numpy.float)

        # allocate memory for W_{lat}(k, i\nu)
        self.__wgw_v = numpy.zeros((self.__nbfrq, self.__nkx, self.__nky), dtype = numpy.complex)

        # allocate memory for G_{loc}(\tau), \Sigma_{GW}(\tau)
        self.__gloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)
        self.__sloc2 = numpy.zeros((self.__nspin, self.__ntime), dtype = numpy.float)

        # allocate memory for W_{loc}(\tau), \Pi_{GW}(\tau)
        self.__wloc2 = numpy.zeros((self.__ntime), dtype = numpy.float)
        self.__ploc2 = numpy.zeros((self.__ntime), dtype = numpy.float)

        # allocate memory for W_{loc}(i\nu)
        self.__wgl_v = numpy.zeros((self.__nbfrq), dtype = numpy.complex)

        # double counting for self-energy function
        self.__sdc   = numpy.zeros((self.__nspin, self.__nffrq), dtype = numpy.complex)

    def cal_glat2_by_fourier(self, p, tmesh, fmesh, glat):
        """ calculate G_{lat}(k, \tau) from G_{lat}(k, i\omega) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            glat  : G_{lat}(k, i\omega)
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : glat // iw -> tau, please waiting ...")

        self.__glat2 = self.__glat2 * 0.0
        Diagram.cal_frq_to_tau_f(p, tmesh, self.__glat2, fmesh, glat)
        pyalps.mpi.world.barrier()
        self.__glat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__glat2, lambda x,y: x + y)

    def cal_plat2_by_gw(self, p):
        """ calculate the k-dependent bosonic self-energy function by GW approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    gw approximation : plat // time mesh, please waiting ...")

        self.__plat2 = self.__plat2 * 0.0
        Diagram.cal_particle_hole_bubble(p, self.__glat2, self.__plat2)
        pyalps.mpi.world.barrier()
        self.__plat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__plat2, lambda x,y: x + y)
        self.__plat2 = self.__plat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_wgw_v_by_gw(self, p, tmesh, bmesh, vk):
        """ calculate W(k, i\nu)
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
            vk    : bare interaction
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : wgwv // tau -> iw, please waiting ...")

        self.__wgw_v = self.__wgw_v * 0.0
        Diagram.cal_tau_to_frq_b(p, tmesh, self.__plat2, bmesh, self.__wgw_v)
        pyalps.mpi.world.barrier()
        self.__wgw_v = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__wgw_v, lambda x,y: x + y)

        col = False
        for n in range(self.__nbfrq):
            if numpy.any( (vk * self.__wgw_v[n,:,:])**2 > 1.0 ):
                col = True
        if col:
            self.__wgw_v = self.__wgw_v * 0.5
            if pyalps.mpi.rank == 0:
                print ("   correct self.__wgw_v ...")
        for b in range(self.__nbfrq):
            self.__wgw_v[b,:,:] = vk * vk * self.__wgw_v[b,:,:] / ( 1.0 - (vk * self.__wgw_v[b,:,:])**2 )

    def cal_wlat2_by_gw(self, p, tmesh, bmesh):
        """ calculate W_{lat}(k, \tau) from W_{lat}(k, i\omega) by using fourier transformation
            p     : param dictionary
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier backward : wlat // iw -> tau, please waiting ...")

        self.__wlat2 = self.__wlat2 * 0.0
        Diagram.cal_frq_to_tau_b(p, tmesh, self.__wlat2, bmesh, self.__wgw_v)
        pyalps.mpi.world.barrier()
        self.__wlat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__wlat2, lambda x,y: x + y)

    def cal_slat2_by_gw(self, p):
        """ calculate the k-dependent self-energy function by GW approximation
            p     : param dictionary
        """
        if pyalps.mpi.rank == 0:
            print ("    gw approximation : slat // time mesh, please waiting ...")

        self.__slat2 = self.__slat2 * 0.0
        Diagram.cal_particle_particle_function(p, self.__glat2, self.__wlat2, self.__slat2)
        pyalps.mpi.world.barrier()
        self.__slat2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__slat2, lambda x,y: x + y)
        self.__slat2 = self.__slat2 / ((self.__nkx - 1) * (self.__nky - 1))

    def cal_sdc_by_gw(self, p, tmesh, fmesh, bmesh, gloc):
        """ calculate double counting term for self-energy function
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            bmesh : bosonic mesh
            gloc  : G_{loc}(i\omega)
        """
        self.__gloc2 = self.__gloc2 * 0.0
        Diagram.cal_frq_to_tau_f_l(p, tmesh, self.__gloc2, fmesh, gloc)

        self.__ploc2 = self.__ploc2 * 0.0
        Diagram.cal_particle_hole_bubble_l(p, self.__gloc2, self.__ploc2)
        pyalps.mpi.world.barrier()
        self.__ploc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__ploc2, lambda x,y: x + y)

        self.__wgl_v = self.__wgl_v * 0.0
        Diagram.cal_tau_to_frq_b_l(p, tmesh, self.__ploc2, bmesh, self.__wgl_v)
        if numpy.any( (self.__U * self.__wgl_v)**2 > 1.0 ):
            self.__wgl_v = self.__wgl_v * 0.5
            if pyalps.mpi.rank == 0:
                print ("   correct self.__wgl_v ...")
        for b in range(self.__nbfrq):
            self.__wgl_v[b] = self.__U * self.__U * self.__wgl_v[b] / ( 1.0 - (self.__U * self.__wgl_v[b])**2 )

        self.__wloc2 = self.__wloc2 * 0.0
        Diagram.cal_frq_to_tau_b_l(p, tmesh, self.__wloc2, bmesh, self.__wgl_v)

        self.__sloc2 = self.__sloc2 * 0.0
        Diagram.cal_particle_particle_function_l(p, self.__gloc2, self.__wloc2, self.__sloc2)
        pyalps.mpi.world.barrier()
        self.__sloc2 = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__sloc2, lambda x,y: x + y)

        self.__sdc = self.__sdc * 0.0
        Diagram.cal_tau_to_frq_f_l(p, tmesh, self.__sloc2, fmesh, self.__sdc)
        for s in range(self.__nspin):
            self.__sdc[s] = self.__sdc[s] + self.__U * (self.__gloc2[s,0] + 1.0)

    def cal_new_slat_by_gw(self, _iter_, p, tmesh, fmesh, bmesh, gloc, sloc, slat):
        """ calculate nonlocal contribution to \Sigma, and then combine nonlocal
            \Sigma and local \Sigma to produce new lattice \Sigma
            _iter_: current iteration number
            p     : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
            sloc  : local self-energy function, \Sigma_{loc}
            slat  : lattice self-energy function, \Sigma_{lat}
        """
        if pyalps.mpi.rank == 0:
            print ("    fourier  forward : s_gw // tau -> iw, please waiting ...")

        # build s_gw: \Sigma_{GW}(k, i\omega)
        s_gw = numpy.zeros_like(slat, dtype = numpy.complex)
        Diagram.cal_tau_to_frq_f(p, tmesh, self.__slat2, fmesh, s_gw)
        pyalps.mpi.world.barrier()
        s_gw = pyalps.mpi.all_reduce(pyalps.mpi.world, s_gw, lambda x,y: x + y)

        if pyalps.mpi.rank == 0:
            print ("    gw approximation : slat // non-local, please waiting ...")

        # calculate the hartree term
        s_hartree = numpy.zeros(self.__nspin, dtype = numpy.float)
        for s in range(self.__nspin):
            s_hartree[s] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, self.__glat2[s,0,:,:] + 1.0 )
            s_gw[s] = s_gw[s] + self.__U * s_hartree[s]

        # calculate double counting term
        if p["dc_ty"][0] == 3:
            for s in range(self.__nspin):
                for f in range(self.__nffrq):
                    gk = 1.0 / ( 1j * fmesh[f] + p["mune"][0] - ek - s_gw[s,f,:,:] )
                    gloc[s,f] = LocVar.cal_ksum_2d( self.__nkx, self.__nky, gk )
        self.cal_sdc_by_gw(p, tmesh, fmesh, bmesh, gloc)

        # determine new self-energy function, s_gw and slat will be updated
        Diagram.cal_final_slat(p, s_gw, sloc, slat, self.__sdc)

        # dump the nonlocal self-energy function to external file
        if pyalps.mpi.rank == 0:
            Dump.dump_kdep_data("GW_S_up.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_gw[0])
            Dump.dump_kdep_data("GW_S_dn.dat." + str(_iter_), self.__nffrq, self.__nkx, self.__nky, s_gw[1])

        # dump the full self-energy function
        if pyalps.mpi.rank == 0:
            slat.tofile("GW_S_lat.dat." + str(_iter_))

        # release memory
        del s_gw
        del s_hartree

if __name__ == '__main__':

    # print the header for the edmft/gw + edmft/flex + edmft code
    gw_dmft_header()

    # setup config parameters
    p = gw_dmft_config()

    # extract some key parameters which are useful for the impurity solver
    s_para = gw_dmft_sconfig(p)

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Create class instance (lat, loc, dmft, gw, flex, solver) ...")
        t_start = time.clock()

    # create lattice instance
    lat = LatVar(p)

    # create local instance
    loc = LocVar(p)

    # create dmft instance
    dmft = DMFT(p)

    # create GW instance if really need
    if p["gw_ty"][0] == 1: gw = GW(p)

    # create FLEX instance if really need
    if p["gw_ty"][0] == 2: flex = FLEX(p)
    if p["gw_ty"][0] == 3: flex = FLEX(p)

    # create TMA instance if really need
    if p["gw_ty"][0] == 4: tma = TMA(p)

    # create SOPT instance if really need
    if p["gw_ty"][0] == 5: sopt = SOPT(p)

    # create quantum impurity solver instance
    if p["alps"][0]: solver = AlpsSolver(s_para)

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Build mesh (tmesh, fmesh, bmesh) ...")
        t_start = time.clock()

    # create imaginary time mesh
    tmesh = Mesh.create_time_mesh(p["beta"][0], p["ntime"][0])

    # create fermionic mesh
    fmesh = Mesh.create_fermion_mesh(p["beta"][0], p["nffrq"][0])

    # create bosonic mesh
    bmesh = Mesh.create_boson_mesh(p["beta"][0], p["nbfrq"][0])

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Build lattice (kx, ky, ek, vk) ...")
        t_start = time.clock()

    # create k-mesh for square lattice
    # note: in fact, we do not need kx and ky
    (kx, ky) = Lattice.create_kxky(p["nkx"][0], p["nky"][0])

    # create band dispersion for square lattice
    ek = Lattice.create_ek(p["t"][0], p["nkx"][0], p["nky"][0])

    # create general interaction for 2d Hubbard model on square lattice
    vk = Lattice.create_vk(p["U"][0], p["V"][0], p["nkx"][0], p["nky"][0])

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Reload iteration variables (sloc, slat) ...")
        t_start = time.clock()

    # if we do not need to start from scratch, just try to initialize
    # loc.sloc from existing data
    if p["start"][0] == False :
        loc.reload_sloc()

    # use loc.sloc to initialize lat.slat
    if p["start"][0] == False :
        lat.sloc_to_slat(loc.sloc)

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    g_conv = False     # convergence flag for gloc, impurity green's function
    s_conv = False     # convergence flag for sloc, self-energy function
    conv_value = 0.002 # convergence condition
    for _iter_ in range (1, p["niter"][0] + 1):

        # print the iteration information
        gw_dmft_iter(p, _iter_)

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build lattice variables (glat) ...")
            t_start = time.clock()

        # build glat from scratch
        lat.cal_glat_by_dyson(fmesh, p["mune"][0], ek)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build local variables (gloc, sloc, bloc) ...")
            t_start = time.clock()

        # build gloc from glat
        loc.cal_gloc_by_ksum(p["nkx"][0], p["nky"][0], lat.glat)

        # build sloc from slat
        loc.cal_sloc_by_ksum(p["nkx"][0], p["nky"][0], lat.slat)

        # build bloc from gloc and sloc
        loc.cal_bloc_by_dyson()

        # dump the data, only the master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data2("c_gloc.dat." + str(_iter_), fmesh, loc.gloc)
            Dump.dump_freq_data2("c_sloc.dat." + str(_iter_), fmesh, loc.sloc)
            Dump.dump_freq_data2("c_bloc.dat." + str(_iter_), fmesh, loc.bloc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build dmft variables (hybf, htau) ...")
            t_start = time.clock()

        if p["sskip"][0]:
            if pyalps.mpi.rank == 0:
                print (">>> SKIP <<<")
        else:
            # build hybf from dmft self-consistent equation
            dmft.cal_hybf_by_dmft(p["mune"][0], fmesh, loc.sloc, loc.gloc)

            # build htau from fourier transformation
            dmft.cal_htau_by_fourier(p, tmesh, fmesh)

            # dump the htau data to disk file, used by ctqmc impurity solver
            if pyalps.mpi.rank == 0:
                if p["alps"][0]: # use cthyb in ALPS as a impurity solver
                    Dump.dump_htau_data("D.dat", tmesh, dmft.htau)

            # dump the data, only master node can do it
            if pyalps.mpi.rank == 0:
                Dump.dump_freq_data2("c_hybf.dat." + str(_iter_), fmesh, dmft.hybf)
                Dump.dump_time_data2("c_htau.dat." + str(_iter_), tmesh, dmft.htau)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Launch quantum impurity solver (alps::cthyb or my::ct-sp) ...")
            t_start = time.clock()

        if p["sskip"][0]:
            if pyalps.mpi.rank == 0:
                print (">>> SKIP <<<")
        else:
            # check the status of quantum impurity solver
            solver.setup()

            # boot the quantum impurity solver
            if p["alps"][0]:
                s_para["MU"] = p["mune"][0]
            solver.start(_iter_, s_para)

            # finalize the quantum impurity solver
            solver.final(_iter_)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Rebuild local variables (gloc, sloc) ...")
            t_start = time.clock()

        if p["sskip"][0]:
            if pyalps.mpi.rank == 0:
                print (">>> SKIP <<<")
        else:
            # extract the measurement results: local green's function
            # only master node can do it
            if p["alps"][0]: loc.gloc = solver.extract(3)
            loc.gloc = pyalps.mpi.broadcast(pyalps.mpi.world, loc.gloc, 0)
            loc.sym_gloc(p)

            # extract the measurement results: self-energy function
            # only master node can do it
            if p["alps"][0]: loc.sloc = solver.extract(4)
            #loc.cal_sloc_by_dmft(fmesh, p["mune"][0], dmft.hybf)
            loc.sloc = pyalps.mpi.broadcast(pyalps.mpi.world, loc.sloc, 0)
            # using the hubbard-I approximation to correct the high-frequency part
            loc.cal_shub_by_hub1(p, fmesh)
            loc.cal_sloc_by_hub1()
            loc.sym_sloc(p)

            # dump the data, only master node can do it
            if pyalps.mpi.rank == 0:
                Dump.dump_freq_data2("s_gloc.dat." + str(_iter_), fmesh, loc.gloc)
                Dump.dump_freq_data2("s_sloc.dat." + str(_iter_), fmesh, loc.sloc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Check and Mix local variables (gloc, sloc) ...")
            t_start = time.clock()

        # mix gloc and check its convergence
        g_err = loc.mix_gloc(_iter_, p["alpha"][0])
        g_conv = False
        if g_err < conv_value : g_conv = True

        # mix sloc and check its convergence
        s_err = loc.mix_sloc(_iter_, p["alpha"][0])
        s_conv = False
        if s_err < conv_value : s_conv = True

        # dump the convergence results
        if pyalps.mpi.rank == 0:
            print ("    curr_iter : %i   g_err : %s   s_err : %s" % (_iter_, g_err, s_err) )

        # dump the data, only master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data2("m_gloc.dat." + str(_iter_), fmesh, loc.gloc)
            Dump.dump_freq_data2("m_sloc.dat." + str(_iter_), fmesh, loc.sloc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Apply self-consistent equation (slat) ...")
            t_start = time.clock()

        if p["gw_ty"][0] == 5 and _iter_ > p["gw_it"][0]: # SOPT + DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    SOPT + DMFT scheme is actived")

            # calculate G(k, \tau)
            sopt.cal_glat2_by_fourier(p, tmesh, fmesh, lat.glat)

            # calculate \Pi_{SOPT}(k, \tau)
            sopt.cal_plat2_by_sopt(p)

            # calculate \Sigma_{SOPT}(k, \tau)
            sopt.cal_slat2_by_sopt(p)

            # calculate new \Sigma(k, i\omega)
            sopt.cal_new_slat_by_sopt(_iter_, p, tmesh, fmesh, loc.gloc, loc.sloc, lat.slat)

        elif p["gw_ty"][0] == 4 and _iter_ > p["gw_it"][0]: # TMA + DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    TMA + DMFT scheme is actived")

            # calculate G(k, \tau)
            tma.cal_glat2_by_fourier(p, tmesh, fmesh, lat.glat)

            # calculate \Pi_{TMA}(k, \tau)
            tma.cal_plat2_by_tma(p)

            # calculate T(k, i\nu)
            tma.cal_vpp_v_by_tma(p, tmesh, bmesh)

            # calculate T(k, \tau)
            tma.cal_vpp_t_by_tma(p, tmesh, bmesh)

            # calculate \Sigma_{TMA}(k, \tau)
            tma.cal_slat2_by_tma(p)

            # calculate new \Sigma(k, i\omega)
            tma.cal_new_slat_by_tma(_iter_, p, tmesh, fmesh, bmesh, loc.gloc, loc.sloc, lat.slat)

        elif p["gw_ty"][0] == 3 and _iter_ > p["gw_it"][0]: # FLEX_pp + DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    FLEX_pp + DMFT scheme is actived")

            # calculate G(k, \tau)
            flex.cal_glat2_by_fourier(p, tmesh, fmesh, lat.glat)

            # calculate \chi_{ph}(k, \tau)
            flex.cal_chiph_by_flex(p)

            # calculate \chi_{pp}(k, \tau)
            flex.cal_chipp_by_flex(p)

            # calculate V_{ph}(k, i\nu)
            flex.cal_vph_v_by_flex(p, tmesh, bmesh)

            # calculate V_{pp}(k, i\nu)
            flex.cal_vpp_v_by_flex(p, tmesh, bmesh)

            # calculate V_{ph}(k, \tau)
            flex.cal_vph_t_by_flex(p, tmesh, bmesh)

            # calculate V_{pp}(k, \tau)
            flex.cal_vpp_t_by_flex(p, tmesh, bmesh)

            # calculate \Sigma_{FLEX}(k, \tau)
            flex.cal_slat2_by_flex(p)

            # calculate new \Sigma(k, i\omega)
            flex.cal_new_slat_by_flex(_iter_, p, tmesh, fmesh, bmesh, loc.gloc, loc.sloc, lat.slat)

        elif p["gw_ty"][0] == 2 and _iter_ > p["gw_it"][0]: # FLEX_ph + DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    FLEX_ph + DMFT scheme is actived")

            # calculate G(k, \tau)
            flex.cal_glat2_by_fourier(p, tmesh, fmesh, lat.glat)

            # calculate \chi_{ph}(k, \tau)
            flex.cal_chiph_by_flex(p)

            # calculate V_{ph}(k, i\nu)
            flex.cal_vph_v_by_flex(p, tmesh, bmesh)

            # calculate V_{ph}(k, \tau)
            flex.cal_vph_t_by_flex(p, tmesh, bmesh)

            # calculate \Sigma_{FLEX}(k, \tau)
            flex.cal_slat2_by_flex(p)

            # calculate new \Sigma(k, i\omega)
            flex.cal_new_slat_by_flex(_iter_, p, tmesh, fmesh, bmesh, loc.gloc, loc.sloc, lat.slat)

        elif p["gw_ty"][0] == 1 and _iter_ > p["gw_it"][0]: # SGW + DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    SGW + DMFT scheme is actived")

            # calculate G(k, \tau)
            gw.cal_glat2_by_fourier(p, tmesh, fmesh, lat.glat)

            # calculate \Pi_{GW}(k, \tau)
            gw.cal_plat2_by_gw(p)

            # calculate W(k, i\nu)
            gw.cal_wgw_v_by_gw(p, tmesh, bmesh, vk)

            # calculate W(k, \tau)
            gw.cal_wlat2_by_gw(p, tmesh, bmesh)

            # calculate \Sigma_{GW}(k, \tau)
            gw.cal_slat2_by_gw(p)

            # calculate new \Sigma(k, i\omega)
            gw.cal_new_slat_by_gw(_iter_, p, tmesh, fmesh, bmesh, loc.gloc, loc.sloc, lat.slat)

        else: # DMFT scheme
            if pyalps.mpi.rank == 0:
                print ("    DMFT scheme is actived : slat = sloc")

            # update \Sigma(k, i\omega) with \Sigma_{loc}(i\omega)
            lat.sloc_to_slat(loc.sloc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Adjust chemical potential (mune) ...")
            t_start = time.clock()

        lat.cal_glat_by_dyson(fmesh, p["mune"][0], ek)
        loc.cal_gloc_by_ksum(p["nkx"][0], p["nky"][0], lat.glat)
        loc.cal_sloc_by_ksum(p["nkx"][0], p["nky"][0], lat.slat)
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data2("c_gloc.data", fmesh, loc.gloc)
            Dump.dump_freq_data2("c_sloc.data", fmesh, loc.sloc)

        # determine the chemical potential
        fourier = Fourier(p, 1)
        gtau = numpy.zeros((p["nspin"][0],p["ntime"][0]), dtype = numpy.float)
        myocc = 0.0
        for s in range(p["nspin"][0]):
            fourier.freq_to_time_F(tmesh, gtau[s,:], fmesh, loc.gloc[s,:])
            myocc = myocc + gtau[s,0] + 1.0
        #p["mune"][0] = p["mune"][0] - (myocc - p["occup"][0])
        if pyalps.mpi.rank == 0:
            print "    curr_occ:", myocc, " curr_mune:", p["mune"][0]

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        # check whether the convergence is obtained
        # note : 20 is minimum iteration loop
        if p["gw_ty"][0] > 0 :
            if ( g_conv or s_conv ) and _iter_ > p["gw_it"][0] + 20 :
                if pyalps.mpi.rank == 0:
                    print(">>> Congratulations! The convergence is obtained.")
                    print ("")
                break
        else :
            if ( g_conv or s_conv ) and _iter_ > 20 :
                if pyalps.mpi.rank == 0:
                    print(">>> Congratulations! The convergence is obtained.")
                    print ("")
                break

    #---------------------------------------------------------------------

    gw_dmft_footer(p["gw_ty"][0])
