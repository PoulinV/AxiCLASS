"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

# JL 14.06.2017: TODO: check whether we should free somewhere the allocated fc.filename and titles, data (4 times)

"""
from math import exp,log, log10
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython
from scipy.interpolate import CubicSpline

# Nils : Added for python 3.x and python 2.x compatibility
import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i



# Import the .pxd containing definitions
from cclassy cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

__version__ = _VERSION_.decode("utf-8")

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass


cdef class Class:
    """
    Class wrapping, creates the glue between C and python

    The actual Class wrapping, the only class we will call from MontePython
    (indeed the only one we will import, with the command:
    from classy import Class

    """
    # List of used structures, defined in the header file. They have to be
    # "cdefined", because they correspond to C structures
    cdef precision pr
    cdef background ba
    cdef thermodynamics th
    cdef perturbations pt
    cdef primordial pm
    cdef fourier fo
    cdef transfer tr
    cdef harmonic hr
    cdef output op
    cdef lensing le
    cdef distortions sd
    cdef file_content fc

    cpdef int computed # Flag to see if classy has already computed with the given pars
    cpdef int allocated # Flag to see if classy structs are allocated already
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    # Defining two new properties to recover, respectively, the parameters used
    # or the age (set after computation). Follow this syntax if you want to
    # access other quantities. Alternatively, you can also define a method, and
    # call it (see _T_cmb method, at the very bottom).
    property pars:
        def __get__(self):
            return self._pars
    property state:
        def __get__(self):
            return True
    property Omega_nu:
        def __get__(self):
            return self.ba.Omega0_ncdm_tot
    property nonlinear_method:
        def __get__(self):
            return self.fo.method

    def set_default(self):
        _pars = {
            "output":"tCl mPk",}
        self.set(**_pars)

    def __cinit__(self, default=False):
        cpdef char* dumc
        self.allocated = False
        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

    def __dealloc__(self):
        if self.allocated:
          self.struct_cleanup()
        self.empty()
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
            free(self.fc.filename)

    # Set up the dictionary
    def set(self,*pars,**kars):
        oldpars = self._pars.copy()
        if len(pars)==1:
            self._pars.update(dict(pars[0]))
        elif len(pars)!=0:
            raise CosmoSevereError("bad call")
        self._pars.update(kars)
        if viewdictitems(self._pars) <= viewdictitems(oldpars):
          return # Don't change the computed states, if the new dict was already contained in the previous dict
        self.computed=False
        return True

    def empty(self):
        self._pars = {}
        self.computed = False

    # Create an equivalent of the parameter file. Non specified values will be
    # taken at their default (in Class)
    def _fillparfile(self):
        cdef char* dumc

        if self.fc.size!=0:
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
        self.fc.size = len(self._pars)
        self.fc.name = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.name!=NULL)

        self.fc.value = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.value!=NULL)

        self.fc.read = <short*> malloc(sizeof(short)*len(self._pars))
        assert(self.fc.read!=NULL)

        # fill parameter file
        i = 0
        for kk in self._pars:

            dumcp = kk.encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk]).strip().encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        if(self.allocated != True):
          return
        if "distortions" in self.ncp:
            distortions_free(&self.sd)
        if "lensing" in self.ncp:
            lensing_free(&self.le)
        if "harmonic" in self.ncp:
            harmonic_free(&self.hr)
        if "transfer" in self.ncp:
            transfer_free(&self.tr)
        if "fourier" in self.ncp:
            fourier_free(&self.fo)
        if "primordial" in self.ncp:
            primordial_free(&self.pm)
        if "perturb" in self.ncp:
            perturbations_free(&self.pt)
        if "thermodynamics" in self.ncp:
            thermodynamics_free(&self.th)
        if "background" in self.ncp:
            background_free(&self.ba)
        self.allocated = False
        self.computed = False

    def _check_task_dependency(self, level):
        """
        Fill the level list with all the needed modules

        .. warning::

            the ordering of modules is obviously dependent on CLASS module order
            in the main.c file. This has to be updated in case of a change to
            this file.

        Parameters
        ----------

        level : list
            list of strings, containing initially only the last module required.
            For instance, to recover all the modules, the input should be
            ['lensing']

        """
        if "distortions" in level:
            if "lensing" not in level:
                level.append("lensing")
        if "lensing" in level:
            if "harmonic" not in level:
                level.append("harmonic")
        if "harmonic" in level:
            if "transfer" not in level:
                level.append("transfer")
        if "transfer" in level:
            if "fourier" not in level:
                level.append("fourier")
        if "fourier" in level:
            if "primordial" not in level:
                level.append("primordial")
        if "primordial" in level:
            if "perturb" not in level:
                level.append("perturb")
        if "perturb" in level:
            if "thermodynamics" not in level:
                level.append("thermodynamics")
        if "thermodynamics" in level:
            if "background" not in level:
                level.append("background")
        if len(level)!=0 :
            if "input" not in level:
                level.append("input")
        return level

    def _pars_check(self, key, value, contains=False, add=""):
        val = ""
        if key in self._pars:
            val = self._pars[key]
            if contains:
                if value in val:
                    return True
            else:
                if value==val:
                    return True
        if add:
            sep = " "
            if isinstance(add,str):
                sep = add

            if contains and val:
                    self.set({key:val+sep+value})
            else:
                self.set({key:value})
            return True
        return False

    def compute(self, level=["distortions"]):
        """
        compute(level=["distortions"])

        Main function, execute all the _init methods for all desired modules.
        This is called in MontePython, and this ensures that the Class instance
        of this class contains all the relevant quantities. Then, one can deduce
        Pk, Cl, etc...

        Parameters
        ----------
        level : list
                list of the last module desired. The internal function
                _check_task_dependency will then add to this list all the
                necessary modules to compute in order to initialize this last
                one. The default last module is "lensing".

        .. warning::

            level default value should be left as an array (it was creating
            problem when casting as a set later on, in _check_task_dependency)

        """
        cdef ErrorMsg errmsg

        # Append to the list level all the modules necessary to compute.
        level = self._check_task_dependency(level)

        # Check if this function ran before (self.computed should be true), and
        # if no other modules were requested, i.e. if self.ncp contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.computed and self.ncp.issuperset(level):
            return

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.computed = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.ncp will contain the list of computed modules (under the form of
        # a set, instead of a python list)
        self.ncp=set()
        # Up until the empty set, all modules are allocated
        # (And then we successively keep track of the ones we allocate additionally)
        self.allocated = True

        # --------------------------------------------------------------------
        # Check the presence for all CLASS modules in the list 'level'. If a
        # module is found in level, executure its "_init" method.
        # --------------------------------------------------------------------
        # The input module should raise a CosmoSevereError, because
        # non-understood parameters asked to the wrapper is a problematic
        # situation.
        if "input" in level:
            if input_read_from_file(&self.fc, &self.pr, &self.ba, &self.th,
                                    &self.pt, &self.tr, &self.pm, &self.hr,
                                    &self.fo, &self.le, &self.sd, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i].decode())
            if problem_flag:
                raise CosmoSevereError(
                    "Class did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a CosmoComputationError
        # with the error message from the faulty module of CLASS.
        if "background" in level:
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.ba.error_message)
            self.ncp.add("background")

        if "thermodynamics" in level:
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.ncp.add("thermodynamics")

        if "perturb" in level:
            if perturbations_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturb")

        if "primordial" in level:
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")

        if "fourier" in level:
            if fourier_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.fo) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.fo.error_message)
            self.ncp.add("fourier")

        if "transfer" in level:
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.fo), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")

        if "harmonic" in level:
            if harmonic_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.fo), &(self.tr),
                            &(self.hr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.hr.error_message)
            self.ncp.add("harmonic")

        if "lensing" in level:
            if lensing_init(&(self.pr), &(self.pt), &(self.hr),
                            &(self.fo), &(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)
            self.ncp.add("lensing")

        if "distortions" in level:
            if distortions_init(&(self.pr), &(self.ba), &(self.th),
                                &(self.pt), &(self.pm), &(self.sd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sd.error_message)
            self.ncp.add("distortions")

        self.computed = True

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return

    def raw_cl(self, lmax=-1, nofail=False):
        """
        raw_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned
                (inclusively). This number will be checked against the maximum l
                at which they were actually computed by CLASS, and an error will
                be raised if the desired lmax is bigger than what CLASS can
                give.
        nofail: bool, optional
                Check and enforce the computation of the harmonic module
                beforehand, with the desired lmax.

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view. It also returns now the
                ell array.
        """
        cdef int lmaxR
        cdef double *rcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.hr.has_tt, self.hr.index_ct_tt, 'tt'),
            (self.hr.has_ee, self.hr.index_ct_ee, 'ee'),
            (self.hr.has_te, self.hr.index_ct_te, 'te'),
            (self.hr.has_bb, self.hr.index_ct_bb, 'bb'),
            (self.hr.has_pp, self.hr.index_ct_pp, 'pp'),
            (self.hr.has_tp, self.hr.index_ct_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No Cl computed")
        lmaxR = self.hr.l_max_tot
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Initialise all the needed Cls arrays
        cl = {}
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)

        # Recover for each ell the information from CLASS
        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, rcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.hr.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = rcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(rcl)
        for index_md in range(self.hr.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def lensed_cl(self, lmax=-1,nofail=False):
        """
        lensed_cl(lmax=-1, nofail=False)

        Return a dictionary of the lensed C_l, computed by CLASS, without the
        density C_ls. They must be asked separately with the function aptly
        named density_cl

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
                Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view.
        """
        cdef int lmaxR
        cdef double *lcl = <double*> calloc(self.le.lt_size,sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.le.has_tt, self.le.index_lt_tt, 'tt'),
            (self.le.has_ee, self.le.index_lt_ee, 'ee'),
            (self.le.has_te, self.le.index_lt_te, 'te'),
            (self.le.has_bb, self.le.index_lt_bb, 'bb'),
            (self.le.has_pp, self.le.index_lt_pp, 'pp'),
            (self.le.has_tp, self.le.index_lt_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No lensed Cl computed")
        lmaxR = self.le.l_lensed_max

        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}
        # Simple Cls, for temperature and polarisation, are not so big in size
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)
        for ell from 2<=ell<lmax+1:
            if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
                raise CosmoSevereError(self.le.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = lcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(lcl)
        return cl

    def density_cl(self, lmax=-1, nofail=False):
        """
        density_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l for the matter

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
            Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : numpy array of numpy.ndarrays
            Array that contains the list (in this order) of self correlation of
            1st bin, then successive correlations (set by non_diagonal) to the
            following bins, then self correlation of 2nd bin, etc. The array
            starts at index_ct_dd.
        """
        cdef int lmaxR
        cdef double *dcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        lmaxR = self.pt.l_lss_max
        has_flags = [
            (self.hr.has_dd, self.hr.index_ct_dd, 'dd'),
            (self.hr.has_td, self.hr.index_ct_td, 'td'),
            (self.hr.has_ll, self.hr.index_ct_ll, 'll'),
            (self.hr.has_dl, self.hr.index_ct_dl, 'dl'),
            (self.hr.has_tl, self.hr.index_ct_tl, 'tl')]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)
                l_max_flag = self.hr.l_max_ct[self.hr.index_md_scalars][index]
                if l_max_flag < lmax and lmax > 0:
                    raise CosmoSevereError(
                        "the %s spectrum was computed until l=%i " % (
                            name.upper(), l_max_flag) +
                        "but you asked a l=%i" % lmax)

        if not spectra:
            raise CosmoSevereError("No density Cl computed")
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_lss",lmax)
                self._pars_check("output",'nCl')
                self.compute()
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}

        # For density Cls, the size is bigger (different redshfit bins)
        # computes the size, given the number of correlations needed to be computed
        size = int((self.hr.d_size*(self.hr.d_size+1)-(self.hr.d_size-self.hr.non_diag)*
                (self.hr.d_size-1-self.hr.non_diag))/2);
        for elem in ['dd', 'll', 'dl']:
            if elem in spectra:
                cl[elem] = {}
                for index in range(size):
                    cl[elem][index] = np.zeros(
                        lmax+1, dtype=np.double)
        for elem in ['td', 'tl']:
            if elem in spectra:
                cl[elem] = np.zeros(lmax+1, dtype=np.double)

        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, dcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.hr.error_message)
            if 'dd' in spectra:
                for index in range(size):
                    cl['dd'][index][ell] = dcl[self.hr.index_ct_dd+index]
            if 'll' in spectra:
                for index in range(size):
                    cl['ll'][index][ell] = dcl[self.hr.index_ct_ll+index]
            if 'dl' in spectra:
                for index in range(size):
                    cl['dl'][index][ell] = dcl[self.hr.index_ct_dl+index]
            if 'td' in spectra:
                cl['td'][ell] = dcl[self.hr.index_ct_td]
            if 'tl' in spectra:
                cl['tl'][ell] = dcl[self.hr.index_ct_tl]
        cl['ell'] = np.arange(lmax+1)

        free(dcl)
        for index_md in range(self.hr.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def z_of_r (self,z_array):
        cdef int last_index=0 #junk
        cdef double * pvecback
        r = np.zeros(len(z_array),'float64')
        dzdr = np.zeros(len(z_array),'float64')

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        i = 0
        for redshift in z_array:

            if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # store r
            r[i] = pvecback[self.ba.index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[self.ba.index_bg_H]

            i += 1

        free(pvecback)
        return r[:],dzdr[:]

    def luminosity_distance(self, z):
        """
        luminosity_distance(z)
        """
        cdef int last_index = 0  # junk
        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba, z, long_info,
                inter_normal, &last_index, pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        lum_distance = pvecback[self.ba.index_bg_lum_distance]
        free(pvecback)
        return lum_distance

    # Gives the total matter pk for a given (k,z)
    def pk(self,double k,double z):
        """
        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.fo.method == nl_none):
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_nonlinear,k,z,self.fo.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return pk

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb(self,double k,double z):
        """
        Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")
        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        if (self.fo.method == nl_none):
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_nonlinear,k,z,self.fo.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return pk_cb

    # Gives the total matter pk for a given (k,z)
    def pk_lin(self,double k,double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_m,&pk_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_lin

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb_lin(self,double k,double z):
        """
        Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_cb,&pk_cb_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_cb_lin

    def get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    # [NS] :: TODO :: check optimization
    def get_pk_all(self, k, z, nonlinear = True, cdmbar = False, z_axis_in_k_arr = 0):
        """ General function to get the P(k,z) for ARBITRARY shapes of k,z
            Additionally, it includes the functionality of selecting wether to use the non-linear parts or not,
            and wether to use the cdm baryon power spectrum only
            For Multi-Dimensional k-arrays, it assumes that one of the dimensions is the z-axis
            This is handled by the z_axis_in_k_arr integer, as described in the source code """
        # z_axis_in_k_arr specifies the integer position of the z_axis wihtin the n-dimensional k_arr
        # Example: 1-d k_array -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the first axis -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the last axis  -> z_axis_in_k_arr = 2


        # 1) Select the correct function
        if nonlinear:
            if cdmbar and not (self.ba.Omega0_ncdm_tot == 0.):
                pk_function = self.pk_cb
            else:
                pk_function = self.pk
        else:
            if cdmbar and not (self.ba.Omega0_ncdm_tot == 0.):
                pk_function = self.pk_cb_lin
            else:
                pk_function = self.pk_lin

        # 2) Check if z array, or z value
        if not isinstance(z,(list,np.ndarray)):
            # Only single z value was passed -> k could still be an array of arbitrary dimension
            if not isinstance(k,(list,np.ndarray)):
                # Only single z value AND only single k value -> just return a value
                # This iterates over ALL remaining dimensions
                return pk_function(k,z)
            else:
                k_arr = np.array(k)
                out_pk = np.empty(np.shape(k_arr))
                iterator = np.nditer(k_arr,flags=['multi_index'])
                while not iterator.finished:
                    out_pk[iterator.multi_index] = pk_function(iterator[0],z)
                    iterator.iternext()
                # This iterates over ALL remaining dimensions
                #for index_k in range(k_arr.shape[-1]):
                #    out_pk[...,index_k] = pk_function(k_arr[...,index_k],z)
                return out_pk

        # 3) An array of z values was passed
        k_arr = np.array(k)
        z_arr = np.array(z)
        if( z_arr.ndim != 1 ):
            raise CosmoSevereError("Can only parse one-dimensional z-arrays, not multi-dimensional")

        if( k_arr.ndim > 1 ):
            # 3.1) If there is a multi-dimensional k-array of EQUAL lenghts
            out_pk = np.empty(np.shape(k_arr))
            # Bring the z_axis to the front
            k_arr = np.moveaxis(k_arr, z_axis_in_k_arr, 0)
            out_pk = np.moveaxis(out_pk, z_axis_in_k_arr, 0)
            if( len(k_arr) != len(z_arr) ):
                raise CosmoSevereError("Mismatching array lengths of the z-array")
            for index_z in range(len(z_arr)):
                iterator = np.nditer(k_arr[index_z],flags=['multi_index'])
                while not iterator.finished:
                    out_pk[index_z][iterator.multi_index] = pk_function(iterator[0],z[index_z])
                    iterator.iternext()
                # This iterates over ALL remaining dimensions
                #for index_k in range(k_arr[index_z].shape[-1]):
                #    out_pk[index_z][...,index_k] = pk_function(k_arr[index_z][...,index_k],z_arr[index_z])
            # Move the z_axis back into position
            k_arr = np.moveaxis(k_arr, 0, z_axis_in_k_arr)
            out_pk = np.moveaxis(out_pk, 0, z_axis_in_k_arr)
            return out_pk
        else:
            # 3.2) If there is a multi-dimensional k-array of UNEQUAL lenghts
            if isinstance(k_arr[0],(list,np.ndarray)):
                # A very special thing happened: The user passed a k array with UNEQUAL lengths of k arrays for each z
                out_pk = []
                for index_z in range(len(z_arr)):
                    k_arr_at_z = np.array(k_arr[index_z])
                    out_pk_at_z = np.empty(np.shape(k_arr_at_z))
                    iterator = np.nditer(k_arr_at_z,flags=['multi_index'])
                    while not iterator.finished:
                        out_pk_at_z[iterator.multi_index] = pk_function(iterator[0],z[index_z])
                        iterator.iternext()
                    #for index_k in range(k_arr_at_z.shape[-1]):
                    #   out_pk_at_z[...,index_k] = pk_function(k_arr_at_z[...,index_k],z_arr[index_z])
                    out_pk.append(out_pk_at_z)
                return out_pk

            # 3.3) If there is a single-dimensional k-array
            # The user passed a z-array, but only a 1-d k array
            # Assume thus, that the k array should be reproduced for all z
            out_pk = np.empty((len(z_arr),len(k_arr)))
            for index_z in range(len(z_arr)):
                for index_k in range(len(k_arr)):
                    out_pk[index_z,index_k] = pk_function(k_arr[index_k],z_arr[index_z])
            return out_pk

    #################################
    # Gives a grid of values of matter and/or cb power spectrum, together with the vectors of corresponding k and z values
    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species = False, h_units=False):
        """
        Returns a grid of matter power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        only_clustering_species : bool
                Whether the returned power spectrum is for galaxy clustering and excludes massive neutrinos, or always includes everything (default)
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        pk : grid of power spectrum values, pk[index_k,index_z]
        k : vector of k values, k[index_k] (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values, z[index_z]
        """

        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.fo.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.fo.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size-self.fo.index_ln_tau_pk),'float64')
        cdef int index_k, index_tau, index_pk
        cdef double z_max_nonlinear, z_max_requested

        # consistency checks

        if self.fo.has_pk_matter == False:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations; add 'mPk' in 'output'")

        if nonlinear == True and self.fo.method == nl_none:
            raise CosmoSevereError("You ask classy to return an array of nonlinear P(k,z) values, but the input parameters sent to CLASS did not require any non-linear P(k,z) calculations; add e.g. 'halofit' or 'HMcode' in 'nonlinear'")

        # check wich type of P(k) to return (total or clustering only, i.e. without massive neutrino contribution)
        if (only_clustering_species == True):
            index_pk = self.fo.index_pk_cluster
        else:
            index_pk = self.fo.index_pk_total

        # get list of redshifts
        # the ln(times) of interest are stored in self.fo.ln_tau[index_tau]
        # with index_tau in the range:
        # self.fo.index_ln_tau_pk <= index_tau < self.fo.ln_tau_size-self.fo.index_ln_tau_pk
        #
        # If the user forgot to ask for z_max_pk>0, then self.fo.ln_tau_size == 1 and self.fo.index_ln_tau_pk==0

        if self.fo.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(self.fo.ln_tau_size-self.fo.index_ln_tau_pk):
                if index_tau == self.fo.ln_tau_size-self.fo.index_ln_tau_pk-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.fo.ln_tau[index_tau+self.fo.index_ln_tau_pk]))

        # check consitency of the list of redshifts

        if nonlinear == True:
            # Check highest value of z at which nl corrections could be computed.
            # In the table tau_sampling it corresponds to index: self.fo.index_tau_min_n
            z_max_nonlinear = self.z_of_tau(self.fo.tau[self.fo.index_tau_min_nl])

            # Check highest value of z in the requested output.
            # In the table tau_sampling it corresponds to index: self.fo.tau_size - self.fo.ln_tau_size + self.fo.index_ln_tau_pk
            z_max_requested = z[0]

            # The first z must be larger or equal to the second one, that is,
            # the first index must be smaller or equal to the second one.
            # If not, raise and error.

            if (self.fo.index_tau_min_nl > self.fo.tau_size - self.fo.ln_tau_size + self.fo.index_ln_tau_pk):
                raise CosmoSevereError("get_pk_and_k_and_z() is trying to return P(k,z) up to z_max=%e (to encompass your requested maximum value of z); but the input parameters sent to CLASS (in particular ppr->nonlinear_min_k_max=%e) were such that the non-linear P(k,z) could only be consistently computed up to z=%e; increase the precision parameter 'nonlinear_min_k_max', or decrease your requested z_max"%(z_max_requested,self.pr.nonlinear_min_k_max,z_max_nonlinear))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        for index_k in xrange(self.fo.k_size_pk):
            k[index_k] = self.fo.k[index_k]*units

        # get P(k,z) array

        for index_tau in xrange(self.fo.ln_tau_size-self.fo.index_ln_tau_pk):
            for index_k in xrange(self.fo.k_size_pk):
                if nonlinear == True:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_nl[index_pk][(index_tau+self.fo.index_ln_tau_pk) * self.fo.k_size + index_k])
                else:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_l[index_pk][(index_tau+self.fo.index_ln_tau_pk) * self.fo.k_size + index_k])

        return pk, k, z

    #################################
    # Gives a grid of each transfer functions arranged in a dictionary, together with the vectors of corresponding k and z values
    def get_transfer_and_k_and_z(self, output_format='class', h_units=False):
        """
        Returns a dictionary of grids of density and/or velocity transfer function values and the z and k at which it has been fully computed.
        Useful for creating interpolators.
        When setting CLASS input parameters, include at least one of 'dTk' (for density transfer functions) or 'vTk' (for velocity transfer functions).
        Following the default output_format='class', all transfer functions will be normalised to 'curvature R=1' at initial time
        (and not 'curvature R = -1/k^2' like in CAMB).
        You may switch to output_format='camb' for the CAMB definition and normalisation of transfer functions.
        (Then, 'dTk' must be in the input: the CAMB format only outputs density transfer functions).
        When sticking to output_format='class', you also get the newtonian metric fluctuations phi and psi.
        If you set the CLASS input parameter 'extra_metric_transfer_functions' to 'yes',
        you get additional metric fluctuations in the synchronous and N-body gauges.

        Parameters
        ----------
        output_format  : ('class' or 'camb')
                Format transfer functions according to CLASS (default) or CAMB
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        tk : dictionary containing all transfer functions.
                For instance, the grid of values of 'd_c' (= delta_cdm) is available in tk['d_c']
                All these grids have indices [index_k,index,z], for instance tk['d_c'][index_k,index,z]
        k : vector of k values (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values
        """
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.pt.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.pt.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef int index_k, index_tau
        cdef char * titles
        cdef double * data
        cdef file_format outf

        # consistency checks
        if (self.pt.has_density_transfers == False) and (self.pt.has_velocity_transfers == False):
            raise CosmoSevereError("You ask classy to return transfer functions, but the input parameters sent to CLASS did not require any T(k,z) calculations; add 'dTk' and/or 'vTk' in 'output'")

        index_md = self.pt.index_md_scalars;

        if (self.pt.ic_size[index_md] > 1):
            raise CosmoSevereError("For simplicity, get_transfer_and_k_and_z() has been written assuming only adiabatic initial conditions. You need to write the generalisation to cases with multiple initial conditions.")

        # check out put format
        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        # check name and number of trnasfer functions computed ghy CLASS

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if perturbations_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)

        # get list of redshifts
        # the ln(times) of interest are stored in self.fo.ln_tau[index_tau]
        # with index_tau in the range:
        # self.fo.index_ln_tau_pk <= index_tau < self.fo.ln_tau_size-self.fo.index_ln_tau_pk
        #
        # If the user forgot to ask for z_max_pk>0, then self.fo.ln_tau_size == 1 and self.fo.index_ln_tau_pk==0

        z_size = self.pt.ln_tau_size-self.pt.index_ln_tau_pk
        if z_size == 1:
            raise CosmoSevereError("You ask classy to return an array of T_x(k,z) values, but the input parameters sent to CLASS did not require any transfer function calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(z_size):
                if index_tau == z_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.pt.ln_tau[index_tau+self.pt.index_ln_tau_pk]))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        k_size = self.pt.k_size_pk
        for index_k in xrange(k_size):
            k[index_k] = self.pt.k[index_md][index_k]*units

        # create output dictionary

        tk = {}
        for index_type,name in enumerate(names):
            if index_type > 0:
                tk[name] = np.zeros((k_size, z_size),'float64')

        # allocate the vector in wich the transfer functions will be stored temporarily for all k and types at a given z
        data = <double*>malloc(sizeof(double)*number_of_titles*self.pt.k_size[index_md])

        # get T(k,z) array

        for index_tau in xrange(z_size):

            if perturbations_output_data_at_index_tau(&self.ba, &self.pt, outf, index_tau+self.pt.index_ln_tau_pk, number_of_titles, data)==_FAILURE_:
                raise CosmoSevereError(self.pt.error_message)

            for index_type,name in enumerate(names):
                if index_type > 0:
                    for index_k in xrange(k_size):
                        tk[name][index_k, index_tau] = data[index_k*number_of_titles+index_type]

        free(data)
        return tk, k, z

    #################################
    # Gives a grid of values of the power spectrum of the quantity [k^2*(phi+psi)/2], where (phi+psi)/2 is the Weyl potential, together with the vectors of corresponding k and z values
    def get_Weyl_pk_and_k_and_z(self, nonlinear=False, h_units=False):
        """
        Returns a grid of Weyl potential (phi+psi) power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.
        Note that this function just calls get_pk_and_k_and_z and corrects the output
        by the ratio of transfer functions [(phi+psi)/d_m]^2.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        Weyl_pk : grid of Weyl potential (phi+psi) spectrum values, Weyl_pk[index_k,index_z]
        k : vector of k values, k[index_k] (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values, z[index_z]
        """
        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk,self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] k4 = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] phi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] psi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] d_m = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] Weyl_pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')

        cdef bint input_nonlinear = nonlinear
        cdef bint input_h_units = h_units

        cdef int index_z

        # get total matter power spectrum
        pk, k, z = self.get_pk_and_k_and_z(nonlinear=input_nonlinear, only_clustering_species = False, h_units=input_h_units)

        # get transfer functions
        tk_and_k_and_z = {}
        tk_and_k_and_z, k, z = self.get_transfer_and_k_and_z(output_format='class',h_units=input_h_units)
        phi = tk_and_k_and_z['phi']
        psi = tk_and_k_and_z['psi']
        d_m = tk_and_k_and_z['d_m']

        # get an array containing k**4 (same for all redshifts)
        for index_z in range(self.fo.ln_tau_size-self.pt.index_ln_tau_pk):
            k4[:,index_z] = k**4

        # rescale total matter power spectrum to get the Weyl power spectrum times k**4
        # (the latter factor is just a convention. Since there is a factor k**2 in the Poisson equation,
        # this rescaled Weyl spectrum has a shape similar to the matter power spectrum).
        Weyl_pk = pk * ((phi+psi)/2./d_m)**2 * k4

        return Weyl_pk, k, z

    #################################
    # Gives sigma(R,z) for a given (R,z)
    def sigma(self,double R,double z, h_units = False):
        """
        Gives sigma (total matter) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z).
         This is unless h_units is set to true, in which case R is the radius in units of Mpc/h,
         and R=8 corresponds to sigma8(z))

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma

        R_in_Mpc = (R if not h_units else R/self.ba.h)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,R_in_Mpc,z,self.fo.index_pk_m,out_sigma,&sigma)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return sigma

    # Gives sigma_cb(R,z) for a given (R,z)
    def sigma_cb(self,double R,double z, h_units = False):
        """
        Gives sigma (cdm+b) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)
         This is unless h_units is set to true, in which case R is the radius in units of Mpc/h,
         and R=8 corresponds to sigma8(z))

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma_cb

        R_in_Mpc = (R if not h_units else R/self.ba.h)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("sigma_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        # If necessary, convert R to units of Mpc
        if h_units:
            R /= self.ba.h

        if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,R,z,self.fo.index_pk_cb,out_sigma,&sigma_cb)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return sigma_cb

    # Gives effective logarithmic slope of P_L(k,z) (total matter) for a given (k,z)
    def pk_tilt(self,double k,double z):
        """
        Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
        (k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

        .. note::

            there is an additional check to verify whether output contains `mPk` and whether k is in the right range

        """
        cdef double pk_tilt

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get pk_tilt(k,z) you must add mPk to the list of outputs.")

        if (k < self.fo.k[1] or k > self.fo.k[self.fo.k_size-2]):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)

        if fourier_pk_tilt_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_total,&pk_tilt)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_tilt

    #calculates the hmcode window_function of the Navarrow Frenk White Profile
    def fourier_hmcode_window_nfw(self,double k,double rv,double c):
        """
        Gives window_nfw for a given wavevector k, virial radius rv and concentration c

        """
        cdef double window_nfw


        if fourier_hmcode_window_nfw(&self.fo,k,rv,c,&window_nfw)==_FAILURE_:
                 raise CosmoSevereError(self.hr.error_message)

        return window_nfw

    def age(self):
        self.compute(["background"])
        return self.ba.age

    def h(self):
        return self.ba.h

    def n_s(self):
        return self.pm.n_s

    def tau_reio(self):
        return self.th.tau_reio


    def fEDE(self):
        return self.ba.f_ede

    def omega_axion(self):
        return self.ba.omega_axion

    def zc(self):
        return 10**self.ba.log10_z_c

    def log10_zc(self):
        return -1*self.ba.log10_axion_ac

    def thetai_scf(self):
        return self.ba.scf_parameters[0]

    def log10_f_axion(self):
        return self.ba.log10_f_axion

    def log10_m_axion(self):
        return self.ba.log10_m_axion

    def phase_shift(self):
        return self.pt.phase_shift

    def amplitude(self):
        return self.pt.amplitude
    def Omega_m(self):
        return self.ba.Omega0_m

    def Omega_r(self):
        return self.ba.Omega0_r

    def theta_s_100(self):
        return 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)

    def theta_star_100(self):
        return 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)

    def Omega_Lambda(self):
        return self.ba.Omega0_lambda

    def Omega_g(self):
        return self.ba.Omega0_g

    def Omega_b(self):
        return self.ba.Omega0_b

    def omega_b(self):
        return self.ba.Omega0_b * self.ba.h * self.ba.h

    def Neff(self):
        return self.ba.Neff


    def rs_rec(self):
        self.compute(["thermodynamics"])
        return self.th.rs_rec

    def rd_rec(self):
        self.compute(["thermodynamics"])
        return self.th.rd_rec

    def da_rec(self):
        self.compute(["thermodynamics"])
        return self.th.da_rec

    def z_rec(self):
        self.compute(["thermodynamics"])
        return self.th.z_rec


    def k_eq(self):
        self.compute(["background"])
        return self.ba.a_eq*self.ba.H_eq

    def z_eq(self):
        self.compute(["background"])
        return 1./self.ba.a_eq-1.

    def sigma8(self):
        self.compute(["fourier"])
        return self.fo.sigma8[self.fo.index_pk_m]

    #def neff(self):
    #    self.compute(["harmonic"])
    #    return self.hr.neff

    def sigma8_cb(self):
        self.compute(["fourier"])
        return self.fo.sigma8[self.fo.index_pk_cb]

    def rs_drag(self):
        self.compute(["thermodynamics"])
        return self.th.rs_d

    def z_reio(self):
        self.compute(["thermodynamics"])
        return self.th.z_reio

    def angular_distance(self, z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A

    #################################
    # Get angular diameter distance of object at z2 as seen by observer at z1,
    def angular_distance_from_to(self, z1, z2):
        """
        angular_distance_from_to(z)

        Return the angular diameter distance of object at z2 as seen by observer at z1,
        that is, sin_K((chi2-chi1)*np.sqrt(|k|))/np.sqrt(|k|)/(1+z2).
        If z1>z2 returns zero.

        Parameters
        ----------
        z1 : float
                Observer redshift
        z2 : float
                Source redshift

        Returns
        -------
        d_A(z1,z2) in Mpc
        """
        cdef int last_index #junk
        cdef double * pvecback

        if z1>=z2:
            return 0.

        else:
            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z1,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # This is the comoving distance to object at z1
            chi1 = pvecback[self.ba.index_bg_conf_distance]

            if background_at_z(&self.ba,z2,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # This is the comoving distance to object at z2
            chi2 = pvecback[self.ba.index_bg_conf_distance]

            free(pvecback)

            if self.ba.K == 0:
                return (chi2-chi1)/(1+z2)
            elif self.ba.K > 0:
                return np.sin(np.sqrt(self.ba.K)*(chi2-chi1))/np.sqrt(self.ba.K)/(1+z2)
            elif self.ba.K < 0:
                return np.sinh(np.sqrt(-self.ba.K)*(chi2-chi1))/np.sqrt(-self.ba.K)/(1+z2)

    def comoving_distance(self, z):
        """
        comoving_distance(z)

        Return the comoving distance

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        r = pvecback[self.ba.index_bg_conf_distance]

        free(pvecback)

        return r

    def scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)

        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return D

    def scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        f = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return f

    #################################
    # gives f(z)*sigma8(z) where f9z) is the scale-inbdependent growth factor
    def scale_independent_f_sigma8(self, z):
        """
        scale_independent_f_sigma8(z)

        Return the scale invariant growth factor f(z) multiplied by sigma8(z)

        Parameters
        ----------
        z : float
                Desired redshift

        Returns
        -------
        f(z)*sigma8(z) (dimensionless)
        """
        return self.scale_independent_growth_factor_f(z)*self.sigma(8,z,h_units=True)

    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    def effective_f_sigma8(self, z, z_step=0.1):
        """
        effective_f_sigma8(z)

        Returns the time derivative of sigma8(z) computed as (d sigma8/d ln a)

        Parameters
        ----------
        z : float
                Desired redshift
        z_step : float
                Default step used for the numerical two-sided derivative. For z < z_step the step is reduced progressively down to z_step/10 while sticking to a double-sided derivative. For z< z_step/10 a single-sided derivative is used instead.

        Returns
        -------
        (d ln sigma8/d ln a)(z) (dimensionless)
        """

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)

        # if possible, use two-sided derivative with default value of z_step
        if z >= z_step:
            return (self.sigma(8,z-z_step,h_units=True)-self.sigma(8,z+z_step,h_units=True))/(2.*z_step)*(1+z)
        else:
            # if z is between z_step/10 and z_step, reduce z_step to z, and then stick to two-sided derivative
            if (z > z_step/10.):
                z_step = z
                return (self.sigma(8,z-z_step,h_units=True)-self.sigma(8,z+z_step,h_units=True))/(2.*z_step)*(1+z)
            # if z is between 0 and z_step/10, use single-sided derivative with z_step/10
            else:
                z_step /=10
                return (self.sigma(8,z,h_units=True)-self.sigma(8,z+z_step,h_units=True))/z_step*(1+z)

    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    def effective_f_sigma8_spline(self, z, Nz=20):
        """
        effective_f_sigma8_spline(z)

        Returns the time derivative of sigma8(z) computed as (d sigma8/d ln a)

        Parameters
        ----------
        z : float
                Desired redshift
        Nz : integer
                Number of values used to spline sigma8(z) in the range [z-0.1,z+0.1]

        Returns
        -------
        (d ln sigma8/d ln a)(z) (dimensionless)
        """

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)
        z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))

        if (z<0) or (z>z_max):
            raise CosmoSevereError("You asked for effective_f_sigma8 at a redshift %e outside of the computed range [0,%e"%(z,z_max))

        if (z<0.1):
            z_array = np.linspace(0, 0.2, num = Nz)
        elif (z<z_max-0.1):
            z_array = np.linspace(z-0.1, z+0.1, num = Nz)
        else:
            z_array = np.linspace(z_max-0.2, z_max, num = Nz)

        sig8_array = np.empty_like(z_array)
        for iz, zval in enumerate(z_array):
            sig8_array[iz] = self.sigma(8,zval,h_units=True)
        return -CubicSpline(z_array,sig8_array).derivative()(z)*(1+z)

   #################################
    def z_of_tau(self, tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        cdef double z
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        z = 1./pvecback[self.ba.index_bg_a]-1.

        free(pvecback)

        return z

    def Hubble(self, z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return H

    def Om_m(self, z):
        """
        Omega_m(z)

        Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_m = pvecback[self.ba.index_bg_Omega_m]

        free(pvecback)

        return Om_m

    def Om_b(self, z):
        """
        Omega_b(z)

        Return the baryon density fraction (exactly, the ratio of quantities defined by Class as
        index_bg_rho_b and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_b = pvecback[self.ba.index_bg_rho_b]/pvecback[self.ba.index_bg_rho_crit]

        free(pvecback)

        return Om_b

    def Om_cdm(self, z):
        """
        Omega_cdm(z)

        Return the cdm density fraction (exactly, the ratio of quantities defined by Class as
        index_bg_rho_cdm and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        if self.ba.has_cdm == True:

            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            Om_cdm = pvecback[self.ba.index_bg_rho_cdm]/pvecback[self.ba.index_bg_rho_crit]

            free(pvecback)

        else:

            Om_cdm = 0.

        return Om_cdm

    def Om_ncdm(self, z):
        """
        Omega_ncdm(z)

        Return the ncdm density fraction (exactly, the ratio of quantities defined by Class as
        Sum_m [ index_bg_rho_ncdm1 + n ], with n=0...N_ncdm-1, and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        if self.ba.has_ncdm == True:

            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            rho_ncdm = 0.
            for n in xrange(self.ba.N_ncdm):
                rho_ncdm += pvecback[self.ba.index_bg_rho_ncdm1+n]
            Om_ncdm = rho_ncdm/pvecback[self.ba.index_bg_rho_crit]

            free(pvecback)

        else:

            Om_ncdm = 0.

        return Om_ncdm

    def ionization_fraction(self, z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        xe = pvecthermo[self.th.index_th_xe]

        free(pvecback)
        free(pvecthermo)

        return xe

    def baryon_temperature(self, z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        Tb = pvecthermo[self.th.index_th_Tb]

        free(pvecback)
        free(pvecthermo)

        return Tb

    def T_cmb(self):
        """
        Return the CMB temperature
        """
        return self.ba.T_cmb

    # redundent with a previous Omega_m() funciton,
    # but we leave it not to break compatibility
    def Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.ba.Omega0_m

    def get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        cdef char *titles
        cdef double* data
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if background_output_titles(&self.ba, titles)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.ba.bt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if background_output_data(&self.ba, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        background = {}

        for i in range(number_of_titles):
            background[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                background[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return background

    def get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if thermodynamics_output_titles(&self.ba, &self.th, titles)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.th.tt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if thermodynamics_output_data(&self.ba, &self.th, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        thermodynamics = {}

        for i in range(number_of_titles):
            thermodynamics[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                thermodynamics[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return thermodynamics

    def get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if primordial_output_titles(&self.pt, &self.pm, titles)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pm.lnk_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if primordial_output_data(&self.pt, &self.pm, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        primordial = {}

        for i in range(number_of_titles):
            primordial[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                primordial[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return primordial

    def get_perturbations(self):
        """
        Return scalar, vector and/or tensor perturbations as arrays for requested
        k-values.

        .. note::

            you need to specify both 'k_output_values', and have some
            perturbations computed, for instance by setting 'output' to 'tCl'.

        Returns
        -------
        perturbations : dict of array of dicts
                perturbations['scalar'] is an array of length 'k_output_values' of
                dictionary containing scalar perturbations.
                Similar for perturbations['vector'] and perturbations['tensor'].
        """

        perturbations = {}

        if self.pt.k_output_values_num<1:
            return perturbations

        cdef:
            Py_ssize_t j
            Py_ssize_t i
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            list names
            list tmparray
            dict tmpdict
            double[:,::1] data_mv
            double ** thedata
            int * thesizes

        # Doing the exact same thing 3 times, for scalar, vector and tensor. Sorry
        # for copy-and-paste here, but I don't know what else to do.
        for mode in ['scalar','vector','tensor']:
            if mode=='scalar' and self.pt.has_scalars:
                thetitles = <bytes> self.pt.scalar_titles
                thedata = self.pt.scalar_perturbations_data
                thesizes = self.pt.size_scalar_perturbation_data
            elif mode=='vector' and self.pt.has_vectors:
                thetitles = <bytes> self.pt.vector_titles
                thedata = self.pt.vector_perturbations_data
                thesizes = self.pt.size_vector_perturbation_data
            elif mode=='tensor' and self.pt.has_tensors:
                thetitles = <bytes> self.pt.tensor_titles
                thedata = self.pt.tensor_perturbations_data
                thesizes = self.pt.size_tensor_perturbation_data
            else:
                continue
            thetitles = str(thetitles.decode())
            names = thetitles.split("\t")[:-1]
            number_of_titles = len(names)
            tmparray = []
            if number_of_titles != 0:
                for j in range(self.pt.k_output_values_num):
                    timesteps = thesizes[j]//number_of_titles
                    tmpdict={}
                    data_mv = <double[:timesteps,:number_of_titles]> thedata[j]
                    for i in range(number_of_titles):
                        tmpdict[names[i]] = np.asarray(data_mv[:,i])
                    tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    def get_transfer(self, z=0., output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'mTk' and/or 'vCTk' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that 0<z<z_pk.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : dictionary containing transfer functions.
        """
        cdef char *titles
        cdef double* data
        cdef char ic_info[1024]
        cdef FileName ic_suffix
        cdef file_format outf

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        index_md = self.pt.index_md_scalars;
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if perturbations_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pt.k_size[index_md]

        size_ic_data = timesteps*number_of_titles;
        ic_num = self.pt.ic_size[index_md];

        data = <double*>malloc(sizeof(double)*size_ic_data*ic_num)

        if perturbations_output_data_at_z(&self.ba, &self.pt, outf, <double> z, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        transfers = {}

        for index_ic in range(ic_num):
            if perturbations_output_firstline_and_ic_suffix(&self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                raise CosmoSevereError(self.pt.error_message)
            ic_key = <bytes> ic_suffix

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.zeros(timesteps, dtype=np.double)
                for index in range(timesteps):
                    tmpdict[names[i]][index] = data[index_ic*size_ic_data+index*number_of_titles+i]

            if ic_num==1:
                transfers = tmpdict
            else:
                transfers[ic_key] = tmpdict

        free(titles)
        free(data)

        return transfers

    def get_current_derived_parameters(self, names):
        """
        get_current_derived_parameters(names)

        Return a dictionary containing an entry for all the names defined in the
        input list.

        Parameters
        ----------
        names : list
                Derived parameters that can be asked from Monte Python, or
                elsewhere.

        Returns
        -------
        derived : dict

        .. warning::

            This method used to take as an argument directly the data class from
            Monte Python. To maintain compatibility with this old feature, a
            check is performed to verify that names is indeed a list. If not, it
            returns a TypeError. The old version of this function, when asked
            with the new argument, will raise an AttributeError.

        """
        if type(names) != type([]):
            raise TypeError("Deprecated")

        derived = {}
        for name in names:
            if name == 'h':
                value = self.ba.h
            elif name == 'H0':
                value = self.ba.h*100
            elif name == 'Omega0_lambda' or name == 'Omega_Lambda':
                value = self.ba.Omega0_lambda
            elif name == 'Omega0_fld':
                value = self.ba.Omega0_fld
            elif name == 'age':
                value = self.ba.age
            elif name == 'conformal_age':
                value = self.ba.conformal_age
            elif name == 'm_ncdm_in_eV':
                value = self.ba.m_ncdm_in_eV[0]
            elif name == 'm_ncdm_tot':
                value = self.ba.Omega0_ncdm_tot*self.ba.h*self.ba.h*93.14
            elif name == 'Omega0_ncdm_tot':
                value = self.ba.Omega0_ncdm_tot               
            elif name == 'Neff':
                value = self.ba.Neff
            elif name == 'Omega_axion_ac':
                value = self.ba.Omega_axion_ac
            elif name == 'log10_axion_ac':
                value = self.ba.log10_axion_ac
            elif name == 'f_axion':
                value = self.ba.f_axion
            elif name == 'omega_axion':
                value = self.ba.omega_axion
            elif name == 'm_axion':
                value = self.ba.m_scf
            elif name == 'log10_m_axion':
                value = self.ba.log10_m_axion
            elif name == 'log10_f_axion':
                value = self.ba.log10_f_axion
            elif name == 'f_ede':
                value = self.ba.f_ede
            elif name == 'log10_z_c':
                value = self.ba.log10_z_c
            elif name == 'phi_scf_c':
                value = self.ba.phi_scf_c
            elif name == 'phi_ini_scf':
                value = self.ba.phi_ini_scf
            elif name == 'V0_phi2n':
                value = self.ba.V0_phi2n
            elif name == 'f_ede_peak':
                value = self.ba.f_ede_peak
            elif name == 'a_peak':
                value = self.ba.a_peak
            elif name == 'Omega0_axion':
                value = self.ba.Omega0_axion
            elif name == 'Omega_m':
                value = self.ba.Omega0_m
            elif name == 'omega_m':
                value = self.ba.Omega0_m*self.ba.h**2
            elif name == 'xi_idr':
                value = self.ba.T_idr/self.ba.T_cmb
            elif name == 'N_dg':
                value = self.ba.Omega0_idr/self.ba.Omega0_g*8./7.*pow(11./4.,4./3.)
            elif name == 'Gamma_0_nadm':
                value = self.th.a_idm_dr*(4./3.)*(self.ba.h*self.ba.h*self.ba.Omega0_idr)
            elif name == 'a_dark':
                value = self.th.a_idm_dr
            elif name == 'tau_reio':
                value = self.th.tau_reio
            elif name == 'z_reio':
                value = self.th.z_reio
            elif name == 'z_rec':
                value = self.th.z_rec
            elif name == 'tau_rec':
                value = self.th.tau_rec
            elif name == 'rs_rec':
                value = self.th.rs_rec
            elif name == 'rd_rec':
                value = self.th.rd_rec
            elif name == 'rs_rec_h':
                value = self.th.rs_rec*self.ba.h
            elif name == 'ds_rec':
                value = self.th.ds_rec
            elif name == 'ds_rec_h':
                value = self.th.ds_rec*self.ba.h
            elif name == 'ra_rec':
                value = self.th.da_rec*(1.+self.th.z_rec)
            elif name == 'ra_rec_h':
                value = self.th.da_rec*(1.+self.th.z_rec)*self.ba.h
            elif name == 'da_rec':
                value = self.th.da_rec
            elif name == 'da_rec_h':
                value = self.th.da_rec*self.ba.h
            elif name == 'z_star':
                value = self.th.z_star
            elif name == 'tau_star':
                value = self.th.tau_star
            elif name == 'rs_star':
                value = self.th.rs_star
            elif name == 'ds_star':
                value = self.th.ds_star
            elif name == 'ra_star':
                value = self.th.ra_star
            elif name == 'da_star':
                value = self.th.da_star
            elif name == 'rd_star':
                value = self.th.rd_star
            elif name == 'z_d':
                value = self.th.z_d
            elif name == 'tau_d':
                value = self.th.tau_d
            elif name == 'ds_d':
                value = self.th.ds_d
            elif name == 'ds_d_h':
                value = self.th.ds_d*self.ba.h
            elif name == 'rs_d':
                value = self.th.rs_d
            elif name == 'rs_d_h':
                value = self.th.rs_d*self.ba.h
            elif name == '100*theta_s':
                value = 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == '100*theta_d':
                value = 100.*self.th.rd_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == '100*theta_star':
                value = 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)
            elif name == 'YHe':
                value = self.th.YHe
            elif name == 'n_e':
                value = self.th.n_e
            elif name == 'A_s':
                value = self.pm.A_s
            elif name == 'ln10^{10}A_s':
                value = log(1.e10*self.pm.A_s)
            elif name == 'n_s':
                value = self.pm.n_s
            elif name == 'alpha_s':
                value = self.pm.alpha_s
            elif name == 'beta_s':
                value = self.pm.beta_s
            elif name == 'r':
                # This is at the pivot scale
                value = self.pm.r
            elif name == 'r_0002':
                # at k_pivot = 0.002/Mpc
                value = self.pm.r*(0.002/self.pm.k_pivot)**(
                    self.pm.n_t-self.pm.n_s-1+0.5*self.pm.alpha_s*log(
                        0.002/self.pm.k_pivot))
            elif name == 'n_t':
                value = self.pm.n_t
            elif name == 'alpha_t':
                value = self.pm.alpha_t
            elif name == 'V_0':
                value = self.pm.V0
            elif name == 'V_1':
                value = self.pm.V1
            elif name == 'V_2':
                value = self.pm.V2
            elif name == 'V_3':
                value = self.pm.V3
            elif name == 'V_4':
                value = self.pm.V4
            elif name == 'epsilon_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                value = eps1*((1.-eps1/3.+eps2/6.)/(1.-eps1/3.))**2
            elif name == 'eta_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = (2.*eps1-eps2/2.-2./3.*eps1**2+5./6.*eps1*eps2-eps2**2/12.-eps23/6.)/(1.-eps1/3.)
            elif name == 'ksi_V^2':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = 2.*(1.-eps1/3.+eps2/6.)*(2.*eps1**2-3./2.*eps1*eps2+eps23/4.)/(1.-eps1/3.)**2
            elif name == 'exp_m_2_tau_As':
                value = exp(-2.*self.th.tau_reio)*self.pm.A_s
            elif name == 'phi_min':
                value = self.pm.phi_min
            elif name == 'phi_max':
                value = self.pm.phi_max
            elif name == 'sigma8':
                value = self.fo.sigma8[self.fo.index_pk_m]
            elif name == 'sigma8_cb':
                value = self.fo.sigma8[self.fo.index_pk_cb]
            elif name == 'sigma12':
                value = self.fo.sigma12
            elif name == 'k_eq':
                value = self.ba.a_eq*self.ba.H_eq
            elif name == 'n_L_lya':
                value = self.fo.n_L_lya
            elif name == 'Delta_Lsquared_lya':
                value = self.fo.Delta_Lsquared_lya
            elif name == 'kp_lya':
                value = self.fo.kp_lya
                value = self.fo.sigma8[self.fo.index_pk_cb]
            elif name == 'k_eq':
                value = self.ba.a_eq*self.ba.H_eq
            elif name == 'a_eq':
                value = self.ba.a_eq
            elif name == 'z_eq':
                value = 1./self.ba.a_eq-1.
            elif name == 'H_eq':
                value = self.ba.H_eq
            elif name == 'tau_eq':
                value = self.ba.tau_eq
            elif name == 'g_sd':
                value = self.sd.sd_parameter_table[0]
            elif name == 'y_sd':
                value = self.sd.sd_parameter_table[1]
            elif name == 'mu_sd':
                value = self.sd.sd_parameter_table[2]
            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value
        return derived

    def nonlinear_scale(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_scale(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        for index_z in range(z_size):
            if fourier_k_nl_at_z(&self.ba,&self.fo,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return k_nl

    def nonlinear_scale_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """

make        nonlinear_scale_cb(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size

        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use pk, rather than pk_cb."
                )
        for index_z in range(z_size):
            if fourier_k_nl_at_z(&self.ba,&self.fo,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return k_nl_cb

    def fourier_hmcode_sigma8(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigma8_at_z(&self.ba,&self.fo,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_8

    def fourier_hmcode_sigma8_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigma8_at_z(&self.ba,&self.fo,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_8_cb

    def fourier_hmcode_sigmadisp(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp

    def fourier_hmcode_sigmadisp_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_cb

    def fourier_hmcode_sigmadisp100(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp100_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_100

    def fourier_hmcode_sigmadisp100_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp100_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_100_cb

    def fourier_hmcode_sigmaprime(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmaprime_at_z(&self.ba,&self.fo,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_prime

    def fourier_hmcode_sigmaprime_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmaprime_at_z(&self.ba,&self.fo,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_prime_cb

    def __call__(self, ctx):
        """
        Function to interface with CosmoHammer

        Parameters
        ----------
        ctx : context
                Contains several dictionaries storing data and cosmological
                information

        """
        data = ctx.get('data')  # recover data from the context

        # If the module has already been called once, clean-up
        if self.state:
            self.struct_cleanup()

        # Set the module to the current values
        self.set(data.cosmo_arguments)
        self.compute(["lensing"])

        # Compute the derived paramter value and store them
        params = ctx.getData()
        self.get_current_derived_parameters(
            data.get_mcmc_parameters(['derived']))
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)

    def get_pk_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk_cb

    def Omega0_k(self):
        """ Curvature contribution """
        return self.ba.Omega0_k

    def Omega0_cdm(self):
        return self.ba.Omega0_cdm

    def spectral_distortion_amplitudes(self):
        if self.sd.type_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_type_amps = np.zeros(self.sd.type_size,'float64')
        for i in range(self.sd.type_size):
          sd_type_amps[i] = self.sd.sd_parameter_table[i]
        return sd_type_amps

    def spectral_distortion(self):
        if self.sd.x_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_amp = np.zeros(self.sd.x_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sd_nu = np.zeros(self.sd.x_size,'float64')
        for i in range(self.sd.x_size):
          sd_amp[i] = self.sd.DI[i]*self.sd.DI_units*1.e26
          sd_nu[i] = self.sd.x[i]*self.sd.x_to_nu
        return sd_nu,sd_amp


    def get_sources(self):
        """
        Return the source functions for all k, tau in the grid.

        Returns
        -------
        sources : dictionary containing source functions.
        k_array : numpy array containing k values.
        tau_array: numpy array containing tau values.
        """
        sources = {}

        cdef:
            int index_k, index_tau, i_index_type;
            int index_type;
            int index_md = self.pt.index_md_scalars;
            double * k = self.pt.k[index_md];
            double * tau = self.pt.tau_sampling;
            int index_ic = self.pt.index_ic_ad;
            int k_size = self.pt.k_size[index_md];
            int tau_size = self.pt.tau_size;
            int tp_size = self.pt.tp_size[index_md];
            double *** sources_ptr = self.pt.sources;
            double [:,:] tmparray = np.zeros((k_size, tau_size)) ;
            double [:] k_array = np.zeros(k_size);
            double [:] tau_array = np.zeros(tau_size);

        names = []

        for index_k in range(k_size):
            k_array[index_k] = k[index_k]
        for index_tau in range(tau_size):
            tau_array[index_tau] = tau[index_tau]

        indices = []

        if self.pt.has_source_t:
            indices.extend([
                self.pt.index_tp_t0,
                self.pt.index_tp_t1,
                self.pt.index_tp_t2
                ])
            names.extend([
                "t0",
                "t1",
                "t2"
                ])
        if self.pt.has_source_p:
            indices.append(self.pt.index_tp_p)
            names.append("p")
        if self.pt.has_source_phi:
            indices.append(self.pt.index_tp_phi)
            names.append("phi")
        if self.pt.has_source_phi_plus_psi:
            indices.append(self.pt.index_tp_phi_plus_psi)
            names.append("phi_plus_psi")
        if self.pt.has_source_phi_prime:
            indices.append(self.pt.index_tp_phi_prime)
            names.append("phi_prime")
        if self.pt.has_source_psi:
            indices.append(self.pt.index_tp_psi)
            names.append("psi")
        if self.pt.has_source_H_T_Nb_prime:
            indices.append(self.pt.index_tp_H_T_Nb_prime)
            names.append("H_T_Nb_prime")
        if self.pt.index_tp_k2gamma_Nb:
            indices.append(self.pt.index_tp_k2gamma_Nb)
            names.append("k2gamma_Nb")
        if self.pt.has_source_h:
            indices.append(self.pt.index_tp_h)
            names.append("h")
        if self.pt.has_source_h_prime:
            indices.append(self.pt.index_tp_h_prime)
            names.append("h_prime")
        if self.pt.has_source_eta:
            indices.append(self.pt.index_tp_eta)
            names.append("eta")
        if self.pt.has_source_eta_prime:
            indices.append(self.pt.index_tp_eta_prime)
            names.append("eta_prime")
        if self.pt.has_source_delta_tot:
            indices.append(self.pt.index_tp_delta_tot)
            names.append("delta_tot")
        if self.pt.has_source_delta_m:
            indices.append(self.pt.index_tp_delta_m)
            names.append("delta_m")
        if self.pt.has_source_delta_cb:
            indices.append(self.pt.index_tp_delta_cb)
            names.append("delta_cb")
        if self.pt.has_source_delta_g:
            indices.append(self.pt.index_tp_delta_g)
            names.append("delta_g")
        if self.pt.has_source_delta_b:
            indices.append(self.pt.index_tp_delta_b)
            names.append("delta_b")
        if self.pt.has_source_delta_cdm:
            indices.append(self.pt.index_tp_delta_cdm)
            names.append("delta_cdm")
        if self.pt.has_source_delta_idm:
            indices.append(self.pt.index_tp_delta_idm)
            names.append("delta_idm")
        if self.pt.has_source_delta_dcdm:
            indices.append(self.pt.index_tp_delta_dcdm)
            names.append("delta_dcdm")
        if self.pt.has_source_delta_fld:
            indices.append(self.pt.index_tp_delta_fld)
            names.append("delta_fld")
        if self.pt.has_source_delta_scf:
            indices.append(self.pt.index_tp_delta_scf)
            names.append("delta_scf")
        if self.pt.has_source_delta_dr:
            indices.append(self.pt.index_tp_delta_dr)
            names.append("delta_dr")
        if self.pt.has_source_delta_ur:
            indices.append(self.pt.index_tp_delta_ur)
            names.append("delta_ur")
        if self.pt.has_source_delta_idr:
            indices.append(self.pt.index_tp_delta_idr)
            names.append("delta_idr")
        if self.pt.has_source_delta_ncdm:
            for incdm in range(self.ba.N_ncdm):
              indices.append(self.pt.index_tp_delta_ncdm1+incdm)
              names.append("delta_ncdm[{}]".format(incdm))
        if self.pt.has_source_theta_tot:
            indices.append(self.pt.index_tp_theta_tot)
            names.append("theta_tot")
        if self.pt.has_source_theta_m:
            indices.append(self.pt.index_tp_theta_m)
            names.append("theta_m")
        if self.pt.has_source_theta_cb:
            indices.append(self.pt.index_tp_theta_cb)
            names.append("theta_cb")
        if self.pt.has_source_theta_g:
            indices.append(self.pt.index_tp_theta_g)
            names.append("theta_g")
        if self.pt.has_source_theta_b:
            indices.append(self.pt.index_tp_theta_b)
            names.append("theta_b")
        if self.pt.has_source_theta_cdm:
            indices.append(self.pt.index_tp_theta_cdm)
            names.append("theta_cdm")
        if self.pt.has_source_theta_idm:
            indices.append(self.pt.index_tp_theta_idm)
            names.append("theta_idm")
        if self.pt.has_source_theta_dcdm:
            indices.append(self.pt.index_tp_theta_dcdm)
            names.append("theta_dcdm")
        if self.pt.has_source_theta_fld:
            indices.append(self.pt.index_tp_theta_fld)
            names.append("theta_fld")
        if self.pt.has_source_theta_scf:
            indices.append(self.pt.index_tp_theta_scf)
            names.append("theta_scf")
        if self.pt.has_source_theta_dr:
            indices.append(self.pt.index_tp_theta_dr)
            names.append("theta_dr")
        if self.pt.has_source_theta_ur:
            indices.append(self.pt.index_tp_theta_ur)
            names.append("theta_ur")
        if self.pt.has_source_theta_idr:
            indices.append(self.pt.index_tp_theta_idr)
            names.append("theta_idr")
        if self.pt.has_source_theta_ncdm:
            for incdm in range(self.ba.N_ncdm):
              indices.append(self.pt.index_tp_theta_ncdm1+incdm)
              names.append("theta_ncdm[{}]".format(incdm))

        for index_type, name in zip(indices, names):
            tmparray = np.empty((k_size,tau_size))
            for index_k in range(k_size):
                for index_tau in range(tau_size):
                    tmparray[index_k][index_tau] = sources_ptr[index_md][index_ic*tp_size+index_type][index_tau*k_size + index_k];

            sources[name] = np.asarray(tmparray)

        return (sources, np.asarray(k_array), np.asarray(tau_array))
