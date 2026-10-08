try:
    from XyceObjects import DeviceOptions, SolverState
except:
    pass
from unittest import TestCase
import inspect
import numpy as np

# specifies functions that must be defined
class BaseDevice(object):

    # ------------------------------------------------------------------
    # Backwards compatibility for pre-nanobind device files.
    #
    # The nanobind refactor changed the calling convention: Xyce-PyMi now
    # calls each device method with NO positional arguments and instead
    # provisions data as attributes on the device object (self.b_params,
    # self.i_params, self.deviceOptions, self.solverState, self.jacStamp,
    # self.F/Q/B/dFdX/... etc.) before the call.  Older device files still
    # declare the legacy positional signatures, e.g.
    #     def getArraySizes(self, b_params, d_params, i_params, s_params):
    # A zero-arg call into those raises TypeError.
    #
    # _pymi_bind(name) resolves this ONCE per device method and returns a
    # zero-argument callable that Xyce-PyMi can invoke (and cache, for the
    # hot computeXyceVectors path):
    #   * new-style method (no positional params) -> the bound method itself
    #   * legacy-style method -> a closure that feeds the arguments from
    #     self.* using _LEGACY_ARGS, in the canonical Xyce order.
    # The attribute names below match exactly what Xyce-PyMi sets on the
    # device prior to the call, so both conventions read the same data.
    # ------------------------------------------------------------------
    _LEGACY_ARGS = {
        "processPythonParams":                ["b_params", "d_params", "i_params", "s_params"],
        "processNumInternalVars":             ["b_params", "d_params", "i_params", "s_params"],
        "processNumStateVars":                ["b_params", "d_params", "i_params", "s_params"],
        "processNumStoreVars":                ["b_params", "d_params", "i_params", "s_params"],
        "processNumBranchDataVars":           ["b_params", "d_params", "i_params", "s_params"],
        "processNumBranchDataVarsIfAllocated":["b_params", "d_params", "i_params", "s_params"],
        "getArraySizes":                      ["b_params", "d_params", "i_params", "s_params"],
        "getJacStampSize":                    ["b_params", "d_params", "i_params", "s_params"],
        "setJacStamp":                        ["jacStamp", "b_params", "d_params", "i_params", "s_params"],
        "setErrorWeightMask":                 ["b_params", "d_params", "i_params", "s_params"],
        "initialize":                         ["deviceOptions", "solverState",
                                               "b_params", "d_params", "i_params", "s_params"],
        "finalize":                           ["b_params", "d_params", "i_params", "s_params"],
        "computeXyceVectors":                 ["fSV", "solV", "stoV", "staV",
                                               "deviceOptions", "solverState",
                                               "origFlag", "F", "Q", "B", "dFdX", "dQdX",
                                               "dFdXdVp", "dQdXdVp",
                                               "b_params", "d_params", "i_params", "s_params"],
    }

    def _pymi_bind(self, name):
        # Return a zero-arg callable for method `name`, adapting a legacy
        # positional signature to the new attribute-based convention.
        # Resolution happens once; the returned callable does no further
        # signature inspection, so the cached computeXyceVectors path stays
        # fast.
        fn = getattr(self, name)
        params = [p for p in inspect.signature(fn).parameters.values()
                  if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
        if not params:
            # new-style, zero-arg method: call directly
            return fn
        # legacy-style method: feed positional args from self.* attributes
        try:
            legacy = self._LEGACY_ARGS[name]
        except KeyError:
            raise TypeError(
                "Xyce-PyMi device method '%s' declares positional parameters "
                "but is not a known device callback; cannot adapt its legacy "
                "signature." % name)
        if len(params) > len(legacy):
            raise TypeError(
                "Xyce-PyMi device method '%s' declares more required "
                "arguments (%d) than the legacy convention provides (%d); it "
                "is not backwards compatible with this Xyce version."
                % (name, len(params), len(legacy)))
        return lambda: fn(*[getattr(self, a) for a in legacy])


    # intended if no p_params to store in b,d,i,s params
    # otherwise, override and be sure last call inside
    # processPythonParams is to call pythonParamsMerge
    def processPythonParams(self, *_legacy_args):
        # has access to:
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        pass

    # default merging into b,d,i,s params
    def pythonParamsMerge(self, *args):
        # merges p_params into b,d,i,s params
        # gives priority to existing values in b,d,i, and s
        # since these were specified in the netlist by the user
        # has access to:
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        if len(args) == 1:
            (p_params,) = args
            b_params, d_params = self.b_params, self.d_params
            i_params, s_params = self.i_params, self.s_params
        elif len(args) == 5:
            b_params, d_params, i_params, s_params, p_params = args
        else:
            raise TypeError(
                "pythonParamsMerge() expects either (p_params) [new style] "
                "or (b_params, d_params, i_params, s_params, p_params) "
                "[legacy], got %d positional arguments." % len(args))
        for item in p_params.items():
            if isinstance(item[1], bool):
                if (item[0] not in b_params.keys()):
                    b_params[item[0]] = 1 if item[1] else 0
            elif isinstance(item[1], float):
                if (item[0] not in d_params.keys()):
                    d_params[item[0]] = item[1]
            elif isinstance(item[1], int):
                if (item[0] not in i_params.keys()):
                    i_params[item[0]] = item[1]
            elif isinstance(item[1], str):
                if (item[0] not in s_params.keys()):
                    s_params[item[0]] = item[1]

    def processNumInternalVars(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        self.i_params["numInternalVars"]=0
        return self.i_params["numInternalVars"]

    def processNumStateVars(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        self.i_params["numStateVars"]=0
        return self.i_params["numStateVars"]

    def processNumStoreVars(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        self.i_params["numStoreVars"]=0
        return self.i_params["numStoreVars"]

    def processNumBranchDataVars(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        self.i_params["numBranchDataVars"]=0
        return self.i_params["numBranchDataVars"]

    def processNumBranchDataVarsIfAllocated(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        try:
            self.i_params["numBranchDataVarsIfAllocated"] = \
                self.i_params["numExternalVars"]
        except:
            assert False, "numExternalVars not yet populated in i_params \
                    (set by the Xyce-PyMi harness) when \
                    processNumBranchDataVarsIfAllocated(...) ran"
        return self.i_params["numBranchDataVarsIfAllocated"]

    def getArraySizes(self, *_legacy_args):
        # has access to:
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        try:
            num_vars = self.i_params["numVars"]
        except:
            assert False, "getArraySizes(...) called before \
                    processTotalVars(...)"
        size_dict = {}
        # 'F' must have size num_vars
        size_dict['F']=[num_vars,]
        # 'Q' can have size num_vars or 0
        size_dict['Q']=[num_vars,]
        # 'B' can have size num_vars or 0
        size_dict['B']=[num_vars,]

        # limiting variables will have same size
        # as variable they shadow

        # gradients of variables F and Q (dFdX and dQdX) will have
        # shape (size F x size F) and (size Q x size Q), respectively
        return size_dict
    
    # tell Xyce-PyMi Jacobian stamp (nonzeros per row) size
    # if not overridden, a dense Jacobian will be created
    def getJacStampSize(self, *_legacy_args):
        # has access to:
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        return np.zeros(shape=(0,),dtype='i4')

    # only needs overridden if getJacStampSize is overridden
    def setJacStamp(self, *_legacy_args):
        # has access to:
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        #   self.jacStamp
        pass

    # can mask intern solution variables by setting to 0.0
    def setErrorWeightMask(self, *_legacy_args):
        # has access to:
        #   self.b_params / self.d_params / self.i_params / self.s_params
        return np.empty((0,), dtype=np.float64)

    # called prior to computeXyceVectors, computed results should be
    # store in self.* so that they persist for future 
    # computeXyceVectors calls
    def initialize(self, *_legacy_args):
        # has access to:
        #   self.deviceOptions
        #   self.solverState
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        pass

    # called at the end of the simulation
    def finalize(self, *_legacy_args):
        # has access to:
        #   self.deviceOptions
        #   self.solverState
        #   self.b_params
        #   self.d_params
        #   self.i_params
        #   self.s_params
        #
        # NOTE: subclasses that override finalize() MUST call
        # super().finalize() so base-class teardown (e.g. Kokkos finalize in
        # KokkosDevice) still runs. This replaces the old __del__ chain.
        pass

    # this function must be overridden to provide device definition
    def computeXyceVectors(self, *_legacy_args):
        # has access to:
        #   self.deviceOptions
        #   self.solverState
        #   self.b_params - boolean parameters
        #   self.d_params - double parameters
        #   self.i_params - integer parameters
        #   self.s_params - string parameters
        #   self.fSV      - flagSolV
        #   self.F        - current (A) vector
        #   self.Q        - charge vector
        #   self.B        - source vector
        #   self.dFdX     - Jacobian
        #   self.dQdX     - Jacobian
        #   self.dFdXdVp  - variable for limiting
        #   self.dQdXdVp  - variable for limiting
        #   self.origFlag - variable for limiting
        #   self.solV     - solution vector of voltages
        #   self.stoV     - storage variables (updated each time step)
        #   self.staV     - state variables
        raise NotImplementedError

import numpy as np
class TestBaseDevice(TestCase):

    def setUp(self):
        self.fSV  = np.array([1,2,3],dtype='f8')
        self.solV = np.array([[1,2,3],[1,2,3],[1,2,3]],dtype='f8')
        self.stoV = np.array([[1,2,3],[1,2,3]],dtype='f8')
        self.staV = np.array([[1,2,3],[1,2,3],[1,2,3]],dtype='f8')
        self.deviceOptions = None
        self.solverState = None
        self.origFlag = True
        self.F = np.array([1,2,3],dtype='f8')
        self.Q = np.array([0,0,0],dtype='f8')
        self.B = np.array([1,1,1],dtype='f8')
        self.dFdX = np.array([[2,3,4],[5,6,7]],dtype='f8')
        self.dQdX = np.array([[8,9,10],[11,12,13]],dtype='f8')
        self.dFdXdVp = np.array([20,21,22],dtype='f8')
        self.dQdXdVp = np.array([23,24,25],dtype='f8')
        self.b_params = {'bool_set_first' : True }
        self.d_params = {'double_set_first' : 12.01 }
        self.i_params = {}
        self.s_params = {}

    def tearDown(self):
        pass

    def test_computeXyceVectors_missing(self):
        BD = BaseDevice()
        self.assertRaises(NotImplementedError, BD.computeXyceVectors, 
                                self.fSV, self.solV, self.stoV, self.staV, 
                                self.deviceOptions, self.solverState,
                                self.origFlag, self.F, self.Q, self.B, self.dFdX, self.dQdX, self.dFdXdVp, self.dQdXdVp, 
                                self.b_params, self.d_params, self.i_params, self.s_params)

    def test_pythonParamsMerge(self):
        BD = BaseDevice()
        p_params = {'string_name':'Test'}
        BD.pythonParamsMerge(self.b_params, self.d_params, self.i_params, self.s_params, p_params)
        assert self.b_params['bool_set_first']==True
        assert self.d_params['double_set_first']==12.01
        assert self.s_params['string_name']=='Test'
