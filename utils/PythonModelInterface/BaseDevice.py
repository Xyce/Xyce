try:
    from XyceObjects import DeviceOptions, SolverState
except:
    pass

# specifies functions that must be defined
class BaseDevice(object):

    # intended if no p_params to store in b,d,i,s params
    # otherwise, override and be sure last call inside
    # processPythonParams is to call pythonParamsMerge
    def processPythonParams(self, b_params, d_params, 
            i_params, s_params):
        pass

    # default merging into b,d,i,s params
    def pythonParamsMerge(self, b_params, d_params, 
            i_params, s_params, p_params):
        # merges p_params into b,d,i,s params
        # gives priority to existing values in b,d,i, and s 
        # since these were specified in the netlist by the user
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

    # do not override, only called from Xyce-PyMi
    def processTotalVars(self, i_params):
        try:
            num_external_vars = i_params["numExternalVars"]
        except:
            assert False, "processTotalVars() called before \
                    setNumExternalVars()"
        try:
            num_internal_vars = i_params["numInternalVars"]
        except:
            assert False, "processTotalVars() called before \
                    processNumInternalVars()"
        i_params["numVars"] = num_external_vars + num_internal_vars

    # do not override, only called from Xyce-PyMi
    # and determined by topology
    def setNumExternalVars(self, i_params, num_external_vars):
        i_params["numExternalVars"]=num_external_vars

    def processNumInternalVars(self, b_params, d_params, 
            i_params, s_params):
        i_params["numInternalVars"]=0
        return i_params["numInternalVars"]

    def processNumStateVars(self, b_params, d_params, 
            i_params, s_params):
        i_params["numStateVars"]=0
        return i_params["numStateVars"]

    def processNumStoreVars(self, b_params, d_params, 
            i_params, s_params):
        i_params["numStoreVars"]=0
        return i_params["numStoreVars"]

    def processNumBranchDataVars(self, b_params, d_params, 
            i_params, s_params):
        i_params["numBranchDataVars"]=0
        return i_params["numBranchDataVars"]

    def processNumBranchDataVarsIfAllocated(self, b_params, d_params, 
            i_params, s_params):
        try:
            i_params["numBranchDataVarsIfAllocated"] = \
                i_params["numExternalVars"]
        except:
            assert False, "processNumBranchDataVarsIfAllocated(...) \
                    called before setNumExternalVars(...)"
        return i_params["numBranchDataVarsIfAllocated"]

    def getArraySizes(self, b_params, d_params, 
            i_params, s_params):
        try:
            num_vars = i_params["numVars"]
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
    def getJacStampSize(self, b_params, d_params, 
            i_params, s_params):
        return np.zeros(shape=(0,),dtype='i4')
    
    # only needs overridden if getJacStampSize is overridden
    def setJacStamp(self, jacStamp, b_params, d_params, 
            i_params, s_params):
        pass

    # called prior to computeXyceVectors, computed results should be
    # store in self.* so that they persist for future 
    # computeXyceVectors calls
    def initialize(self, deviceOptions, solverState, 
            b_params, d_params, i_params, s_params):
        pass

    # called at the end of the simulation
    def finalize(self, b_params, d_params, i_params, s_params):
        pass

    # this function must be overridden to provide device definition
    def computeXyceVectors(self, fSV, solV, stoV, staV, 
            deviceOptions, solverState,
            origFlag, F, Q, B, dFdX, dQdX, dFdXdVp, dQdXdVp, 
            b_params, d_params, i_params, s_params):
        raise NotImplementedError

try:
    from unittest import TestCase
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


except:
    pass
