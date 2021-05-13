from XyceObjects import DeviceOptions, SolverState
# specifies functions that must be defined
class BaseDevice(object):

    # intended if no p_params to store in b,d,i,s params
    # otherwise, override and be sure last call inside
    # processPythonParams is to call pythonParamsMerge
    def processPythonParams(self, b_params, d_params, i_params, s_params):
        pass

    # default merging into b,d,i,s params
    def pythonParamsMerge(self, b_params, d_params, i_params, s_params, p_params):
        # merges p_params into b,d,i,s params
        # gives priority to existing values in b,d,i, and s since these were
        # specified in the netlist by the user, default back to p_params
        for item in p_params.items():
            if (item[0] not in d_params.keys()):
                if isinstance(item[1], bool):
                    b_params[item[0]] = 1 if item[1] else 0
                elif isinstance(item[1], float):
                    d_params[item[0]] = item[1]
                elif isinstance(item[1], int): 
                    i_params[item[0]] = item[1]
                elif isinstance(item[1], str): 
                    s_params[item[0]] = item[1]

    # do not override, only called from Xyce-PyMi
    def processTotalVars(self, i_params):
        try:
            num_external_vars = i_params["numExternalVars"]
        except:
            assert False, "processTotalVars() called before setNumExternalVars()"
        try:
            num_internal_vars = i_params["numInternalVars"]
        except:
            assert False, "processTotalVars() called before processNumInternalVars()"
        i_params["numVars"] = num_external_vars + num_internal_vars

    # do not override, only called from Xyce-PyMi
    # and determined by topology
    def setNumExternalVars(self, i_params, num_external_vars):
        i_params["numExternalVars"]=num_external_vars

    def processNumInternalVars(self, b_params, d_params, i_params, s_params):
        i_params["numInternalVars"]=0
        return i_params["numInternalVars"]

    def processNumStateVars(self, b_params, d_params, i_params, s_params):
        i_params["numStateVars"]=0
        return i_params["numStateVars"]

    def processNumStoreVars(self, b_params, d_params, i_params, s_params):
        i_params["numStoreVars"]=0
        return i_params["numStoreVars"]

    def processNumBranchDataVars(self, b_params, d_params, i_params, s_params):
        i_params["numBranchDataVars"]=0
        return i_params["numBranchDataVars"]

    def processNumBranchDataVarsIfAllocated(self, b_params, d_params, i_params, s_params):
        try:
            i_params["numBranchDataVarsIfAllocated"]=i_params["numExternalVars"]
        except:
            assert False, "processNumBranchDataVarsIfAllocated(...) called before setNumExternalVars(...)"
        return i_params["numBranchDataVarsIfAllocated"]

    def getArraySizes(self, b_params, d_params, i_params, s_params):
        try:
            num_vars = i_params["numVars"]
        except:
            assert False, "getArraySizes(...) called before processTotalVars(...)"
        size_dict = {}
        size_dict['F']=[num_vars,]
        size_dict['Q']=[num_vars,]
        size_dict['B']=[num_vars,]
        size_dict['dFdX']=[num_vars,num_vars]
        size_dict['dQdX']=[num_vars,num_vars]
        return size_dict
    
    # tell Xyce-PyMi Jacobian stamp (nonzeros per row) size
    # if not overridden, a dense Jacobian will be created
    def getJacStampSize(self, b_params, d_params, i_params, s_params):
        return np.zeros(shape=(0,),dtype='i4')
    
    # only needs overridden if getJacStampSize is overridden
    def setJacStamp(self, jacStamp, b_params, d_params, i_params, s_params):
        pass

    # called prior to computeXyceVectors, computed results should be
    # store in self.* so that they persist for future 
    # computeXyceVectors calls
    def initialize(self, deviceOptions, solverState, b_params, d_params, i_params, s_params):
        pass

    # this function must be overridden to provide device definition
    def computeXyceVectors(self, solV, fSV, stoV, t, deviceOptions, solverState,
            origFlag, F, Q, B, dFdX, dQdX, dFdXdVp, dQdXdVp, 
            b_params, d_params, i_params, s_params):
        raise NotImplementedError
