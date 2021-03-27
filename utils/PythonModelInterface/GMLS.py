import pycompadre 
import numpy as np

class GMLS(object):

    def __init__(self, source_sites, polynomial_order, weighting_power = 2, epsilon_multiplier = 1.5):
        self.last_target_site = np.zeros((1,))

        self.kokkos_obj=pycompadre.KokkosParser()

        assert len(source_sites.shape)==2, "2D array must be given to GMLS for source_sites (#sites x spatial dimension)"
        self.input_dimensions = source_sites.shape[1]
        self.polynomial_order = polynomial_order

        # initialize 3rd order reconstruction using 2nd order basis in 3D (GMLS)
        self.gmls_obj=pycompadre.GMLS(polynomial_order, self.input_dimensions, "QR", "STANDARD")
        self.gmls_obj.setWeightingPower(weighting_power)
        self.weighting_power = weighting_power
        self.gmls_obj.setWeightingType("power")

        # neighbor search
        self.epsilon_multiplier = epsilon_multiplier
        self.gmls_helper = pycompadre.ParticleHelper(self.gmls_obj)
        self.gmls_helper.generateKDTree(source_sites)

        self.gmls_obj.addTargets(pycompadre.TargetOperation.ScalarPointEvaluation)
        self.gmls_obj.addTargets(pycompadre.TargetOperation.PartialXOfScalarPointEvaluation)
        self.gmls_obj.addTargets(pycompadre.TargetOperation.PartialYOfScalarPointEvaluation)

    def __del__(self):
        del self.gmls_obj
        del self.gmls_helper
        del self.kokkos_obj

    def predict(self, target_site, data_vector):
        assert target_site.shape[1]==self.input_dimensions, "Wrong dimensions for target site (%d vs %d)" % (target_site.shape[1], self.input_dimensions)

        # set new target_site
        #self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)

        # generate stencil with number of batches set to 1, and keeping coefficients (not necessary)
        if (not(np.array_equal(target_site, self.last_target_site))):
            self.last_target_site = target_site
            self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)
            self.gmls_obj.generateAlphas(1, False)

        #if (self.last_target_site is not None):
        #    if (not(np.array_equal(target_site, self.last_target_site))):
        #        self.last_target_site = target_site
        #        # set new target_site
        #        self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)
        #        self.gmls_obj.generateAlphas(1, False)
        #else:
        #    self.last_target_site = target_site
        #    # set new target_site
        #    self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)
        #    self.gmls_obj.generateAlphas(1, False)

        # apply stencil to sample data for all targets
        #out = np.array([0],dtype=np.float64)
        return self.gmls_helper.applyStencilSingleTarget(data_vector, pycompadre.TargetOperation.ScalarPointEvaluation)
        #return out[0]
        #output = self.gmls_helper.applyStencil(data_vector, pycompadre.TargetOperation.ScalarPointEvaluation)
        #return output
        #return np.array([1,])#self.gmls_helper.applyStencil(data_vector, pycompadre.TargetOperation.ScalarPointEvaluation)

    def gradient(self, target_site, data_vector, derivative_direction):
        assert target_site.shape[1]==self.input_dimensions, "Wrong dimensions for target site (%d vs %d)" % (target_site.shape[1], self.input_dimensions)


        #self.last_target_site = target_site
        # generate stencil with number of batches set to 1, and keeping coefficients (not necessary)
        #if (self.last_target_site is None):

        if (not(np.array_equal(target_site, self.last_target_site))):
            self.last_target_site = target_site
            self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)
            self.gmls_obj.generateAlphas(1, False)

        #if (self.last_target_site is not None):
        #    if (not(np.array_equal(target_site, self.last_target_site))):
        #        #self.gmls_obj.generateAlphas(1, False)
        #        self.last_target_site = target_site
        #        # set new target_site
        #        self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)
        #else:
        #    #self.gmls_obj.generateAlphas(1, False)
        #    self.last_target_site = target_site
        #    # set new target_site
        #    self.gmls_helper.generateNeighborListsFromKNNSearchAndSet(target_site, self.polynomial_order, self.input_dimensions, self.epsilon_multiplier)

        # apply stencil to sample data for all targets
        #out = np.reshape(np.array([0,0],dtype=np.float64), newshape=(1,2))
        #out = np.array([0],dtype=np.float64)
        if (derivative_direction==0):
            return self.gmls_helper.applyStencilSingleTarget(data_vector, pycompadre.TargetOperation.PartialXOfScalarPointEvaluation)
        elif (derivative_direction==1):
            return self.gmls_helper.applyStencilSingleTarget(data_vector, pycompadre.TargetOperation.PartialYOfScalarPointEvaluation)
        else:
            return 0
        #return out[0]

        #return self.gmls_helper.applyStencil(data_vector, pycompadre.TargetOperation.GradientOfScalarPointEvaluation)[:,derivative_direction]
        #return np.array([1,])#self.gmls_helper.applyStencil(data_vector, pycompadre.TargetOperation.GradientOfScalarPointEvaluation)[:,derivative_direction]
