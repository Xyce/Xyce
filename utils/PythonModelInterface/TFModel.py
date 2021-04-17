import os, sys
if not hasattr(sys, 'argv'):
    sys.argv  = ['']
if len(sys.argv)==0:
    sys.argv.append('') 


import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np
import tensorflow.keras.backend as K
import h5py
import importlib

import numpy as np

class TFModel(object):

    def __init__(self, fname):

        assert os.path.isabs(fname), "Relative path given to TFModel construction."

        EAGER = False
        # disable eager execution for TF2 and initialize with a seed
        # to get reproducible results
        if not EAGER:
            tf.compat.v1.disable_eager_execution()
            tf.keras.backend.clear_session()
            sess = tf.compat.v1.Session()
            tf.random.set_seed(100234)
            tf.compat.v1.keras.backend.set_session(sess)
            tf.compat.v1.global_variables_initializer()
        else:
            tf.keras.backend.clear_session()
            tf.random.set_seed(100234)

        with h5py.File(fname, 'r') as f:
            for group in f.keys():
                if (group=="parameters"):
                    this_parameter_data = f[group][()]
                    parameter_data = this_parameter_data.decode("unicode_escape")
                if (group=="custom_layers"):
                    this_custom_layers_data = f[group][()]
                    custom_layers_data = this_custom_layers_data.decode("unicode_escape")
                if (group=="custom_losses"):
                    this_custom_losses_data = f[group][()]
                    custom_losses_data = this_custom_losses_data.decode("unicode_escape")
                if (group=="data_transform"):
                    this_data_transform_data = f[group][()]
                    data_transform_data = this_data_transform_data.decode("unicode_escape")

        # look for DATA_* groups to populate DATA with
        DATA={}
        DATA_names_list=[]
        # this data file name could be changed if DATA_* vars were stored elsewhere
        with h5py.File(fname, 'r') as f:
            for group in f.keys():
                if group.startswith("DATA_"):
                    DATA[group[len("DATA_"):]]=f[group][()]
        # get each * from DATA_* and make it a variable (i.e. DATA_ipn -> DATA.ipn)
        globals().update({'DATA':DATA})
        
        global calling_from_h5
        calling_from_h5 = True

        custom_layers_spec = importlib.util.spec_from_loader('CustomLayers', loader=None)
        custom_losses_spec = importlib.util.spec_from_loader('CustomLosses', loader=None)
        data_transform_spec = importlib.util.spec_from_loader('DataTransform', loader=None)
        custom_layers_module = importlib.util.module_from_spec(custom_layers_spec)
        custom_losses_module = importlib.util.module_from_spec(custom_losses_spec)
        data_transform_module = importlib.util.module_from_spec(data_transform_spec)
        if 'custom_layers_data' in locals() and 'custom_losses_data' in locals() \
                and 'data_transform_data' in locals() and 'parameter_data' in locals():
            exec(custom_layers_data, custom_layers_module.__dict__)
            exec(custom_losses_data, custom_losses_module.__dict__)
            exec(data_transform_data, data_transform_module.__dict__)

        prms_spec = importlib.util.spec_from_loader('prms', loader=None)
        prms = importlib.util.module_from_spec(prms_spec)
        prms.__dict__.update(custom_layers_module.__dict__)
        prms.__dict__.update(custom_losses_module.__dict__)
        prms.__dict__.update(data_transform_module.__dict__)
        if 'custom_layers_data' in locals() and 'custom_losses_data' in locals() \
                and 'data_transform_data' in locals() and 'parameter_data' in locals():
            exec(parameter_data, globals(), prms.__dict__)
            if (prms.PRECISION.lower()=="double"):
                K.set_floatx('float64')
        else:
            K.set_floatx('float64')


        customLayersClasses = dict([(name, cls) for name, cls in custom_layers_module.__dict__.items() if isinstance(cls, type)])
        customLossesClasses = dict([(name, cls) for name, cls in custom_losses_module.__dict__.items() if isinstance(cls, type)])
        customObjects = {}
        customObjects.update(customLayersClasses)
        customObjects.update(customLossesClasses)

        self.tf_model = keras.models.load_model(fname, custom_objects=customObjects, compile=True)
        self.sess = tf.compat.v1.keras.backend.get_session()
        self.loss = self.tf_model.output
        self.inputs = self.tf_model.input
        self.grads = tf.keras.backend.gradients(self.loss, self.inputs)
        if 'custom_layers_data' in locals() and 'custom_losses_data' in locals() \
                and 'data_transform_data' in locals() and 'parameter_data' in locals():
            self.transform = prms.DATA_TRANSFORM
        else:
            self.transform = None

        # for storing last evaluated value
        self.last_input_value = sys.float_info.max
        self.last_output_value = None
        self.last_output_gradient = None

    def eval_if_needed(self, input_value):
        if (input_value!=self.last_input_value):
            input_array = np.ndarray([1,1],dtype='f8')
            if self.transform is not None:
                input_array[0][0] = self.transform.transformInput(input_value)
            else:
                input_array[0][0] = input_value
            self.last_input_value = input_value
            self.last_output_value = self.tf_model.predict(input_array)
            self.last_output_gradient = self.sess.run(self.grads, {self.inputs:input_array})
            #(self.last_output_value, self.last_output_gradient) = f_and_df(net=self.tf_model, v=input_array)
        return (self.last_output_value, self.last_output_gradient)

    def predict(self, input_value):
        (output_array, blah) = self.eval_if_needed(input_value)
        if self.transform is not None:
            return self.transform.invertOutput(output_array[0][0], input_value)
        else:
            return output_array[0][0]

    def gradient(self, input_value, output_comp=0):
        #return 0
        (prediction, out) = self.eval_if_needed(input_value)
        #out_val = self.transform.invertOutput(prediction[0][0], input_value)
        if self.transform is not None:
            deriv = self.transform.derivativeOfRawOutputWithRespectToRawInput(self.last_input_value, prediction[0][0], out[0][0][0], self.last_input_value)
        else:
            deriv = out[0][0][0]
        return deriv

