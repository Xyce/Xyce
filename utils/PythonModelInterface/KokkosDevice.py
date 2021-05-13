import pycompadre 
import numpy as np
from BaseDevice import BaseDevice

class KokkosDevice(BaseDevice):

    instance_count = 0

    @staticmethod
    def safe_instantiate_kokkos():
        if KokkosDevice.instance_count == 0:
            KokkosDevice.kokkos_obj=pycompadre.KokkosParser()
        KokkosDevice.instance_count += 1

    @staticmethod
    def safe_finalize_kokkos():
        if KokkosDevice.instance_count > 1:
            KokkosDevice.instance_count -= 1
        elif KokkosDevice.instance_count == 1:
            del KokkosDevice.kokkos_obj
            KokkosDevice.instance_count -= 1

    def __init__(self):
        KokkosDevice.safe_instantiate_kokkos()

    def __del__(self):
        KokkosDevice.safe_finalize_kokkos()

