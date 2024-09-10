# This file is used to demonstrate the Python initialize(),
# runSimulation() and close() methods.

import sys

from xyce_interface import xyce_interface

# this calls the xyce_interface.open() method to
# make a xyce object

xyceObj = xyce_interface()

argv= ['runACircuit.cir']

result = xyceObj.initialize(argv)

result = xyceObj.runSimulation()

xyceObj.close()

