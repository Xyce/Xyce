# This file is used to demonstrate the Python initialize(),
# runSimulation() and close() methods.

import sys

from xyce_interface import xyce_interface

# this calls the xyce_interface.open() method to
# make a xyce object

# if the user specifies the Xyce library directory, pass that
# to the xyce_interface constructor.  Otherwise use the default.
xyceObj = None
if( len(sys.argv) > 1):
  libDirectory = sys.argv[1]
  xyceObj = xyce_interface(libdir=libDirectory)
else:
  xyceObj = xyce_interface()

if( xyceObj == None):
  print("Error:  Could not create a Xyce object.  Exiting")
  exit(-1)

argv= ['runACircuit.cir']
print( "calling initialize with netlist %s" % argv[0] )

result = xyceObj.initialize(argv)
print( "return value from initialize is %d" % result )

print( "Calling runSimulation..." )
result = xyceObj.runSimulation()
print( "return value from runSimulation is %d" % result )

print( "calling close")
xyceObj.close()

