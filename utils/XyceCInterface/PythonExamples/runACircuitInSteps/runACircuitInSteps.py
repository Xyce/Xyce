# This file is used to demonstrate the Python simulateUntil() method

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

argv= ['runACircuitInSteps.cir']
print( "calling initialize with netlist %s" % argv[0] )

result = xyceObj.initialize(argv)
print( "return value from initialize is %d" % result )

steps = range(0,10)
for i in steps:
  requested_time = 0.0 + (i+1) * 0.1
  print( "Calling simulateUntil for requested_time = %f" % requested_time )
  actual_time = 0.0
  (result, actual_time) = xyceObj.simulateUntil( requested_time )
  print( "simulateUntil status = %d and actual_time = %f" % (result, actual_time) )

print( "calling close")
xyceObj.close()

