# This file is used to demonstrate the Python initialize(),
# runSimulation() and close() methods.

from xyce_interface import xyce_interface

# this calls the xyce_interface.open() method to
# make a xyce object

xyceObj = xyce_interface()
print( xyceObj )

argv= ['runACircuit.cir']
print( "calling initialize with netlist %s" % argv[0] )

result = xyceObj.initialize(argv)
print( "return value from initialize is %d" % result )

print( "Calling runSimulation..." )
result = xyceObj.runSimulation()
print( "return value from runSimulation is %d" % result )

print( "calling close")
xyceObj.close()

