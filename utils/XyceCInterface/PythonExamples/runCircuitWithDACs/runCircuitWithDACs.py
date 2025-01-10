# This file is used to demonstrate the Python getDeviceNames(),
# getDACDeviceNames(), updateTimeVoltagePairs() and 
# obtainResponse() methods

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

argv= ['runCircuitWithDACs.cir']
print( "calling initialize with netlist %s" % argv[0] )

result = xyceObj.initialize(argv)
print( "return value from initialize is %d" % result )

# get ADC and DAC names;
(result, names) = xyceObj.getDeviceNames("YADC")
print( "return value from getDeviceNames is %d" % result )
print( names )

(result, DACnames) = xyceObj.getDACDeviceNames()
print( "return value from getDACDeviceNames is %d" % result )
print( DACnames )

#
# A bug in the DAC device (put there for Habinero support) only takes
# the last time in the time voltage pairs list if the current sim time is 0.
# So simulate a bit first.
requested_time = 1.0e-10
(result, actual_time) = xyceObj.simulateUntil( requested_time )
print( "simulationUntil status = %d and actual_time = %f" % (result, actual_time) )

# try setting up the DAC to pulse twice
timeArrayBase = [ 0.0, 0.1e-4, 0.2e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.8e-4, 1.0e-4, 1.1e-4 ]
timeArray = timeArrayBase[:]
voltageArray= [ 0.0, 0.0,    3.0,    3.0,    0.0,    0.0,    3.0,    3.0,    0.0    ] 

steps = range(0,10)
total_sim_time = 20.0e-4
for i in steps:
  result = xyceObj.updateTimeVoltagePairs( DACnames[0], timeArray, voltageArray )
  requested_time = 0.0 + (i+1) * 0.1 * total_sim_time
  print( "Calling simulateUntil for requested_time %f" % requested_time )
  actual_time = 0.0
  (result, actual_time) = xyceObj.simulateUntil( requested_time )
  print( "simulateUntil status = %d and actual_time = %f" % (result, actual_time) )
  
  # get some result from the ciruit
  (result,value) = xyceObj.obtainResponse('YMEMRISTORRES')
  print( "return value from obtainResponse = %d" % result)
  print( "R= %f " % value )

  # update timeArray to repeat pulse 
  for j in range(0,len(timeArray)):
    timeArray[j] = timeArrayBase[j] + requested_time
  
print( "calling close")
xyceObj.close()

