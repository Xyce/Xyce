# This file is used to demonstrate the Python getDeviceNames(),
# getDACDeviceNames(), updateTimeVoltagePairs() and 
# obtainResponse() methods

import sys
import math

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

argv= ['MultiPort2.cir']
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
timeArrayBase = [ 0.0, 0.1e-5, 0.2e-5, 0.4e-5, 0.5e-5, 0.7e-5, 0.8e-5, 1.0e-5, 1.1e-5 ]
timeArray = timeArrayBase[:]
voltageArray= [ 0.0, 0.0,    3.0,    3.0,    0.0,    0.0,    3.0,    3.0,    0.0    ] 

numSteps = 250
steps = range(0,numSteps)
total_sim_time = 5.0e-3
simResult = 0
for i in steps:
  #result = xyceObj.updateTimeVoltagePairs( DACnames[0], timeArray, voltageArray )
  requested_time = 0.0 + (i+1) * (1.0/numSteps) * total_sim_time
  print( "Calling simulateUntil for requested_time %f" % requested_time )
  actual_time = 0.0
  (result, actual_time) = xyceObj.simulateUntil( requested_time )
  print( "simulateUntil status = %d and actual_time = %f" % (result, actual_time) )
  (status, ADCnames, numADCnames, numPoints, timeArray, voltageArray) = xyceObj.getTimeVoltagePairsADC()
  #(status, ADCnames, numADCnames, numPoints, timeArray, voltageArray) = xyceObj.getTimeVoltagePairsADCLimitData()
  for i in range(len(ADCnames)):
    for j in range(numPoints):
      expectedValue = 0.0
      if( ADCnames[i] == 'YADC!ADCHIGH' ):
        expectedValue = 5.0*math.sin( 2*math.pi*5e3*timeArray[i][j]) 
      if( ADCnames[i] == 'YADC!ADCLOW' ):
        expectedValue = 3.0*math.cos( 2*math.pi*2e3*timeArray[i][j])
      if( math.fabs( voltageArray[i][j] -  expectedValue) < 1.0e-7):
        print("Compare test PASS" )
      else:
        print( "Compare test FAIL  -> %s time=%g value=%g != expected value = %g" % (ADCnames[i], timeArray[i][j], voltageArray[i][j], expectedValue ) )
        simResult = -1
  # update timeArray to repeat pulse 
  #for j in range(0,len(timeArray)):
  #  timeArray[j] = timeArrayBase[j] + requested_time
  
print( "calling close")
xyceObj.close()
exit( simResult )
