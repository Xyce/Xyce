# This file is used to demonstrate the Python 
# setCircuitParameter() getCircuitValue() methods
import sys
from xyce_interface import xyce_interface

xyceObj = xyce_interface()
result = xyceObj.initialize(['TestNetlist2.cir'])

# get circuit values before running 
r1res = 0.0
cirTemp = 0.0
r1res = xyceObj.getCircuitValue("R1:R")
print("==> R1 R=%g" % (r1res))
cirTemp = xyceObj.getCircuitValue("TEMP")
print("==> Circuit Temp = %f" %(cirTemp))

numSteps=100
deltaTime = 1.0 / numSteps
steps = range(0,10)
for i in range(0,numSteps):
  if( i==(numSteps/2)) :
    # in the middle of the run change the resistance on resistor R1
    setParamResult = xyceObj.setCircuitParameter( "R1:R", 1000.0)
    setParamResult = xyceObj.setCircuitParameter( "TEMP", -50.0)

  if( i==((numSteps/2)+1)):
    # in the middle of the run change the resistance on resistor R1
    circuitParamValue = 0.0
    circuitParamValue = xyceObj.getCircuitValue( "R1:R")
    print( "==> R1 R = %f" % (circuitParamValue))
    circuitParamValue = xyceObj.getCircuitValue( "TEMP")
    print( "==> Circuit temperature is now = %f" % (circuitParamValue))
 
  print("==> Taking step %d of %d" % (i, numSteps))
  requested_time = 0.0 + (i+1) * deltaTime
  actual_time = 0.0
  (result, actual_time) = xyceObj.simulateUntil( requested_time )

# get some measure results at the end of the run

maxV1=xyceObj.getCircuitValue( "MAXV1")
minV1=xyceObj.getCircuitValue( "MINV1")
maxV2=xyceObj.getCircuitValue( "MAXV2")
minV2=xyceObj.getCircuitValue( "MINV2")

print("==> Measure results from Xyce V1 = %f, %f  V2 = %f, %f" % (maxV1, minV1, maxV2, minV2))
print( "calling close")
xyceObj.close()

