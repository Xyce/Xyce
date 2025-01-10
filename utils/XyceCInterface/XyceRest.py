
#
# To run this from flask type:
#
# flask --app XyceRest.py run
#
# or on MacOS with MacPorts flask installed:
#
# flask-3.10 --app XyceRest.py run
#
from xyce_interface import xyce_interface

from flask import Flask
from flask import request

import os
import json 

app = Flask(__name__)

XyceObjectsDict = {}

# REST HTTP Status codes used on return
# 1xx - informational
#       100 client continue
#       101 client should upgrade to returned protocol
#       102 Processing (WebDAV) server is processing the request
#       103 Early Hints for preloading
# 2xx - Success
#       200 request succeeded
#       201 request successful and new resource has been created
#       202 request received but not completed 
#       203 returned info is local and not authoritative
#       204 success, but no content to return
#       205 reset content 
#       206 partial content
#       207 multi-status (WebDAV)
#       208 already reported (WebDAV)
#       226 IM Used
# 3xx - Redirection, Client must take some further action
#       300--308 reasons for redirecting to another server
# 4xx - Client error.  The client did something wrong
#       400 Bad request
#       401 Unauthorized 
#       402 Payment required 
#       403 Forbidden
#       404 Not found
#       405 Method not allowed
#       406 Not acceptable
#       407 Proxy authentication required
#       408 Request timeout 
#       409 Conflict
#       410 Gone 
#       411-499 obscure errors not needed to list
# 5xx - Server error. The server hit an error
#       500 Internal server error
#       501 Not implemented 
#       502 Bad gateway 
#       503 Service unavailable
#       504-511 More obscure server errors.



@app.route("/test")
def test_response():
  return "Active and responding."
  
@app.route("/status")
def status():
  print(XyceObjectsDict)
  # just return the number of allocated Xyce objects.
  return dict({'numInstance':len(XyceObjectsDict.keys())})
  
@app.route("/xyce_open", methods=['POST'])  
def open():
  args=json.loads(request.get_json())
  
  # use a random 16 digit ID as a session ID
  # convert it to a string for easier rep. in json
  uuid = str(os.urandom(16))[2:-1].replace('\\','').replace('"','').replace('\'','')
  xycdObj = None
  libDirectory = None
  if ('libdir' not in args):
    # use the interface's default lib directory.
  	xyceObj = xyce_interface()
  else:
    libDirectory=args['libdir']
    xyceObj = xyce_interface(libdir=libDirectory)
  
  XyceObjectsDict[uuid] = {'libdir': libDirectory, 'xyceObj': xyceObj}
  return dict({'uuid': uuid})
  
@app.route("/xyce_initialize", methods=['POST'])
def initialize():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  if( 'simfile' not in args):
    return 'simfile not specified.', 400
  simFile = [ args['simfile'] ]
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  xyceObj.initialize(simFile)
  return 'success',200
  
@app.route("/xyce_getsimtime", methods=['POST'])
def getTime():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  time = xyceObj.getSimTime()
  return dict({'time': time})
  
@app.route("/xyce_getfinaltime", methods=['POST'])
def getFinalTime():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  time = xyceObj.getFinalTime()
  return  dict({'time': time})

@app.route("/xyce_getdacnames", methods=['POST'])
def getDACDeviceNames():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status,dacNames) = xyceObj.getDACDeviceNames()
  return dict({'status':status, 'dacNames': dacNames})
  
@app.route("/xyce_getadcmap", methods=['POST'])
def getADCMap():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status, ADCnames, widths, resistances, upperVLimits, lowerVLimits, settlingTimes) = xyceObj.getADCMap()
  return dict({'status':status, 'ADCnames':ADCnames, 'widths':widths, 'resistances':resistances, 'upperVLimits':upperVLimits, 'lowerVLimits':lowerVLimits, 'settlingTimes':settlingTimes})

  
@app.route("/xyce_checkcircuitparamexists", methods=['POST'])
def checkCircuitParameterExists():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in args):
    return 'paramname not supplied., 400'
  paramName = args['paramname']
  result = xyceObj.checkCircuitParameterExists(paramName)
  return dict({'result': result})
  
@app.route("/xyce_getcircuitvalue", methods=['POST'])
def getCircuitValue():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in args):
    return 'paramname not supplied., 400'
  paramName = args['paramname']
  result = xyceObj.getCircuitValue(paramName)
  return dict({'value':result})

@app.route("/xyce_setcircuitparameter", methods=['POST'])
def setCircuitParameter():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in args):
    return 'paramname not supplied., 400'
  paramName = args['paramname']
  if( 'paramval' not in args):
    return 'paramval not supplied., 400'
  paramValue = args['paramval']
  result = xyceObj.setCircuitParameter(paramName, paramValue)
  return dict({'status':result})
  
@app.route("/xyce_gettimevoltagepairs", methods=['POST'])
def getTimeVoltagePairs():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status, ADCnames, numADCnames, numPointsArray, timeArray, voltageArray)  = xyceObj.getTimeVoltagePairsADCLimitData()
  timeArrayLinearized = []
  voltageArrayLinearized = []
  for i in range(0,numADCnames):
    for j in range(0, numPointsArray[i]):
      timeArrayLinearized.append( timeArray[i][j])
      voltageArrayLinearized.append( voltageArray[i][j])
  return dict({'status':status, 'ADCnames':ADCnames, 'numADCnames':numADCnames, 'numPointsInArray':numPointsArray, 'timeArray':timeArrayLinearized, 'voltageArray':voltageArrayLinearized})

@app.route("/xyce_updatetimevoltagepairs", methods=['POST'])
def updateTimeVoltagePairs():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'devname' not in args):
    return 'devname not supplied., 400'
  devName = args['devname']  
  if( 'timearray' not in args):
    return 'timearray not supplied., 400'
  timeArray = args['timearray']
  if( 'voltarray' not in args):
    return 'voltarray not supplied., 400'
  voltArray = args['voltarray']
  result = xyceObj.updateTimeVoltagePairs( devName, timeArray, voltArray)
  return dict({'status':result})

  
@app.route("/xyce_simulateuntil", methods=['POST'])
def simulateUntil():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'simtime' not in args):
    return 'simtime not supplied., 400'
  simTime = args['simtime']
  (status, simulatedTime ) = xyceObj.simulateUntil( simTime )
  return  dict({'status': status, 'simulatedTime': simulatedTime})


@app.route("/xyce_run", methods=['POST'])
def run():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj != None):
    xyceObj.runSimulation()
  return 'success',200
    
@app.route("/xyce_close", methods=['POST'])
def close():
  args=json.loads(request.get_json())
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj != None):
    xyceObj.close()
  del XyceObjectsDict[theKey]
  return 'success',200
  
@app.route("/xyce_closeall")
def closeall():
  for aKey in XyceObjectsDict.keys():
    if (XyceObjectsDict[aKey]['xyceObj'] != None):
      XyceObjectsDict[aKey]['xyceObj'].close()
  XyceObjectsDict.clear()
  return 'success',200
  

