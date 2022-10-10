from xyce_interface import xyce_interface

from flask import Flask
from flask import request

import os

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
  return "Active and responding.\n"
  
@app.route("/status")
def status():
  print(XyceObjectsDict)
  # probably need to return this as a json object
  return list(XyceObjectsDict.keys())
  
@app.route("/xyce_open", methods=['POST'])  
def open():
  args=request.get_json()
  if ('libdir' not in args):
  	return 'libdir not specified.',400
  print(args)
  # use a random 16 digit ID as a session ID
  # convert it to a string for easier rep. in json
  uuid = str(os.urandom(16))[2:-1].replace('\\','').replace('"','').replace('\'','')
  libDirectory=args['libdir']
  #print(libDirectory)
  xyceObj = xyce_interface(libdir=libDirectory)
  print(xyceObj)
  XyceObjectsDict[uuid] = {'libdir': libDirectory, 'xyceObj': xyceObj}
  return uuid
  
@app.route("/xyce_initialize", methods=['POST'])
def initialize():
  args=request.get_json()
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
  return '',200
  
@app.route("/xyce_getsimtime", methods=['POST'])
def getTime():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  time = xyceObj.getSimTime()
  return time
  
@app.route("/xyce_getfinaltime", methods=['POST'])
def getFinalTime():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  time = xyceObj.getFinalTime()
  return time

@app.route("/xyce_getdacnames", methods=['POST'])
def getDACDeviceNames():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status,dacNames) = xyceObj.getDACDeviceNames()
  return (status,dacNames)
  
  
@app.route("/xyce_checkcircuitparamexists", methods=['POST'])
def checkCircuitParameterExists():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in XyceObjectsDict):
    return 'paramname not supplied., 400'
  paramName = [ args['paramname'] ]
  result = xyceObj.checkCircuitParameterExists(paramName)
  return result
  
@app.route("/xyce_getadcmap", methods=['POST'])
def getADCMap():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status, ADCnames, widths, resistances, upperVLimits, lowerVLimits, settlingTimes) = xyceObj.getADCMap()
  return (status, ADCnames, widths, resistances, upperVLimits, lowerVLimits, settlingTimes)

@app.route("/xyce_getcircuitvalue", methods=['POST'])
def getCircuitValue():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in XyceObjectsDict):
    return 'paramname not supplied., 400'
  paramName = [ args['paramname'] ]
  result = xyceObj.getCircuitValue(paramName)
  return result

@app.route("/xyce_setcircuitparameter", methods=['POST'])
def setCircuitParameter():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'paramname' not in XyceObjectsDict):
    return 'paramname not supplied., 400'
  paramName = [ args['paramname'] ]
  if( 'paramval' not in XyceObjectsDict):
    return 'paramval not supplied., 400'
  paramValue = [ args['paramval'] ]
  result = xyceObj.setCircuitParameter(paramName, paramValue)
  return result
  
@app.route("/xyce_gettimevoltagepairs", methods=['POST'])
def getTimeVoltagePairs():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  (status, ADCnames, numADCnames, numPointsArray, timeArray, voltageArray)  = xyceObj.getTimeVoltagePairsADCLimitData
  return (status, ADCnames, numADCnames, numPointsArray, timeArray, voltageArray)

@app.route("/xyce_updatetimevoltagepairs", methods=['POST'])
def updateTimeVoltagePairs():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'devname' not in XyceObjectsDict):
    return 'devname not supplied., 400'
  devName = [ args['devname'] ]  
  if( 'timearray' not in XyceObjectsDict):
    return 'timearray not supplied., 400'
  timeArray = [ args['timearray'] ]
  if( 'voltarray' not in XyceObjectsDict):
    return 'voltarray not supplied., 400'
  voltArray = [ args['voltarray'] ]
  result = xyceObj.updateTimeVoltagePairs( devName, timeArray, voltArray)

  
@app.route("/xyce_simulateuntil", methods=['POST'])
def simulateUntil():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj == None):
    return 'No existing Xyce Object', 200
  if( 'simtime' not in XyceObjectsDict):
    return 'simtime not supplied., 400'
  simTime = [ args['simtime'] ]
  (status, simulatedTime ) = xyceObj.simulateUntil( simTime )
  return (status, simulatedTime )


@app.route("/xyce_run", methods=['POST'])
def run():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj != None):
    xyceObj.runSimulation()
  return '',200
    
@app.route("/xyce_close", methods=['POST'])
def close():
  args=request.get_json()
  if( 'uuid' not in args):
    return 'Session ID not specified.', 400
  theKey = args['uuid']
  if( theKey not in XyceObjectsDict):
    return 'Invalid Session ID., 400'
  
  xyceObj = XyceObjectsDict[theKey]['xyceObj']
  if( xyceObj != None):
    xyceObj.close()
  del XyceObjectsDict[theKey]
  return '',200
  
@app.route("/xyce_closeall")
def closeall():
  for aKey in XyceObjectsDict.keys():
    if (XyceObjectsDict[aKey]['xyceObj'] != None):
      XyceObjectsDict[aKey]['xyceObj'].close()
  XyceObjectsDict.clear()
  return '',200
  

