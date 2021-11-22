

# Python wrapper on Xyce via Xyce library mode via ctypes

import sys,traceback,types
#import ctypes 
from ctypes import *
from ctypes.util import *

class xyce_interface:
  def __init__(self,libdir="",name="",cmdargs=None):
    try:
      libName=find_library('xycecinterface')
      if( libName != None ):
        self.lib = cdll.LoadLibrary(libName)
        #self.lib = CDLL(libName)
        print( self.lib )
      else:
        # library wasn't found on normal system paths so 
        # try appending libdir to the name
        if sys.platform.startswith('darwin'):
          libName=libdir + "/" + "libxycecinterface.dylib"
          print("Trying to load " + libName )
          self.lib = CDLL(libName,RTLD_GLOBAL)
        else:
          libName=libdir + "/" + "libxycecinterface.so"
          print("Trying to load " + libName )
          self.lib = CDLL(libName,RTLD_GLOBAL)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise OSError("Could not load libxyce dynamic library")
    if cmdargs:
      self.xycePtr = c_void_p()
      self.lib.xyce_open(byref(self.xycePtr))
    else:
      self.xycePtr = c_void_p()
      self.lib.xyce_open(byref(self.xycePtr))

  def __del__(self):
    if self.xycePtr: 
      self.lib.xyce_close(byref(self.xycePtr))

  def close(self):
    self.lib.xyce_close(byref(self.xycePtr))
    self.xycePtr = None

  def initialize(self, args):
    # convert args to a c-array that can be passed to the
    # underlying code
    args.insert(0, "xyce_interface.py")
    narg = len(args)
    # allocate new C array 
    cargs = (c_char_p*narg)()
    # initialize the array
    for i in range(0,narg):
      cargs[i] = args[i].encode('utf-8')
    print( cargs )
    status = self.lib.xyce_initialize( byref(self.xycePtr), narg, cargs )
    return status

  def runSimulation( self ):
    status = self.lib.xyce_runSimulation( byref(self.xycePtr) )
    return status

  def simulateUntil( self, requestedTime ):
    csimTime = c_double(0)
    status = self.lib.xyce_simulateUntil( byref(self.xycePtr), c_double(requestedTime), byref(csimTime) )
    simulatedTime = csimTime.value
    return (status, simulatedTime )

  def getNumDevices( self, basename):
    cBaseName = c_char_p(basename.encode('utf-8'))
    cNumDevices = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    status = self.lib.xyce_getNumDevices( byref(self.xycePtr), cBaseName, byref(cNumDevices), byref(cMaxDeviceNameLength) )
    numDevices = cNumDevices.value
    maxDeviceNameLength = cMaxDeviceNameLength.value
    return (status, numDevices, maxDeviceNameLength)

  def getTotalNumDevices( self):
    cNumDevices = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    status = self.lib.xyce_getTotalNumDevices( byref(self.xycePtr), byref(cNumDevices), byref(cMaxDeviceNameLength) )
    numDevices = cNumDevices.value
    maxDeviceNameLength = cMaxDeviceNameLength.value
    return (status, numDevices, maxDeviceNameLength)

  def getDeviceNames( self, basename):
    # calling xyce_getDeviceNames(void ** ptr, char * modelGroupName, int & numDevNames, char ** deviceNames)
    cBaseName = c_char_p(basename.encode('utf-8'))
    cNumDeviceNames = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    names = []

    status = self.lib.xyce_getNumDevices( byref(self.xycePtr), cBaseName, byref(cNumDeviceNames), byref(cMaxDeviceNameLength) )
    # if the call to xyce_getNumDevices() fails then return an empty array
    if status != 1:
      return (status, names)

    deviceNameBuff = [create_string_buffer(cMaxDeviceNameLength.value) for i in range(cNumDeviceNames.value)]
    cDeviceNameArray = (c_char_p*cNumDeviceNames.value)(*map(addressof, deviceNameBuff))
    #print( cDeviceNameArray )
    status = self.lib.xyce_getDeviceNames( byref(self.xycePtr), cBaseName, byref(cNumDeviceNames), cDeviceNameArray)
    #print( cNumDeviceNames.value, cDeviceNameArray[0],  cDeviceNameArray[1])
    for i in range(0, cNumDeviceNames.value):
      names.insert(i, cDeviceNameArray[i] )
    return (status, names)

  def getAllDeviceNames( self):
    cNumDeviceNames = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    names = []

    status = self.lib.xyce_getTotalNumDevices( byref(self.xycePtr),  byref(cNumDeviceNames), byref(cMaxDeviceNameLength) )
    # if the call to xyce_getTotalNumDevices() fails then return an empty array
    if status != 1:
      return (status, names)

    deviceNameBuff = [create_string_buffer(cMaxDeviceNameLength.value) for i in range(cNumDeviceNames.value)]
    cDeviceNameArray = (c_char_p*cNumDeviceNames.value)(*map(addressof, deviceNameBuff))
    status = self.lib.xyce_getAllDeviceNames( byref(self.xycePtr), byref(cNumDeviceNames), cDeviceNameArray)
    for i in range(0, cNumDeviceNames.value):
      names.insert(i, cDeviceNameArray[i] )
    return (status, names)

  def getDACDeviceNames( self ):
    basename = "YADC"
    cBaseName = c_char_p(basename.encode('utf-8'))
    cNumDeviceNames = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    DACnames = []

    # if the call to xyce_getNumDevices() fails then return an empty array
    status = self.lib.xyce_getNumDevices( byref(self.xycePtr), cBaseName, byref(cNumDeviceNames), byref(cMaxDeviceNameLength) )
    if status != 1:
      return (status, DACnames)

    deviceNameBuff = [create_string_buffer(cMaxDeviceNameLength.value) for i in range(cNumDeviceNames.value)]
    cDeviceNameArray = (c_char_p*cNumDeviceNames.value)(*map(addressof, deviceNameBuff))
    #print( cDeviceNameArray )
    status = self.lib.xyce_getDACDeviceNames( byref(self.xycePtr), byref(cNumDeviceNames), cDeviceNameArray)
    #print( cNumDeviceNames.value, cDeviceNameArray[0],  cDeviceNameArray[1])
    for i in range(0, cNumDeviceNames.value):
      DACnames.insert(i, cDeviceNameArray[i] )
    return (status, DACnames)
    return status

  def checkDeviceParamName(self , paramName):
    cvarName = c_char_p(paramName)
    status = self.lib.xyce_checkDeviceParamName( byref(self.xycePtr), cvarName )
    return status

  def getDeviceParamVal(self, paramName):
    cparamName = c_char_p(paramName.encode('utf-8'))
    cValue = c_double(0.0)
    status = self.lib.xyce_getDeviceParamVal( byref(self.xycePtr), cparamName, byref(cValue) )
    return (status, (cValue.value))

  def getNumAdjNodesForDevice( self, deviceName):
    cdeviceName = c_char_p(deviceName.encode('utf-8'))
    cNumAdjNodes = c_int(0)
    status = self.lib.xyce_getNumAdjNodesForDevice( byref(self.xycePtr), cdeviceName, byref(cNumAdjNodes) )
    return (status, (cNumAdjNodes.value))

  def getAdjGIDsForDevice( self, deviceName):
    cdeviceName = c_char_p(deviceName.encode('utf-8'))
    cNumAdjNodes = c_int(0)
    GIDs=[]

    # if the call to xyce_getNumAdjNodesForDevice() fails then return an empty array
    status = self.lib.xyce_getNumAdjNodesForDevice( byref(self.xycePtr), cdeviceName, byref(cNumAdjNodes) )
    if status != 1:
      return (status, GIDs)

    cGIDs = (c_int*cNumAdjNodes.value)()
    status = self.lib.xyce_getAdjGIDsForDevice( byref(self.xycePtr), cdeviceName, byref(cNumAdjNodes), byref(cGIDs) )

    for i in range(0,cNumAdjNodes.value):
      GIDs.insert(i,cGIDs[i])
    return (status, GIDs)

  def getADCMap( self ):
    basename = "YADC"
    cBaseName = c_char_p(basename.encode('utf-8'))
    cNumDeviceNames = c_int( 0 )
    cMaxDeviceNameLength = c_int( 0 )
    ADCnames=[]
    widths=[]
    resistances=[]
    upperVLimits=[]
    lowerVLimits=[]
    settlingTimes=[]

    status = self.lib.xyce_getNumDevices( byref(self.xycePtr), cBaseName, byref(cNumDeviceNames), byref(cMaxDeviceNameLength) )
    # if the call to xyce_getNumDevices() fails then return empty arrays
    if status != 1:
      return (status, ADCnames, widths, resistances, upperVLimits, lowerVLimits, settlingTimes)

    deviceNameBuff = [create_string_buffer(cMaxDeviceNameLength.value) for i in range(cNumDeviceNames.value)]
    cADCNames = (c_char_p*cNumDeviceNames.value)(*map(addressof, deviceNameBuff))
    cWidths = (c_int*cNumDeviceNames.value)()
    cResistances = (c_double*cNumDeviceNames.value)()
    cUpperVLimits = (c_double*cNumDeviceNames.value)()
    cLowerVLimits = (c_double*cNumDeviceNames.value)()
    cSettlingTimes = (c_double*cNumDeviceNames.value)()

    status = self.lib.xyce_getADCMap( byref(self.xycePtr), byref(cNumDeviceNames), cADCNames,
                                      byref(cWidths), byref(cResistances),
                                      byref(cUpperVLimits), byref(cLowerVLimits),
                                      byref(cSettlingTimes) )

    for i in range(0, cNumDeviceNames.value):
      ADCnames.insert(i,cADCNames[i])
      widths.insert(i,cWidths[i])
      resistances.insert(i,cResistances[i])
      upperVLimits.insert(i,cUpperVLimits[i])
      lowerVLimits.insert(i,cLowerVLimits[i])
      settlingTimes.insert(i,cSettlingTimes[i])

    return (status, ADCnames, widths, resistances, upperVLimits, lowerVLimits, settlingTimes)

  def updateTimeVoltagePairs( self, basename,  time, voltage):
    cBaseName = c_char_p(basename)
    if( len( time ) != len( voltage ) ):
      print( "Time and Voltage arrays passed to updateTimeVoltagePairs are not of the same length.")
      return -1
    cNumPts = c_int(len(time))
    cArray = c_double*len(time)
    cTimeArray = (c_double*len(time))(*time)
    cVoltageArray = (c_double*len(time))(*voltage)

    status = self.lib.xyce_updateTimeVoltagePairs( byref(self.xycePtr), cBaseName, cNumPts, byref(cTimeArray), byref(cVoltageArray) )
    return status

  def checkResponseVarName(self , varName):
    cvarName = c_char_p(varName)
    status = self.lib.xyce_checkResponseVar( byref(self.xycePtr), cvarName )
    return status

  def obtainResponse(self, varName):
    cvarName = c_char_p(varName.encode('utf-8'))
    cValue = c_double(0.0)
    status = self.lib.xyce_obtainResponse( byref(self.xycePtr), cvarName, byref(cValue) )
    return (status, (cValue.value))

  def getTimeVoltagePairsADCsz( self ):
    cNumPoints = c_int(0)
    status = self.lib.xyce_getTimeVoltagePairsADCsz( byref(self.xycePtr), byref(cNumPoints) )
    numPoints = (cNumPoints.value)
    return (status, numPoints)

  def getTimeVoltagePairsADC( self ):
    cNumADCnames = c_int(0)
    cNumPoints = c_int(0)
    maxNumDevices=1000
    maxDeviceNameLength=1000
    maxNumTimeVoltagePairs=1000
    ADCDeviceNameBuff = [create_string_buffer(maxDeviceNameLength) for i in range(maxNumDevices)]
    cADCDeviceNameArray = (c_char_p*maxNumDevices)(*map(addressof, ADCDeviceNameBuff))

    #make the two double** arrays
    double_ctime = POINTER(c_double)
    inner_ctime_array = (c_double * maxNumTimeVoltagePairs)
    cTimeArray = (double_ctime * maxNumDevices) ()
    for i in range(maxNumDevices):
       cTimeArray[i] = inner_ctime_array()

    double_cvoltage = POINTER(c_double)
    inner_cvoltage_array = (c_double * maxNumTimeVoltagePairs)
    cVoltageArray = (double_cvoltage * maxNumDevices) ()
    for i in range(maxNumDevices):
       cVoltageArray[i] = inner_cvoltage_array()
    
    status = self.lib.xyce_getTimeVoltagePairsADC( byref(self.xycePtr), byref(cNumADCnames), cADCDeviceNameArray, byref(cNumPoints), byref(cTimeArray), byref(cVoltageArray) )

    # make the integer return values
    numADCnames = (cNumADCnames.value)
    numPoints = (cNumPoints.value)

    # make and populate the arrays that will be used in
    # the return statement
    ADCnames = []
    for i in range(0, numADCnames):
      ADCnames.insert(i, cADCDeviceNameArray[i] ) 

    # make the returned timeArray
    double_time = POINTER(c_double)
    inner_time_array = (c_double * numPoints)
    timeArray = (double_time * numADCnames) ()
    for i in range(numADCnames):
       timeArray[i] = inner_time_array()

    # copy over the elements from cTimeArray to timeArray
    for i in range(numADCnames):
      for j in range(numPoints):
       timeArray[i][j] = cTimeArray[i][j]

    # make the returned voltageArray
    double_voltage = POINTER(c_double)
    inner_voltage_array = (c_double *  numPoints)
    voltageArray = (double_voltage * numADCnames) ()
    for i in range((cNumADCnames.value)):
       voltageArray[i] = inner_voltage_array()

    # copy over the elements from cVoltageArray to voltageArray
    for i in range(numADCnames):
      for j in range(numPoints):
       voltageArray[i][j] = cVoltageArray[i][j]

    return (status, ADCnames, numADCnames, numPoints, timeArray, voltageArray)

  def getTimeStatePairsADC( self ):
    cNumADCnames = c_int(0)
    cNumPoints = c_int(0)
    maxNumDevices=1000
    maxDeviceNameLength=1000
    maxNumTimeStatePairs=1000
    ADCDeviceNameBuff = [create_string_buffer(maxDeviceNameLength) for i in range(maxNumDevices)]
    cADCDeviceNameArray = (c_char_p*maxNumDevices)(*map(addressof, ADCDeviceNameBuff))

    #make the double** and int** arrays
    double_ctime = POINTER(c_double)
    inner_ctime_array = (c_double * maxNumTimeStatePairs)
    cTimeArray = (double_ctime * maxNumDevices) ()
    for i in range(maxNumDevices):
       cTimeArray[i] = inner_ctime_array()

    int_cstate = POINTER(c_int)
    inner_cstate_array = (c_int * maxNumTimeStatePairs)
    cStateArray = (int_cstate * maxNumDevices) ()
    for i in range(maxNumDevices):
       cStateArray[i] = inner_cstate_array()

    status = self.lib.xyce_getTimeStatePairsADC( byref(self.xycePtr), byref(cNumADCnames), cADCDeviceNameArray, byref(cNumPoints), byref(cTimeArray), byref(cStateArray) )

    # make the integer return values
    numADCnames = (cNumADCnames.value)
    numPoints = (cNumPoints.value)

    # make and populate the arrays that will be used in
    # the return statement
    ADCnames = []
    for i in range(0, numADCnames):
      ADCnames.insert(i, cADCDeviceNameArray[i] ) 

    # make the returned timeArray
    double_time = POINTER(c_double)
    inner_time_array = (c_double * numPoints)
    timeArray = (double_time * numADCnames) ()
    for i in range(numADCnames):
       timeArray[i] = inner_time_array()

    # copy over the elements from cTimeArray to timeArray
    for i in range(numADCnames):
      for j in range(numPoints):
       timeArray[i][j] = cTimeArray[i][j]

    # make the returned stateArray
    double_state = POINTER(c_int)
    inner_state_array = (c_int *  numPoints)
    stateArray = (double_state * numADCnames) ()
    for i in range((cNumADCnames.value)):
       stateArray[i] = inner_state_array()

    # copy over the elements from cStateArray to stateArray
    for i in range(numADCnames):
      for j in range(numPoints):
       stateArray[i][j] = cStateArray[i][j]

    return (status, ADCnames, numADCnames, numPoints, timeArray, stateArray)

  def setADCWidths( self, ADCnames, ADCwidths ):
    # make a char ** structure for ADCnames and an int[] for ADCwidths
    if( len(ADCnames) != len(ADCwidths) ):
      print( "Names and widths arrays passed to setADCWidths are not of the same length.")
      return -1
    nADCs = len(ADCnames)
    cADCnames = (c_char_p*nADCs)()
    for i in range(0,nADCs):
      cADCnames[i] = ADCnames[i].encode('utf-8')

    cADCwidths = (c_int*nADCs)()
    for i in range(0,nADCs):
      cADCwidths[i]=ADCwidths[i]

    status = self.lib.xyce_setADCWidths( byref(self.xycePtr), nADCs, byref(cADCnames), byref(cADCwidths) )
    return status

  def getADCWidths( self, ADCnames ):
    nADCs = len(ADCnames)
    cADCnames = (c_char_p*nADCs)()
    for i in range(0,nADCs):
      cADCnames[i] = ADCnames[i].encode('utf-8')

    cADCwidths = (c_int*nADCs)()

    status = self.lib.xyce_getADCWidths( byref(self.xycePtr), nADCs, byref(cADCnames), byref(cADCwidths) )

    width = []
    for i in range(0,nADCs):
      width.insert(i,cADCwidths[i])
    return (status,width)     


#   def file(self,file):
#     self.lib.spparks_file(self.xycePtr,file)
# 
#   def command(self,cmd):
#     self.lib.spparks_command(self.xycePtr,cmd)
# 
#   def extract(self,name,type):
#     if type == 0:
#       self.lib.spparks_extract.restype = POINTER(c_int)
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr[0]
#     if type == 1:
#       self.lib.spparks_extract.restype = POINTER(c_int)
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr
#     if type == 2:
#       self.lib.spparks_extract.restype = POINTER(POINTER(c_int))
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr
#     if type == 3:
#       self.lib.spparks_extract.restype = POINTER(c_double)
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr[0]
#     if type == 4:
#       self.lib.spparks_extract.restype = POINTER(c_double)
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr
#     if type == 5:
#       self.lib.spparks_extract.restype = POINTER(POINTER(c_double))
#       ptr = self.lib.spparks_extract(self.xycePtr,name)
#       return ptr
#     return None
# 
#   def energy(self):
#     self.lib.spparks_energy.restype = c_double
#     return self.lib.spparks_energy(self.xycePtr)
