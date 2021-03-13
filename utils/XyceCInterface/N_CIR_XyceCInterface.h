//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        : C callable functions to link to Xyce's C++ interface
//
// Special Notes  : These functions expose Xyce's C++ interface to 
//                  other programs (i.e. Python) that can dynamically 
//                  link to a C library (via ctypes in Python)
//
//
//-----------------------------------------------------------------------------


/*
 simple class file for Xyce interface
*/

#ifdef __cplusplus
extern "C" {
#endif

void xyce_open( void** ptr );
void xyce_close( void** ptr);
int xyce_initialize( void** ptr, int narg, char ** argv);

int xyce_runSimulation(void ** ptr);

int xyce_simulateUntil(  void ** ptr, double requestedUntilTime, double* completedUntilTime );

int xyce_getNumDevices(void **ptr, char * modelGroupName, int* numDevNames, int* maxDevNameLength);
int xyce_getDeviceNames(void ** ptr, char * modelGroupName, int* numDevNames, char ** deviceNames);
int xyce_getTotalNumDevices(void **ptr, int* numDevNames, int* maxDevNameLength);
int xyce_getAllDeviceNames(void ** ptr, int* numDevNames, char ** deviceNames);
int xyce_getDACDeviceNames(void ** ptr, int* numDevNames, char ** deviceNames);

int xyce_checkDeviceParamName(void **ptr, char* full_param_name);
int xyce_getDeviceParamVal(void **ptr, char* full_param_name, double* value);

int xyce_getNumAdjNodesForDevice(void **ptr, char* deviceName, int* numAdjNodes);
int xyce_getAdjGIDsForDevice(void **ptr, char* deviceName, int* numAdjNodes, int* adjGIDs);

int xyce_updateTimeVoltagePairs(void ** ptr, char * DACname, int numPoints, double * timeArray, double * voltageArray);
                
int xyce_getTimeVoltagePairsADC( void** ptr, int * numADCnames, char ** ADCnames, int * numPoints, double ** timeArray, double ** voltageArray );
/* 
   Note: ADCnames, numPoints, timeArray and voltageArray must be preallocated!  this function cannot 
   allocate that storage or the ctypes layer breaks 
*/
int xyce_getTimeStatePairsADC( void** ptr, int * numADCnames, char ** ADCnames, int * numPoints, double ** timeArray, int ** stateArray ); 
/* 
   Note: ADCnames, numPoints, timeArray and stateArray must be preallocated!  this function cannot 
   allocate that storage or the ctypes layer breaks 
*/

int xyce_getADCMap(void ** ptr, int *numADCnames, char ** ADCnames, int *  widths, double * resistances,
	           double * upperVLimits, double * lowerVLimits, double * settlingTimes);

int xyce_setADCWidths(void ** ptr, int numADCnames, char ** ADCnames, int *  widths);
int xyce_getADCWidths(void ** ptr, int numADCnames, char ** ADCnames, int *  widths);

int xyce_checkResponseVar(void ** ptr, char * variable_name);

int xyce_obtainResponse(void ** ptr, char *  variable_name, double *value);


#ifdef __cplusplus
}
#endif

