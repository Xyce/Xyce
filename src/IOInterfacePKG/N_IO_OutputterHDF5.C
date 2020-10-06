//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputMgr.h>

#include <N_PDS_ParMap.h>
#include <N_UTL_ExtendedString.h>

#ifdef Xyce_USE_HDF5

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

#endif // Xyce_USE_HDF5

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputMgr::prepareHDF5Output
// Purpose       : creates HDF5 output object and prepares for writing
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool OutputMgr::prepareHDF5Output(Parallel::Machine comm)
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  hdf5PlistId_ = H5Pcreate(H5P_FILE_ACCESS);

  // based on if we're running in serial or parallel, add parallel output
  // option to hdf5PlistId_
  if (Parallel::size(comm) > 1)
  {
#ifdef Xyce_PARALLEL_MPI
    H5Pset_fapl_mpio(hdf5PlistId_, comm, MPI_INFO_NULL);
    H5Pset_fclose_degree(hdf5PlistId_, H5F_CLOSE_STRONG);
#endif
  }
  std::string hdf5ExtendedFileName = hdf5FileName_ + filenameSuffix_ + ".h5d";
  hdf5FileId_ = H5Fcreate(hdf5ExtendedFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, hdf5PlistId_);
#endif  // Xyce_USE_HDF5

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::updateHDF5Output
// Purpose       : writes curent solution data to the HDF5 output object
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool
OutputMgr::updateHDF5Output(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVecPtr)
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  int status=0;
  if (!hdf5HeaderWritten_)
  {
    hdf5HeaderWritten_ = true;

    // in parallel we can't write an array of variale length node names.
    // Thus we need to make an array of fixed length strings the longest
    // of which is the maximum length of a node name.  This wastes a bit of
    // output space, but we only do this once in the solution var map.

    // first on each processor, find the max length of a node name.
    int globalMaxNodeNameLength = 0;
    for (NodeNameMap::const_iterator nameIter = getSolutionNodeMap().begin(), nameEnd = getSolutionNodeMap().end(); nameIter != nameEnd ; ++nameIter)
    {
      int nodeNameLength = nameIter->first.size();
      globalMaxNodeNameLength = std::max(globalMaxNodeNameLength, nodeNameLength);
    }

    Parallel::AllReduce(comm, MPI_MAX, &globalMaxNodeNameLength, 1);

    globalMaxNodeNameLength++;  // add one for string terminator
    // now make a char[][] array of size [numLocalNodes][maxNodeNameLength]
    int numLocalNodes = solnVecPtr.localLength();
    char* nodeNameSpace = new char[ numLocalNodes * globalMaxNodeNameLength ];
    // could simplify this with a placement new operator(i.e. new nodeNameSpace nodeNameArray).
    char** nodeNameArray = new char * [ numLocalNodes ];
    for (int i=0; i< numLocalNodes; i++)
    {
      nodeNameArray[i] = nodeNameSpace + i * sizeof(char) * globalMaxNodeNameLength;//new char(globalMaxNodeNameLength);
    }
    // also need an array for global ID's
    int* nodeGIDArray = new int[ numLocalNodes ];

    // to fill the arrays in parallel, each processor needs to know how many
    // unknowns are on prior processors.  We could infer this from the GID values,
    // however that would not work in block analysis modes like MPDE, HB and AC
    // Thus, we will have the individual processes communicate this info .
    // This duplicates what could be done with an epetra multivector, if we could
    // construct one that held strings or char*.  We can't so for now this is what
    // we'll do.
    std::vector<int> unknownsPerProc(Parallel::size(comm));

    Parallel::AllGather(comm, numLocalNodes, unknownsPerProc);
    std::vector<int> totalUnknownsPriorToProc;
    totalUnknownsPriorToProc.resize(Parallel::size(comm));
    int sum = 0;
    for (int i=0; i<Parallel::size(comm); i++)
    {
      totalUnknownsPriorToProc[i] = sum;
      sum += unknownsPerProc[i];
    }

    // fill up the arrays
    int index = 0;
    for (NodeNameMap::const_iterator nameIter = getSolutionNodeMap().begin(), nameEnd = getSolutionNodeMap().end(); nameIter != nameEnd ; ++nameIter, index++)
    {
      strncpy(nodeNameArray[index], (*nameIter).first.c_str(), globalMaxNodeNameLength);
      nodeGIDArray[index] = solnVecPtr.pmap()->localToGlobalIndex((*nameIter).second);
    }

    // make the group for writing
    // in parallel, all processors must create groups and datasets.
    hid_t solVarMapGroup = H5Gcreate(hdf5FileId_, "SolVarMap", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (solVarMapGroup < 0)
      Xyce::dout() << "Error in making solVarMapGroup on procID_ " << Parallel::rank(comm) << std::endl;

    // make the file data space
    hsize_t nodeNamesDim[1] = {solnVecPtr.globalLength()};
    hid_t fileDataspace = H5Screate_simple(1,  nodeNamesDim, NULL);

    // make the memorySpace
    hsize_t nodeNamesLocalDim[1] = {solnVecPtr.localLength()};
    hid_t memorySpace = H5Screate_simple(1, nodeNamesLocalDim, NULL);

    // make the dataset
    hid_t fixedLenString = H5Tcreate(H5T_STRING, globalMaxNodeNameLength);
    hid_t nodeNameDataSet = H5Dcreate(solVarMapGroup, "VarNames", fixedLenString, fileDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // select portion for writing
    hid_t fileDataspaceSelection = H5Dget_space(nodeNameDataSet);
    hsize_t offset[1] = {totalUnknownsPriorToProc[Parallel::rank(comm)]};
    hsize_t stride[1] = {1};
    hsize_t count[1] = {unknownsPerProc[Parallel::rank(comm)]};
    hsize_t block[1] = {1};
    H5Sselect_hyperslab(fileDataspaceSelection, H5S_SELECT_SET, offset, stride, count, block);

    // writing property list
    hid_t writePropertyList = H5Pcreate(H5P_DATASET_XFER);

#ifdef Xyce_PARALLEL_MPI
    if (Parallel::size(comm) > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
#endif

    status = H5Dwrite(nodeNameDataSet, fixedLenString, memorySpace, fileDataspaceSelection, writePropertyList, nodeNameSpace);

    // write the GID array
    hid_t gidDataSet = H5Dcreate(solVarMapGroup, "VarGID", H5T_NATIVE_INT, fileDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(fileDataspaceSelection, H5S_SELECT_SET, offset, stride, count, block);
    status = H5Dwrite(gidDataSet, H5T_NATIVE_INT, memorySpace, fileDataspaceSelection, writePropertyList, nodeGIDArray);
    status = H5Dclose(gidDataSet);

    delete [] nodeNameArray;
    delete [] nodeGIDArray;

    status = H5Tclose(fixedLenString);
    if (status < 0)
      Xyce::dout() << "Error closing fixedLenString. on " << Parallel::rank(comm) << std::endl;
    status = H5Dclose(nodeNameDataSet);
    if (status < 0)
      Xyce::dout() << "Error closing nodeNameDataSet. on " << Parallel::rank(comm) << std::endl;
    status = H5Sclose(fileDataspace);
    if (status < 0)
      Xyce::dout() << "Error closing solVarMapDataspace. on " << Parallel::rank(comm) << std::endl;
    status = H5Gclose(solVarMapGroup);
    if (status < 0)
      Xyce::dout() << "Error closing solVarMapGroup. on " << Parallel::rank(comm) << std::endl;

    // now output the independent variables
    hid_t independentVarGroup = H5Gcreate(hdf5FileId_, "IndepVars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (independentVarGroup < 0)
      Xyce::dout() << "Error in making independentVarGroup on Parallel::rank(comm) " << Parallel::rank(comm) << std::endl;

    // filespace & memspace
    hsize_t independentDataDim[1] = {1};
    hsize_t independentMaxDataDim[1] = {H5S_UNLIMITED};
    hid_t independentFileSpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);
    hid_t independentMemorySpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);

    // array will be unlimited in length because we don't know how many steps Xyce will take
    // So we've set the max dimension to H5S_UNLIMITED above. we need to set the size of the
    // chunk by which it will grow as well
    // now the dataset
    hid_t independentDataProperty_ = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunkSize[1] = {1};
    H5Pset_chunk(independentDataProperty_, 1, chunkSize);

    // make the dataset
    hid_t independentVarDataSet = H5Dcreate(independentVarGroup, "Time", H5T_NATIVE_DOUBLE, independentFileSpace, H5P_DEFAULT, independentDataProperty_, H5P_DEFAULT);

    // write
    writePropertyList=H5Pcreate(H5P_DATASET_XFER);

#ifdef Xyce_PARALLEL_MPI
    if (Parallel::size(comm) > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_INDEPENDENT);
    }
#endif

    if (Parallel::rank(comm) == 0)
    {
      H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentMemorySpace, independentFileSpace, writePropertyList, &outputState_.circuitTime_);
    }
    H5Sclose(independentFileSpace);
    H5Sclose(independentMemorySpace);
    H5Dclose(independentVarDataSet);
    H5Gclose(independentVarGroup);
    H5Pclose(independentDataProperty_);

    // now output the first solution vector
    hid_t solutionVecGroup = H5Gcreate(hdf5FileId_, "SolutionVec", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (solutionVecGroup < 0)
      Xyce::dout() << "Error in making solutionVecGroup on Parallel::rank(comm) " << Parallel::rank(comm) << std::endl;

    // filespace & memspace
    hsize_t dependentDataLocalDim[2] = {1, solnVecPtr.localLength()};
    hsize_t dependentDataGlobalDim[2] = {1, solnVecPtr.globalLength()};
    hsize_t dependentMaxDataDim[2] = {H5S_UNLIMITED, solnVecPtr.globalLength()};
    // remember, filespace is the size of the global data on the disk.
    hid_t dependentFileSpace = H5Screate_simple(2, dependentDataGlobalDim, dependentMaxDataDim);
    // memspace defines the space used locally on a given processor.
    hid_t dependentMemorySpace = H5Screate_simple(2, dependentDataLocalDim, dependentMaxDataDim);

    // array will be unlimited in length because we don't know how many steps Xyce will take
    // So we've set the max dimension to H5S_UNLIMITED above. we need to set the size of the
    // chunk by which it will grow as well
    // now the dataset
    hid_t dependentDataProperty_ = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t solChunkSize[2] = {1, solnVecPtr.globalLength()};
    H5Pset_chunk(dependentDataProperty_, 2, solChunkSize);

    // make the dataset
    hid_t dependentVarDataSet = H5Dcreate(solutionVecGroup, "Solution", H5T_NATIVE_DOUBLE, dependentFileSpace, H5P_DEFAULT, dependentDataProperty_, H5P_DEFAULT);

    // set up hyperslab to define the relationship between memspace and filespace.
    hid_t dependentVarFilespace = H5Dget_space(dependentVarDataSet);
    hsize_t solVecStart[2] = {0, solnVecPtr.pmap()->localToGlobalIndex(0) - solnVecPtr.pmap()->indexBase()};
    hsize_t solVecStride[2] = {1, 1};
    hsize_t solVecCount[2] = {1, solnVecPtr.localLength()};
    hsize_t solVecBlock[2] = {1, 1};
    H5Sselect_hyperslab(dependentFileSpace, H5S_SELECT_SET, solVecStart, solVecStride, solVecCount, solVecBlock);

    writePropertyList=H5Pcreate(H5P_DATASET_XFER);

#ifdef Xyce_PARALLEL_MPI
    if (Parallel::size(comm) > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
#endif

    // write the data
    H5Dwrite(dependentVarDataSet, H5T_NATIVE_DOUBLE, dependentMemorySpace, dependentFileSpace, writePropertyList, &solnVecPtr[0]);

    H5Pclose(writePropertyList);
    H5Sclose(dependentFileSpace);
    H5Sclose(dependentMemorySpace);
    H5Dclose(dependentVarDataSet);
    H5Gclose(solutionVecGroup);
  }
  else
  {
    // update  independent variables
    hid_t independentVarGroup = H5Gopen(hdf5FileId_, "IndepVars", H5P_DEFAULT);
    hid_t independentVarDataSet = H5Dopen(independentVarGroup, "Time", H5P_DEFAULT);

    // data sets have been written to once.  Now extend them with new data
    hsize_t newDim[1] = {hdf5IndexValue_+1};
    hsize_t coords[1] = {hdf5IndexValue_};

    hid_t writePropertyList=H5Pcreate(H5P_DATASET_XFER);

#ifdef Xyce_PARALLEL_MPI
    if (Parallel::size(comm) > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_INDEPENDENT);
    }
#endif

    H5Dset_extent(independentVarDataSet, newDim);

    // make the memspace
    hsize_t independentDataDim[1] = {1};
    hsize_t independentMaxDataDim[1] = {H5S_UNLIMITED};
    hid_t independentVarMemSpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);
    //hid_t independentVarMemSpace = H5Dget_space(independentVarDataSet);
    // get filespace
    hid_t independentVarFilespace = H5Dget_space(independentVarDataSet);

    // indicate that we have only updated the last element in this space
    //H5Sselect_elements(independentVarFilespace, H5S_SELECT_SET, 1, coords);
    hsize_t start[1] = {hdf5IndexValue_};
    hsize_t stride[1] = {1};
    hsize_t count[1] = {1};
    hsize_t block[1] = {1};
    H5Sselect_hyperslab(independentVarFilespace, H5S_SELECT_SET, start, stride, count, block);
    // write the update to the data space

    if (Parallel::rank(comm) == 0)
    {
      //H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentVarMemSpace, independentVarFilespace, writePropertyList, &circuitTime_-hdf5IndexValue_);
      H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentVarMemSpace, independentVarFilespace, writePropertyList, &outputState_.circuitTime_);
    }

    H5Dclose(independentVarDataSet);
    H5Gclose(independentVarGroup);
    H5Sclose(independentVarMemSpace);
    H5Sclose(independentVarFilespace);

    // update solution
    hid_t solutionVecGroup = H5Gopen(hdf5FileId_, "SolutionVec", H5P_DEFAULT);
    hid_t dependentVarDataSet = H5Dopen(solutionVecGroup, "Solution", H5P_DEFAULT);

    // data sets have been written to once.  Now extend them with new data
    hsize_t solutionNewDim[2] = {hdf5IndexValue_+1, solnVecPtr.globalLength() };
    hsize_t solutionCoords[2] = {hdf5IndexValue_, 0};
    H5Dset_extent(dependentVarDataSet, solutionNewDim);

    // filespace & memspace
    hsize_t dependentDataLocalDim[2] = {1, solnVecPtr.localLength()};
    hsize_t dependentDataGlobalDim[2] = {1, solnVecPtr.globalLength()};
    hsize_t dependentMaxDataDim[2] = {H5S_UNLIMITED, solnVecPtr.globalLength()};
    // remember, filespace is the size of the global data on the disk.
    // hid_t dependentFileSpace = H5Screate_simple(2, dependentDataGlobalDim, dependentMaxDataDim);
    // memspace defines the space used locally on a given processor.
    hid_t dependentMemorySpace = H5Screate_simple(2, dependentDataLocalDim, dependentMaxDataDim);
    // get filespace
    hid_t dependentFileSpace = H5Dget_space(dependentVarDataSet);

    // set up hyperslab to define the relationship between memspace and filespace.
    hsize_t solVecStart[2] = {hdf5IndexValue_, solnVecPtr.pmap()->localToGlobalIndex(0) - solnVecPtr.pmap()->indexBase()};
    hsize_t solVecStride[2] = {1, 1};
    hsize_t solVecCount[2] = {1, solnVecPtr.localLength()};
    hsize_t solVecBlock[2] = {1, 1};
    H5Sselect_hyperslab(dependentFileSpace, H5S_SELECT_SET, solVecStart, solVecStride, solVecCount, solVecBlock);

    writePropertyList=H5Pcreate(H5P_DATASET_XFER);

#ifdef Xyce_PARALLEL_MPI
    if (Parallel::size(comm) > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
#endif

    // write the data
    H5Dwrite(dependentVarDataSet, H5T_NATIVE_DOUBLE, dependentMemorySpace, dependentFileSpace, writePropertyList, &solnVecPtr[0]);

    H5Pclose(writePropertyList);
    H5Sclose(dependentFileSpace);
    H5Sclose(dependentMemorySpace);
    H5Dclose(dependentVarDataSet);
    H5Gclose(solutionVecGroup);
  }

  // update index before return
  hdf5IndexValue_++;
#endif  // Xyce_USE_HDF5

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::closeHDF5Output
// Purpose       : closes HDF5 ouput object
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool OutputMgr::closeHDF5Output()
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  int status = 0;
  H5Pclose(hdf5PlistId_);
  status = H5Fclose(hdf5FileId_);
  if (status < 0)
      Xyce::dout() << "Error closing hdf5FileId_." << std::endl;
#endif  // Xyce_USE_HDF5

  return result;
}

} // namespace IO
} // namespace Xyce
