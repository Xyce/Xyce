//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 12/11/09
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <algorithm>
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_ROM.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Util.h>

#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_IO_mmio.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>

// Keep redundant (assuming Trilinos and Xyce are configured in the same environment) definitions from whining
#undef HAVE_INTTYPES_H
#undef HAVE_STDINT_H
#include <Trilinos_Util.h>

namespace Xyce {
namespace Device {


namespace ROM {


void Traits::loadInstanceParameters(ParametricData<ROM::Instance> &p)
{
p.addPar("BASE_FILENAME", std::string("rom_input"), &ROM::Instance::baseFileName);
  p.addPar("MASK_VARS", false, &ROM::Instance::maskROMVars);
  p.addPar("USE_PORT_DESCRIPTION", 0, &ROM::Instance::usePortDesc);
}

void Traits::loadModelParameters(ParametricData<ROM::Model> &p)
{}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & ROMiter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(ROMiter),
    isCSparse(false),
    isGSparse(false),
    maskROMVars(false),
    numROMVars(0),
    baseFileName(""),
    dt(0),
    dt_last(0),
    alph(0),
    alph_last(0),
    coef(0),
    coefLast(0),
    currentOrder(0),
    usedOrder(0),
    lastTimeStepNumber(0)
{
  numExtVars   = IB.numExtVars;  // we have as many as were specified on the
                                 // instance line

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Read in reduced-order model.
  if (given("BASE_FILENAME"))
  {
    FILE *c_file, *g_file, *b_file, *l_file;
    Xyce::IO::MMIO::MM_typecode mat_code;
    int M=0, N=0, nz=0;
    std::string cfile = Instance::baseFileName + ".Chat";
    std::string gfile = Instance::baseFileName + ".Ghat";
    std::string bfile = Instance::baseFileName + ".Bhat";
    std::string lfile = Instance::baseFileName + ".Lhat";
    c_file = fopen(cfile.c_str(), "r");
    g_file = fopen(gfile.c_str(), "r");
    b_file = fopen(bfile.c_str(), "r");
    l_file = fopen(lfile.c_str(), "r");
    if (c_file == NULL || g_file == NULL || b_file == NULL || l_file == NULL)
    {
      UserFatal(*this) << "Cannot open one of the ROM files: " << cfile << "," << gfile << "," << bfile << "," << lfile;
    }

    // Get the input-output numbers from the projected B matrix
    Xyce::IO::MMIO::mm_read_banner( b_file, &mat_code );
    Xyce::IO::MMIO::mm_read_mtx_array_size( b_file, &M, &N );

    // Set number of internal variables to dimension of reduced-order model
    if(Instance::usePortDesc)
    {
      numIntVars= 0;
      numROMVars = 0;
      numStateVars = M+N;
    }
    else
    {
      numIntVars   = M+N;  // Number of ROM variables + current variables for ports
      numROMVars   = M;
      numStateVars = 0;
    }

    numExtVars   = N;

    // Read in B and L matrices.
    // NOTE: B and L are assumed to be dense arrays at this time, so read them in.
    Bhat.resize( M*N );
    Lhat.resize( M*N );

    // NOTE: The banners and array sizes need to be read in
    //       before the rest of the file can be read.
    int tmpM=0, tmpN=0;
    Xyce::IO::MMIO::mm_read_banner( l_file, &mat_code );   // Already read b_file banner
    Xyce::IO::MMIO::mm_read_mtx_array_size( l_file, &tmpM, &tmpN );  // Already read b_file array size

    // Read in Bhat and Lhat multivectors
    for (int i=0; i<M*N; i++)
    {
      fscanf(b_file, "%lg\n", &Bhat[i]);
      fscanf(l_file, "%lg\n", &Lhat[i]);
    }

    // Read in C matrix.
    // NOTE:  This matrix may have been sparsified or stored in a symmetric format
    Xyce::IO::MMIO::mm_read_banner( c_file, &mat_code );
    isCSparse = mm_is_sparse(mat_code);

    // If the input matrix C is dense, read in the dense data
    if (!isCSparse)
    {
      // Read in array size.
      Xyce::IO::MMIO::mm_read_mtx_array_size( c_file, &tmpM, &tmpN );
      // TODO:  Add check of tmpM and tmpN

      Chat.resize( M*M );
      if (mm_is_general(mat_code))
      {
        // Read in Chat matrix, stored column-wise
        for (int i=0; i<M*M; i++)
        {
          fscanf(c_file, "%lg\n", &Chat[i]);
        }
      }
      else if (mm_is_symmetric(mat_code) || mm_is_skew(mat_code))
      {
        int arraySize = M*(M+1)/2;  // Only one triangle is stored
        std::vector<double> Chat_tmp( arraySize );  // Only one triangle is stored
        // Read in Chat matrix, stored column-wise
        for (int i=0; i<arraySize; i++)
        {
          fscanf(c_file, "%lg\n", &Chat_tmp[i]);
        }
        // Copy Chat matrix into full dense matrix
        for (int j=0; j<M; j++)
        {
          for (int i=j; i<M; i++)
          {
            double val = Chat_tmp[j*M - j*(j-1)/2 + i - j];
            Chat[j*M+i] = val;
            if (i!=j)
            {
              if (mm_is_symmetric(mat_code))
                Chat[i*M+j] = val;  // Symmetric
              else
                Chat[i*M+j] = -val; // Skew-symmetric
            }
          }
        }
      }
      else
      {
        UserFatal(*this) << "Do not recognize the Matrix Market format for Chat (matrix is not general or symmetric)";
      }
    }
    else
    {
      int nnz=0;
      Xyce::IO::MMIO::mm_read_mtx_crd_size( c_file, &tmpM, &tmpN, &nnz );
      if (nnz==0)
      {
        UserFatal(*this) << "Chat has zero entries according to the Matrix Market file " << cfile;
      }

      // Temporary storage for Chat, to read in coordinate format
      std::vector<double> Chat_tmp( nnz );
      std::vector<int> rowIdx( nnz ), colIdx( nnz );

      if (nnz > 0)
        Xyce::IO::MMIO::mm_read_mtx_crd_data( c_file, tmpM, tmpN, nnz, &rowIdx[0], &colIdx[0], &Chat_tmp[0], mat_code );

      // Reduce the row and column indices to 0-based indexing and add entries if the matrix is (skew) symmetric.
      for (int i=0; i<nnz; ++i)
      {
        rowIdx[i]--;
        colIdx[i]--;

        // If matrices are symmetric or skew symmetric, add entries to the coordinate storage format.
        if (mm_is_symmetric(mat_code) || mm_is_skew(mat_code))
        {
          if (rowIdx[i]!=colIdx[i])
          {
            if (mm_is_symmetric(mat_code))
            {
              rowIdx.push_back(colIdx[i]);
              colIdx.push_back(rowIdx[i]);
              Chat_tmp.push_back(Chat_tmp[i]);  // Symmetric
            }
            else
            {
              rowIdx.push_back(colIdx[i]);
              colIdx.push_back(rowIdx[i]);
              Chat_tmp.push_back(-Chat_tmp[i]);  // Skew-symmetric
            }
          }
        }
      }

      // Use nnz2, the number of nonzeros for the symmetrized matrix, to compute the new CSR format.
      int nnz2 = rowIdx.size();

      // Allocate storage for Chat in CSR format
      Chat.resize( nnz2 );
      Chat_rowPtr.resize( tmpM+1 );
      Chat_colIdx.resize( nnz2 );

      // Convert coordinate (Matrix Market) storage to CSR format
      if (nnz2 > 0)
        Trilinos_Util_coocsr( tmpM, nnz2, &Chat_tmp[0], &rowIdx[0], &colIdx[0], &Chat[0], &Chat_colIdx[0], &Chat_rowPtr[0] );
    }

    // Read in G matrix.
    // NOTE:  This matrix may have been sparsified or stored in a symmetric format
    Xyce::IO::MMIO::mm_read_banner( g_file, &mat_code );
    isGSparse = mm_is_sparse(mat_code);

    // If the input matrix G is dense, read in the dense data
    if (!isGSparse)
    {
      // Read in array size.
      Xyce::IO::MMIO::mm_read_mtx_array_size( g_file, &tmpM, &tmpN );
      // TODO:  Add check of tmpM and tmpN

      Ghat.resize( M*M );
      if (mm_is_general(mat_code))
      {
        // Read in Ghat matrix, stored column-wise
        for (int i=0; i<M*M; i++)
        {
          fscanf(g_file, "%lg\n", &Ghat[i]);
        }
      }
      else if (mm_is_symmetric(mat_code) || mm_is_skew(mat_code))
      {
        int arraySize = M*(M+1)/2;
        std::vector<double> Ghat_tmp( arraySize );  // Only one triangle is stored
        // Read in Ghat matrix, stored column-wise
        for (int i=0; i<arraySize; i++)
        {
          fscanf(g_file, "%lg\n", &Ghat_tmp[i]);
        }
        // Copy Ghat matrix into full dense matrix
        for (int j=0; j<M; j++)
        {
          for (int i=j; i<M; i++)
          {
            double val = Ghat_tmp[j*M - j*(j-1)/2 + i - j];
            Ghat[j*M+i] = val;
            if (i!=j)
            {
              if (mm_is_symmetric(mat_code))
                Ghat[i*M+j] = val;  // Symmetric
              else
                Ghat[i*M+j] = -val; // Skew-symmetric
            }
          }
        }
      }
      else
      {
        UserFatal(*this) << "Do not recognize the Matrix Market format for Ghat (matrix is not general or symmetric)";
      }
    }
    else
    {
      int nnz=0;
      Xyce::IO::MMIO::mm_read_mtx_crd_size( g_file, &tmpM, &tmpN, &nnz );
      if (nnz==0)
      {
        UserFatal(*this) << "Ghat has zero entries according to the Matrix Market file " << gfile;
      }

      // Temporary storage for Ghat, to read in coordinate format
      std::vector<double> Ghat_tmp( nnz );
      std::vector<int> rowIdx( nnz ), colIdx( nnz );

      if (nnz > 0)
        Xyce::IO::MMIO::mm_read_mtx_crd_data( g_file, tmpM, tmpN, nnz, &rowIdx[0], &colIdx[0], &Ghat_tmp[0], mat_code );

      // Reduce the row and column indices to 0-based indexing
      for (int i=0; i<nnz; ++i)
      {
        rowIdx[i]--;
        colIdx[i]--;

        // If matrices are symmetric or skew symmetric, add entries to the coordinate storage format.
        if (mm_is_symmetric(mat_code) || mm_is_skew(mat_code))
        {
          if (rowIdx[i]!=colIdx[i])
          {
            if (mm_is_symmetric(mat_code))
            {
              rowIdx.push_back(colIdx[i]);
              colIdx.push_back(rowIdx[i]);
              Ghat_tmp.push_back(Ghat_tmp[i]);  // Symmetric
            }
            else
            {
              rowIdx.push_back(colIdx[i]);
              colIdx.push_back(rowIdx[i]);
              Ghat_tmp.push_back(-Ghat_tmp[i]);  // Skew-symmetric
            }
          }
        }
      }

      // Use nnz2, the number of nonzeros for the symmetrized matrix, to compute the new CSR format.
      int nnz2 = rowIdx.size();

      // Allocate storage for Ghat in CSR format
      Ghat.resize( nnz2 );
      Ghat_rowPtr.resize( tmpM+1 );
      Ghat_colIdx.resize( nnz2 );

      // Convert coordinate (Matrix Market) storage to CSR format
      if (nnz2 > 0)
        Trilinos_Util_coocsr( tmpM, nnz2, &Ghat_tmp[0], &rowIdx[0], &colIdx[0], &Ghat[0], &Ghat_colIdx[0], &Ghat_rowPtr[0] );
    }

    // Create a union map for Chat and Ghat if both are sparse, otherwise the union is a dense block.
    if (isCSparse && isGSparse)
    {
      CG_rowPtr.resize( M+1 );
      CG_colIdx.resize( Chat_colIdx.size() + Ghat_colIdx.size() );  // Maximum number of nonzeros in the union

      std::vector<int>::iterator it;
      std::vector<int>::iterator itCG = CG_colIdx.begin();
      std::vector<int>::iterator itChat = Chat_colIdx.begin();
      std::vector<int>::iterator itGhat = Ghat_colIdx.begin();
      CG_rowPtr[0] = 0;
      for (int i=0; i<M; ++i)
      {
        // Get the number of nonzero entries for this row
        int numEntriesChat = Chat_rowPtr[i+1] - Chat_rowPtr[i];
        int numEntriesGhat = Ghat_rowPtr[i+1] - Ghat_rowPtr[i];
        // Compute the union of the entries
        it = set_union( itChat, itChat + numEntriesChat, itGhat, itGhat + numEntriesGhat, itCG );
        CG_rowPtr[i+1] = CG_rowPtr[i] + (int)(it - itCG);
        // Check if we need to sort the entries HERE!!!!

        // Increment the iterators
        itCG = it;
        itChat += numEntriesChat;
        itGhat += numEntriesGhat;
      }
    }

    // Allocate memory for projected matrices
    Qhat.resize( M );
    Fhat.resize( M+N );
    i_ip.resize( N );

    if(Instance::usePortDesc>0)
    {
      Jstamp.resize(N*N);
      Fstamp.resize(N);
      G2.resize((M+N)*(M+N));
      C2.resize((M+N)*(M+N));
      A2.resize((M+N)*(M+N));
      A2last.resize((M+N)*(M+N));
      G2p.resize((M+N)*N);
      Gp2.resize((M+N)*N);
      A2sol.resize((M+N)*N);

      ipiv_A2.resize(M+N);
    }

    fclose(c_file);
    fclose(g_file);
    fclose(b_file);
    fclose(l_file);

    if(Instance::usePortDesc>0)
    {
      if (!isGSparse)
      {
        // Construct constant matrices for two-level stamp
        // G2 = [eye(N),-Lhat'; zeros(M,N), Ghat];
        for(int iy=0; iy<M; iy++) // add Ghat term
        {
          for(int ix=0; ix<M; ix++)
            G2[((N+M+1)*N)+ix+(M+N)*iy] = Ghat[ix+iy*M];
        }
      }
      else
      {
        for (int ix=0; ix<M; ix++)
        {
          for (int j=Ghat_rowPtr[ix]; j<Ghat_rowPtr[ix+1]; ++j)
            G2[((N+M+1)*N)+ix+(M+N)*Ghat_colIdx[j]] = Ghat[j];
        }
      }

      for(int ix=0; ix<N; ix++)  // add eye(N) term
        G2[ix+(M+N)*ix] = 1;
      for(int ix=0; ix<M; ix++) // add Lhat' term
      {
        for(int iy=0; iy<N; iy++)
          G2[((M+N)*N)+iy+(M+N)*ix] = -Lhat[ix+iy*M];  // L is transposed
      }
      // C2 = [zeros(N),zeros(N,M); zeros(M,N), Chat];
      if (!isCSparse)
      {
        for(int iy=0; iy<M; iy++) // add Chat term
        {
          for(int ix=0; ix<M; ix++)
            C2[((N+M+1)*N)+ix+(M+N)*iy] = Chat[ix+iy*M];
        }
      }
      else
      {
        for (int ix=0; ix<M; ix++)
        {
          for (int j=Chat_rowPtr[ix]; j<Chat_rowPtr[ix+1]; ++j)
            C2[((N+M+1)*N)+ix+(M+N)*Chat_colIdx[j]] = Chat[j];
        }
      }

      // G2p = -[zeros(N); Bhat];
      for(int iy=0; iy<N; iy++) // add Bhat term
      {
        for(int ix=0; ix<M; ix++)
          G2p[N+ix+(M+N)*iy] = -Bhat[ix+iy*M];
      }

      // Gp2 = [eye(N), zeros(M,N)];
      for(int iy=0; iy<N; iy++)  // add eye(N) term
      {
        Gp2[iy+(N)*iy] = 1;
      }

      // Construct jacStamp
      if( jacStamp.empty() )
      {
        jacStamp.resize(numExtVars);

        // Put in external variables (node voltages)
        for (int i=0; i<numExtVars; i++)
        {
          jacStamp[i].resize(numExtVars);
          for(int j=0; j<numExtVars; j++)
            jacStamp[i][j] = j;
        }
      }
    }
    else
    {
      // Create Jacobian stamp for direct stamping of ROM into full system
      //
      // [ Stamps for     0          0          ] [ x_NL ]   [ v_NL ]
      // [ f(x_NL,u_p)   I_N         0          ] [ u_p  ] = [ v_p  ]   (a)
      // [  0      0     I_N      -Lhat^T       ] [ i_p  ]   [  0   ]   (b)
      // [  0    -Bhat    0  (Ghat + Chat d/dt) ] [xhat_q]   [  0   ]   (c)
      //
      if( jacStamp.empty() )
      {
        // Resize Jacobian stamp to size of reduced system + 2 x #ports
        jacStamp.resize(numIntVars+numExtVars);

        // Equations (a): Put in external variables (port voltages)
        for (int i=0; i<numExtVars; i++)
        {
          jacStamp[i].resize(1);
          jacStamp[i][0] = numExtVars+i;
        }

        // Equations (b): Put in internal variables (size of reduced model + port currents)
        for (int i=numExtVars; i<2*numExtVars; i++)
        {
          jacStamp[i].resize(numROMVars+1);
          jacStamp[i][0] = i;
          for (int j=0; j<numROMVars; j++)
            jacStamp[i][j+1] = 2*numExtVars+j;
        }

        // Equations (c): Put in projected system
        for (int i=2*numExtVars; i<numIntVars+numExtVars; i++)
        {
          int numEntries = numIntVars;
          if (isCSparse && isGSparse)
            numEntries = numExtVars+(CG_rowPtr[i-2*numExtVars+1]-CG_rowPtr[i-2*numExtVars]);
          jacStamp[i].resize( numEntries );

          for (int j=0; j<numExtVars; j++)
            jacStamp[i][j] = j;

          // Insert entries for (Ghat + Chat d/dt)
          if (isCSparse && isGSparse)
          {
            for (int j=numExtVars; j<numEntries; j++)
              jacStamp[i][j] = 2*numExtVars + CG_colIdx[CG_rowPtr[i-2*numExtVars] + j - numExtVars];
          }
          else
          {
            for (int j=numExtVars; j<numEntries; j++)
              jacStamp[i][j] = numExtVars + j;
          }
        }
      }
    }
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params:

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.  For the matrix indices, first do the
  // rows.

  // Obtain indices for internal variables
  li_ROM.resize(numROMVars);
  if(usePortDesc==0)
  {
    for (int i=0; i<numROMVars; i++)
      li_ROM[i] = intLIDVec[i+numExtVars];
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << section_divider << std::endl;

    Xyce::dout() << "::registerLIDs:\n";
    Xyce::dout() << "  name = " << getName() << std::endl;

    Xyce::dout() << "\nsolution indices:\n";
    for (int i=0; i<numExtVars; ++i)
      Xyce::dout() << "  li_up[" << i << "] = " << extLIDVec[i] << std::endl;
    if (usePortDesc==0)
    {
      for (int i=0; i<numExtVars; ++i)
        Xyce::dout() << "  li_ip[" << i << "] = " << intLIDVec[i] << std::endl;
    }
    for (int i=0; i<numROMVars; i++)
      Xyce::dout() << "  li_ROM[" << i << "] = " << li_ROM[i] << std::endl;

    Xyce::dout() << section_divider << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  li_state.resize(numStateVars, 0);
  for (int ix=0; ix<numStateVars; ix++) {
    li_state[ix] = staLIDVec[ix];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  if (numIntVars > 0)  // Then add the current variables at the ports
    for (int i = 0; i < numExtVars; i++)
      addInternalNode(symbol_table, intLIDVec[i], getName(), "ip_Node" + Teuchos::Utils::toString(i + 1));

  for (int i = 0; i < numROMVars; i++)
    addInternalNode(symbol_table, li_ROM[i], getName(), "ROM_Node"+Teuchos::Utils::toString(i+1));
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{

  DeviceInstance::registerJacLIDs( jacLIDVec );

  if(usePortDesc>0)
  {
    AEqu_NodeOffset.resize(numExtVars);
    for (int i=0; i<numExtVars; i++)
    {
      AEqu_NodeOffset[i].resize(numExtVars);
      for (int j=0; j<numExtVars; j++)
      {
        AEqu_NodeOffset[i][j] = jacLIDVec[i][j];
      }
    }
  }
  else
  {
    AEqu_up_NodeOffset.resize(numExtVars);
    AEqu_ip_NodeOffset.resize(numExtVars);
    for (int i=0; i<numExtVars; i++)
    {
      AEqu_up_NodeOffset[i] = jacLIDVec[i][0];
      AEqu_ip_NodeOffset[i] = jacLIDVec[numExtVars+i][0];
    }

    // Get the offsets for -L^T
    ROMEqu_Lt_NodeOffset.resize(numROMVars);
    for (int i=0; i<numROMVars; i++)
    {
      ROMEqu_Lt_NodeOffset[i] = jacLIDVec[numExtVars][i+1];
    }

    // Get the offsets for -B
    ROMEqu_B_NodeOffset.resize(numExtVars*numROMVars);
    for (int i=0; i<numROMVars; i++)
    {
      for (int j=0; j<numExtVars; j++)
        ROMEqu_B_NodeOffset[i*numExtVars+j] = jacLIDVec[2*numExtVars+i][j];
    }

    // Get the offsets for C d/dt
    // If the matrix is sparse, then an offset will be stored for each nonzero.
    if (isCSparse)
    {
      ROMEqu_C_NodeOffset.resize(Chat_rowPtr[numROMVars]);
      for (int i=0; i<numROMVars; i++)
      {
        int nnz = CG_rowPtr[i+1]-CG_rowPtr[i];
        int CrowPtr = Chat_rowPtr[i], Cnnz = Chat_rowPtr[i+1]-Chat_rowPtr[i];
        for (int j=0, Cidx=0; j<nnz && Cidx<Cnnz; j++)
        {
          int colIdx = CG_colIdx[CG_rowPtr[i]+j];
          if (colIdx == Chat_colIdx[CrowPtr])
          {
            ROMEqu_C_NodeOffset[CrowPtr++] = jacLIDVec[2*numExtVars+i][numExtVars+j];
            Cidx++;
          }
        }
      }
    }

    // Get the offsets for G
    // If the matrix is sparse, then an offset will be stored for each nonzero.
    if (isGSparse)
    {
      ROMEqu_G_NodeOffset.resize(Ghat_rowPtr[numROMVars]);
      for (int i=0; i<numROMVars; i++)
      {
        int nnz = CG_rowPtr[i+1]-CG_rowPtr[i];
        int GrowPtr = Ghat_rowPtr[i], Gnnz = Ghat_rowPtr[i+1]-Ghat_rowPtr[i];
        for (int j=0, Gidx=0; j<nnz && Gidx<Gnnz; j++)
        {
          int colIdx = CG_colIdx[CG_rowPtr[i]+j];
          if (colIdx == Ghat_colIdx[GrowPtr])
          {
            ROMEqu_G_NodeOffset[GrowPtr++] = jacLIDVec[2*numExtVars+i][numExtVars+j];
            Gidx++;
          }
        }
      }
    }

    // If either C or G is dense, we need these node offsets, which assume the storage is contiguous.
    if (!isCSparse || !isGSparse)
    {
      ROMEqu_GpC_NodeOffset.resize(numROMVars);
      for (int i=0; i<numROMVars; i++)
      {
        ROMEqu_GpC_NodeOffset[i] = jacLIDVec[2*numExtVars][numExtVars+i];
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << Xyce::section_divider << std::endl;
      Xyce::dout() << "Instance::registerJacLIDs\n";
      if (usePortDesc==0)
      {
        Xyce::dout() << " AEqu_up_NodeOffset: ";
        for (int i=0; i<numExtVars; i++)
          Xyce::dout() << AEqu_up_NodeOffset[i] << "  ";
        Xyce::dout() << std::endl;
        Xyce::dout() << " AEqu_ip_NodeOffset: ";
        for (int i=0; i<numExtVars; i++)
          Xyce::dout() << AEqu_ip_NodeOffset[i] << "  ";
        Xyce::dout() << std::endl;
        Xyce::dout() << " AROMEqu_Lt_NodeOffset: ";
        for (int i=0; i<numROMVars; i++)
          Xyce::dout() << ROMEqu_Lt_NodeOffset[i] << "  ";
        Xyce::dout() << std::endl;
        Xyce::dout() << " AROMEqu_B_NodeOffset: " << std::endl;
        for (int i=0; i<numROMVars; i++)
        {
          for (int j=0; j<numExtVars; j++)
            Xyce::dout() << ROMEqu_B_NodeOffset[i*numExtVars+j] << "  ";
          Xyce::dout() << std::endl;
        }
        Xyce::dout() << " AROMEqu_GpC_NodeOffset: ";
        for (int i=0; i<numROMVars; i++)
          Xyce::dout() << ROMEqu_GpC_NodeOffset[i] << "  ";
        Xyce::dout() << std::endl;
      }
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // Get pointers for KCL equations and Ghat, Lhat, and Bhat
  if(usePortDesc==0)
  {
    fEqu_up_NodePtr.resize(numExtVars);
    for (int i=0; i<numExtVars; ++i)
    {
      fEqu_up_NodePtr[i] = &(dFdx[extLIDVec[i]][AEqu_up_NodeOffset[i]]);
    }

    fEqu_ip_NodePtr.resize(numExtVars);
    for (int i=0; i<numExtVars; ++i)
    {
      fEqu_ip_NodePtr[i] = &(dFdx[intLIDVec[i]][AEqu_ip_NodeOffset[i]]);
    }

    // Get pointers for Ghat (assuming contiguous), only if Ghat is not sparse.
    if (isGSparse)
    {
      fROMEqu_Ghat_VarsPtrs.resize(Ghat_rowPtr[numROMVars]);
      for (int i=0; i<numROMVars; ++i)
      {
        for (int j=Ghat_rowPtr[i]; j<Ghat_rowPtr[i+1]; j++)
        {
          fROMEqu_Ghat_VarsPtrs[j]=&(dFdx[li_ROM[i]][ROMEqu_G_NodeOffset[j]]);
        }
      }
    }
    else
    {
      fROMEqu_Ghat_VarsPtrs.resize(numROMVars);
      for (int i=0; i<numROMVars; ++i)
      {
        fROMEqu_Ghat_VarsPtrs[i]=&(dFdx[li_ROM[i]][ROMEqu_GpC_NodeOffset[0]]);
      }
    }

    // Lhat (assuming contiguous)
    fROMEqu_Lhat_VarsPtrs.resize(numExtVars);
    for (int i=0; i<numExtVars; ++i)
    {
      fROMEqu_Lhat_VarsPtrs[i]=&(dFdx[intLIDVec[i]][ROMEqu_Lt_NodeOffset[0]]);
    }

    // Bhat (these are not guaranteed to be contiguous)
    fROMEqu_Bhat_VarsPtrs.resize(numExtVars*numROMVars);
    for (int i=0; i<numROMVars; ++i)
    {
      for (int j=0; j<numExtVars; ++j)
      {
        fROMEqu_Bhat_VarsPtrs[numExtVars*i+j]=&(dFdx[li_ROM[i]][ROMEqu_B_NodeOffset[i*numExtVars+j]]);
      }
    }

    // Get pointers for Chat (assuming contiguous), only if Chat is not sparse.
    if (isCSparse)
    {
      qROMEqu_Chat_VarsPtrs.resize(Chat_rowPtr[numROMVars]);
      for (int i=0; i<numROMVars; ++i)
      {
        for (int j=Chat_rowPtr[i]; j<Chat_rowPtr[i+1]; j++)
        {
          qROMEqu_Chat_VarsPtrs[j]=&(dQdx[li_ROM[i]][ROMEqu_C_NodeOffset[j]]);
        }
      }
    }
    else
    {
      qROMEqu_Chat_VarsPtrs.resize(numROMVars);
      for (int i=0; i<numROMVars; ++i)
        qROMEqu_Chat_VarsPtrs[i]=&(dQdx[li_ROM[i]][ROMEqu_GpC_NodeOffset[0]]);
    }
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadErrorWeightMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Keith Santarelli, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/12/08
//-----------------------------------------------------------------------------
void Instance::loadErrorWeightMask()
{
  if (maskROMVars)
  {
    Linear::Vector * maskVectorPtr = extData.deviceErrorWeightMask_;

    for (int i=0; i<numROMVars; i++)
    {
      (*maskVectorPtr)[li_ROM[i]] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 ROM instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;
  double * xhat = &solVec[li_ROM[0]];

  Teuchos::BLAS<int, double> blas;

  // Compute Qhat = Chat * xhat
  if (isCSparse)
    Linear::crsAxpy( numROMVars, 1.0, &Chat[0], &Chat_rowPtr[0], &Chat_colIdx[0], xhat, 0.0, &Qhat[0] );
  else
    blas.GEMV( Teuchos::NO_TRANS, numROMVars, numROMVars, 1.0, &Chat[0], numROMVars, xhat, 1, 0.0, &Qhat[0], 1 );

  for (int i=0; i<numROMVars; i++)
  {
    qVec[li_ROM[i]] += Qhat[i];
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 ROM instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  std::vector<double> v_up(numExtVars);
  for (int i=0; i<numExtVars; ++i)
  {
    v_up[i] = solVec[extLIDVec[i]];
    Fhat[i] = solVec[intLIDVec[i]];
    i_ip[i] = solVec[intLIDVec[i]];
  }
  double * xhat = &solVec[li_ROM[0]];

  Teuchos::BLAS<int, double> blas;

  // Compute Fhat[0:numExtVars-1] = i_ip - Lhat'*xhat
  blas.GEMV( Teuchos::TRANS, numROMVars, numExtVars, -1.0, &Lhat[0], numROMVars, xhat, 1, 1.0, &Fhat[0], 1 );

  // Compute Fhat[numExtVars:numIntVars] = Ghat*xhat - Bhat*v_up
  if (isGSparse)
    Linear::crsAxpy( numROMVars, 1.0, &Ghat[0], &Ghat_rowPtr[0], &Ghat_colIdx[0], xhat, 0.0, &Fhat[numExtVars] );
  else
    blas.GEMV( Teuchos::NO_TRANS, numROMVars, numROMVars, 1.0, &Ghat[0], numROMVars, xhat, 1, 0.0, &Fhat[numExtVars], 1 );
  blas.GEMV( Teuchos::NO_TRANS, numROMVars, numExtVars, -1.0, &Bhat[0], numROMVars, &v_up[0], 1, 1.0, &Fhat[numExtVars], 1 );

  // Load F vector
  for (int i=0; i<numExtVars; ++i)
  {
    fVec[extLIDVec[i]] += i_ip[i];
    fVec[intLIDVec[i]] += Fhat[i];
  }
  for (int i=0; i<numROMVars; i++)
  {
    fVec[li_ROM[i]] += Fhat[numExtVars+i];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 ROM instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // Load Chat portion of dQ/dx
  for (int i=0; i<numROMVars; i++)
  {
    if (isCSparse)
    {
      for (int j=Chat_rowPtr[i]; j<Chat_rowPtr[i+1]; j++)
      {

        dQdx[li_ROM[i]][ROMEqu_C_NodeOffset[j]] += Chat[j];
      }
    }
    else
    {
      for (int j=0; j<numROMVars; j++)
      {

        dQdx[li_ROM[i]][ROMEqu_GpC_NodeOffset[j]] += Chat[j*numROMVars + i];
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 ROM instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  for (int i=0; i<numExtVars; ++i)
  {

    dFdx[extLIDVec[i]][AEqu_up_NodeOffset[i]] += 1.0;

    dFdx[intLIDVec[i]][AEqu_ip_NodeOffset[i]] += 1.0;
  }

  // Load -Lhat portion of dF/dx
  for (int j=0; j<numROMVars; j++)
  {
    for (int i=0; i<numExtVars; i++)
    {

      dFdx[intLIDVec[i]][ROMEqu_Lt_NodeOffset[j]] -= Lhat[j];
    }
  }

  // Load -Bhat portion of dF/dx
  for (int i=0; i<numROMVars; i++)
  {
    for (int j=0; j<numExtVars; j++)
    {

      dFdx[li_ROM[i]][ROMEqu_B_NodeOffset[i*numExtVars+j]] -= Bhat[j*numROMVars + i];
    }
  }

  // Load Ghat portion of dF/dx
  for (int i=0; i<numROMVars; i++)
  {
    if (isGSparse)
    {
      for (int j=Ghat_rowPtr[i]; j<Ghat_rowPtr[i+1]; j++)
      {

        dFdx[li_ROM[i]][ROMEqu_G_NodeOffset[j]] += Ghat[j];
      }
    }
    else
    {
      for (int j=0; j<numROMVars; j++)
      {

        dFdx[li_ROM[i]][ROMEqu_GpC_NodeOffset[j]] += Ghat[j*numROMVars + i];
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
}


// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 5/17/00
//-----------------------------------------------------------------------------

Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

Model::~Model()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }

}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------

std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of ROM instances: " << isize << std::endl;
  os << "    name\t\tmodelName\tParameters" << std::endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


// ROM Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec, int loadType)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " ----------------------------------" << std::endl;
    Xyce::dout() << " Master::updateState: " << std::endl;
  }

  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }
 
  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);

    Teuchos::BLAS<int, double> blas;
    if (ci.usePortDesc>0) // TWO-LEVEL DEVICE STAMPS
    {
      const char test = 'N';
      int MN = ci.numStateVars;
      int N = ci.numExtVars;
      int M = MN-N;

      // Check if this is the first device call after a successful timestep
      int sameTimeStep=0;
      if (ci.getSolverState().newtonIter==1) { ci.lastTimeStepNumber = ci.getSolverState().timeStepNumber_; }
      if (ci.lastTimeStepNumber==ci.getSolverState().timeStepNumber_) {sameTimeStep=1;}

      //********************************************************************
      // Set up the coefficients for discretization for the INNER solve
      //********************************************************************
      if (ci.usePortDesc == 1)  // DYNAMIC ORDER SELECTION
      {
        if (sameTimeStep == 0)
          ci.coefLast = ci.coef;

        if (ci.getSolverState().currentOrder_ == 0)
          ci.coef = 1.0;     // Force BE
        else
          ci.coef = (1.0/ci.getSolverState().currentOrder_);

        if (ci.getSolverState().usedOrder_ == 0)
          ci.coefLast = 1.0; // Force BE
      }
      else if (ci.usePortDesc == 2) { // BACKWARD EULER
        ci.coef = 1.0;
        ci.coefLast = 1.0;
      }
      else if (ci.usePortDesc == 3) { // TRAP
        ci.coef=0.5;
        ci.coefLast = 0.5;
      }
      else
        Xyce::dout() << "Bad 'USE_PORT_DESCRIPTION' flag" << std::endl;


      // Get alph=1/dt for current and last time steps
      if (sameTimeStep==0) {ci.alph_last = ci.alph;}
      ci.alph=1/ci.getSolverState().currTimeStep_;
      if (ci.getSolverState().dcopFlag==1) { ci.alph=0.0; ci.alph_last=0.0; }
      if (ci.getSolverState().timeStepNumber_ == 0) { ci.alph_last=0.0; }


      //********************************************************************
      // Load port and internal state variables
      //********************************************************************
      std::vector<double> lastStaVec, currStaVec, nextStaVec; // Internal states
      lastStaVec.resize(M+N,0);
      currStaVec.resize(M+N,0);
      nextStaVec.resize(M+N,0);
      //for(int ix=0; ix<(M+N); ++ix) { lastStaVec[ix] = (*ci.extData.lastStaVectorPtr)[ix]; }  // WRONG
      for(int ix=0; ix<(M+N); ++ix) { lastStaVec[ix] = (*ci.extData.lastStaVectorPtr)[ ci.li_state[ix] ]; }
      std::vector<double> lastPortVec, currPortVec, nextPortVec; // Port voltages
      lastPortVec.resize(N,0);
      currPortVec.resize(N,0);
      nextPortVec.resize(N,0);
      for(int ix=0; ix<N; ++ix)
      {
        lastPortVec[ix] = (*ci.extData.lastSolVectorPtr)[ci.extLIDVec[ix]];
        currPortVec[ix] = (*ci.extData.currSolVectorPtr)[ci.extLIDVec[ix]];
        nextPortVec[ix] = solVec[ci.extLIDVec[ix]];
      }

      /*
        std::vector<int> lidState = set this up in registerStateLIDs
        localState[ix] = ci.extradata.currState[ ci.lidState[ix] ];
        ci.extData.currState[ ci.lidState[ix] ] = localState[ix];
      */

      // Begin math stuff
      Teuchos::LAPACK<int,double> lapack;
      std::vector<int> ipiv_A2last (M+N, 0);
      int info_A2, info_A2last, info2_A2, info2_A2last;


      //****************************************************************
      // Compute, or load, old internal state (currStaVec) [INNER SOLVE]
      // x2[t] = A2\(-coef*G2p*xp[t] - (1-coef)*G2p*xp[t-1] - [(1-coef)*G2-C2/dt]*x2[t-1])
      //******************************************************
      int updateCurrStaVec=0;
      //      if (ci.getSolverState().newtonIter==0) { updateCurrStaVec=1; }
      if (sameTimeStep==0) {updateCurrStaVec=1; }
      if (updateCurrStaVec==1) {  // compute new `last' state
        blas.GEMV( Teuchos::NO_TRANS,N+M,N+M,-(1-ci.coefLast),&ci.G2[0],M+N,&lastStaVec[0],1,0.0,&currStaVec[0],1 ); // -(1-coefLast)*G2*lastStaVec
        blas.GEMV( Teuchos::NO_TRANS,N+M,N+M,ci.alph_last,&ci.C2[0],M+N,&lastStaVec[0],1,1.0,&currStaVec[0],1 ); // C2/dtlast*lastStaVec
        blas.GEMV( Teuchos::NO_TRANS,N+M,N,-ci.coefLast,&ci.G2p[0],N+M,&currPortVec[0],1, 1.0,&currStaVec[0],1 ); // -coefLast*G2p*currPortVec
        blas.GEMV( Teuchos::NO_TRANS,N+M,N,-(1-ci.coefLast),&ci.G2p[0],N+M,&lastPortVec[0],1,1.0,&currStaVec[0],1 ); // -(1-coefLast)*G2p*lastPortVec
        int useOldA2=1;
        if (ci.getSolverState().timeStepNumber_ < 2) { useOldA2=0; }
        if (ci.lastTimeStepNumber==ci.getSolverState().timeStepNumber_) { useOldA2=0; }
        if (useOldA2==1)  // use [previously factored] A2 as A2last
          lapack.GETRS(test, M+N, 1, &ci.A2[0], M+N, &ci.ipiv_A2[0], &currStaVec[0], M+N, &info2_A2); // A2\(xtmp)
        else {  // compute new A2last
          for(int ix=0; ix<(M+N)*(M+N); ix++) {  ci.A2last[ix]= (ci.alph_last*ci.C2[ix]) + (ci.coefLast*ci.G2[ix]); } // A2 using previous dt
          lapack.GETRF( M+N, M+N, &ci.A2last[0], M+N, &ipiv_A2last[0], &info_A2last);  // factor A2last
          lapack.GETRS(test, M+N, 1, &ci.A2last[0], M+N, &ipiv_A2last[0], &currStaVec[0], M+N, &info2_A2last); // A2last\(xtmp)
        }
        for(int ix=0; ix<(M+N); ++ix) { (*ci.extData.currStaVectorPtr)[ ci.li_state[ix] ] = currStaVec[ix]; } 	// Store answer
      }
      else
      {  // Load currStaVec
        for(int ix=0; ix<(M+N); ++ix) { currStaVec[ix] = (*ci.extData.currStaVectorPtr)[ ci.li_state[ix] ]; }
      }


      //**************************************************
      // construct A2 and Jacobian (Jstamp), if necessary
      // Jstamp = -Gp2*(A2\(coef*G2p)) = -Gp2 * A2sol
      //**************************************************
      int updateA2 = 1;
      if (ci.getSolverState().timeStepNumber_ > 2){
        if (ci.getSolverState().newtonIter>0) { updateA2=0; }
        if (ci.alph==ci.alph_last) { if (ci.coef==ci.coefLast) { updateA2=0; }}}
      if (updateA2==1) {
        for(int ix=0; ix<(M+N)*(M+N); ix++) {	  ci.A2[ix]= (ci.alph*ci.C2[ix]) + (ci.coef*ci.G2[ix]); } // A2 = 1/dt*C2 + coef*G2;
        lapack.GETRF( M+N, M+N, &ci.A2[0], M+N, &ci.ipiv_A2[0], &info_A2);  // Factor A2
        for(int ix=0; ix<(M+N)*N; ix++) { ci.A2sol[ix] = ci.coef * ci.G2p[ix]; }  // copy to RHS vector
        lapack.GETRS(test, M+N, N, &ci.A2[0], M+N, &ci.ipiv_A2[0], &ci.A2sol[0], M+N, &info2_A2);
        blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, N, N, M+N, -1.0,&ci.Gp2[0],N,&ci.A2sol[0],M+N,0.0,&ci.Jstamp[0], N);
      }


      //************************************************
      // Compute and save new internal state (nextStaVec)
      // x2[t] = A2\(-coef*G2p*xp[t] - (1-coef)*G2p*xp[t-1] - [(1-coef)*G2-C2/dt]*x2[t-1])
      //************************************************
      blas.GEMV( Teuchos::NO_TRANS,N+M,N+M,-(1-ci.coef),&ci.G2[0],M+N,&currStaVec[0],1,0.0,&nextStaVec[0],1 ); // -(1-coef)*G2*currStaVec
      blas.GEMV( Teuchos::NO_TRANS,N+M,N+M,ci.alph,&ci.C2[0],M+N,&currStaVec[0],1,1.0,&nextStaVec[0],1 ); // C2/dt*currStaVec
      blas.GEMV( Teuchos::NO_TRANS,N+M,N,-ci.coef,&ci.G2p[0],N+M,&nextPortVec[0],1, 1.0,&nextStaVec[0],1 ); // -coef*G2p*nextPortVec
      blas.GEMV( Teuchos::NO_TRANS,N+M,N,-(1-ci.coef),&ci.G2p[0],N+M,&currPortVec[0],1,1.0,&nextStaVec[0],1 ); // -(1-coef)*G2p*currPortVec
      lapack.GETRS(test, M+N, 1, &ci.A2[0], M+N, &ci.ipiv_A2[0], &nextStaVec[0], M+N, &info2_A2);               // A2\(prev)
      //      for(int ix=0; ix<(M+N); ++ix) { 	(*ci.extData.nextStaVectorPtr)[ix] = nextStaVec[ix]; }  // save internal state
      for(int ix=0; ix<(M+N); ++ix) { 	(*ci.extData.nextStaVectorPtr)[ ci.li_state[ix] ] = nextStaVec[ix]; }  // save internal state


      //********************************************************
      // Compute residual stamp (Fstamp) using new internal state
      // Fstamp = Gp2 * nextStaVec
      //********************************************************
      blas.GEMV( Teuchos::NO_TRANS,N,N+M,1.0,&ci.Gp2[0],N,&nextStaVec[0],1,0.0,&ci.Fstamp[0],1 ); // Gp2*nextStaVec


    }  // END 2-LEVEL STAMP

    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------

    if(ci.usePortDesc==0) // REGULAR DEVICE STAMPS
    {
      int N = ci.numExtVars;
      std::vector<double> v_up(N);
      for (int i=0; i<N; ++i)
      {
        v_up[i] = solVec[ci.extLIDVec[i]];
        ci.i_ip[i] = solVec[ci.intLIDVec[i]];
        ci.Fhat[i] = solVec[ci.intLIDVec[i]];
      }
      double * xhat = &solVec[ci.li_ROM[0]];

      // Compute Fhat[0:numExtVars-1] = i_ip - Lhat'*xhat
      blas.GEMV( Teuchos::TRANS, ci.numROMVars, ci.numExtVars, -1.0, &ci.Lhat[0], ci.numROMVars, xhat, 1, 1.0, &ci.Fhat[0], 1 );

      //    Xyce::dout() << "Fhat (should just change first two elements): " << std::endl;
      //    for (int i=0; i<ci.numIntVars; i++)
      //      Xyce::dout() << "Fhat [" << i << "] = " << ci.Fhat[i] << std::endl;

      // Compute Fhat[numExtVars:numIntVars] = Ghat*xhat - Bhat*v_up
      if (ci.isGSparse)
        Linear::crsAxpy( ci.numROMVars, 1.0, &ci.Ghat[0], &ci.Ghat_rowPtr[0], &ci.Ghat_colIdx[0], xhat, 0.0, &ci.Fhat[ci.numExtVars] );
      else
        blas.GEMV( Teuchos::NO_TRANS, ci.numROMVars, ci.numROMVars, 1.0, &ci.Ghat[0], ci.numROMVars, xhat, 1, 0.0, &ci.Fhat[ci.numExtVars], 1 );
      blas.GEMV( Teuchos::NO_TRANS, ci.numROMVars, ci.numExtVars, -1.0, &ci.Bhat[0], ci.numROMVars, &v_up[0], 1, 1.0, &ci.Fhat[ci.numExtVars], 1 );

      //    Xyce::dout() << "Fhat (should just change last elements representing ROM): " << std::endl;
      //    for (int i=0; i<ci.numIntVars; i++)
      //      Xyce::dout() << "Fhat [" << i << "] = " << ci.Fhat[i] << std::endl;

      // Compute Qhat = Chat * xhat
      if (ci.isCSparse)
        Linear::crsAxpy( ci.numROMVars, 1.0, &ci.Chat[0], &ci.Chat_rowPtr[0], &ci.Chat_colIdx[0], xhat, 0.0, &ci.Qhat[0] );
      else
        blas.GEMV( Teuchos::NO_TRANS, ci.numROMVars, ci.numROMVars, 1.0, &ci.Chat[0], ci.numROMVars, xhat, 1, 0.0, &ci.Qhat[0], 1 );

      //    Xyce::dout() << "Qhat (should change all elements): " << std::endl;
      //    for (int i=0; i<ci.numROMVars; i++)
      //      Xyce::dout() << "Qhat [" << i << "] = " << ci.Qhat[i] << std::endl;
    }
  }

  return true;
}



//-----------------------------------------------------------------------------
// Function      : Master::printMatrix
// Purpose       :
// Special Notes : For debugging 2-level device stamp
// Scope         : public
// Creator       : Brad Bond, SNL
// Creation Date : 2011-03-11
//-----------------------------------------------------------------------------
void Master::printMatrix ( std::string vname, double * Matrix, int Nrows, int Ncols )
{
  Xyce::dout() << std::endl << vname << ": " << std::endl;
  for(int ix=0; ix < Nrows; ix++)
  {
    for(int iy=0; iy < Ncols; iy++)
    {
      Xyce::dout() << Matrix[iy*Nrows+ix] << " ";
    }
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, 
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " ----------------------------------" << std::endl;
    Xyce::dout() << " Master::loadDAEVectors: " << std::endl;
  }

  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);
    if(ci.usePortDesc==0)
    {
      // Load F vector
      for (int i=0; i<ci.numExtVars; ++i)
      {
        fVec[ci.extLIDVec[i]] += ci.i_ip[i];
        fVec[ci.intLIDVec[i]] += ci.Fhat[i];
      }

      // Load ROM part of F and Q vector
      for (int i=0; i<ci.numROMVars; i++)
      {
        // F vector

        fVec[ci.li_ROM[i]] += ci.Fhat[ci.numExtVars+i];
        // Q vector

        qVec[ci.li_ROM[i]] += ci.Qhat[i];
      }
    }

    // TWO-LEVEL STAMP
    else
    {
      for (int i=0; i<ci.numExtVars; ++i)
      {
        fVec[ci.extLIDVec[i]] += ci.Fstamp[i];
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << " ----------------------------------" << std::endl;
    Xyce::dout() << " Master::loadDAEMatrices: " << std::endl;
  }

  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ci = *(*it);

    if(ci.usePortDesc==0)
    {

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
      {
        Xyce::dout() << " loads for ROM " << ci.getName() << std::endl;
      }

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      // Load dF/dx
      for (int i=0; i<ci.numExtVars; ++i)
      {

        *(ci.fEqu_up_NodePtr[i]) += 1.0;

        *(ci.fEqu_ip_NodePtr[i]) += 1.0;
      }

      // Load -Lhat portion of dF/dx
      for (int i=0; i<ci.numExtVars; i++)
      {
        double * LhatPtr = ci.fROMEqu_Lhat_VarsPtrs[i];
        for (int j=0; j<ci.numROMVars; j++)
        {

          LhatPtr[j] -= ci.Lhat[i*ci.numROMVars + j];
        }
      }

      // Load -Bhat portion of dF/dx
      for (int i=0; i<ci.numROMVars; i++)
      {
        for (int j=0; j<ci.numExtVars; j++)
        {

          *(ci.fROMEqu_Bhat_VarsPtrs[ci.numExtVars*i+j]) -= ci.Bhat[j*ci.numROMVars + i];
        }
      }

      // Load Ghat portion of dF/dx
      if (ci.isGSparse)
      {
        int nnz = ci.Ghat_rowPtr[ci.numROMVars];
        for (int i=0; i<nnz; ++i)
          *(ci.fROMEqu_Ghat_VarsPtrs[i]) += ci.Ghat[i];
      }
      else
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          double * GhatPtr = ci.fROMEqu_Ghat_VarsPtrs[i];
          for (int j=0; j<ci.numROMVars; j++)
          {

            GhatPtr[j] += ci.Ghat[j*ci.numROMVars + i];
          }
        }
      }

      // Load dQ/dx
      // Load Chat portion of dQ/dx
      if (ci.isCSparse)
      {
        int nnz=ci.Chat_rowPtr[ci.numROMVars];
        for(int i=0; i<nnz; i++)
          *(ci.qROMEqu_Chat_VarsPtrs[i]) += ci.Chat[i];
      }
      else
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          double * ChatPtr = ci.qROMEqu_Chat_VarsPtrs[i];
          for (int j=0; j<ci.numROMVars; j++)
          {

            ChatPtr[j] += ci.Chat[j*ci.numROMVars + i];
          }
        }
      }

#else

      for (int i=0; i<ci.numExtVars; i++)
      {

        dFdx[ci.extLIDVec[i]][ci.AEqu_up_NodeOffset[i]] += 1.0;

        dFdx[ci.intLIDVec[i]][ci.AEqu_ip_NodeOffset[i]] += 1.0;
      }

      // Load -Lhat portion of dF/dx
      for (int j=0; j<ci.numROMVars; j++)
      {
        for (int i=0; i<ci.numExtVars; i++)
        {

          dFdx[ci.intLIDVec[i]][ci.ROMEqu_Lt_NodeOffset[j]] -= ci.Lhat[i*ci.numROMVars + j];
        }
      }

      // Load -Bhat portion of dF/dx
      for (int i=0; i<ci.numROMVars; i++)
      {
        for (int j=0; j<ci.numExtVars; j++)
        {

          dFdx[ci.li_ROM[i]][ci.ROMEqu_B_NodeOffset[i*ci.numExtVars+j]] -= ci.Bhat[j*ci.numROMVars + i];
        }
      }

      // Load Ghat portion of dF/dx
      if (ci.isGSparse)
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          for (int j=ci.Ghat_rowPtr[i]; j<ci.Ghat_rowPtr[i+1]; j++)
          {
            
            dFdx[ci.li_ROM[i]][ci.ROMEqu_G_NodeOffset[j]] += ci.Ghat[j];
          }
        }
      }
      else
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          for (int j=0; j<ci.numROMVars; j++)
          {

            dFdx[ci.li_ROM[i]][ci.ROMEqu_GpC_NodeOffset[j]] += ci.Ghat[j*ci.numROMVars + i];
          }
        }
      }

      // Load Chat portion of dQ/dx
      if (ci.isCSparse)
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          for (int j=ci.Chat_rowPtr[i]; j<ci.Chat_rowPtr[i+1]; j++)
          {
            
            dQdx[ci.li_ROM[i]][ci.ROMEqu_C_NodeOffset[j]] += ci.Chat[j];
          }
        }
      }
      else
      {
        for (int i=0; i<ci.numROMVars; i++)
        {
          for (int j=0; j<ci.numROMVars; j++)
          {

            dQdx[ci.li_ROM[i]][ci.ROMEqu_GpC_NodeOffset[j]] += ci.Chat[j*ci.numROMVars + i];
          }
        }
      }
#endif
    }
    else
    {
      for (int i=0; i<ci.numExtVars; i++)
      {
        for (int j=0; j<ci.numExtVars; j++)
        {
          dFdx[ci.extLIDVec[i]][ci.AEqu_NodeOffset[i][j]] += ci.Jstamp[i*ci.numExtVars + j];
        }
      }
    }
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("ROM")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("rom", 1)
      .registerModelType("rom", 1);
  }
}

} // namespace ROM
} // namespace Device
} // namespace Xyce
