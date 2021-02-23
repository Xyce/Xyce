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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoeksra, SNL, Parallel Computational Sciences
//
// Creation Date  : 4/2/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <EpetraExt_LPTrans_From_GraphTrans.h>
#include <EpetraExt_LPTrans_From_MatrixTrans.h>

#include <EpetraExt_CrsSingletonFilter_LinearProblem.h>

#include <EpetraExt_SolverMap_LinearProblem.h>

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_USE_ISORROPIA
#include <EpetraExt_AmesosBTFGlobal_LinearProblem.h>
#include <EpetraExt_Isorropia_CrsGraph.h>
#endif
#endif

#ifdef Xyce_AMD
#include <EpetraExt_AMD_CrsGraph.h>
#ifdef Xyce_PARALLEL_MPI
#include <EpetraExt_AmesosAMDGlobal_CrsGraph.h>
#endif
#endif

#include <EpetraExt_AmesosBTF_CrsGraph.h>

#include <EpetraExt_Scale_LinearProblem.h>

#include <EpetraExt_Reindex_LinearProblem.h>

#include <Teuchos_ScalarTraits.hpp>

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>

#include <N_LAS_TransformTool.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_FeatureTest.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : TransformTool::operator()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/2/03
//-----------------------------------------------------------------------------
Teuchos::RCP<Transform>
TransformTool::
operator()( const Util::OptionBlock & options )
{
  int    rcm       = 0;
  int    rcmTestWidth = 0;

  int    sFilter   = 0;

  bool   scale     = false;
  int    lScale    = 0;
  int    rScale    = 0;
  double expScale  = 1.0;
  int    iterScale = 1;
  bool   btf       = false;
  int    globalbtf = 0;
  double globalbtftol = Teuchos::ScalarTraits<double>::eps();
  bool   globalbtfverb = false;
  int    globalamd = 0;
  bool   globalamdverb = false;

  int    partition = 1;   // Isorropia partitioning
  std::string partition_type="HYPERGRAPH";

  // AMD should not be used in serial if KLU is being used.
#ifdef Xyce_PARALLEL_MPI
  bool   amd       = true;
#else
  bool   amd       = false;
#endif
  bool   amd_verbose = false;

  bool   reindex   = true;
  bool   solverMap = true;

  for (Util::ParamList::const_iterator it = options.begin(), end = options.end(); it != end; ++it)
  {
    std::string tag = (*it).uTag();

    if     ( tag == "TR_SINGLETON_FILTER" )  sFilter = (*it).getImmutableValue<int>();
    else if( tag == "TR_SCALE" )             scale = (*it).getImmutableValue<int>();
    else if( tag == "TR_SCALE_LEFT" )        lScale = (*it).getImmutableValue<int>();
    else if( tag == "TR_SCALE_RIGHT" )       rScale = (*it).getImmutableValue<int>();
    else if( tag == "TR_SCALE_EXP" )         expScale = (*it).getImmutableValue<double>();
    else if( tag == "TR_SCALE_ITER" )        iterScale = (*it).getImmutableValue<int>();
    else if( tag == "TR_BTF" )               btf = (*it).getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_BTF" )        globalbtf = (*it).getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_BTF_DROPTOL" )globalbtftol = (*it).getImmutableValue<double>();
    else if( tag == "TR_GLOBAL_BTF_VERBOSE" )globalbtfverb = (*it).getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_AMD" )        globalamd = (*it).getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_AMD_VERBOSE" )globalamdverb = (*it).getImmutableValue<int>();
    else if( tag == "TR_PARTITION" )         partition = (*it).getImmutableValue<int>();
    else if( tag == "TR_PARTITION_TYPE" )    partition_type = (*it).usVal();
    else if( tag == "TR_AMD" )               amd = (*it).getImmutableValue<int>();
    else if( tag == "TR_AMD_VERBOSE" )       amd_verbose = (*it).getImmutableValue<int>();
    else if( tag == "TR_REINDEX" )           reindex = (*it).getImmutableValue<int>();
    else if( tag == "TR_SOLVERMAP" )         solverMap = (*it).getImmutableValue<int>();
  }

// Sort out which partitioning to use.
#ifndef Xyce_PARALLEL_MPI

  // If we are not running in parallel partitioning is not necessary.
  partition = 0;

#else

  // Turn off partitioning if Isorropia is not enabled.
#ifndef Xyce_USE_ISORROPIA
  if ( partition ) {
    partition = 0;
  }
#endif

  // Change graph to hypergraph partitioning if ParMETIS is not enabled.
#ifndef HAVE_LIBPARMETIS
  if ( partition && ( partition_type == "GRAPH" ) ) {
    Report::UserWarning0()
        << "TransformTool::operator():  ParMETIS not enabled, changing partitioning to HYPERGRAPH.";
    partition_type = "HYPERGRAPH";
  }
#endif

#endif // Xyce_PARALLEL_MPI

#ifndef Xyce_AMD
  amd = false;
#endif

  if (VERBOSE_LINEAR)
  {
    Xyce::lout() << "Linear Transforms" << std::endl
                 << "-----------------" << std::endl;
    if( sFilter )
      Xyce::lout() << "Singleton Filter" << std::endl;
    if( globalamd )
      Xyce::lout() << "Global AMD" << std::endl;
    if( globalbtf )
      Xyce::lout() <<  "Global BTF" << std::endl;
    if( scale )
      Xyce::lout() << "Scaling" << std::endl;
#ifdef Xyce_USE_ISORROPIA
    if( partition ) {
      Xyce::lout() << "Isorropia Partitioning (" << partition_type << ") " << std::endl;
    }
#endif
    if( amd )
      Xyce::lout() << "AMD" << std::endl;
    if( btf )
      Xyce::lout() << "BTF" << std::endl;
    if( reindex )
      Xyce::lout() << "Reindexing" << std::endl;
    if( solverMap )
      Xyce::lout() << "Column Remapping" << std::endl;
    Xyce::lout() << "-----------------" << std::endl;
  }

  Teuchos::RCP<Transform> CompTrans = Teuchos::rcp( new Transform() );

  bool TransFlag = false;

  // First remove any dense rows or columns using singleton filtering
  if( sFilter > 0 )
  {
    CompTrans->addTransform(
      dynamic_cast<EpetraExt::SameTypeTransform<Epetra_LinearProblem>*>
      (new EpetraExt::LinearProblem_CrsSingletonFilter(true)) );
    TransFlag = true;

#ifdef Xyce_PARALLEL_MPI
    if ( partition ) 
    {
      // Need to reindex because Isorropia does not do that at this time,
      // and symmetrization using different GID indices for row and column maps
      // can cause an issue.
      CompTrans->addTransform( new EpetraExt::LinearProblem_SolverMap() );

      CompTrans->addTransform( new EpetraExt::LinearProblem_Reindex( 0 ) );
    }
#endif

  }

  // Permute the global graph using AMD ordering (partitioning should be used in parallel)
  if ( globalamd )
  {
#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_AMD
    EpetraExt::AmesosAMDGlobal_CrsGraph * globalAMDTrans = new EpetraExt::AmesosAMDGlobal_CrsGraph( globalamdverb );
    EpetraExt::LinearProblem_GraphTrans * globalAMD_LPTrans =
      new EpetraExt::LinearProblem_GraphTrans(
        *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(globalAMDTrans)) );
    CompTrans->addTransform( globalAMD_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( globalAMDTrans ), "globalAMDTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
#endif
#endif
  }

  // Permute the global graph using BTF ordering and a block partitioning (partitioning should *not* be used in parallel)
  if( globalbtf )
  {
#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_USE_ISORROPIA
    std::string balanceType = "linear";
    if (globalbtf == 2)
      balanceType = "isorropia";
    EpetraExt::AmesosBTFGlobal_LinearProblem * GlobalBTF_LPTrans =
      new EpetraExt::AmesosBTFGlobal_LinearProblem( globalbtftol, balanceType, false, globalbtfverb );
    CompTrans->addTransform( GlobalBTF_LPTrans );
    TransFlag = true;
#endif // Xyce_USE_ISORROPIA
#endif
  }

  // Partition the global graph using Isorropia
#ifdef Xyce_PARALLEL_MPI
  if( partition )
  {
    // Create the parameter list and fill it with values.
    Teuchos::ParameterList paramlist;
    Teuchos::ParameterList& sublist = paramlist.sublist("ZOLTAN");
    if (partition_type == "HYPERGRAPH") {
      paramlist.set("PARTITIONING METHOD", "HYPERGRAPH");
      sublist.set("CHECK_GRAPH", "0");
      sublist.set("LB_APPROACH", "PARTITION");
    }
    else if (partition_type == "GRAPH") {
      paramlist.set("PARTITIONING METHOD", "GRAPH");
      sublist.set("CHECK_GRAPH", "0");
      sublist.set("LB_APPROACH", "PARTITION");
      sublist.set("GRAPH_PACKAGE", "PARMETIS");
      sublist.set("PARMETIS_METHOD", "PARTKWAY");
    }
    else {
      Report::DevelFatal() << "Selected Partitioning Type is Invalid!";
    }

    if (VERBOSE_LINEAR)
      sublist.set("DEBUG_LEVEL", "2" );

    EpetraExt::Isorropia_CrsGraph * ITrans = new EpetraExt::Isorropia_CrsGraph( paramlist );
    EpetraExt::LinearProblem_GraphTrans * I_LPTrans =
      new EpetraExt::LinearProblem_GraphTrans(
        *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(ITrans)) );
    CompTrans->addTransform( I_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( ITrans ), "ITrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
  }
#endif

  // Perform AMD ordering on the local graph
  if( amd )
  {
#ifdef Xyce_AMD
    EpetraExt::CrsGraph_AMD * AMDTrans = 0;
    AMDTrans = new EpetraExt::CrsGraph_AMD( amd_verbose );
    EpetraExt::LinearProblem_GraphTrans * AMD_LPTrans =
      new EpetraExt::LinearProblem_GraphTrans(
	*(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(AMDTrans)) );
    CompTrans->addTransform( AMD_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( AMDTrans ), "AMDTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
#else
    Report::DevelFatal() << "TR_AMD not Enabled!";
#endif
  }

  // Perform BTF ordering on the local graph
  if( btf )
  {
    EpetraExt::AmesosBTF_CrsGraph * BTFTrans = new EpetraExt::AmesosBTF_CrsGraph();
    EpetraExt::LinearProblem_GraphTrans * BTF_LPTrans =
      new EpetraExt::LinearProblem_GraphTrans(
	*(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(BTFTrans)) );
    CompTrans->addTransform( BTF_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( BTFTrans ), "BTFTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
  }

  // Scale the linear problem
  if ( scale )
  {
    CompTrans->addTransform(
      new EpetraExt::LinearProblem_Scale(
        static_cast<EpetraExt::LinearProblem_Scale::ScaleType>(lScale),
        static_cast<EpetraExt::LinearProblem_Scale::ScaleType>(rScale),
        expScale,
        iterScale ) );
    TransFlag = true;
  }

  if( solverMap && TransFlag )
  {
    CompTrans->addTransform( new EpetraExt::LinearProblem_SolverMap() );
    TransFlag = true;
  }

  if( reindex && TransFlag )
  {
    CompTrans->addTransform( new EpetraExt::LinearProblem_Reindex( 0 ) );
    TransFlag = true;
  }

  if( TransFlag ) return CompTrans;
  else
  {
    return Teuchos::null;
  }

}

} // namespace Linear
} // namespace Xyce
