# -----------------------------------------------------------------------------
# Default Xyce Config/Build/Test Options 
# -----------------------------------------------------------------------------

# -- configure the build environment ------------------------------------------

# set CMake project name
project ( Xyce )

# set Xyce release version string
set ( Xyce_VERSION_MAJOR           7                CACHE STRING " " )
set ( Xyce_VERSION_MINOR           1                CACHE STRING " " )
set ( Xyce_VERSION_PATCH           0                CACHE STRING " " )
set ( Xyce_VERSION_EXTRA           "Q"    CACHE STRING " " )

# DEBUG:  rewrite the version code; be sure to use quotes
set ( VERSION                      "R:Q:7.1"         CACHE INTERNAL " " )

# set the top-level arch directory (eg. /Net/Proj/Xyce/arch/linux)
# set ( Xyce_ARCHDIR                 "/Net/Proj/Xyce/arch/linux" )

# config options
option ( Xyce_FORTRAN_TEST         "test for fortran compiler"                     ON  )


# target options
option ( Xyce_ENABLE_SHARED        "build shared libraries"                        OFF )
option ( Xyce_ENABLE_STATIC        "build static libraries"                        ON  )

option ( Xyce_DISABLE_BINARY       "disable building of Xyce binary"               OFF )
option ( Xyce_STATIC_BINARY        "create statically linked binary"               OFF )


# build options
option ( Xyce_AMD                  "enable AMD reordering"                         ON )
option ( Xyce_BELOS                "enable Belos iterative solver support"         ON )
option ( Xyce_CHARON               "Charon device support"                         OFF )
option ( Xyce_Dakota               "Dakota direct linkage support"                 OFF )
option ( Xyce_DEPENDENCY_TRACKING  "do not reject slow dependency extractors"      OFF )
option ( Xyce_EXTDEV               "Ext device extensions"                         ON  )
option ( Xyce_ML                   "enable ML preconditioner support"              OFF )
option ( Xyce_NEW_DAE_FORMULATION  "new DAE formulation"                           OFF )
option ( Xyce_OP_START             "allow startup from previous operating point"   ON  )
option ( Xyce_PARALLEL_MPI         "enable parallel build with MPI"                OFF )
option ( Xyce_PARDISO_MKL          "Pardiso direct solver"                         OFF )
option ( Xyce_RAD_MODELS           "Radiation Model Library"                       ON  )
option ( Xyce_NONFREE_MODELS       "Non-Free Models Library"                       ON  )
option ( Xyce_ADMS_MODELS          "ADMS Models Library"                           ON  )
option ( Xyce_NEURON_MODELS        "Neuron Models Library"                         ON  )
option ( Xyce_REACTION_PARSER      "enable reaction parser"                        ON  )
option ( Xyce_SENSITIVITY_ENABLE   "Use sensitivity capability"                    OFF )
option ( Xyce_SHYLU                "ShyLU hybrid solver support"                   OFF )
option ( Xyce_SPICE_NORMS          "using SPICE type norms"                        ON  )
option ( Xyce_SUPERLU              "SuperLU direct solver"                         OFF )
option ( Xyce_SUPERLUDIST          "SuperLU_Dist linear solver"                    OFF )
option ( Xyce_TRILINOS_DEV         "Trilinos development support"                  OFF )
option ( Xyce_UMFPACK              "UmfPack direct solver"                         ON  )
option ( Xyce_USE_BSIM3_CONST      "using constants from BSIM3"                    OFF )
option ( Xyce_USE_ISORROPIA        "Use Isorropia for partitioning"                ON  )
option ( Xyce_USE_ZOLTAN           "Use Zoltan for partitioning"                   ON  )
option ( Xyce_USE_PARMETIS         "Use ParMETIS for partitioning"                 ON  )
option ( Xyce_ATHENA               "Build the ATHENA device"                       ON )
option ( Xyce_ADMS_SENSITIVITIES               "Enable analytic sensitivities in ADMS-generated devices"                       ON )
option ( Xyce_REPO_ACCESS       "access to xyce source git repository"             OFF )

# verbosity options
option ( Xyce_VERBOSE_LINEAR       "verbosity in linear solver"                    OFF )
option ( Xyce_VERBOSE_NONLINEAR    "verbosity in nonlinear solver"                 OFF )
option ( Xyce_VERBOSE_NOX          "verbosity in NOX nonlinear solver library"     OFF )
option ( Xyce_VERBOSE_TIME         "verbosity in time integrator"                  OFF )

# debugging options
option ( Xyce_DEBUG                "enable debugging symbols"                      OFF )
option ( Xyce_DEBUG_ANALYSIS       "analysis package"                              OFF )
option ( Xyce_DEBUG_CIRCUIT        "circuit package"                               OFF )
option ( Xyce_DEBUG_DEVICE         "device package"                                OFF )
option ( Xyce_DEBUG_DIRECTSOLVE    "direct solver package"                         OFF )
option ( Xyce_DEBUG_DISTRIBUTION   "distribution package"                          OFF )
option ( Xyce_DEBUG_EXPRESSION     "expression package"                            OFF )
option ( Xyce_DEBUG_IO             "I/O interface package"                         OFF )
option ( Xyce_DEBUG_LINEAR         "linear solver package"                         OFF )
option ( Xyce_DEBUG_NONLINEAR      "nonlinear solver package"                      OFF )
option ( Xyce_DEBUG_PARALLEL       "parallel distribution package"                 OFF )
option ( Xyce_DEBUG_RESTART        "restart"                                       OFF )
option ( Xyce_DEBUG_TIME           "time integrator package"                       OFF )
option ( Xyce_DEBUG_TOPOLOGY       "topology package"                              OFF )
option ( Xyce_TEST_SOLN_VAR_MAP    "mappping of internal variables to external"    OFF )
