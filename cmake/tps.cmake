# This is to define virtual target out of some of the Third Party Softwares
# that provide clumps of variables but not nice little virtual target 
# keywords, like 'fftw3'

# The below text is copied directly from Bryan Hughes' tps.cmake file on another project. It's art.
# ---
# include third party info

# ---------------------------------------------------------------------
# ***** READ THE TEXT BELOW BEFORE YOU CHANGE ANYTHING IN THIS FILE *****

# You are modifying a file that affects the ENTIRE BUILD. Before making changes, understand what you are doing!

# DO NOT add "include_directories" or "link_directories" statements to this file. Yes, I am talking to YOU.

# Try the following instead:
# * If there is a Find<yourlibrary>.cmake file, add find_package(yourlibrary REQUIRED) to this file.
# * If there is NOT a Find<yourlibrary>.cmake file, and you are are unable to find one on the internet,
#   then "set" the following variables:
#     a variable named <yourlibrary>_INCLUDE_DIR (or something similiar)
#     a variable named <yourlibrary>_LIBRARIES (or something similar)

# After you have done one of the things above, change the leaf CMakeLists.txt file inside your target directory to
# include the directories or libraries, local to your target.

#   THE WHOLE WORLD DOES NOT NEED YOUR THIRD PARTY LIBRARY'S INCLUDE AND LINKER PATHS.

# If you do not understand these instructions or are otherwise unable to follow them, ASK FOR HELP.

# If you disregard this message, I reserve the right to spray you with a hose, because you are a bad person.

# ***** READ THE TEXT ABOVE BEFORE YOU CHANGE ANYTHING IN THIS FILE *****
# ---------------------------------------------------------------------

###################
## Find libraries.
###################

# JCV NOTE: BLAS, LAPACK, and AMD should be probed via Trilinos!!!!
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(SuiteSparse OPTIONAL_COMPONENTS AMD)
add_library(SuiteSparse::AMD INTERFACE IMPORTED GLOBAL)
set_target_properties(SuiteSparse::AMD PROPERTIES 
	INTERFACE_INCLUDE_DIRECTORIES "${SUITESPARSE_INCLUDE_DIR}"
	INTERFACE_LINK_LIBRARIES "${SuiteSparse_LIBRARIES}")

# From Ross Bartlett "You just do find_package(Trilinos REQUIRED) and it reads
# in all the vars defined in the installed file TrilinosConfig.cmake.  If you
# look at that file, you should see it has all of the info you need."
# Serial:
find_package(Trilinos REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos
     Sacado Triutils AztecOO Belos TrilinosCouplings OPTIONAL_COMPONENTS
     ShyLU Basker Amesos2 Stokhos ROL)
# Parallel:
#   find_package(Trilinos REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos
#        Sacado Triutils AztecOO Belos TrilinosCouplings Isorropia Zoltan
#        OPTIONAL_COMPONENTS ShyLU Basker Amesos2 Stokhos ROL)

message(" Your Trilinos_DIR is " ${Trilinos_DIR})
add_library(trilinos INTERFACE IMPORTED GLOBAL)
set_target_properties(trilinos PROPERTIES
   	INTERFACE_INCLUDE_DIRECTORIES "${Trilinos_INCLUDE_DIRS}" 
)
# find the right fftw library.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   message("Using Intel provided MKL, the Math Kernel Library")
   find_package(MKL)
else ()
   message("Using FFTW3")
   find_package(FFTW REQUIRED)
   add_library(FFTW::FFTW INTERFACE IMPORTED GLOBAL)
   set_target_properties(FFTW::FFTW PROPERTIES
   	INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}" 
	INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}")

endif ()

if (Xyce_REACTION_PARSER)
   find_package(BISON 2.3 REQUIRED)
   find_package(FLEX REQUIRED)
   #[[ "Define YYTEXT_POINTER if yytext defaults to 'char *' instead of 'char
   []'"]]
endif ()
