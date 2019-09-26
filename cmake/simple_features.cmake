
#This is being set here, because the check for iota in "check_include_file_cxx.cmake" is not working
set (HAVE_IOTA true)

#### Open Source ONLY!!!!!!! ########
#add_definitions( -DXyce_USE_FFTW)
add_definitions( -DXyce_USE_FFT) 
set(Xyce_USE_FFTW true)

#### Only on the systems that have bison and flex for now #######
# Also this is both a C++ and CMAKE variable set. There are both ifs 
# in the CmakeLists.txt in src/DeviceModelPKG/Core/CMakeLists.txt but 
# also #ifndefs inside the source code on the variable. 
#add_definitions( -DXyce_REACTION_PARSER) # This is a #cmakedefine in Xyce_config.h.cmake
set(Xyce_REACTION_PARSER true CACHE BOOL "The reaction parser depends on Flex and Bison. Turn this off if you don't need
this feature. This is a super sepcialized feature that is for chemical reactions." ) 
add_definitions( -DYYTEXT_POINTER) 
#
# If everything is in place, define `Xyce_REACTION_PARSER`
#* If everything is _not_ in place, give a warning (and make sure
#  `Xyce_REACTION_PARSER` is disabled)
#* Allow the developer to turn it off
#
# If Xyce_REACTION_PARSER is both a CMAKE variable and a C++ variable 
# the code in `src/DeviceModelPKG/Core` will
# through Bison and Flex. The `CMakeLists.txt` has an if on the CMAKE 
# variable Xyce_REACTION_PARSER which executes the Bison and Flex Find_package 
# macros names BISON_TARGET and FLEX_TARGET respectively. They're detailed in the
# CMAKE documentation for each find_package. 

