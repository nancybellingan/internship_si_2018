#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(CMAKE_C_COMPILER "gcc-7")
set(CMAKE_CXX_COMPILER "g++-7")

#-fsanitize=address
add_compile_options(-g -O0  -fno-omit-frame-pointer -Wall -Wextra -Werror -pedantic -Wshadow -Wno-unused-variable)

SET(ASAN_FLAGS    "-fsanitize=address")

#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${ASAN_FLAGS}" )

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
#find_package(Sanitizers)


project(SiDetPhant LANGUAGES C CXX)
set(CMAKE_VERBOSE_MAKEFILE ON)



set(CMAKE_CXX_STANDARD 11)




#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#-----------------------------------------------------------------------------
# User defines for Gemetry settings
#if one of the following option is ON, then the correspondent definition is enabled and back in the SRC files the ifdef with such definition is enabled and executed


option(WITH_SPECSHOOT "Using ProbVec.in for energy distribution" OFF)
if(WITH_SPECSHOOT)
	add_definitions(-DSPECSCHUSS)
endif()

#FLUKA, fluctuating cascade in german, simulation package for particles and nuclei in matter
option(WITH_FLUKA_PSF "Using Fluka particle wise phase space file" OFF)
if(WITH_FLUKA_PSF)
	add_definitions(-DFLUKA_PSF)
endif()

#range shiftler as in a special kind of proton therapy with a stopping material between beam and patient
option(WITH_RANGE_SHIFTER "Using Range Shifter" OFF)
if(WITH_RANGE_SHIFTER)
	add_definitions(-DRSHIFTER)
endif()

option(WITH_MFD "Using MultiFunctionalDetector" OFF)
if(WITH_MFD)
	add_definitions(-DMFD)
endif()

option(WITH_DEFMATERIAL "Using G4Materials" ON)
if(WITH_DEFMATERIAL)
	add_definitions(-DDEFMATERIAL)
endif()

# DUMP Voxel only together with VoxelPhantom
option(WITH_DUMPVOXEL "Dumping voxel values on screen" ON)
if(WITH_DUMPVOXEL)
	add_definitions(-DDUMPVOXEL)
endif()

option(WITH_DEBUG "User defined Debug information to standard output" ON)
if(WITH_DEBUG)
	add_definitions(-DDGB)
endif()


option(WITH_FAST_DET "Using Fast/Albedo detector on measurement position" ON)
if(WITH_FAST_DET)
	add_definitions(-DFAST_DET)
endif()

option(WITH_ISOT_SOURCE "Using isotropic source instead of parallel square source" OFF)
if(WITH_ISOT_SOURCE)
	add_definitions(-DISOT_SOURCE)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(SiDetPhant SiDetPhant.cc ${sources} ${headers})
target_link_libraries(SiDetPhant ${Geant4_LIBRARIES} )
#target_link_libraries(SiDetPhant ${ASAN_FLAGS} )



#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build SiDetPhant. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SiDetPhant_SCRIPTS
	macros/run1.mac macros/run.mac macros/plot.mac macros/vis.mac
	inputs/Binning.ebin scripts/singleSOBP.sh
	outputs/totED.txt outputs/Ergfile.dat outputs/ErgfileVoxels.dat
  )

foreach(_script ${SiDetPhant_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SiDetPhant DESTINATION bin)
