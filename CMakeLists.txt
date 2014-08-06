project(hermes-examples)
	cmake_minimum_required(VERSION 2.8)
    
###########################################################################
#############//////---------------------------------\\\\\\#################
#############!!!!!! NEVER CHANGE THIS FILE DIRECTLY !!!!!!#################
################!!!!!! USE the "CMake.vars" file !!!!!!####################
#############\\\\\\---------------------------------//////#################
###########################################################################
#   __      __      __      __      __      __      __      __      __    #
# /    \  /    \  /    \  /    \  /    \  /    \  /    \  /    \  /    \  #
#| STOP || STOP || STOP || STOP || STOP || STOP || STOP || STOP || STOP | #
# \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  \ __ /  #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#   ||      ||      ||      ||      ||      ||      ||      ||      ||    #
#  /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \    /  \   #
# ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ ~~~~~~~ #
###########################################################################
#############//////---------------------------------\\\\\\#################
#############!!!!!! NEVER CHANGE THIS FILE DIRECTLY !!!!!!#################
################!!!!!! USE the "CMake.vars" file !!!!!!####################
#############\\\\\\---------------------------------//////#################
###########################################################################
  
  SET(WITH_1d YES)
  SET(WITH_2d-advanced YES)
    SET(WITH_acoustics YES)
    SET(WITH_advection-diffusion-reaction YES)
    SET(WITH_heat-transfer YES)
    SET(WITH_helmholtz YES)
    SET(WITH_euler YES)
    SET(WITH_elasticity-linear YES)
    SET(WITH_maxwell YES)
    SET(WITH_navier-stokes YES)
    SET(WITH_nernst-planck YES)
    SET(WITH_neutronics NO)
    SET(WITH_richards YES)
    SET(WITH_schroedinger YES)
    SET(WITH_wave-equation YES)
  SET(WITH_2d-benchmarks-general YES)
  SET(WITH_2d-benchmarks-nist YES)

  set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
  include(CommonTargetProperties)

  # For Win64
  if(${CMAKE_CL_64})
    set(WIN64 YES)
  else(${CMAKE_CL_64})
    set(WIN64 NO)
  endif(${CMAKE_CL_64})

  # Where to look for the static libraries.
  set(HERMES_DIRECTORY /usr/local/lib)
  set(HERMES_INCLUDE_PATH /usr/local/include)
  set(HERMES_SOURCE_DIRECTORY ~/hermes)
  
  # Allow to override the default values in CMake.vars:
  include("CMake.vars.default.Linux" OPTIONAL)
  include(CMake.vars OPTIONAL RESULT_VARIABLE CMAKE_VARS_LOADED)
  if(${CMAKE_VARS_LOADED} MATCHES "NOTFOUND")
    message(FATAL_ERROR "You need to have a CMake.vars file - see the documentation and the example CMake.vars.example.Windows|Linux files for instructions.")
  endif()

  if(DEFINED HERMES_CMAKE_VARS_FILE_LOCATION)
    if(NOT ${HERMES_CMAKE_VARS_FILE_LOCATION} STREQUAL "")
      include(${HERMES_CMAKE_VARS_FILE_LOCATION} OPTIONAL RESULT_VARIABLE HERMES_CMAKE_VARS_LOADED)
      if(${HERMES_CMAKE_VARS_LOADED} MATCHES "NOTFOUND")
        message(FATAL_ERROR "Your specification of HERMES_CMAKE_VARS_FILE_LOCATION in your CMake.vars file is wrong (file does not exist?).")
      endif()
    endif()
  endif()
  
   if(NOT DEFINED HERMES_STATIC_LIBS)
     message(FATAL_ERROR "You must either set HERMES_CMAKE_VARS_FILE_LOCATION in your CMake.vars file, or use a CMake.vars file according to the instructions in the debian package if you installed Hermes using that package.")
   endif()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Disable all warnings and turn on only important ones:
    IF(DEFINED OLD_GCC)
      set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -fPIC")
    ELSE()
      set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -fPIC -std=c++11")
    ENDIF()
    set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

    set(HERMES_RELEASE_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -O3")
    set(HERMES_DEBUG_FLAGS   "${CMAKE_CXX_FLAGS} -g -D_DEBUG")
    
    # look for libbfd
    find_library( BFD_LIBRARY  NAMES bfd libbfd bfd-2.23.2-system PATH /usr/lib /usr/lib64 )

    include(CheckCSourceCompiles)
    check_c_source_compiles(
    "#include <bfd.h>
    int main(void) {
    return 0;
    } " BFD_WORKS)

    if (BFD_LIBRARY AND BFD_WORKS AND WITH_BFD)
      set( HAVE_BFD YES )
    else()
      set(BFD_LIBRARY "")
    endif()
    
    if(H2D_DEBUG)
      set(HERMES_FLAGS ${HERMES_DEBUG_FLAGS})
    else()
      set(HERMES_FLAGS ${HERMES_RELEASE_FLAGS})
    endif()
  endif()

  # This overrides CXX flags for MSVC
  if(MSVC)
    include(MSVC_FLAGS)
  endif(MSVC)
  
 if(WITH_TC_MALLOC)
    find_package(TCMALLOC REQUIRED)
    include_directories(${TCMALLOC_INCLUDE_DIR})
  endif(WITH_TC_MALLOC)
  
  if(WITH_PJLIB)
    find_package(PJLIB REQUIRED)
    include_directories(${PJLIB_INCLUDE_DIR})
  endif(WITH_PJLIB)
  
  IF(DEFINED BSON_LIBRARY)
    IF(DEFINED BSON_INCLUDE_DIR)
      include_directories(${BSON_INCLUDE_DIR})
    ELSE()
      MESSAGE(FATAL_ERROR When BSON_LIBRARY is defined, so must be BSON_INCLUDE_DIR.)
    ENDIF()
  ELSE()
    if(WITH_BSON)
      find_package(BSON REQUIRED)
      include_directories(${BSON_INCLUDE_DIR})
    endif(WITH_BSON)
  ENDIF()

  IF(DEFINED MATIO_LIBRARY)
    IF(DEFINED MATIO_INCLUDE_DIR)
      include_directories(${MATIO_INCLUDE_DIR})
    ELSE()
      MESSAGE(FATAL_ERROR When MATIO_LIBRARY is defined, so must be MATIO_INCLUDE_DIR.)
    ENDIF()
  ELSE()
    if(WITH_MATIO)
      find_package(MATIO REQUIRED)
      include_directories(${MATIO_INCLUDE_DIR})
    endif(WITH_MATIO)
  ENDIF()

  find_package(XSD REQUIRED)
  find_package(XERCES REQUIRED)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(WITH_STACKTRACE)
      find_package(EXECINFO)
    endif(WITH_STACKTRACE)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  # Solvers.
  if(WITH_UMFPACK)
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIRS})
  endif(WITH_UMFPACK)

  if(WITH_MATIO AND MATIO_WITH_HDF5)
    find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIR})
  else()
    set(MATIO_WITH_HDF5 NO)
  endif()

  IF(DEFINED PARALUTION_LIBRARY)
    IF(DEFINED PARALUTION_INCLUDE_DIR)
      include_directories(${PARALUTION_INCLUDE_DIR})
    ELSE()
      MESSAGE(FATAL_ERROR When PARALUTION_LIBRARY is defined, so must be PARALUTION_INCLUDE_DIR.)
    ENDIF()
  ELSE()
    if(WITH_PARALUTION)
    find_package(PARALUTION REQUIRED)
    include_directories(${PARALUTION_INCLUDE_DIR})
    endif(WITH_PARALUTION)
  ENDIF()

  if(MSVC)
    if(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
        # The search for CLAPACK and associated libraries will be performed
        # every time the project is configured in MSVC (TODO: Find out if
        # CLAPACK_LIBRARY, BLAS_LIBRARY, F2C_LIBRARY and CLAPACK_INCLUDE_DIRS
        # could be cached to improve subsequent configuration performance).
        find_package(CLAPACK REQUIRED)
        include_directories(${CLAPACK_INCLUDE_DIRS})
    endif(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
    if(WITH_MUMPS)
      find_package(WINBLAS REQUIRED)
    endif(WITH_MUMPS)
  else(MSVC)
    # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
    # FindLAPACK is Fortran-only so Fortran must be enabled as well.
    if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
      if (NOT LAPACK_FOUND)
          enable_language(Fortran)
          find_package(LAPACK REQUIRED)
            
        # If no error occured, LAPACK library has been found. Save the path to
        # it to cache, so that it will not be searched for during next 'cmake .'
          set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES}
                  CACHE STRING  "Path to LAPACK/BLAS libraries.")
          set(LAPACK_FOUND      YES
                  CACHE STRING  "Have LAPACK/BLAS libraries been found?")
      endif (NOT LAPACK_FOUND)
      add_definitions(-DWITH_BLAS)  
    endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
  endif(MSVC)

  if(WITH_PETSC)
    find_package(PETSC REQUIRED)
    # For some header files, PETSc uses the same name for both real and complex
    # versions (they only reside in different directories), so we can't globally
    # set the architecture-specific include directory here, but rather on a
    # per-target basis.
    include_directories(${PETSC_INCLUDE_DIRS})
    include_directories(${PETSC_CPLX_INCLUDE_DIRS})
    set(PETSC_LIBRARIES ${PETSC_CPLX_LIBRARIES})
  endif(WITH_PETSC)

  if(WITH_MUMPS)
    find_package(MUMPS REQUIRED)
    if(WITH_MPI)
      find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS
    endif(WITH_MPI)
    include_directories(${MUMPS_INCLUDE_DIR})
    set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})
    LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
    if(NOT(MSVC))
      LIST(APPEND MUMPS_LIBRARIES ${GFORTRAN_LIBRARY})
    endif()
  endif(WITH_MUMPS)

  if(WITH_TRILINOS)
    find_package(TRILINOS REQUIRED)
    include_directories(${TRILINOS_INCLUDE_DIR})
  endif(WITH_TRILINOS)
      
  if(WITH_SUPERLU)
    find_package(SUPERLU REQUIRED)
    include_directories(${SUPERLU_INCLUDE_DIR})
    if(SUPERLU_MT)
      add_definitions(-DSLU_MT)
    endif(SUPERLU_MT)
  endif(WITH_SUPERLU)

  # SMP.
  if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(WITH_OPENMP)

  if(WITH_EXODUSII)
    find_package(EXODUSII REQUIRED)
    include_directories(${EXODUSII_INCLUDE_DIR})
  endif(WITH_EXODUSII)

  # If using any package that requires MPI (e.g. parallel versions of MUMPS, PETSC).
  if(WITH_MPI)
    if(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH) # If MPI was not defined by the user
      find_package(MPI REQUIRED)
    else(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
      if(NOT EXISTS MPI_INCLUDE_PATH)
        message(FATAL_ERROR "MPI include directory ${MPI_INCLUDE_PATH} not found")
      endif(NOT EXISTS MPI_INCLUDE_PATH)
    endif(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
    include_directories(${MPI_INCLUDE_PATH})
  endif(WITH_MPI)
    
  if(NOT H2D_WITH_GLUT)
    set(NOGLUT  YES)
  endif(NOT H2D_WITH_GLUT)
  
  if(HAVE_BFD)
    find_package(LIBIBERTY REQUIRED)
    include_directories(${LIBIBERTY_INCLUDE_DIR})
  endif()
  
  if(H2D_WITH_GLUT)
    if(NOT APPLE)   # See hermes/OSX.cmake for APPLE-specific search procedures.
      find_package(GLUT REQUIRED)
      find_package(GLEW REQUIRED)
    endif(NOT APPLE)
    # For views synchronization, pthreads are required.
    find_package(PTHREAD REQUIRED)
    include_directories(${PTHREAD_INCLUDE_DIR})
  else(H2D_WITH_GLUT)
    add_definitions(-DNOGLUT)
  endif(H2D_WITH_GLUT)

  # Mesh format.
  if(WITH_EXODUSII)
    find_package(EXODUSII REQUIRED)
    include_directories(${EXODUSII_INCLUDE_DIR})
  endif(WITH_EXODUSII)
  
  include_directories(${XSD_INCLUDE_DIR})
  include_directories(${XERCES_INCLUDE_DIR})
  
  set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
  set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

  find_package(HERMES_COMMON REQUIRED)

  include_directories(${HERMES_COMMON_INCLUDE_PATH})
  include_directories(${HERMES2D_INCLUDE_PATH})

  IF(WITH_1d)
    add_subdirectory(1d)
	ENDIF(WITH_1d)
  IF(WITH_2d-advanced)
    add_subdirectory(2d-advanced)
  ENDIF(WITH_2d-advanced)
  IF(WITH_2d-benchmarks-general)
    add_subdirectory(2d-benchmarks-general)
  ENDIF(WITH_2d-benchmarks-general)
	IF(WITH_2d-benchmarks-nist)
    add_subdirectory(2d-benchmarks-nist)
  ENDIF(WITH_2d-benchmarks-nist)

  #------------------------------------------------------------------------------
  # Report.
  #~~~~~~~~
  message("\nHermes-examples Configuration results")
  message("---------------------")
  if(WIN32)
    message("Build for Win64: ${WIN64}")
  endif(WIN32)

  message("\n-------Features-------")
  message("Build with OpenMP: ${WITH_OPENMP}")
  message("Build with TCMalloc: ${WITH_TC_MALLOC}")
  message("Build with BSON: ${WITH_BSON}")
  message("Build with MATIO: ${WITH_MATIO}")
  if(${WITH_MATIO})
    message(" MATIO with HDF5: ${MATIO_WITH_HDF5}")
  endif()
  if(${WITH_MPI})
    message("Build with MPI: ${WITH_MPI}")
  endif()
  message("Build with EXODUSII: ${WITH_EXODUSII}")

  message("\n-------Solvers-------")
  message("Build with UMFPACK: ${WITH_UMFPACK}")
  message("Build with PARALUTION: ${WITH_PARALUTION}")
  message("Build with PETSC: !TO_DO!")
  message("Build with MUMPS: ${WITH_MUMPS}")
  message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
  message("Build with TRILINOS: ${WITH_TRILINOS}")

  message("---------------------")
  message("\n")
