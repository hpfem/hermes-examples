project(hermes-examples)
	cmake_minimum_required(VERSION 2.8)
  
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
    SET(WITH_neutronics YES)
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

  # This has to be the same as in the library's CMake.
  # Default "yes"
  set(H2D_WITH_GLUT YES)
  set(WITH_TC_MALLOC YES)
  set(WITH_PARALUTION YES)
  set(WITH_UMFPACK YES)
  set(WITH_OPENMP YES)
  # Default "no"
  set(WITH_EXODUSII NO)
  set(WITH_MUMPS NO)
  set(WITH_MATIO NO)
  set(WITH_BSON NO)
  set(WITH_SUPERLU NO)

  # Where to look for the static libraries.
  set(HERMES_DIRECTORY /usr/local/lib)
  set(HERMES_INCLUDE_PATH /usr/local/include)
  set(DEP_INCLUDE_PATHS /usr/local/include)

  # Allow to override the default values in CMake.vars:
  include(CMake.vars OPTIONAL)

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Disable all warnings and turn on only important ones:
    set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

    set(RELEASE_FLAGS "-DNDEBUG -O3")
    set(DEBUG_FLAGS   "-g")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  # Enabling multiprocessor build on MSVC
  if(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    if(WITH_OPENMP)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    endif(WITH_OPENMP)
  endif(MSVC)

  # This overrides CXX flags for MSVC
  if(MSVC)
    if(WIN64)
      set(MSVC_DEFINES "/DWIN64 /D_WINDOWS /Dfinite=_finite /wd4275 /wd4251")
    else(WIN64)
      set(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dfinite=_finite /wd4275 /wd4251")
    endif(WIN64)
    set(CMAKE_CXX_FLAGS_DEBUG          "/D_DEBUG /Od /Ob2 /MDd /Zi ${MSVC_DEFINES}")
    set(HERMES_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(CMAKE_CXX_FLAGS_RELEASE        "/DNDEBUG /O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(HERMES_RELEASE_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/DNDEBUG /O2 /Ob2 /MD /Zi ${MSVC_DEFINES}")
    
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO")
    set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO")
    set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO")
  endif(MSVC)

  include_directories("${PTHREAD_ROOT}/include")

  if(${WITH_UMFPACK})
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIRS})
  endif()

  find_package(PTHREAD REQUIRED)

  if(${WITH_TC_MALLOC})
    find_package(TCMALLOC REQUIRED)
    include_directories(${TCMALLOC_INCLUDE_DIR})
  endif()

  if(${WITH_PARALUTION})
    find_package(PARALUTION REQUIRED)
    include_directories(${PARALUTION_INCLUDE_DIR})
  endif()

  if(${WITH_MATIO})
    find_package(MATIO REQUIRED)
    include_directories(${MATIO_INCLUDE_DIR})
  endif()

  if(MSVC)
    if(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
        find_package(CLAPACK REQUIRED)
        include_directories(${CLAPACK_INCLUDE_DIRS})
    endif(WITH_PETSC OR WITH_TRILINOS OR WITH_SUPERLU)
    if(WITH_MUMPS)
      find_package(WINBLAS REQUIRED)
    endif(WITH_MUMPS)
  else(MSVC)
    if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
      if (NOT LAPACK_FOUND)
          enable_language(Fortran)
          find_package(LAPACK REQUIRED)
          set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES}
                  CACHE STRING  "Path to LAPACK/BLAS libraries.")
          set(LAPACK_FOUND      YES
                  CACHE STRING  "Have LAPACK/BLAS libraries been found?")
      endif (NOT LAPACK_FOUND)
      add_definitions(-DWITH_BLAS)  
    endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
  endif(MSVC)

  if(WITH_MUMPS)
    find_package(MUMPS REQUIRED)
    if(WITH_MPI)
      find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS
    endif(WITH_MPI)
    include_directories(${MUMPS_INCLUDE_DIR})
    set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})
    LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
  endif(WITH_MUMPS)

  if(WITH_TRILINOS)
    find_package(TRILINOS REQUIRED)
    include_directories(${TRILINOS_INCLUDE_DIR})
  endif(WITH_TRILINOS)
      
  if(${WITH_BSON})
    find_package(BSON REQUIRED)
    include_directories(${BSON_INCLUDE_DIR})
  endif()
        
  if(WITH_SUPERLU)
    find_package(SUPERLU REQUIRED)
    include_directories(${SUPERLU_INCLUDE_DIR})
    if(SUPERLU_MT)
      add_definitions(-DSLU_MT)
    endif(SUPERLU_MT)
  endif(WITH_SUPERLU)
  set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
  set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

  find_package(HERMES_COMMON REQUIRED)

  include_directories(${HERMES_COMMON_INCLUDE_PATH})
  include_directories(${HERMES2D_INCLUDE_PATH})
  include_directories(${DEP_INCLUDE_PATHS})

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
  message("Build with TRILINOS: !TO_DO!")

  message("---------------------")
  message("\n")
