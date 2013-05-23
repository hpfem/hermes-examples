set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_HOME_DIRECTORY}/cmake/cxx_flag_overrides.cmake)
project(hermes_examples)
	cmake_minimum_required(VERSION 2.6)
  
  set(HERMES_EXAMPLES_DEBUG YES)

	set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
        include(CommonTargetProperties)
	
	# Where to look for the static libraries.
	set(HERMES_DIRECTORY /usr/local/lib)
	set(HERMES_INCLUDE_PATH /usr/local/include)
	set(DEP_INCLUDE_PATHS /usr/local/include)

	# Trilinos
	# Enable support for Trilinos solvers.
	set(WITH_TRILINOS           NO)
				
	# Experimental
	set(WITH_ZOLTAN             NO)
	# If MPI is enabled, the MPI library installed on the system should be found by 
  # CMake automatically. If the found library doesn't match the one used to compile the
  # particular MPI-dependent package, the other two options should be used to specify it.
  #
  set(WITH_MPI                NO)
  # set(MPI_LIBRARIES         -lmpi)
  # set(MPI_INCLUDE_PATH      /usr/include/openmpi
  

SET(WITH_1d NO)
  SET(WITH_2d-advanced YES)
    SET(WITH_acoustics NO)
    SET(WITH_advection-diffusion-reaction NO)
    SET(WITH_heat-transfer NO)
    SET(WITH_helmholtz YES)
    SET(WITH_euler YES)
    SET(WITH_elasticity-linear YES)
    SET(WITH_maxwell YES)
    SET(WITH_navier-stokes NO)
    SET(WITH_nernst-planck NO)
    SET(WITH_neutronics NO)
    SET(WITH_richards NO)
    SET(WITH_schroedinger NO)
    SET(WITH_wave-equation NO)
  SET(WITH_2d-benchmarks-general NO)
  SET(WITH_2d-benchmarks-nist NO)
  
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)

	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
			# Disable all warnings and turn on only important ones:
			set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS}")
			set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

			set(RELEASE_FLAGS "-DNDEBUG -O3")
			set(DEBUG_FLAGS   "-g")
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

	# Enabling multiprocessor build on MSVC
	if(MSVC)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") 
	endif(MSVC)
	
  include_directories("${PTHREAD_ROOT}/include")
	find_package(PTHREAD REQUIRED)
  find_package(TCMALLOC REQUIRED)
  include_directories(${TCMALLOC_INCLUDE_DIR})

	if(WITH_TRILINOS)
		find_package(TRILINOS REQUIRED)
		include_directories(${TRILINOS_INCLUDE_DIR})
	endif(WITH_TRILINOS)
	
	if(WITH_MPI)
		if(NOT MPI_INCLUDE_PATH)
			find_package(MPI REQUIRED)
		endif(NOT MPI_INCLUDE_PATH)
		include_directories(${MPI_INCLUDE_PATH})
	endif(WITH_MPI)
	
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
  
  
  
  
  
  