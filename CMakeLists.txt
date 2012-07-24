project(hermes_examples)
	cmake_minimum_required(VERSION 2.6)

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
  
  # By default all examples are turned on.
  # This can be changed in your CMake.vars file.
  SET(WITH_1d YES)
  SET(WITH_2d-advanced YES)
    SET(WITH_acoustics YES)
    SET(WITH_advection-diffusion-reaction YES)
    SET(WITH_euler YES)
    SET(WITH_heat-transfer YES)
    SET(WITH_helmholtz YES)
    SET(WITH_elasticity-linear YES)
    SET(WITH_maxwell YES)
    SET(WITH_navier-stokes YES)
    SET(WITH_nernst-planck YES)
    SET(WITH_neutronics YES)
    SET(WITH_richards YES)
    SET(WITH_schroedinger YES)
    SET(WITH_wave-equation YES)
    SET(WITH_miscellaneous YES)
  SET(WITH_2d-benchmarks-general YES)
  SET(WITH_2d-benchmarks-nist YES)

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
	
	find_package(PTHREAD REQUIRED)

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
  
  
  
  
  
  