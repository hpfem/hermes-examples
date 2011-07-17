project(hermes_examples)
	cmake_minimum_required(VERSION 2.6)

	set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
        include(CommonTargetProperties)
	
	# Where to look for the static libraries.
	set(HERMES_DIRECTORY usr/local/lib)
	set(HERMES_INCLUDE_PATH usr/local/include)
	set(DEP_INCLUDE_PATHS usr/local/include)

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
	
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)

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

	find_package(HERMES_COMMON REQUIRED)
	
	set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
	set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

	include_directories(${HERMES_COMMON_INCLUDE_PATH})
	include_directories(${HERMES2D_INCLUDE_PATH})
	include_directories(${DEP_INCLUDE_PATHS})

	add_subdirectory(benchmarks-general)
	add_subdirectory(benchmarks-nist)
	add_subdirectory(tutorials)
	add_subdirectory(advanced_examples)
