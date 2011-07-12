project(hermes_examples)
	cmake_minimum_required(VERSION 2.6)

	set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
	
	# Where to look for the static libraries.
	set(HERMES_DIRECTORY /lib)
	set(HERMES_INCLUDE_PATH "/include")
	set(DEP_INCLUDE_PATHS /include)

	macro(SET_COMMON_TARGET_PROPERTIES TRGT HERMES_VERSION)
		
		set(HERMES_VERSION ${HERMES_VERSION})
		
		find_package(HERMES_COMMON REQUIRED)
		find_package(HERMES REQUIRED)
		
		target_link_libraries(${TRGT} ${HERMES_COMMON_LIBRARY})
		target_link_libraries(${TRGT} ${HERMES_LIBRARY})
				
	endmacro(SET_COMMON_TARGET_PROPERTIES)
		
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)

	find_package(PTHREAD REQUIRED)
	
	set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
	set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

	include_directories(${HERMES_COMMON_INCLUDE_PATH})
	include_directories(${HERMES2D_INCLUDE_PATH})
	include_directories(${DEP_INCLUDE_PATHS})

	add_subdirectory(benchmarks-general)
	add_subdirectory(benchmarks-nist)
	add_subdirectory(tutorials)
	add_subdirectory(advanced_examples)
