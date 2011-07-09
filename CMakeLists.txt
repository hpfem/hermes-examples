project(hermes_examples)
	cmake_minimum_required(VERSION 2.6)

	set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
	include(FindHermes)
	include(FindPTHREAD)
	
	# Where to look for the static libraries.
	set(HERMES_DIRECTORY "/usr/lib")
	set(HERMES_INCLUDE_PATHS "/usr/include")
	set(DEP_INCLUDE_PATHS "/usr/include")

	macro(SET_COMMON_TARGET_PROPERTIES TRGT HERMES_LIBRARY)
		if(HERMES_LIBRARY STREQUAL "HERMES2D")
			FIND_PATH(INCLUDE_DIRS hermes2d.h ${HERMES_INCLUDE_PATHS})
			FIND_HERMES(HERMES_LIBRARY HERMES_DIRECTORY)
		endif(HERMES_LIBRARY STREQUAL "HERMES2D")
			
		# Link to Hermes.
		target_link_libraries(${TRGT} ${HERMES})
	endmacro(SET_COMMON_TARGET_PROPERTIES)
		
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)

	find_package(PTHREAD REQUIRED)
	
	include_directories(${HERMES_INCLUDE_PATHS})
	include_directories(${DEP_INCLUDE_PATHS})
	
	add_subdirectory(benchmarks-general)
	add_subdirectory(benchmarks-nist)
	add_subdirectory(tutorials)
	add_subdirectory(advanced_examples)