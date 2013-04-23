macro(SET_COMMON_TARGET_PROPERTIES TRGT HERMES_VERSION)
	
	set(HERMES_VERSION ${HERMES_VERSION})
		
  find_package(HERMES REQUIRED)

  if (COMPILE_FLAGS)
    set_target_properties(${TRGT} PROPERTIES COMPILE_FLAGS ${COMPILE_FLAGS})
  endif (COMPILE_FLAGS)
  
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  if(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  endif(MSVC)

	target_link_libraries(${TRGT} ${HERMES_COMMON_LIBRARY})
	target_link_libraries(${TRGT} ${HERMES_LIBRARY})
	# Is empty if WITH_TRILINOS = NO
	target_link_libraries(${TRGT} ${TRILINOS_LIBRARIES})
			
endmacro(SET_COMMON_TARGET_PROPERTIES)
