macro(SET_COMMON_TARGET_PROPERTIES TRGT HERMES_VERSION)
	set(HERMES_VERSION ${HERMES_VERSION})
  find_package(HERMES REQUIRED)

	target_link_libraries(${TRGT} ${HERMES_COMMON_LIBRARY})
	target_link_libraries(${TRGT} ${HERMES_LIBRARY} ${MATIO_LIBRARY} ${BSON_LIBRARY})
  target_link_libraries(${TRGT} ${TESTING_CORE_LIBRARY})

	# Is empty if WITH_TRILINOS = NO
	target_link_libraries(${TRGT} ${TRILINOS_LIBRARIES})
endmacro(SET_COMMON_TARGET_PROPERTIES)
