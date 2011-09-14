project(test-example-euler-gamm-explicit)

add_executable(${PROJECT_NAME} main.cpp ../../euler_util.cpp ../../numerical_flux.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D") 
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-euler-gamm-explicit ${BIN})
set_tests_properties(test-example-euler-gamm-explicit PROPERTIES LABELS slow)

