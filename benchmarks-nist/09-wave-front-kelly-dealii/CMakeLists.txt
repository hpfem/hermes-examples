project(nist-09-kelly-dealii) 
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set(COMPILE_FLAGS   "-O3")
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
