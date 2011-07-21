project(1d-system)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
