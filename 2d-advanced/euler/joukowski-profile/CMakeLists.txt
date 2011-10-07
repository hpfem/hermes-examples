project(joukowski-profile)

add_executable(${PROJECT_NAME} main.cpp ../euler_util.cpp ../numerical_flux.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")