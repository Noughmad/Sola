project(integrator)

add_executable(integrator main.cpp)
target_link_libraries(integrator gsl gslcblas)
