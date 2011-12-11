project(metropolis)

add_executable(metropolis main.cpp)

target_link_libraries(metropolis m gsl gslcblas)