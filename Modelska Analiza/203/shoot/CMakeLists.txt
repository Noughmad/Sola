project(shoot)

add_executable(shoot main.cpp)
target_link_libraries(shoot gsl gslcblas)

