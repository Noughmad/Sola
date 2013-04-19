project(metropolis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(metropolis oscilator.cpp main.cpp)
target_link_libraries(metropolis nutils gsl gslcblas)

install(TARGETS metropolis RUNTIME DESTINATION bin)
