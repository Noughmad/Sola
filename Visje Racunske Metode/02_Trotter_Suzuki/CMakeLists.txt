project(trotter)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

add_executable(trotter main.cpp)
target_link_libraries(trotter gsl gslcblas)

