project(kopeli)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(kopeli veriga.cpp main.cpp)
target_link_libraries(kopeli gsl gslcblas)

install(TARGETS kopeli RUNTIME DESTINATION bin)
