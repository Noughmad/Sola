project(vrtavka)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)
find_package(Qt4 REQUIRED)


add_executable(vrtavka main.cpp)
target_link_libraries(vrtavka gsl gslcblas)

