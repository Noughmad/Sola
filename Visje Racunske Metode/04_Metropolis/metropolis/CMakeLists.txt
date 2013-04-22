project(metropolis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_package(Qt5Core)
include_directories(${Qt5Core_INCLUDE_DIRS})
set(CMAKE_POSITION_INDEPENDENT_CODE ${Qt5_POSITION_INDEPENDENT_CODE})

add_executable(metropolis oscilator.cpp main.cpp)
target_link_libraries(metropolis nutils gsl gslcblas ${Qt5Core_LIBRARIES})

install(TARGETS metropolis RUNTIME DESTINATION bin)
