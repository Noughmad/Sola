project(fem)

find_package(Qt4 REQUIRED)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
add_executable(fem main.cpp)
target_link_libraries(fem ${QT_QTCORE_LIBRARIES} ${QT_QTGUI_LIBRARIES} gsl gslcblas spqr lapack blas cholmod amd colamd)

