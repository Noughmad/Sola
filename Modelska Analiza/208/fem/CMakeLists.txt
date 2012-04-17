project(fem)

find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} /usr/include/arpack++)

add_definitions(-fpermissive)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
add_executable(fem galerkin.cpp delitev.cpp delitev_graf.cpp delitev_arpack.cpp delitev_lastne.cpp main.cpp)
target_link_libraries(fem ${QT_QTCORE_LIBRARIES} ${QT_QTGUI_LIBRARIES} gsl gslcblas spqr lapack blas cholmod amd klu colamd arpack parpack superlu /usr/lib/openmpi/libmpi_f77.so mgl)

