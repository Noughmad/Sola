project(fem)

find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} /usr/include/arpack++)

add_definitions(-fpermissive)

option(WITH_ARPACK "Use ARPACK for finding eigenvalues")
option(WITH_MATHGL "Use MathGL for plotting")
option(WITH_GSL "Use the GNU Scientific Library" ON)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

set (LIBS
    ${QT_QTCORE_LIBRARIES} 
    ${QT_QTGUI_LIBRARIES} 
    spqr 
    lapack 
    blas
    cholmod
    amd
    colamd
)

if (WITH_ARPACK)
    add_definitions(-DWITH_ARPACK)
    set (LIBS ${LIBS}
        arpack
        parpack
        superlu
        /usr/lib/openmpi/libmpi_f77.so
    )
endif()

if (WITH_MATHGL)
    add_definitions(-DWITH_MATHGL)
    set (LIBS ${LIBS}
        mgl
    )
endif()

if (WITH_GSL)
    add_definitions(-DWITH_GSL)
    set (LIBS ${LIBS}
        gsl
        gslcblas
    )
endif()

add_executable(fem galerkin.cpp delitev.cpp delitev_graf.cpp delitev_arpack.cpp delitev_lastne.cpp main.cpp)
target_link_libraries(fem ${LIBS})

