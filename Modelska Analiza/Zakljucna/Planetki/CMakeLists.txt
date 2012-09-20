project(planetki)

find_package(Qt4 COMPONENTS QtCore QtGui QtSVG)

include_directories(${QT_INCLUDES})

add_executable(planetki
    planets.cpp
    main.cpp
    toolbox.cpp
    KEP/src/lambert_problem.cpp
)
target_link_libraries(planetki ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTSVG_LIBRARY} gsl gslcblas)
