
add_library(grate SHARED integrator.cpp)
target_link_libraries(grate gsl gslcblas ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})