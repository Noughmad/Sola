project(bem)

find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)


add_executable(bem main.cpp)
target_link_libraries(bem ${QT_QTCORE_LIBRARIES} ${QT_QTGUI_LIBRARIES} gsl gslcblas)
