project(strobo)
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)

add_executable(strobo main.cpp)
target_link_libraries(strobo ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})

