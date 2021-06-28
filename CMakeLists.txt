cmake_minimum_required(VERSION 3.5)

project(fl4 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

#find_package(Boost 1.71.0)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkFiltersSources
  vtkFiltersGeometry
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsQt
  vtkIOGeometry
  vtkRenderingAnnotation
  vtkRenderingLabel
  vtkIOImage
  QUIET
)
if (NOT VTK_FOUND)
  message("${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)
add_executable(fl4
    #    gui/main_test.cpp
    gui/main.cpp
    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    gui/icons.qrc
#    gui/preferences_gui.h

    #property browser
    gui/propertybrowser/qtvariantproperty.h
    gui/propertybrowser/qtbuttonpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.h
    gui/propertybrowser/qteditorfactory.cpp
    gui/propertybrowser/qteditorfactory.h
    gui/propertybrowser/qtpropertybrowserutils.cpp
    gui/propertybrowser/qtvariantproperty.cpp
    gui/propertybrowser/qtpropertymanager.cpp
    gui/propertybrowser/qtpropertymanager.h
    gui/propertybrowser/qtpropertybrowserutils_p.h
    gui/propertybrowser/qtpropertybrowser.qrc
    gui/propertybrowser/qttreepropertybrowser.h
    gui/propertybrowser/qttreepropertybrowser.cpp
    gui/propertybrowser/qtpropertybrowser.h
    gui/propertybrowser/qtbuttonpropertybrowser.h
    gui/propertybrowser/qtpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.h

    gui/backgroundworker.cpp
    gui/backgroundworker.h

    simulation/modelcontrollerinterface.h

    geometry/node.cpp
    geometry/node.h

    geometry/element.cpp
    geometry/element.h

    geometry/meshfragment.cpp
    geometry/meshfragment.h
    geometry/mesh.cpp
    geometry/mesh.h

    simulation/parameters_sim.h
    simulation/model.cpp
    simulation/model.h
    simulation/interaction.h
    simulation/interaction.cpp

    simulation/equationofmotionsolver.cpp
    simulation/equationofmotionsolver.h

    gui/interactor/SpecialSelector2D.cpp
    gui/interactor/SpecialSelector2D.h

    bvh/bvhn.cpp
    bvh/bvhn.h
    bvh/kdop8.cpp
    bvh/kdop8.h
    bvh/ConcurrentPool.h
)


target_include_directories(fl4 PUBLIC gui geometry simulation gui/propertybrowser)
target_include_directories(fl4 PUBLIC gui/objectpropertybrowser gui/windows)
target_include_directories(fl4 PUBLIC gui/interactor)

# some paths need to be set manually
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/eigen)
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/vtk/Common/Core/)
#target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/gmsh-4.8.4-source/api)

target_link_directories(fl4 PUBLIC /usr/local/lib)
#target_link_directories(fl4 PUBLIC /home/s2/Projects-libraries/gmsh-build)

target_link_libraries(fl4 PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
target_link_libraries(fl4 PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fl4 PUBLIC gmsh)
target_link_libraries(fl4 PUBLIC tbb)
target_link_libraries(fl4 PUBLIC gomp pthread m dl)

#MOSEK
target_include_directories(fl4 PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/h")
target_link_directories(fl4 PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/bin")
target_link_libraries(fl4 PRIVATE mosek64)

target_compile_options(fl4 PUBLIC -fopenmp)


