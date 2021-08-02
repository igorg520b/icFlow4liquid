cmake_minimum_required(VERSION 3.5)

project(fl4liquid LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Boost REQUIRED)

find_package(VTK 9.0 REQUIRED COMPONENTS
  CommonColor
  CommonCore
  FiltersSources
  FiltersGeometry
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  ViewsQt
  IOGeometry
  RenderingAnnotation
  RenderingLabel
  IOImage
  QUIET
)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)
add_executable(fl4liquid
    gui/main.cpp
    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    gui/icons.qrc

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

target_include_directories(fl4liquid PUBLIC ${HDF5_INCLUDE_DIRS})
target_include_directories(fl4liquid PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(fl4liquid PUBLIC gui geometry simulation gui/propertybrowser)
target_include_directories(fl4liquid PUBLIC gui/objectpropertybrowser gui/windows)
target_include_directories(fl4liquid PUBLIC gui/interactor)

target_link_libraries(fl4liquid PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
target_link_libraries(fl4liquid PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fl4liquid PUBLIC gmsh)
target_link_libraries(fl4liquid PUBLIC tbb)
target_link_libraries(fl4liquid PUBLIC gomp pthread m dl)
#target_link_libraries(fl4liquid PUBLIC ${HDF5_C_LIBRARIES})
target_link_libraries(fl4liquid PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(fl4liquid PUBLIC Eigen3::Eigen)
target_link_libraries(fl4liquid PUBLIC ${Boost_LIBRARIES})

#MOSEK
target_include_directories(fl4liquid PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/h")
target_link_directories(fl4liquid PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/bin")
target_link_libraries(fl4liquid PRIVATE mosek64)

target_compile_options(fl4liquid PUBLIC -fopenmp)


