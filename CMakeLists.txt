cmake_minimum_required(VERSION 3.14)

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
find_package(spdlog REQUIRED)

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

    geometry/elements/baseelement.h
    geometry/elements/element.cpp
    geometry/elements/element.h
    geometry/elements/edge.cpp
    geometry/elements/edge.h
    geometry/elements/interaction.h
    geometry/elements/interaction.cpp


    geometry/meshfragment.cpp
    geometry/meshfragment.h
    geometry/mesh.cpp
    geometry/mesh.h
    geometry/elements/boundaryedge.h
    geometry/elements/cohesivezone.cpp
    geometry/elements/cohesivezone.h

    simulation/parameters_sim.h
    simulation/model.cpp
    simulation/model.h

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

target_include_directories(fl4liquid PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(fl4liquid PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(fl4liquid PRIVATE gui geometry geometry/elements simulation gui/propertybrowser)
target_include_directories(fl4liquid PRIVATE gui/objectpropertybrowser gui/windows)
target_include_directories(fl4liquid PRIVATE gui/interactor)

target_link_libraries(fl4liquid PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
target_link_libraries(fl4liquid PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fl4liquid PRIVATE gmsh)
target_link_libraries(fl4liquid PRIVATE tbb)
target_link_libraries(fl4liquid PRIVATE gomp m dl)
#target_link_libraries(fl4liquid PRIVATE pthread)
target_link_libraries(fl4liquid PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(fl4liquid PRIVATE Eigen3::Eigen)
target_link_libraries(fl4liquid PRIVATE ${Boost_LIBRARIES})
target_link_libraries(fl4liquid PRIVATE spdlog::spdlog)

#MOSEK
target_include_directories(fl4liquid PRIVATE "~/mosek/9.3/tools/platform/linux64x86/h")
target_link_directories(fl4liquid PRIVATE "~/mosek/9.3/tools/platform/linux64x86/bin")
target_link_libraries(fl4liquid PRIVATE mosek64)

target_compile_options(fl4liquid PRIVATE -fopenmp -Wno-stringop-overflow)
#target_compile_options(fl4liquid PRIVATE "$<$<CONFIG:RELEASE>:-fopenmp>")
