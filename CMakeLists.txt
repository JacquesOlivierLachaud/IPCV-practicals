CMAKE_MINIMUM_REQUIRED(VERSION 3.20)
project(IPCV-tutorials)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/recipes)
# C++20 (instead of C++17) is necessary for windows environment.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(dgtal)
include(polyscope)

include_directories(${DGTAL_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})

add_executable(helloworld helloworld.cpp)
target_link_libraries(helloworld DGtal::DGtal)

add_executable(exampleBoard2D exampleBoard2D.cpp)
target_link_libraries(exampleBoard2D DGtal::DGtal)

add_executable(examplePolyscope examplePolyscope.cpp)
target_link_libraries(examplePolyscope DGtal::DGtal polyscope)

add_executable(volViewer volViewer.cpp)
target_link_libraries(volViewer DGtal::DGtal polyscope)



## Practicals

add_executable(vol-filter practical-filtering/vol-filter.cpp)
target_link_libraries(vol-filter DGtal::DGtal polyscope)

add_executable(homotopic-thinning practical-homotopic-thinning/homotopic-thinning.cpp)
target_link_libraries(homotopic-thinning DGtal::DGtal polyscope)

add_executable(scaleaxis practical-scaleaxis/scaleaxis.cpp)
target_link_libraries(scaleaxis DGtal::DGtal polyscope)

add_executable(2D-estimation practical-2D-estimation/2D-estimation.cpp)
target_link_libraries(2D-estimation DGtal::DGtal)

add_executable(3D-estimation practical-3D-estimation/3D-estimation.cpp)
target_link_libraries(3D-estimation DGtal::DGtal polyscope)

## Utils

add_executable(calculus utils/calculus.cpp)
target_link_libraries(calculus DGtal::DGtal polyscope)

add_executable(geodesics utils/geodesics.cpp)
target_link_libraries(geodesics DGtal::DGtal polyscope)


