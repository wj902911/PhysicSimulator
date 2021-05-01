cmake_minimum_required(VERSION 3.1)
project(simulator)

set (CMAKE_CXX_STANDARD 17)
#set (CMAKE_POSITION_INDEPENDENT_CODE ON)
#set (CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
#option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
add_compile_options(/bigobj)

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include(cmake/common_functions.cmake)

find_package(LIBIGL REQUIRED QUIET)
#find_package(Threads REQUIRED)
#find_package(Eigen3 REQUIRED NO_MODULE)

#add_library(psim_std_dep INTERFACE)
#target_link_libraries(psim_std_dep INTERFACE Eigen3::Eigen Threads::Threads)
#target_include_directories(psim_std_dep SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}/libigl/include/)

add_subdirectory(lib)

#igl_set_folders()

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui psim)
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_SOURCE_DIR}/lib")
