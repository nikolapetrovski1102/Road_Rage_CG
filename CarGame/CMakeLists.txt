cmake_minimum_required(VERSION 3.24)
project(untitled51)

set(CMAKE_CXX_STANDARD 17)

# Set the path to TDM-GCC
set(TDM_GCC_PATH "C:/TDM-GCC-64")

# Set the path to GLUT
set(GLUT_INCLUDE_DIR "${TDM_GCC_PATH}/x86_64-w64-mingw32/include/GL")

# Link GLUT, OpenGL, and GLU libraries
set(GLUT_LIBRARY "${TDM_GCC_PATH}/x86_64-w64-mingw32/lib/libfreeglut.a")
set(OPENGL_LIBRARY "${TDM_GCC_PATH}/x86_64-w64-mingw32/lib/libopengl32.a")
set(GLU_LIBRARY "${TDM_GCC_PATH}/x86_64-w64-mingw32/lib/libglu32.a")

# Set include directories
include_directories(${GLUT_INCLUDE_DIR})

# Add the executable
add_executable(untitled51 main.cpp BmpLoader.cpp)

# Link the libraries to the executable
target_link_libraries(untitled51 ${GLUT_LIBRARY} ${OPENGL_LIBRARY} ${GLU_LIBRARY})
