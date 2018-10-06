cmake_minimum_required(VERSION 3.0.1)
project(FentonExperiment)

# Set version of the C++ to be C++11
set (CMAKE_CXX_STANDARD 11)

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

# Set the compilation flags
#SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Find and link the OpenMP package
FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

# Find Eigen 3.3 package
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

# Compile the dependencies modules
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/config)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/purkinje)
ADD_SUBDIRECTORY(src/cell_models)

ADD_EXECUTABLE(FentonExperiment src/main.cpp)

ADD_DEPENDENCIES(FentonExperiment graph config solvers purkinje cellmodels ${OPT_DEPS})

TARGET_LINK_LIBRARIES(FentonExperiment graph config solvers purkinje cellmodels ${OPT_DEPS} dl m Eigen3::Eigen)

