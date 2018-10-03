cmake_minimum_required(VERSION 3.0.1)
project(FentonExperiment)

if (CMAKE_VERSION VERSION_LESS "3.9.0")
    # Set only because there are old target_link_libraries() statements in the
    # FindCUDA.cmake file that comes with CMake that error out if the new behavior
    # is used.  In newer versions of CMake we can instead set CUDA_LINK_LIBRARIES_KEYWORD which fixes this issue.
    cmake_policy(SET CMP0023 OLD)
else()
    set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
endif()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

IF (CMAKE_VERSION VERSION_LESS "3.1")
    IF (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        SET(CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
    ENDIF(CMAKE_C_COMPILER_ID)
ELSE()
    SET(CMAKE_C_STANDARD 99)
    SET (CMAKE_CXX_STANDARD 11)
ENDIF()

FIND_PACKAGE(OpenMP)

IF (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

FIND_PACKAGE(CUDA)

IF (CUDA_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")

    ADD_SUBDIRECTORY(src/gpu_utils)
    SET(OPT_DEPS gpu_utils)
ENDIF(CUDA_FOUND)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/grid)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/purkinje)

ADD_EXECUTABLE(FentonExperiment src/main.cpp)

ADD_DEPENDENCIES(FentonExperiment utils grid solvers hashes ini_parser graph purkinje ${OPT_DEPS})

TARGET_LINK_LIBRARIES(FentonExperiment grid solvers utils hashes ini_parser graph purkinje ${OPT_DEPS} Eigen3::Eigen ${CUDA_LIBRARIES} dl m)

ADD_SUBDIRECTORY(src/models_library)

SET(COMMON_SRC ../libraries_common/config_helpers.cpp ../libraries_common/config_helpers.h)

###
