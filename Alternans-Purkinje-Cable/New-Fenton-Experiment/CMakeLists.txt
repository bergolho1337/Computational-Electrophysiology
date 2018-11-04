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

IF (CMAKE_VERSION VERSION_LESS "3.1")
    IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    ENDIF(CMAKE_CXX_COMPILER_ID)
ELSE()
    set(CMAKE_CXX_STANDARD 11)
ENDIF()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

FIND_PACKAGE(CUDA)
IF (CUDA_FOUND)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")
    ADD_SUBDIRECTORY(src/gpu_utils)
    SET(OPT_DEPS gpu_utils)
ENDIF(CUDA_FOUND)

#ADD_SUBDIRECTORY(src/utils)
#ADD_SUBDIRECTORY(src/alg)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
#ADD_SUBDIRECTORY(src/string)
ADD_SUBDIRECTORY(src/grid)
ADD_SUBDIRECTORY(src/purkinje)

ADD_EXECUTABLE(FentonExperiment src/main.cpp)

#ADD_DEPENDENCIES(FentonExperiment utils alg solvers hashes ini_parser string graph ${OPT_DEPS})
ADD_DEPENDENCIES(FentonExperiment solvers ini_parser hashes grid purkinje ${OPT_DEPS})

#TARGET_LINK_LIBRARIES(FentonExperiment alg solvers utils hashes ini_parser string graph ${OPT_DEPS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CUDA_LIBRARIES} dl m)
TARGET_LINK_LIBRARIES(FentonExperiment solvers ini_parser hashes grid purkinje ${OPT_DEPS} ${CUDA_LIBRARIES} dl m)

ADD_SUBDIRECTORY(src/models_library)

#SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.h)

#ADD_SUBDIRECTORY(src/stimuli_library)
#ADD_SUBDIRECTORY(src/domains_library)
#ADD_SUBDIRECTORY(src/purkinje_library)
#ADD_SUBDIRECTORY(src/matrix_assembly_library)
#ADD_SUBDIRECTORY(src/linear_system_solver_library)
#ADD_SUBDIRECTORY(src/extra_data_library)
###

#SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

#ADD_SUBDIRECTORY(src/tests)
