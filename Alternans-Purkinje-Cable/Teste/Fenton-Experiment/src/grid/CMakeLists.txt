set(SOURCE_FILES_GRID grid.h grid.cpp
                      cell.h cell.cpp)

ADD_LIBRARY(grid STATIC ${SOURCE_FILES_GRID})

if(NOT MSVC)
    set_property(TARGET grid PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()