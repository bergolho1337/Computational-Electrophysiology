set(SOURCE_FILES_ALG
        grid/grid.h
        grid/grid.c
        cell/cell.h
        cell/cell.c
        ../vector/stretchy_buffer.h
        )

ADD_LIBRARY(alg STATIC ${SOURCE_FILES_ALG})

if(NOT MSVC)
    set_property(TARGET alg PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
