set(SOURCE_FILES_HASH hash_common.h
        point_hash.cpp point_hash.h
        string_hash.cpp string_hash.h
       )

ADD_LIBRARY(hashes STATIC ${SOURCE_FILES_HASH})

if(NOT MSVC)
    set_property(TARGET hashes PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()