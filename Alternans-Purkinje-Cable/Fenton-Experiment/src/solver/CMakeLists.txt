set(SOURCE_FILES_SOLVER solver.cpp solver.h
    config/config_parser.cpp config/config_parser.h
    config/purkinje_config.cpp config/purkinje_config.h
    config/stim_config.cpp config/stim_config.cpp)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
