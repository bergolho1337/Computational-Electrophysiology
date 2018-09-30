set(SOURCE_FILES_SOLVER
        config/config_common.cpp config/config_common.h
        config/stim_config_hash.cpp config/stim_config_hash.h
        config/stim_config.h config/stim_config.cpp
        config/purkinje_config.h config/purkinje_config.cpp
        ode_solver.cpp ode_solver.h
        output_utils.cpp output_utils.h
        monodomain_solver.h monodomain_solver.cpp
        config/config_parser.cpp config/config_parser.h
        constants.h)

#INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
