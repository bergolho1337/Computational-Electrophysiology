set(SOURCE_FILES_SOLVER
        output_utils.h
        output_utils.c
        config/config_common.c config/config_common.h
        config/stim_config_hash.c config/stim_config_hash.h
        config/stim_config.h config/stim_config.c
        config/purkinje_config.c config/purkinje_config.h
        config/extra_data_config.c config/extra_data_config.h
        ode_solver.c ode_solver.h
        monodomain_solver.h monodomain_solver.c
        config/config_parser.c config/config_parser.h
        constants.h)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
