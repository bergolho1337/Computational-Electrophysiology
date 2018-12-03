set(SOURCE_FILES_SOLVER
        #output_utils.h output_utils.c
        config/config_common.cpp config/config_common.h
        config/config_helper.cpp config/config_helper.h
        config/stim_config_hash.cpp config/stim_config_hash.h
        config/stim_config.h config/stim_config.cpp
        #config/domain_config.c config/domain_config.h
        #config/assembly_matrix_config.c config/assembly_matrix_config.h
        #config/purkinje_config.c config/purkinje_config.h
        #config/extra_data_config.c config/extra_data_config.h
        #config/linear_system_solver_config.c config/linear_system_solver_config.h
        #ode_solver.c ode_solver.h
        monodomain_solver.h monodomain_solver.cpp
        config/config_parser.cpp config/config_parser.h
        constants.h
        timer.h
        )

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
