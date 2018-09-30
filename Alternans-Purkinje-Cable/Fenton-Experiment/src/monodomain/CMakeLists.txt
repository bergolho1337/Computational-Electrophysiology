set(SOURCE_FILES_SOLVER
        monodomain_solver.h monodomain_solver.cpp
        ode_solver.h ode_solver.cpp
        config/config_parser.cpp config/config_parser.h
        config/stim_config.cpp config/stim_config.h
        config/stim_config_hash.cpp config/stim_config_hash.h)

#INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
