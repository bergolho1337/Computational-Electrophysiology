set(SOURCE_FILES_SOLVER monodomain_solver.h monodomain_solver.cpp
                        ode_solver.h ode_solver.cpp)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
