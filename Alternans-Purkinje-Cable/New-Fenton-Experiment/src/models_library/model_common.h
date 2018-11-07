//
// Created by sachetto on 05/10/17.
//

// Every model need to implement the functions described in this model file in order to be loaded correctly from the
// edo solver. This models_library should compile without using any dependency of our codebase

#ifndef MONOALG3D_MODEL_COMMON_H
#define MONOALG3D_MODEL_COMMON_H

#include "../monodomain/constants.h"

#include <cstdio>
#include <cstdbool>
#include <cstdint>
#include <cstddef>

struct cell_model_data 
{
    int number_of_ode_equations;
    double initial_v;
    char *model_library_path;
};

struct control_volume
{
  double *y_old;
  double *y_star;
  double *y_new;
};

#define GET_CELL_MODEL_DATA(name) EXPORT_FN void name (struct cell_model_data *cell_model,\
                                                      struct control_volume *volumes, int num_volumes)
typedef GET_CELL_MODEL_DATA (get_cell_model_data_fn);

// CPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_CPU(name) EXPORT_FN void name (struct control_volume *volumes, int num_volumes)
typedef SET_ODE_INITIAL_CONDITIONS_CPU (set_ode_initial_conditions_cpu_fn);

#define SOLVE_MODEL_ODES_CPU(name)                                                                                     \
EXPORT_FN void name (double dt, double *stim_currents, int num_volumes,\
                     struct control_volume *volumes)
typedef SOLVE_MODEL_ODES_CPU (solve_model_ode_cpu_fn);


#endif // MONOALG3D_MODEL_COMMON_H
