//
// Created by sachetto on 05/10/17.
//

// Every model need to implement the functions described in this model file in order to be loaded correctly from the
// edo solver. This models_library should compile without using any dependency of our codebase

#ifndef _MODEL_COMMON_H_
#define _MODEL_COMMON_H_

#include <cstdbool>
#include <cstdint>
#include <cstddef>
#include <cmath>

#define EXPORT_FN

#define GET_CELL_MODEL_DATA(name) EXPORT_FN void name (uint32_t *num_equations)
typedef GET_CELL_MODEL_DATA (get_cell_model_data_fn);

// CPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_CPU(name) EXPORT_FN void name (double *sv)
typedef SET_ODE_INITIAL_CONDITIONS_CPU (set_ode_initial_conditions_cpu_fn);

#define SOLVE_MODEL_ODES_CPU(name)                                                                                     \
EXPORT_FN void name (const double dt, double *sv, double *stim_currents, const uint32_t num_cells_to_solve)
typedef SOLVE_MODEL_ODES_CPU (solve_model_ode_cpu_fn);

#endif // MONOALG3D_MODEL_COMMON_H
