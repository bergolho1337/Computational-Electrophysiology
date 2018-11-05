//
// Created by sachetto on 14/10/17.
//

#ifndef FENTON_CONFIG_HELPER_H
#define FENTON_CONFIG_HELPER_H

#include "config_parser.h"

#include <cstdbool>
#include <cstdint>
#include <dlfcn.h>

#include "../../monodomain/monodomain_solver.h"
#include "../../grid/grid.h"
#include "../../models_library/model_common.h"

void set_celular_model (struct monodomain_solver *solver, struct user_options *configs);
void set_control_volumes (struct monodomain_solver *solver, struct grid *the_grid);
void set_derivative (struct monodomain_solver *solver, struct grid *the_grid);
void set_velocity_points (struct monodomain_solver *solver, struct grid *the_grid);
void set_plot_points (struct monodomain_solver *solver);

#endif //MONOALG3D_CONFIG_COMMON_H
