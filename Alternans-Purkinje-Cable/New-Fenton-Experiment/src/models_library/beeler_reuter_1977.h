#ifndef MONOALG3D_MODEL_BEELER_REUTER_1977_H
#define MONOALG3D_MODEL_BEELER_REUTER_1977_H

#include <cstdint>
#include "model_common.h"

#define NEQ 8
#define INITIAL_V (-84.624f)

#ifdef __CUDACC__

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(double *sv, int num_volumes);

__global__ void solve_gpu(double dt, double *sv, double* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(double *sv_, double *rDY_, double stim_current, int threadID_);

#endif

void solve_model_ode_cpu(double dt, struct control_volume &volume,\
                         double stim_current);
void RHS_cpu(double *rDY_, const double *y_old, const double *y_star, double stim_current);

#endif // MONOALG3D_MODEL_BEELER_REUTER_1977_H

