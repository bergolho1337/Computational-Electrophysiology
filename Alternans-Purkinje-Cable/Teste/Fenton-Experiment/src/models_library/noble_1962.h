#ifndef MONOALG3D_MODEL_NOBLE_1962_H
#define MONOALG3D_MODEL_NOBLE_1962_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 4
#define INITIAL_V (-75.5344986658f)

#ifdef __CUDACC__

extern "C" {
    #include "../utils/logfile_utils.h"
}

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(float *sv, int num_volumes);

__global__ void solve_gpu(float dt, float *sv, float* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(float *sv_, float *rDY_, float stim_current, int threadID_);

#else
#include "../utils/logfile_utils.h"
#endif


void solve_model_ode_cpu(float dt, float *sv, float stim_current);
void RHS_cpu(const float *sv, float *rDY_, float stim_current);

#endif //MONOALG3D_MODEL_NOBLE_1962_H
