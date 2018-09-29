#include "fitzhugh_1961.h"
#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) 
{

    print_to_stdout_and_file("Using fitzhugh_1961 GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;


    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));


    //the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

         *((real * )((char *) sv + pitch * 0) + threadID) = 0.000000f; //V millivolt 
         *((real * )((char *) sv + pitch * 1) + threadID) = 0.000000f; //h dimensionless 
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id);

            for(int i = 0; i < NEQ; i++) {
                *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
            }            

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_) {

    //State variables
    const real V_old_ =  *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real h_old_ =  *((real*)((char*)sv_ + pitch * 1) + threadID_);

    //Parameters
    const real alpha = -0.100000000000000e+00f;
    const real gamma = 3.000000000000000e+00f;
    const real epsilon = 5.000000000000000e-03f;

    real calc_I_stim = stim_current;

    rDY_[0] = (( V_old_*(V_old_ - alpha)*(1.00000 - V_old_) - h_old_) + calc_I_stim);
    rDY_[1] = epsilon*(V_old_ -  gamma*h_old_);

}
