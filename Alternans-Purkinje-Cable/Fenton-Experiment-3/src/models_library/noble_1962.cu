#include "noble_1962.h"
#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) 
{

    print_to_stdout_and_file("Using noble_1962 GPU model\n");

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

         *((real * )((char *) sv + pitch * 0) + threadID) = -87.0f;    //V millivolt 
         *((real * )((char *) sv + pitch * 1) + threadID) = 0.01f;     //m dimensionless
         *((real * )((char *) sv + pitch * 2) + threadID) = 0.8f;      //h millivolt 
         *((real * )((char *) sv + pitch * 3) + threadID) = 0.01f;     //n dimensionless 
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
    const real V_old_ = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real m_old_ = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real h_old_ = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real n_old_ = *((real*)((char*)sv_ + pitch * 3) + threadID_);

    //___________________________________________________________________________
    //Parameters (miliseconds)
    const real Cm = 12.0f;                                 // (microF)
    const real g_na_max = 400.0f;                       // (microS)
    const real E_na = 40.0f;                               // (millivolt)
    const real g_L = 0.075f;                                // (microS)
    const real E_L = -60.0f;                               // (millivolt)

    real calc_I_stim = stim_current;

    // Algebraics
    real g_na = m_old_*m_old_*m_old_*h_old_*g_na_max;
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    real alpha_h =  ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    real i_na =  (g_na+1.4e-01)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    real beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    real beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    real beta_n =  ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    //real g_K1 =  1.3f*exp((- V_old_ - 90.0000)/50.0000)+ 0.015f*exp((V_old_+90.0000)/60.0000);
    real g_K1 = (((1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)))+(1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)))));
    real g_K2 =  1.2f*n_old_*n_old_*n_old_*n_old_;
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    rDY_[0] = ( - (i_na + i_k + i_leak + calc_I_stim)) / Cm;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  (beta_m*m_old_) );
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  (beta_h*h_old_) );
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  (beta_n*n_old_) );

}
