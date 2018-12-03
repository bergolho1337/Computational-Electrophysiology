#ifndef FENTON_GPU_UTILS_H
#define FENTON_GPU_UTILS_H

#include "cuda_runtime.h"

#define check_cuda_errors(ans) { cuda_assert((ans), __FILE__, __LINE__); }
void cuda_assert(cudaError_t code, const char *file, int line);

#endif // FENTON_GPU_UTILS_H