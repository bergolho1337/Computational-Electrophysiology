set(SOURCE_FILES_GPU_UTILS
        gpu_utils.h gpu_utils.cpp)

CUDA_ADD_LIBRARY(gpu_utils STATIC ${SOURCE_FILES_GPU_UTILS})
