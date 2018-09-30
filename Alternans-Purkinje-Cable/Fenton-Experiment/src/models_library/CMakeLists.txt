SET(LIB_TYPE MODULE) #set the lib type

############## NOBLE 1962 ##############################
#set(MODEL_FILE_CPU "noble_1962.c")
#set(MODEL_FILE_GPU "noble_1962.cu")
#set(COMMON_HEADERS model_common.h noble_1962.h)

#if(CUDA_FOUND)
#    CUDA_ADD_LIBRARY(noble_1962 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
#	 ADD_LIBRARY(noble_1962 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
#endif(CUDA_FOUND)
#########################################################