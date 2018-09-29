SET(LIB_TYPE MODULE) #set the lib type

############## TEN TUSCHER 2006 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2006_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2006_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_3_RS.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(ten_tusscher_2006 PUBLIC utils)	
	ENDIF()
else()
    ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)

#########################################################

############## BONDARENKO ##############################
set(MODEL_FILE_CPU "bondarenko_2004.c")
set(MODEL_FILE_GPU "bondarenko_2004_GPU.cu")
set(COMMON_HEADERS model_common.h bondarenko_2004.h)


if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(bondarenko_2004 PUBLIC utils)
	ENDIF()
else()
    ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 3 ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_RS_GPU.cu")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3 ${LIB_TYPE} model_gpu_utils.h model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(ten_tusscher_3 PUBLIC utils)
	endif()
else()
    ADD_LIBRARY(ten_tusscher_3 ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## FABERRUDY 2000 ##############################
set(MODEL_FILE_CPU "luo_rudy_1991.c")
set(MODEL_FILE_GPU "luo_rudy_1991.cu")
set(COMMON_HEADERS model_common.h luo_rudy_1991.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(luo_rudy_1991 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(luo_rudy_1991 PUBLIC utils)
	endif()
else()
    ADD_LIBRARY(luo_rudy_1991 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## FITZHUGH NAGUMO 1961 ##############################
set(MODEL_FILE_CPU "fitzhugh_1961.c")
set(MODEL_FILE_GPU "fitzhugh_1961.cu")
set(COMMON_HEADERS model_common.h fitzhugh_1961.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(fitzhugh_1961 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(fitzhugh_1961 PUBLIC utils)
	endif()
else()
    ADD_LIBRARY(fitzhugh_1961 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## NOBLE 1962 ##############################
set(MODEL_FILE_CPU "noble_1962.c")
set(MODEL_FILE_GPU "noble_1962.cu")
set(COMMON_HEADERS model_common.h noble_1962.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(noble_1962 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(noble_1962 PUBLIC utils)
	endif()
else()
    ADD_LIBRARY(noble_1962 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################