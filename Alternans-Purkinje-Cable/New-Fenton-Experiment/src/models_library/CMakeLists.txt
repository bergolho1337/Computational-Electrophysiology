SET(LIB_TYPE MODULE) #set the lib type

############## FITZHUGH NAGUMO 1961 ##############################
set(MODEL_FILE_CPU "fitzhugh_1961.cpp")
set(COMMON_HEADERS model_common.h fitzhugh_1961.h)

ADD_LIBRARY(fitzhugh_1961 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
#########################################################

############## BEELER REUTER 1977 ##############################
set(MODEL_FILE_CPU "beeler_reuter_1977.cpp")
set(COMMON_HEADERS model_common.h beeler_reuter_1977.h)

ADD_LIBRARY(beeler_reuter_1977 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
#########################################################

