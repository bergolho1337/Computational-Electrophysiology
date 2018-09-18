SET(LIB_TYPE MODULE) #set the lib type

############## MITCHELL SHAEFFER 2002 ##############################
set(MODEL_FILE_CPU "mitchell_shaeffer.c")
set(COMMON_HEADERS model_common.h mitchell_shaeffer.h)

ADD_LIBRARY(mitchell_shaeffer ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

