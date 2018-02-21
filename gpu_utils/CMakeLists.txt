
cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (gpu_utils)

find_package(CUDA)

set (CUDA_VERBOSE_BUILD OFF)
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

CUDA_ADD_LIBRARY(MWAC_GPU_Utils SHARED gpu_utils.cu reduce_kernel.cu)
if (GALAXY)
    #on a CRAY - probably Galaxy
    set (CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 ")
    set_target_properties(MWAC_GPU_Utils PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA ${CRAY_CUDA_TOOLKIT_INCLUDE_OPTS}")
else ()
    CUDA_INCLUDE_DIRECTORIES (${CUDA_INCLUDE_DIR} "/usr/local/cuda/samples/common/inc/")

    set (CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_52,code=compute_52 ")

endif ()

CUDA_ADD_CUFFT_TO_TARGET(MWAC_GPU_Utils)

#install(FILES gpu_utils.h reduce.h cuda_error.h DESTINATION include)
#install(TARGETS MWAC_GPU_Utils LIBRARY DESTINATION lib)

