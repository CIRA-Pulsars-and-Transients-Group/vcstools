
cmake_minimum_required (VERSION 2.6)
project (gpu_utils)

find_package(CUDA)
set (CUDA_VERBOSE_BUILD OFF)
set (CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_52,code=compute_52 ")
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

CUDA_INCLUDE_DIRECTORIES (${CUDA_INCLUDE_DIR} "/usr/local/cuda/samples/common/inc/")
CUDA_ADD_LIBRARY(MWA_GPUtils SHARED gpu_utils.cu reduce_kernel.cu)
CUDA_ADD_CUFFT_TO_TARGET(MWA_GPUtils)

install(FILES gpu_utils.h reduce.h cuda_error.h DESTINATION include)
install(TARGETS MWA_GPUtils LIBRARY DESTINATION lib)

