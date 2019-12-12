CC=g++
NVCC=nvcc
CXXFLAGS= -O3 -Wextra -std=c++11 -lm 


CUDAFLAGS= -arch=sm_60  -std=c++11 -c
#CUDAFLAGS= -arch=compute_50 -code=sm_50 -std=c++11 -c

EXE=laserCAv6_1024

# Install in Ubuntu (Linux) "libhighgui-dev" and "libopencv-highgui-dev" packages for OpenCV.
# This is optional if it is not used the "#ifdef __VIDEO" option.
OPENCV = `pkg-config opencv --cflags --libs`
LIBS= -lcudart $(OPENCV)


all: $(EXE)

cuda_functions.o : cuda_functions.cu
	$(NVCC) $(CUDAFLAGS) cuda_functions.cu 

kernelV6.o : kernelV6.cu
	$(NVCC) $(CUDAFLAGS) kernelV6.cu

kernel_initv6.o : kernel_initv6.cu
	$(NVCC) $(CUDAFLAGS) kernel_initv6.cu

laserCAv6_1024: cuda_functions.o kernelV6.o kernel_initv6.o
	$(CC) -o laserCAv6_1024 main.cpp kernelV6.o cuda_functions.o kernel_initv6.o $(LIBS) $(CXXFLAGS)


clean:
	rm -rf laserCAv6_1024 *.o


# TO TAKE INTO ACCOUNT:
# 1) The real architecture must be an implementation of the virtual architecture
# 2) The virtual architecture should always be chosen as low as possible, thereby maximizing the actual GPUs to run on. 
# 3) The real architecture should be chosen as high as possible 
#
# GPU REAL Architecture Feature List:
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
# 
# 	                        Basic features
# sm_30 and sm_32           + Kepler support
#                           + Unified memory programming
# sm_35	                    + Dynamic parallelism support
# sm_50, sm_52, and sm_53	+ Maxwell support
# sm_60, sm_61, and sm_62	+ Pascal support
# sm_70 and sm_72	        + Volta support
# sm_75	                    + Turing support
#
#
# GPU VIRTUAL Architecture Feature List:
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
#
# 	                                        Basic features
# compute_30 and compute_32                 + Kepler support
#                                           + Unified memory programming
# compute_35	                            + Dynamic parallelism support
# compute_50, compute_52, and compute_53	+ Maxwell support
# compute_60, compute_61, and compute_62	+ Pascal support
# compute_70 and compute_72	                + Volta support
# compute_75	                            + Turing support

