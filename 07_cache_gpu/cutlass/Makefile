sgemm : gemm.cu
	nvcc -DTEST_SGEMM -arch=sm_60 gemm.cu -O3 -std=c++11 -I. -lcublas

clean :
	rm -f *.o *.out
