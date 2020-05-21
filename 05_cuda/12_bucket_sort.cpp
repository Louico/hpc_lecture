#include <cstdio>
#include <cstdlib>
#include <vector>


__global__ void bucketSort(int* key, int * bucket, int n, int range){

    //__deice__ __managed__ int
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    bucket[i] = 0;
    //extern __shared__ std::vector<int> bucket(range);
    //for(int j=1; j<range; j++){
    //    bucket[j] = 0;
    //}



    for(int j=0; j<n; j++){
        if(key[j]==i){
            bucket[i]++;
        }
    }

    __syncthreads();
    int k = 0;
    for(int t = 0;t<i;t++){
        k+=bucket[t];
    }


    for(;bucket[i]>0;bucket[i]--){
        key[k++]=i;

    }

    __syncthreads();
};
int main() {
  int n = 50;
  int range = 5;
  //std::vector<int> key(n); no stl in device
  //std::vector<int> bucket(range);
  int  *key,*bucket;

  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");



  bucketSort<<<1,range>>>(key, bucket, n,range);
  cudaDeviceSynchronize();

/*
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }

  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
  */




  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(bucket);
  cudaFree(key);
}
