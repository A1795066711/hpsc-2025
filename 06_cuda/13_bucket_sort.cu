#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void initialize_1(int * bucket){
	bucket[threadIdx.x] = 0;
}

__global__ void initialize_2(int * key , int * bucket , int * sum){
        atomicAdd(&bucket[key[threadIdx.x]] , 1);
	sum[threadIdx.x] = 0;
}

__global__ void reduce(int * sum , int * bucket , int n){
        int index = int( threadIdx.x / n );
        int index_ = threadIdx.x % n;
	if( index <= index_ ){
		atomicAdd(&sum[index_] , bucket[index]);
	}
	// printf( "[%i %i %i]" , threadIdx.x , index , index_ );
}

__global__ void sort(int * key , int * sum , int * bucket , int n){
	int index = threadIdx.x;
	int position;
	for( int i = 0 ; i < n ; i++ ){
		if( index < sum[i] ){
			position = i;
			break;
		}
	}
	key[index] = position;
}

int main() {
  int n = 50;
  int range = 5;
  int * key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int * sum;
  cudaMallocManaged(&sum, range*sizeof(int));

  int * bucket;
  cudaMallocManaged(&bucket, range*sizeof(int)); 
  initialize_1<<<1,range>>>(bucket);
  initialize_2<<<1,n>>>(key , bucket , sum);
  reduce<<<1,range*range>>>(sum , bucket , range);
  sort<<<1,n>>>( key , sum , bucket , range );
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
