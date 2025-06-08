#include <cstdio>
#include <cstdlib>
#include <vector>

#include <omp.h>

void merge(std::vector<int>& vec, int begin, int mid, int end) {
  std::vector<int> tmp(end-begin+1);
  int left = begin;
  int right = mid+1; 
  for (int i=0; i<tmp.size(); i++) { 
    if (left > mid)
      tmp[i] = vec[right++];
    else if (right > end)
      tmp[i] = vec[left++];
    else if (vec[left] <= vec[right])
      tmp[i] = vec[left++];
    else
      tmp[i] = vec[right++]; 
  }
#pragma omp parrallel for
  for (int i=0; i<tmp.size(); i++) 
    vec[begin++] = tmp[i];
}

void merge_sort_parameters(std::vector<int>& vec, int begin, int end , std::vector<std::vector<int>> & result) {
  if(begin < end) {
    int mid = (begin + end) / 2;
    merge_sort_parameters(vec, begin, mid , result);
    merge_sort_parameters(vec, mid+1, end , result);
    std::vector <int> indices { begin, mid, end };
    result.push_back( indices );
  }
}

void merge_sort( std::vector<int>& vec, std::vector<std::vector<int>> & indices ) {
#pragma omp parrallel for 
  for ( int i = 0 ; i < indices.size() ; i++ ) {
     std::vector<int> parameters = indices[i];
     merge( vec , parameters[0] , parameters[1] , parameters[2] );
  }
}

int main() {
  int n = 20;
  std::vector<int> vec(n);
#pragma omp parrallel for
  for (int i=0; i<n; i++) {
    vec[i] = rand() % (10 * n);
    printf("%d ",vec[i]);
  }
  printf("\n");
  std::vector<std::vector<int>> parameters{};
  merge_sort_parameters(vec, 0, n-1 , parameters);
  merge_sort( vec, parameters );

  for (int i=0; i<n; i++) {
    printf("%d ",vec[i]);
  }
  printf("\n");
}
