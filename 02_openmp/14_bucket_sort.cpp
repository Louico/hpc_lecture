#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range),index(range);
#pragma omp parallel for
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
#pragma omp parallel for shared(bucket)
  for (int i=0; i<n; i++) {
#pragma omp atomic update
    bucket[key[i]]++;
  }
/*
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
*/

#pragma omp parallel for
for (int i=0; i<range; i++){
    index[i] = 0;
#pragma omp parallel for
    for (int j=0; j<i; j++){
#pragma omp atomic update
        index[i] += bucket[j];
    }
#pragma omp parallel for
    for(int k = index[i]; k<index[i]+bucket[i]; k++){
        key[k] = i;
    }

}

/*
#pragma omp parallel for
    for (int i=0; i<range; i++){
        int index = 0;
#pragma omp parallel for reduction(+:index)
        for (int j=0; j<i; j++){
            index += bucket[j];
        }
#pragma omp parallel for
        for(int k = bucket[i]; k>0; k--){
            key[index++] = i;
        }

    }
*/
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
