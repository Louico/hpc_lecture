#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <chrono>
#include <immintrin.h>
#include "timers.h"

using namespace std;
typedef vector<vector<float>> matrix;

void init_block(float *Ac, int mc, int nc) {
  for (int i=0; i<mc; i++)
    for (int j=0; j<nc; j++)
      Ac[i*nc+j] = 0;
}

void load_block(float *Ac, const matrix &A, int mc, int nc, int ic, int jc) {
  for (int i=0; i<mc; i++)
    for (int j=0; j<nc; j++)
      Ac[i*nc+j] = A[i+ic][j+jc];
}

void store_block(float *Ac, matrix &A, int mc, int nc, int ic, int jc) {
  for (int i=0; i<mc; i++)
    for (int j=0; j<nc; j++)
      A[i+ic][j+jc] += Ac[i*nc+j];
}

void matmult(matrix &A, matrix &B, matrix &C, int N) {
  const int m = N, n = N, k = N;
  const int kc = 512;
  const int nc = 64;
  const int mc = 256;
  const int nr = 64;
  const int mr = 32;
  float Cc[mc*nc];
#pragma omp parallel for collapse(2) private(Cc)
  for (int jc=0; jc<n; jc+=nc) {
    for (int pc=0; pc<k; pc+=kc) {
      float Bc[kc*nc];
      load_block(Bc,B,kc,nc,pc,jc);
      for (int ic=0; ic<m; ic+=mc) {
	float Ac[mc*kc];
	load_block(Ac,A,mc,kc,ic,pc);
	init_block(Cc,mc,nc);
        for (int jr=0; jr<nc; jr+=nr) {
          for (int ir=0; ir<mc; ir+=mr) {
            for (int kr=0; kr<kc; kr++) {
              for (int i=ir; i<ir+mr; i++) {
		__m256 Avec = _mm256_broadcast_ss(Ac+i*kc+kr);
                for (int j=jr; j<jr+nr; j+=8) {
                  __m256 Bvec = _mm256_load_ps(Bc+kr*nc+j);
                  __m256 Cvec = _mm256_load_ps(Cc+i*nc+j);
                  Cvec = _mm256_fmadd_ps(Avec, Bvec, Cvec);
                  _mm256_store_ps(Cc+i*nc+j, Cvec);
		}
              }
            }
          }
        }
	store_block(Cc,C,mc,nc,ic,jc);
      }
    }
  }
}

int main(int argc, char **argv) {
  const int N = 4096;
  matrix A(N,vector<float>(N));
  matrix B(N,vector<float>(N));
  matrix C(N,vector<float>(N));
  matmult(A,B,C,N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      A[i][j] = drand48();
      B[i][j] = drand48();
      C[i][j] = 0;
    }
  }
  startTimer();
  matmult(A,B,C,N);
  stopTimer();
  double time = getTime();
  printf("N=%d: %lf s (%lf GFlops)\n",N,time,2.*N*N*N/time/1e9);
#pragma omp parallel for
  for (int i=0; i<N; i++)
    for (int k=0; k<N; k++)
      for (int j=0; j<N; j++)
        C[i][j] -= A[i][k] * B[k][j];
  double err = 0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      err += fabs(C[i][j]);
  printf("error: %lf\n",err/N/N);
}
