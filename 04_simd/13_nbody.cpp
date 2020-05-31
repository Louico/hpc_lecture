#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], j[N], fxj[N], fyj[N];
  srand(1);
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }

  __m256 zero = _mm256_setzero_ps();
  for(int i=0; i<N; i+=8) {
    __m256 xi = _mm256_load_ps(x+i);
    __m256 yi = _mm256_load_ps(y+i);
    __m256 fxi = zero;
    __m256 fyi = zero;
    for(int j=0; j<N; j++) {
      __m256 dx = _mm256_set1_ps(x[j]);
      __m256 dy = _mm256_set1_ps(y[j]);
      __m256 mj = _mm256_set1_ps(m[j]);
      __m256 r2 = zero;
      dx = _mm256_sub_ps(xi, dx);
      dy = _mm256_sub_ps(yi, dy);
      r2 = _mm256_fmadd_ps(dx, dx, r2);
      r2 = _mm256_fmadd_ps(dy, dy, r2);
      __m256 mask = _mm256_cmp_ps(r2, zero, _CMP_GT_OQ);
      __m256 invR = _mm256_rsqrt_ps(r2);
      invR = _mm256_blendv_ps(zero, invR, mask);
      mj = _mm256_mul_ps(mj, invR);
      invR = _mm256_mul_ps(invR, invR);
      mj = _mm256_mul_ps(mj, invR);
      fxi = _mm256_fmadd_ps(dx, mj, fxi);
      fyi = _mm256_fmadd_ps(dy, mj, fyi);
    }
    _mm256_store_ps(fx+i, fxi);
    _mm256_store_ps(fy+i, fyi);
  }
  for(int i=0; i<N; i++)
    printf("%d %g %g\n",i,fx[i],fy[i]);

}
