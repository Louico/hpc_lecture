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
  for(int k=0; k<N; k++) {
      j[k] = k;
  }
  for(int i=0; i<N; i++) {


      __m256 ivec = _mm256_set1_ps(i);
      __m256 jvec = _mm256_load_ps(j);
      __m256 mask = _mm256_cmp_ps(ivec, jvec, _CMP_NEQ_OQ);



      /*if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }*/
      __m256 xivec = _mm256_set1_ps(x[i]);
      __m256 xjvec = _mm256_load_ps(x);
      __m256 yivec = _mm256_set1_ps(y[i]);
      __m256 yjvec = _mm256_load_ps(y);
      __m256 rxvec,ryvec;
      rxvec = _mm256_sub_ps(xivec, xjvec);
      ryvec = _mm256_sub_ps(yivec, yjvec);



      //float rx = x[i] - x[j];
      //float ry = y[i] - y[j];

      __m256 rx2vec = _mm256_mul_ps(rxvec, rxvec);
      __m256 ry2vec = _mm256_mul_ps(ryvec, ryvec);


      __m256 rr3vec = _mm256_add_ps(rx2vec, ry2vec);
      rr3vec = _mm256_rsqrt_ps(rr3vec);

      __m256 rr2vec = _mm256_mul_ps(rr3vec, rr3vec);
      rr3vec = _mm256_mul_ps(rr2vec, rr3vec);
      //float r = std::sqrt(rx * rx + ry * ry);



      __m256 mvec = _mm256_load_ps(m);
      __m256 zerovec = _mm256_set1_ps(0);
      //mvec = _mm256_blendv_ps(zerovec, mvec, mask);
      //mvec = _mm256_mul_ps(mvec, mask);




      __m256 fxjvec,fyjvec;
      fxjvec = _mm256_mul_ps(rxvec, rr3vec);
      fxjvec = _mm256_mul_ps(fxjvec, mvec);
      fxjvec = _mm256_blendv_ps(zerovec, fxjvec, mask);
      fyjvec = _mm256_mul_ps(ryvec, rr3vec);
      fyjvec = _mm256_mul_ps(fyjvec, mvec);
      fyjvec = _mm256_blendv_ps(zerovec, fyjvec, mask);






      //__m256 fxvec = _mm256_load_ps(fx);
      //__m256 fyvec = _mm256_load_ps(fy);
      _mm256_store_ps(fxj,fxjvec);
      _mm256_store_ps(fyj,fyjvec);
      float fxjv = 0,fyjv = 0;
      for(int i=0; i<N ; i++) {
          fxjv+=fxj[i];
          fyjv+=fyj[i];
      }


      //fx[i] -= rx * m[j] / (r * r * r);
      //fy[i] -= ry * m[j] / (r * r * r);
      //fx[i] -= fxj[i];
      //fy[i] -= fyj[i];

      fx[i] -= fxjv;
      fy[i] -= fyjv;

      printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
