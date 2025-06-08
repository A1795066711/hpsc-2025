#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N] , mask[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    mask[i] = i;
  }
  __m256 M = _mm256_load_ps(mask);
  __m256 mvec = _mm256_load_ps(m);
  for(int i=0; i<N; i++) {
    __m256 xi = _mm256_set1_ps(i);
    __mmask8 mask = _mm256_cmp_ps_mask(xi, M, _MM_CMPINT_EQ);
    
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);

    __m256 xvec = _mm256_load_ps(x);
    __m256 yvec = _mm256_load_ps(y);

    __m256 Fx = _mm256_load_ps(fx);
    __m256 Fy = _mm256_load_ps(fy);

    __m256 rx = _mm256_sub_ps( xi , xvec );
    __m256 ry = _mm256_sub_ps( yi , xvec );

    __m256 m_rx = _mm256_mul_ps( mvec , xvec );
    __m256 m_ry = _mm256_mul_ps( mvec , yvec );

    __m256 r_rx = _mm256_mul_ps( rx , rx );
    __m256 r_ry = _mm256_mul_ps( ry , ry );
    __m256 r = _mm256_add_ps( rx , ry );

    __m256 r3 = _mm256_mul_ps( r , r );
    __m256 r3 = _mm256_mul_ps( r3 , r );

    __m256 r_Fx = _mm256_mul_ps(Fx , mvec);
    __m256 r_Fy = _mm256_mul_ps(Fy , mvec);

    r_Fx = _mm256_div_ps(r_Fx , r3);
    r_Fy = _mm256_div_ps(r_Fy , r3);
    r_Fx = _mm256_sub_ps(Fx , r_Fx);
    r_Fy = _mm256_sub_ps(Fy , r_Fy);

    Fx = _mm512_mask_blend_ps(mask, r_Fx, Fx);
    Fy = _mm512_mask_blend_ps(mask, r_Fy, Fy);
    _mm512_store_ps(fx, Fx);
    _mm512_store_ps(fy, Fy);

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
