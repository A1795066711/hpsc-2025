#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 16;
  float x[N], y[N], m[N], fx[N], fy[N] , mask[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    mask[i] = i;
  }
  __m512 M = _mm512_load_ps(mask);
  __m512 mvec = _mm512_load_ps(m);
  for(int i=0; i<N; i++) {
    __m512 xi = _mm512_set1_ps(i);
    __mmask16 mask = _mm512_cmp_ps_mask(xi, M, _MM_CMPINT_EQ);
    
    __m512 xi = _mm512_set1_ps(x[i]);
    __m512 yi = _mm512_set1_ps(y[i]);

    __m512 xvec = _mm512_load_ps(x);
    __m512 yvec = _mm512_load_ps(y);

    __m512 Fx = _mm512_load_ps(fx);
    __m512 Fy = _mm512_load_ps(fy);

    __m512 rx = _mm512_sub_ps( xi , xvec );
    __m512 ry = _mm512_sub_ps( yi , xvec );

    __m512 m_rx = _mm512_mul_ps( mvec , xvec );
    __m512 m_ry = _mm512_mul_ps( mvec , yvec );

    __m512 r_rx = _mm512_mul_ps( rx , rx );
    __m512 r_ry = _mm512_mul_ps( ry , ry );
    __m512 r = _mm512_add_ps( rx , ry );
    __m512 r = _mm512_rsqrt14_ps( r );

    __m512 r3 = _mm512_mul_ps( r , r );
    __m512 r3 = _mm512_mul_ps( r3 , r );

    __m512 r_Fx = _mm512_mul_ps(Fx , mvec);
    __m512 r_Fy = _mm512_mul_ps(Fy , mvec);

    r_Fx = _mm512_div_ps(r_Fx , r3);
    r_Fy = _mm512_div_ps(r_Fy , r3);
    r_Fx = _mm512_sub_ps(Fx , r_Fx);
    r_Fy = _mm512_sub_ps(Fy , r_Fy);

    Fx = _mm512_mask_blend_ps(mask, r_Fx, Fx);
    Fy = _mm512_mask_blend_ps(mask, r_Fy, Fy);
    _mm512_store_ps(fx, Fx);
    _mm512_store_ps(fy, Fy);

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
