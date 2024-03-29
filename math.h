#define NAN (0.0f / 0.0f)
#define INFINITY 1e5000f
#define HUGE_VALF INFINITY

#define SIGN_MASK 0x80000000u
#define EXPONENT_MASK 0x7f800000u
#define MANTISSA_MASK 0x007fffffu

static float utf(unsigned int i) {
  union {
    unsigned int i;
    float f;
  } u;
  u.i = i;
  return u.f;
}

static unsigned int ftu(float f) {
  union {
    unsigned int i;
    float f;
  } u;
  u.f = f;
  return u.i;
}

float fabsf(float x) { return utf(ftu(x) & 0x7fffffffu); }

#define FP_INFINITE 0
#define FP_NAN 1
#define FP_NORMAL 2
#define FP_SUBNORMAL 3
#define FP_ZERO 4

int fpclassify(float x) {
  unsigned int i = ftu(fabsf(x));
  return x == 0x7f800000u  ? FP_INFINITE
         : x > 0x7f800000u ? FP_NAN
         : x < 0x00800000u ? FP_SUBNORMAL
         : x == 0          ? FP_ZERO
                           : FP_NORMAL;
}

#define isinf(x) (ftu(fabsf(x)) == 0x7f800000u) 

#define isnan(x) (ftu(fabsf(x)) > 0x7f800000u)

int isnormal(float x) { return (ftu(x) & EXPONENT_MASK) > 0; }

int isfinite(float x) { return ftu(fabsf(x)) < 0x7f800000u; }

#define signbit(x) (ftu(x) >> 31)

#define isunordered(x, y) (isnan(x) || isnan(y) ? 1.f : 0.f)

#define isless(x, y) ((x) < (y) ? 1.f : 0.f)
#define islessequal(x, y) ((x) <= (y) ? 1.f : 0.f)
#define islessgreater(x, y) ((x) != (y) ? 1.f : 0.f)
#define isgreater(x, y) ((x) > (y) ? 1.f : 0.f)
#define isgreaterequal(x, y) ((x) >= (y) ? 1.f : 0.f)

// x * sign(y)
static float mulsign(float x, float y) {
  return utf(ftu(x) ^ (ftu(y) & SIGN_MASK));
}

// every mainstream cpu architecture has a sqrt instruction
float sqrtf(float x) { return __builtin_sqrtf(x); }

float fmaf(float a, float b, float c) {
#ifdef __FMA__
  return __builtin_fmaf(x, y, z);
#endif
  float ah = utf(ftu(a) & 0xfffff000);
  float bh = utf(ftu(b) & 0xfffff000);
  float al = a - ah;
  float bl = b - bh;
  return ((ah * bh + c) + ah * bl + al * bh) + al * bl;
}

static float acosf_poly(float x){
    float u = 2.2960134e-3f;
    u = u * x + -1.1146357e-2f;
    u = u * x + 2.6900099e-2f;
    u = u * x + -4.8802612e-2f;
    u = u * x + 8.875567e-2f;
    u = u * x + -2.1458527e-1f;
    return u * x + 1.5707962f;
}

float acosf(float x){
    float pi=3.14159265359f;
    float a = fabsf(x);
    float y = sqrtf(1.f-a)*acosf_poly(a);
    return mulsign(y,x)+(x<0?pi:0);
}

// y = (sqrt(1-((a*a-a)/f(a)+a)
// )-1)*1.57079632679489661923+1.57079632679489661923

static float log2f_mantissa(float x) {
#ifdef __FMA__
  return fmaf(fmaf(utf(0x40153ebb), x, utf(0x413b8af9)), x, utf(0x409c1a68)) *
         (x - 1.f) /
         fmaf(fmaf(fmaf(utf(0x3ecfca47), x, utf(0x409f8156)), x,
                   utf(0x40d76ca4)),
              x, 1.f);
#else
  return 5.7474026229742f * (x - 1.f) * (x + 0.45769302384949f) *
         (x + 4.5715238426130f) /
         ((x + 0.16943400804003f) * (x + 1.3506082324849f) *
          (x + 10.770128254704f));
#endif
}

float log2f(float x) {
  // log2(x*y) == log2(x)+log2(y)
  float mantissa = utf(ftu(1.f) | (ftu(x) & MANTISSA_MASK));
  float log2exponent =
      utf(ftu(256.f) | ((ftu(x) & EXPONENT_MASK) >> 8u)) - 383.f;
  return log2exponent + log2f_mantissa(mantissa);
}

float logf(float x) {
  float rlog2e = utf(0x3f317218u);
  return log2f(x) * rlog2e;
}

float log10f(float x) {
  float rlog210 = utf(0x3e9a209b);
  return log2f(x) * rlog210;
}

float log1pf(float x) { return logf(1.f + x); }

float acoshf(float x) { return logf(x + sqrtf(x * x - 1.)); };

float asinf(float x) {
  const float hpi = 1.57079632679489661923f;
  float a = fabsf(x);
  a = (a * a - a) /
          (-3.66063f - a * (1.75866f - a * (0.179323f - 0.0392588f * a))) +
      a;
  a = mulsign(sqrtf(1.f - a) - 1.f, x) * (-hpi);
  return a;
}

float asinhf(float x) { return logf(x + sqrtf(x * x + 1.)); };

float atanf_poly(float x) {
  float a = utf(0x3d267031);
  float b = utf(0x3f28513c);
  float c = utf(0x3e2f725f);
  float d = utf(0x3f7da425);
  float x2 = x * x;
#ifdef __FMA__
  return (fmaf(fmaf(a, x2, b), x2, 1.f) * x) / fmaf(fmaf(x2, c, d), x2, 1.f);
#endif
  return (((a * x2 + b) * x2 + 1.f) * x) / ((x2 * c + d) * x2 + 1.f);
}

float atanf(float x) {
  float a = fabsf(x);
  const float hpi = 1.57079632679489661923f;
  // silly code needed for gcc to avoid branching and to vectorize
  // float y = (a < 1.0f) ? a : 1.0f / a;
  unsigned int mask = a < 1.f ? 0xffffffffu : 0;
  float y = utf((mask & ftu(a)) | ((~mask) & ftu(1.f / a)));
  y = atanf_poly(y);
  return mulsign(
      //(a < 1.0f) ? y : hpi - y
      utf((mask & ftu(y)) | ((~mask) & ftu(hpi - y))), x);
}

float atan2f(float y, float x) {
  const float hpi = 1.57079632679489661923f;
  unsigned int nonzerox = fabsf(x) != 0.f ? 0xffffffffu : 0;
  unsigned int nonzeroy = fabsf(y) != 0.f ? 0xffffffffu : 0;
  unsigned int bothzero = (~nonzerox) & (~nonzeroy);
  float hpisignx = utf((nonzerox | bothzero) & ftu(mulsign(hpi, x)));
  return utf(nonzerox & ftu(atanf(y / x))) + mulsign(hpi - hpisignx, y);
}

float atanhf(float x) { return .5f * logf((1.f + x) / (1.f - x)); }

float cbrtf(float x) {
  float s = utf(0x2a510000u + ftu(x) / 3);     // wtf
  s = s - ((s * s) * s - x) / (3.f * (s * s)); // newton iteration
  s = s - ((s * s) * s - x) / (3.f * (s * s)); // newton iteration
  return s;
}

float ceilf(float x) {
#ifdef __SSE4_1__
  return __builtin_ceilf(x);
#endif
  unsigned int u = ftu(x);
  unsigned int e = u >> 23u & 0xff;
  unsigned int positive = x > 0.f ? 0xffffffffu : 0u;
  unsigned int zero = e < 127u ? 0xffffffffu : 0u;
  unsigned int m = MANTISSA_MASK >> (e - 127u);
  u += positive & m;
  u &= ~m;
  u = (zero & (positive & ftu(1.0f))) | (~zero & u);
  return utf(u);
}

float copysignf(float x, float y) {
  return utf((ftu(x) & 0x7fffffffu) | ftu(y) & 0x80000000u);
}

float roundf(float x) {
#ifdef __SSE4_1__
  return __builtin_roundf(x);
#endif
    float big = mulsign(utf((0x7f+23)<<23),x);
    return x+big-big;
}

float floorf(float x) {
#ifdef __SSE4_1__
  return __builtin_floorf(x);
#endif
  unsigned int u = ftu(x);
  unsigned int e = u >> 23u & 0xff;
  unsigned int negative = x < 0.f ? 0xffffffffu : 0u;
  unsigned int zero = e < 127u ? 0xffffffffu : 0u;
  unsigned int m = MANTISSA_MASK >> (e - 127u);
  u += negative & m;
  u &= ~m;
  u = (zero & (negative & ftu(-1.0f))) | (~zero & u);
  return utf(u);
}

static float sinf_poly(float x) {
  float a = utf(0xb2cc0ff1);
  float b = utf(0x3638a80e);
  float c = utf(0xb9500b44);
  float d = utf(0x3c088883);
  float e = utf(0xbe2aaaaa);
  float x2 = x * x;
#ifdef __FMA__
  return fmaf(fmaf(fmaf(fmaf(fmaf(a, x2, b), x2, c), x2, d), x2, e), x2, 1.f) *
         x;
#endif
  return (((((a * x2 + b) * x2 + c) * x2 + d) * x2 + e) * x2 + 1.f) * x;
}

float cosf(float x) {
  double tau =
      6.2831853071795864769252867665590057683943387987502116419498891846;
  float tauh = tau;
  float taul = tau - ((double)tauh);
  double rtau =
      0.1591549430918953357688837633725143620344596457404564487476673440;
  float rtauh = rtau;
  float rtaul = rtau - ((double)rtauh);
#ifdef __FMA__
  float m = floorf(fmaf(x, (rtauh * 2.f), x * (rtaul * 2.f)));
  float s = ((int)m & 1) * 2 - 1;

  float high = fmaf(m, tauh / 2.f, tauh / 4.f);
  float errorhigh = tauh / 4.f + fmaf(m, tauh / 2.f, -high);
  float low = fmaf(m, taul / 2.f, taul / 4.f);
#else
  float m = floorf(x * (rtauh * 2.f) + x * (rtaul * 2.f));
  float s = ((int)m & 1) * 2 - 1;

  float high = m * tauh / 2.f + tauh / 4.f;
  float errorhigh = tauh / 4.f + fmaf(m, tauh / 2.f, -high);
  float low = m * taul / 2.f + taul / 4.f;
#endif

  x -= high;
  x -= errorhigh;
  x -= low;

  return sinf_poly(x * s);
}
float exp2f_fract(float x) {
  float a = utf(0x3af71c15);
  float b = utf(0x3c130514);
  float c = utf(0x3d64b437);
  float d = utf(0x3e75ea9e);
  float e = utf(0x3f317271);
#ifdef __FMA__
  return fmaf(fmaf(fmaf(fmaf(fmaf(a, x, b), x, c), x, d), x, e), x, 1.f);
#endif
  return ((((a * x + b) * x + c) * x + d) * x + e) * x + 1.f;
}

float exp2f(float x) {
  // exp2(floor(x))*exp2(fract(x)) == exp2(x)
  unsigned int i = ftu(x + 383.f) << 8u;
  float exp2int = utf(i & EXPONENT_MASK);
#ifdef __SSE4_1__
  float fract = x - __builtin_floorf(x);
#else
  float fract = utf(ftu(1.f) | (i & MANTISSA_MASK)) - 1.f;
  fract += x - ((x + 383.f) - 383.f);
#endif
  return exp2int * exp2f_fract(fract);
}

float expf(float x) { return exp2f(x * utf(0x3fb8aa3bu)); }

float expm1f(float x) {
  unsigned int mask = fabsf(x) < .5f ? 0xffffffffu : 0;
  // pade approximant
  float a = x * (-2.f * x * x - 120.f) / (x * ((x - 12.f) * x + 60.f) - 120.f);
  float b = expf(x) - 1.f;
  return utf((mask & ftu(a)) | ((~mask) & ftu(b)));
}

static float erff_poly(float x) {
  float u = 3.118769e-4f;
  u = u * x + -4.67225e-3f;
  u = u * x + 3.3162573e-2f;
  u = u * x + -1.5214339e-1f;
  u = u * x + -9.1684705e-1f;
  u = u * x + -1.6282598f;
  return u * x + 3.1332566e-5f;
}

float erff(float x) {
  // pade approximant
  float a = x * (utf(0x3f906ebbu) + utf(0x3f174f6eu) * (x * x)) /
            (1.f + (x * x) * (utf(0x3f5b6db7u) + utf(0x3e3e2be3u) * (x * x)));
  float b = mulsign(1.f - exp2f(erff_poly(fabsf(x))), x);
  return fabsf(x) < 0.28f ? a : b;
}

float erfcf(float x) {
    float z = x<0.f?-1.f:1.f;
    float w = x<0.f?2.f:0.f;
    x = fabsf(x);
    x = x>10.f?10.f:x;
  float y = expf(-(x*x))*
    fmaf(fmaf(fmaf(fmaf(
        utf(0x35c42f59),x,
        utf(0x3daf42dc)),x,
        utf(0x3ee32e3c)),x,
        utf(0x3f7a7520)),x,1.f)/
    fmaf(fmaf(fmaf(fmaf(
        utf(0x3e1b69eb),x,
        utf(0x3f48fde0)),x,
        utf(0x3fe918cc)),x,
        utf(0x4006d465)),x,1.f);
    return fmaf(y,z,w);
}

float coshf(float x) { return .5f * (expf(x) + expf(-x)); }

float fdimf(float x, float y) {
  unsigned int greater = x > y ? 0xffffffffu : 0u;
  return utf(greater & ftu(x - y));
}



float fmaxf(float x, float y) { return x < y ? y : x; }

float fminf(float x, float y) { return x > y ? y : x; }

float frexpf(float x, int *e) {
  *e = (int)(ftu(x) >> 23u & 0xffu) - 126;
  return utf(ftu(0.5) | (ftu(x) & MANTISSA_MASK));
}

float hypotf(float x, float y) { return sqrtf(x * x + y * y); }

int ilogbf(float x) { return (int)(ftu(x) >> 23u & 0xffu) - 126; }

float ldexpf(float x, int n) { return utf(ftu(x) + (n << 23)); }

// float       lgammaf(float);

long long llrintf(float x) { return __builtin_llrintf(x); }

long long llroundf(float x) { return __builtin_llrintf(roundf(x)); }

float logbf(float x) {
  return utf(ftu(256.f) | ((ftu(x) & EXPONENT_MASK) >> 8u)) - 383.f;
}

long lrintf(float x) { return __builtin_lrintf(x); }

long lroundf(float x) { return __builtin_lrintf(roundf(x)); }

// float       nanf(const char *);
// float       nearbyintf(float);

float nextafterf(float x, float y) {
  int i = y < x ? -1 : 0;
  i |= 1;
  return utf(ftu(x) + i);
}

// float nexttowardf(float, long double);

float powf(float x, float y) { return exp2f(log2f(x) * y); }

float remainderf(float x, float y) { return x - roundf(x / y) * y; }

// float       remquof(float, float, int *);
// float       rintf(float);

float scalblnf(float x, long n) { return utf(ftu(x) + (n << 23)); }

float scalbnf(float x, int n) { return utf(ftu(x) + (n << 23)); }

float sinf(float x) {
  const float tau = 6.28318530717958647692f;
  const float taulo = utf(0x343bbd2e);
  const float rtau = 0.15915494309189533576f;
  const float pi = 3.14159265358979323846f;
  const float pilo = utf(0x33bbbd2e);
  const float rpi = 0.3183098861837907f;
  float z = x - roundf(x * rtau) * tau + roundf(x * rtau) * taulo;
  float y = fabsf(x - roundf(x * rpi) * pi + roundf(x * rpi) * pilo);
  return sinf_poly(mulsign(y, z));
}

float sinhf(float x) { return .5f * (expf(x) - expf(-x)); }

float tanf(float x) { return sinf(x) / cosf(x); }

float tanhf(float x) { return 1.f - 2.f / (expf(2.f * x) + 1.f); }

// float       tgammaf(float);

float truncf(float x) {
#ifdef __SSE4_1__
  return __builtin_truncf(x);
#endif
  int e = (int)(ftu(x) >> 23 & 0xff) - 127;
  unsigned int m = MANTISSA_MASK >> e;
  unsigned int nonzero = e < 0 ? 0u : 0xffffffffu;
  nonzero |= SIGN_MASK;
  return utf((ftu(x) & ~m) & nonzero);
}

float fmodf(float x, float y) { return x - truncf(x / y) * y; }

float modff(float x, float *y) {
  *y = truncf(x);
  return x - (*y);
}

#undef SIGN_MASK
#undef EXPONENT_MASK
#undef MANTISSA_MASK
