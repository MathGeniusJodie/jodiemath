#define NAN       (0.0f/0.0f)
#define INFINITY  1e5000f
#define HUGE_VALF INFINITY

#define SIGN_MASK 0x80000000u
#define EXPONENT_MASK 0x7f800000u
#define MANTISSA_MASK 0x007fffffu

static float utf(unsigned int i){
	union {unsigned int i;float f;} u;
	u.i = i;
	return u.f;
}

static unsigned int ftu(float f){
	union {unsigned int i;float f;} u;
	u.f = f;
	return u.i;
}

float fabsf(float x){
	return utf(ftu(x) & 0x7fffffffu);
}

#define FP_INFINITE 0
#define FP_NAN 1
#define FP_NORMAL 2
#define FP_SUBNORMAL 3
#define FP_ZERO 4

int fpclassify(float x){
	unsigned int i = ftu(fabsf(x));
	return x == 0x7f800000u?FP_INFINITE:
	x > 0x7f800000u?FP_NAN:
	x < 0x00800000u?FP_SUBNORMAL:
	x == 0?FP_ZERO:FP_NORMAL;
}

int isinf(float x){
	return ftu(fabsf(x))==0x7f800000u;
}

int isnan(float x){
	return ftu(fabsf(x))>0x7f800000u;
}

int isnormal(float x){
	return (ftu(x) & EXPONENT_MASK) > 0;
}

int isfinite(float x){
	return ftu(fabsf(x))<0x7f800000u;
}

int signbit(float x){
	return ftu(x)>>31;
}

#define isunordered(x,y) (isnan(x)||isnan(y)?1.f:0.f)

#define isless(x, y)           ((x)<(y)?1.f:0.f)
#define islessequal(x, y)      ((x)<=(y)?1.f:0.f)
#define islessgreater(x, y)    ((x)!=(y)?1.f:0.f)
#define isgreater(x, y)        ((x)>(y)?1.f:0.f)
#define isgreaterequal(x, y)   ((x)>=(y)?1.f:0.f)

// x * sign(y)
static float mulsign(float x, float y){
    return utf(ftu(x)^(ftu(y)&SIGN_MASK));
}

// every mainstream cpu architecture has a sqrt instruction
float sqrtf(float x){
	return __builtin_sqrtf(x);
}

float acosf(float x){
    const float hpi = 1.57079632679489661923f; 
    float a = fabsf(x);
    a = (a*a-a)/(
        -3.66063f-a*(1.75866f-a*(0.179323f-0.0392588f*a))
    )+a;
    a = mulsign(sqrtf(1.f-a)-1.f,x)*hpi+hpi;
    return a;
}

// y = (sqrt(1-((a*a-a)/f(a)+a)  )-1)*1.57079632679489661923+1.57079632679489661923

static float log2f_mantissa(float x){
    return 5.7474026229742f*(x - 1.f)*(x + 0.45769302384949f)*(x + 4.5715238426130f)/
    ((x + 0.16943400804003f)*(x + 1.3506082324849f)*(x + 10.770128254704f));
}

float log2f(float x){
  //log2(x*y) == log2(x)+log2(y)
  float mantissa = utf(ftu(1.f) | (ftu(x) & MANTISSA_MASK));
  float log2exponent = utf(ftu(256.f) | ((ftu(x)&EXPONENT_MASK) >> 8u))-383.f;
  return log2exponent+log2f_mantissa(mantissa);
}

float logf(float x){
	float rlog2e = utf(0x3f317218u);
	return log2f(x)*rlog2e;
}

float log10f(float x){
    float rlog210 = utf(0x3e9a209b);
    return log2f(x)*rlog210;
}

float log1pf(float x){
	return logf(1.f+x);
}

float acoshf(float x) {
	return logf(x + sqrtf(x * x - 1.));
};

float asinf(float x){
	const float hpi = 1.57079632679489661923f; 
	float a = fabsf(x);
	a = (a*a-a)/(
		-3.66063f-a*(1.75866f-a*(0.179323f-0.0392588f*a))
	)+a;
	a = mulsign(sqrtf(1.f-a)-1.f,x)*(-hpi);
	return a;
}

float asinhf(float x) {
	return logf(x + sqrtf(x * x + 1.));
};

float atanf(float x) {
	float a = fabsf(x);
	const float hpi = 1.57079632679489661923f;
	// silly code needed for gcc to avoid branching and to vectorize
	//float y = (a < 1.0f) ? a : 1.0f / a;
	unsigned int mask = a<1.f?0xffffffffu:0;
	float y = utf((mask&ftu(a))|((~mask)&ftu(1.f/a)));
	y = y+y*y*y/((0.203984f*y-1.86607f)*y*y-2.99766f);
	return mulsign(
		//(a < 1.0f) ? y : hpi - y
		utf((mask&ftu(y)) | ((~mask)&ftu(hpi - y)))
	,x);
}

float atan2f(float y, float x){
	const float hpi = 1.57079632679489661923f;
	unsigned int nonzerox = fabsf(x)!=0.f?0xffffffffu:0;
	unsigned int nonzeroy = fabsf(y)!=0.f?0xffffffffu:0;
	unsigned int bothzero = (~nonzerox)&(~nonzeroy);
	float hpisignx = utf((nonzerox|bothzero) & ftu(mulsign(hpi,x)));
	return utf(nonzerox & ftu(atanf(y/x)))
		+ mulsign(hpi-hpisignx,y);
}

float atanhf(float x){
	return .5f * logf((1.f+x)/(1.f-x));
}

float cbrtf(float x){
	float s = utf(0x2a510000u+ftu(x)/3); //wtf
	s=s-((s*s)*s-x)/(3.f*(s*s)); // newton iteration
	s=s-((s*s)*s-x)/(3.f*(s*s)); // newton iteration
	return s;
}

float ceilf(float x){
	#ifdef __SSE4_1__
	return __builtin_ceilf(x);
	#endif
	unsigned int u = ftu(x);
	unsigned int e = u >> 23u & 0xff;
	unsigned int positive = x>0.f?0xffffffffu:0u;
	unsigned int zero = e < 127u?0xffffffffu:0u;
	unsigned int m = MANTISSA_MASK >> (e-127u);
	u+=positive&m;
	u &= ~m;
	u = (zero & (positive & ftu(1.0f))) | (~zero&u);
	return utf(u);
}

float copysignf( float x, float y ){
	return utf((ftu(x) & 0x7fffffffu)| ftu(y) & 0x80000000u);
}

float roundf(float x){
	#ifdef __SSE4_1__
	return __builtin_roundf(x);
	#endif
	const float big=utf((0x7f+23)<<23);
	return mulsign(fabsf(x)+big-big, x);
}

static float cosf_poly(float x){
    x*=x;
    return x*(x*((
        utf(0x37bf85f0u)*x+
        utf(0xbab559edu))*x+
        utf(0x3d2aa41eu))+
        utf(0xbeffff9eu))+
        utf(0x3f7ffffbu);
}

float cosf(float x){
	const float tau = 6.28318530717958647692f;
	const float rtau = 0.15915494309189533576f;
	const float hpi = 1.57079632679489661923f;
	float z = (x+hpi)*rtau;
	float b = z-roundf(z);
	return mulsign(cosf_poly(fabsf(b)*tau-hpi),b);
}

static float exp2f_fract(float x) {
    float a = utf(0xc0d7a065u);
    float b = utf(0x412a97f3u);
    float c = utf(0x42573661u);
    return (a-x)*((x+b)*x+c)/
          ((a+x)*((x-b)*x+c));
}

float exp2f(float x){
	//exp2(floor(x))*exp2(fract(x)) == exp2(x)
	unsigned int i = ftu(x + 383.f)<<8u;
	float exp2int = utf(i & EXPONENT_MASK);
	#ifdef __SSE4_1__
	float fract = x-__builtin_floorf(x);
	#else
	float fract = utf(ftu(1.f) | (i & MANTISSA_MASK)) - 1.f;
	fract += x-((x+383.f)-383.f);
	#endif
	return exp2int * exp2f_fract(fract);
}

float expf(float x) {
	return exp2f(x*utf(0x3fb8aa3bu));
}

float expm1f(float x) {
	return expf(x)-1.f;
}

static float erff_poly(float x){
	return ((((-0.00289436*x+0.0292923)*x-0.148118)*x-0.918879)*x-1.62781)*x;
}

float erff(float x){
	return mulsign(
		1.f-exp2f(erff_poly(fabsf(x)))
	,x);
}

float erfcf(float x) {
	return 1.f-erff(x);
}

float coshf(float x){
	return .5f*(expf(x)+expf(-x));
}

float fdimf(float x, float y){
	unsigned int greater = x > y ? 0xffffffffu : 0u;
	return utf(greater & ftu(x-y));
}

float floorf(float x){
	#ifdef __SSE4_1__
	return __builtin_floorf(x);
	#endif
	unsigned int u = ftu(x);
	unsigned int e = u >> 23u & 0xff;
	unsigned int negative = x<0.f?0xffffffffu:0u;
	unsigned int zero = e < 127u?0xffffffffu:0u;
	unsigned int m = MANTISSA_MASK >> (e-127u);
	u+=negative&m;
	u &= ~m;
	u = (zero & (negative & ftu(-1.0f))) | (~zero&u);
	return utf(u);
}

float fmaf(float x,float y,float z){
	#ifdef __FMA__
	return __builtin_fmaf(x,y,z);
	#endif
	return x*y+z;
}

float fmaxf(float x,float y){
	return x<y?y:x;
}

float fminf(float x,float y){
	return x>y?y:x;
}

float frexpf(float x, int *e){
	*e = (int)(ftu(x)>>23u & 0xffu) - 126;
    return utf(ftu(0.5)|(ftu(x)&MANTISSA_MASK));
}

float hypotf(float x, float y){
	return sqrtf(x*x+y*y);
}

int ilogbf(float x){
	return (int)(ftu(x)>>23u & 0xffu) - 126;
}

float ldexpf(float x, int n){
	return utf(ftu(x)+(n<<23));
}

//float       lgammaf(float);

long long llrintf(float x){
	return __builtin_llrintf(x);
}

long long llroundf(float x){
	return __builtin_llrintf(roundf(x));
}

float logbf(float x){
  return utf(ftu(256.f) | ((ftu(x)&EXPONENT_MASK) >> 8u))-383.f;
}

long lrintf(float x){
	return __builtin_lrintf(x);
}

long lroundf(float x){
	return __builtin_lrintf(roundf(x));
}

//float       nanf(const char *);
//float       nearbyintf(float);

float nextafterf(float x, float y){
	int i = y<x?-1:0;
	i |= 1;
	return utf(ftu(x)+i);
}

//float nexttowardf(float, long double);

float powf(float x, float y){
	return exp2f(log2f(x)*y);
}

float remainderf(float x, float y){
    return x-roundf(x/y)*y;
}

//float       remquof(float, float, int *);
//float       rintf(float);

float scalblnf(float x, long n){
	return utf(ftu(x)+(n<<23));
}

float scalbnf(float x, int n){
	return utf(ftu(x)+(n<<23));
}

float sinf(float x){
	const float tau = 6.28318530717958647692f;
	const float rtau = 0.15915494309189533576f;
	const float hpi = 1.57079632679489661923f;
	const float pi = 3.14159265358979323846f;

	float z = (pi-x)*rtau;
	float b = z-roundf(z);
	return mulsign(cosf_poly(fabsf(b)*tau-hpi),b);
}

float sinhf(float x){
	return .5f*(expf(x)-expf(-x));
}

float tanf(float x){
	return sinf(x)/cosf(x);
}

float tanhf(float x){
	return 1.f - 2.f / (expf(2.f*x)+1.f);
}

//float       tgammaf(float);

float truncf(float x){
	#ifdef __SSE4_1__
	return __builtin_truncf(x);
	#endif
	int e = (int)(ftu(x) >> 23 & 0xff) - 127;
	unsigned int m = MANTISSA_MASK >> e;
	unsigned int nonzero = e < 0?0u:0xffffffffu;
	nonzero |= SIGN_MASK;
	return utf((ftu(x) & ~m) & nonzero);
}

float fmodf(float x,float y){
	return x-truncf(x/y)*y;
}

float modff(float x,float * y){
    *y = truncf(x);
    return x-(*y);
}

#undef SIGN_MASK
#undef EXPONENT_MASK
#undef MANTISSA_MASK

