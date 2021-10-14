# jodiemath
a small, fast, vectorizeable libm

## goals:
* all functions have absolute or relative error ~= 0.000010 or less, this isn't a fast math library, it's meant to be general purpose
* at least attempt to be readable
* try to keep source and assembly output small first and fast second
* only for 32 bit floats, if you need the accuracy of doubles, this library is not for you
* all functions are vectorizeable on x86
* non-correct NaN and infinity behavior in most cases
* no handling of subnormal numbers

## caveats
* needs -fno-math-errno flag to for proper vectorization
* needs -msse4 for functions that rely on rounding to be fast
* gcc and clang might have bugs that prevent proper vectorization

## todo
* better expm1f
* better tanf
* implement gamma functions
* implement remquof
* write tests
