#include "fft.h"

#include <assert.h>
#include <limits.h>
#include <math.h>

#if defined (__GNUC__) || defined (__clang__)
#define fft_clz(n)                                \
  _Generic(n,                                     \
        unsigned int: __builtin_clz(n),           \
        unsigned long: __builtin_clzl(n),         \
        unsigned long long: __builtin_clzll(n))

#define likely(expr)    __builtin_expect(!!(expr), 1)
#define unlikely(expr)  __builtin_expect(!!(expr), 0)
#else

static inline size_t fft_clz(size_t n) {
  size_t count = sizeof (size_t) * CHAR_BIT / 2;
  size_t half = count / 2;
  while (half) {
    if (n >> count)
      count -= half;
    else
      count += half;
    half >> 1;
  }
}

#endif

static inline void rader(const fft_complex_t *restrict array, fft_complex_t *restrict target, size_t logsize) {
  size_t size = (size_t)1 << logsize;
  size_t shift = sizeof (size_t) * CHAR_BIT - logsize;
  for (size_t n = 0, reversed_n = 0; n < size; ++n) {
    FFT_COMPLEX_COPY(target[reversed_n], array[n]);

    /* get next reversed_n */
    reversed_n <<= shift;
    size_t count_leading_ones = fft_clz(~reversed_n);
    reversed_n <<= count_leading_ones;  /* remove leading ones */
    reversed_n |= (size_t)1 << (sizeof (size_t) * CHAR_BIT - 1);
    reversed_n >>= (shift + count_leading_ones);
  }
}

static inline void rader_inplace(fft_complex_t *array, size_t logsize) {
  size_t size = (size_t)1 << logsize;
  size_t shift = sizeof (size_t) * CHAR_BIT - logsize;
  /* nothing should be done for 0 and 0b111...11(size - 1). */
  for (size_t n = 1, reversed_n = size >> 1; n < size - 1; ++n) {
    if (n < reversed_n)
      FFT_COMPLEX_SWAP(array[n], array[reversed_n]);

    /* get next reversed_n */
    reversed_n <<= shift;
    size_t count_leading_ones = fft_clz(~reversed_n);
    reversed_n <<= count_leading_ones;  /* remove leading ones */
    reversed_n |= (size_t)1 << (sizeof (size_t) * CHAR_BIT - 1);
    reversed_n >>= (shift + count_leading_ones);
  }
}

#define DO_BUTTERFLY(begin, end, step) do {                       \
  fft_complex_t unit;                                             \
  FFT_COMPLEX_UNITROOT_RECIP(unit, step);                         \
  size_t half_step = (step) / 2;                                  \
  for (fft_complex_t *p = (begin); p != (end); p += (step)) {     \
    fft_complex_t root;                                           \
    FFT_COMPLEX_SETONE(root);                                     \
    for (size_t i = 0, j = half_step; i < half_step; ++i, ++j) {  \
      fft_complex_t t;                                            \
      FFT_COMPLEX_MUL(t, root, p[j]);                             \
      FFT_COMPLEX_SELFMUL(root, unit);                            \
      fft_complex_t u;                                            \
      FFT_COMPLEX_COPY(u, p[i]);                                  \
      FFT_COMPLEX_ADD(p[i], u, t);                                \
      FFT_COMPLEX_SUB(p[j], u, t);                                \
    }                                                             \
  }                                                               \
} while (0)

void fft_raw(fft_complex_t *x, size_t logsize) {
  if (unlikely(logsize == 0))
    return;

  fft_complex_t *begin = x;
  fft_complex_t *end = begin + ((size_t)1 << logsize);

  /* handle small step separately, compiler will possibly do loop unrolling */
  DO_BUTTERFLY(begin, end, 2);

  if (unlikely(logsize == 1)) /* size == 2 ? */
    return;

  DO_BUTTERFLY(begin, end, 4);

  if (unlikely(logsize == 2)) /* size == 4 ? */
    return;

  /* do generic butterfly in a loop */
  for (size_t step = 8; step <= (size_t)1 << logsize; step *= 2)
    DO_BUTTERFLY(begin, end, step);

}

void fft(const fft_complex_t *restrict x, fft_complex_t *restrict X, size_t logsize) {
  rader(x, X, logsize);
  fft_raw(X, logsize);
}

void fft_inplace(fft_complex_t *x, size_t logsize) {
  rader_inplace(x, logsize);
  fft_raw(x, logsize);
}
