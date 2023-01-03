#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <stdio.h>

#define s21_LN2 0.693147180559945309417
#define S21_PI 3.14159265358979323846
#define S21_DBL_MAX 1.7976931348623158e+308
#define S21_DBL_MIN 2.2250738585072014e-308
#define S21_FLT_MAX 3.40282346638528859811704183484516925e+38F
#define S21_LONG_MAX 0x7fffffffffffffffL
#define S21_EXP 2.7182818284590452353602874713526624
#define S21_EPS 1e-9
#define S21_FLT_MAX 3.40282346638528859811704183484516925e+38F
#define S21_LONG_MAX 0x7fffffffffffffffL

#define S21_INF 1.0 / 0.0
#define S21_NAN 0.0 / 0.0

#define is_fin(x) __builtin_isfinite(x)
#define is_nan(x) __builtin_isnan(x)
#define is_inf(x) __builtin_isinf(x)

int s21_abs(int x);
long double s21_fabs(double x);
long double s21_sqrt(double x);
long double s21_log(double x);
long double s21_exp(double x);
long double s21_floor(double x);
long double pow_fmod(double x, double y);
long double s21_ceil(double x);
long double s21_fmod(double x, double y);
long double s21_pow(double x, double y);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
double bi_pow(double x, unsigned long long i);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_acos(double x);

#endif  // SRC_S21_MATH_H_
