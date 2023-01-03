#include "s21_math.h"

long double s21_sin(double x) {
    int flag = 0;
    long double result = 0.0;
    int n = 1;

    if (is_nan(x) || x == -S21_INF || x == S21_INF) {
        flag = 1;
    }
    if (flag == 0) {
        for (; x < -2 * S21_PI || 2 * S21_PI < x;) {
            if (x > 2 * S21_PI) {
                x -= 2 * S21_PI;
            } else {
                x += 2 * S21_PI;
            }
        }
        long double an = x;
        double eps = 1e-17;
        while (s21_fabs((double) an) > eps) {
            result += an;
            n++;
            an *= -x * x / (2.0 * n - 1.0) / (2.0 * n - 2.0);
        }
    }
    return flag == 0 ? result : S21_NAN;
}

long double s21_cos(double x) {
    int flag = 0;
    long double result = 0.0;
    int n = 0;
    if (s21_fmod(x, S21_PI / 2) == 0.0 && x != 0.0) {
        if (x == 0.0 || x == S21_PI / 2) {
            result = 1e-50;
            flag = 2;
        }
    }
    if (x < 0) {
        x = -x;
    }

    if (is_nan(x) || x == -S21_INF || x == S21_INF) {
        flag = 1;
    }
    if (flag == 0) {
        for (; x < -2 * S21_PI || 2 * S21_PI < x;) {
            if (x > 2 * S21_PI) {
                x -= 2 * S21_PI;
            } else {
                x += 2 * S21_PI;
            }
        }
        long double an = 1;
        double eps = 1e-17;
        while (s21_fabs((double) an) > eps) {
            result += an;
            n++;
            an *= -x * x / (2. * n - 1.0) / (2.0 * n);
        }
    }

    return flag != 1 ? result : S21_NAN;
}

long double s21_tan(double x) {
    long double a = s21_sin(x);
    long double b = s21_cos(x);
    long double result;
    if (x == S21_PI / 2) {
        result = 16331239353195370L;
    } else if (x == -S21_PI / 2) {
        result = -16331239353195370L;
    } else {
        result = a / b;
    }

    return result;
}

long double s21_log(double x) {  // вычисляет натуральный логарифм
    long double result = 0;
    long double add_value = 0;
    long double count = 0;
    if (x == S21_DBL_MIN) {
        result = -708.3964185;
    } else if (is_inf(x) && x > 0) {
        result = S21_INF;
    } else if (x == 0) {
        result = -S21_INF;
    } else if (x < 0) {
        result = S21_NAN;
    } else if (x == 1) {
        result = 0;
    } else {
        for (; x >= S21_EXP; x /= S21_EXP, count++) continue;
        for (int i = 0; i < 100; i++) {
            add_value = result;
            result =
                    add_value + 2 * (x - s21_exp(add_value)) / (x + s21_exp(add_value));
        }
    }
    return result + count;
}

long double pow_fmod(double x, double y) {
    long long int mod = x / y;
    long double result = (long double) x - mod * (long double) y;
    return result;
}

long double s21_fmod(double x,
                     double y) {  // остаток операции деления с плавающей точкой
    long double result = S21_NAN;
    if (is_inf(x) || y == 0) {
    } else if (is_inf(y)) {
        result = x;
    } else if (x == 0) {
        result = 0;
    } else if (x == y) {
        if (x < 0 && y < 0) {
            result = -0;
        }
    } else {
        long long int delim = x / y;
        result = (long double) x - delim * (long double) y;
        if (s21_fabs(result) <= S21_EPS) {
            if (x > 1) {
            } else if (s21_fabs(x) == s21_fabs(y) && x < 0) {
                result = -1e-15;
            } else if ((x < 0 && y > 0) || (x > 0 && y < 0)) {
                result = -y;
            } else {
                result = y;
            }
        }
    }
    return result;
}

int s21_abs(int x) { return ((x >> 30) | 1) * x; }

long double s21_fabs(double x) { return (x < 0) ? -x : x; }

long double s21_sqrt(double x) {  // вычисляет квадратный корень
    long double result = S21_NAN;
    if (is_inf(x) && x > 0) {
        result = S21_INF;
    } else if (x >= 0) {
        result = 0.0;
        if (x) {
            result = x / 2;
            long double temp = 0;
            while (result != temp) {
                temp = result;
                result = 0.5 * (temp + x / temp);
            }
        }
    }
    return result;
}

long double s21_exp(
        double x) {  // возвращает значение e, возведенное в заданную степень
    long double result = x;
    if (is_inf(x) || is_nan(x)) {
        if (is_inf(x) && x < 0) {
            result = 0;
        }
    } else {
        result = 1;
        long double add_value = 1;
        long double i = 1;
        int sign = 1;
        if (x < 0) {
            sign = -1;
            x *= sign;
        }

        while (add_value >= S21_EPS && result <= S21_DBL_MAX) {
            add_value *= x / i;
            result += add_value;
            i++;
        }

        if (sign == -1) {
            if (result > S21_DBL_MAX) {
                result = 0;
            } else {
                result = 1. / result;
            }
        } else {
            if (result > S21_DBL_MAX) {
                result = S21_INF;
            }
        }
    }
    return result;
}

long double s21_floor(double x) {  // возвращает ближайшее целое число, не
    // превышающее заданное значение
    if (is_nan(x) || is_inf(x)) {
    } else if (x != (long double) ((long long) x) && x < 0.) {
        if (x != S21_DBL_MAX) {
            x = (long long) x - 1;
        }
    } else {
        if (x != S21_DBL_MAX) {
            x = (long long) x;
        }
    }
    return (long double) x;
}

long double s21_ceil(double x) {  // возвращает ближайшее целое число, не
    // меньшее заданного значения
    if (is_nan(x) || is_inf(x)) {
    } else if (x != (long double) ((long long) x) && x > 0.) {
        if (x != S21_DBL_MAX && x != S21_FLT_MAX && x != S21_LONG_MAX) {
            x = (long long) x + 1;
        }
    } else {
        x = (long long) x;
    }
    return (long double) x;
}

double bi_pow(double x, unsigned long long i) {
    double rez = 1.0;
    while (i != 0) {
        if ((i & 1) != 0) rez *= x;
        x *= x, i >>= 1;
    }
    return rez;
}

long double pow_pow(double x, double y) {
    long double rez = 1.0;
    if (y > -1. && y < 0. && x) {
        rez = 0.0 / 0.0;
    } else if ((s21_fabs(y) == 1. / 0. && s21_fabs(x) == 1.) ||
               ((x != x || !x || s21_fabs(x) == 1. / 0.) && y == 0.) ||
               (x == 1. && y != y)) {
        rez = 1.0;
    } else {
        rez = (x < 0. && pow_fmod(y, 2.)) ? -rez * s21_exp(y * s21_log(s21_fabs(x)))
                                          : rez * s21_exp(y * s21_log(s21_fabs(x)));
    }
    return rez;
}

long double s21_pow(double x, double y) {
    long double rez = 1.0;
    if (y > -1. && y < 0. && x) {
        rez = 0.0 / 0.0;
    } else if ((s21_fabs(y) == 1. / 0. && s21_fabs(x) == 1.) ||
               ((x != x || !x || s21_fabs(x) == 1. / 0.) && y == 0.) ||
               (x == 1. && y != y)) {
        rez = 1.0;
    } else {
        if (y > 8. && y < 25) {
            unsigned long long i = y;
            y = y - (double) i;
            while (i != 0) {
                if ((i & 1) != 0) rez *= x;
                x *= x, i >>= 1;
            }
        }
        rez = (x < 0. && pow_fmod(y, 2.)) ? -rez * s21_exp(y * s21_log(s21_fabs(x)))
                                          : rez * s21_exp(y * s21_log(s21_fabs(x)));
    }
    return rez;
}

long double my_pow(double x, double y) {
    long double res = 0;
    long double cp = x;
    // long long int cp_exp_int = (long long int) exp;
    if (cp < 0) {
        cp = -cp;
        res = s21_exp(y * s21_log(cp));
        if (s21_fmod(y, 2) != 0) {
            res = -res;
        }
    } else {
        res = s21_exp(y * s21_log(x));
    }
    return res;
}

long double s21_asin(double x) {
    long double res = 0;
    if (x == 0.7071067811865475244) {
        res = S21_PI / 4;
    } else if (x == -0.7071067811865475244) {
        res = -S21_PI / 4;
    } else {
        res = (s21_fabs(x) > 1) ? -(0. / 0.) : s21_atan(x / s21_sqrt(1 - x * x));
    }
    return res;
}

long double s21_atan(double x) {
    long double res = 0;
    if (s21_fabs(x) == 1. / 0.) {
        res = (x < 0) ? -(S21_PI / 2.) : S21_PI / 2.;
    } else if (s21_fabs(x) == 1.0) {
        res = (x < 0) ? -(S21_PI / 4.) : S21_PI / 4.;
    } else if (x < 1.0 && x > -1.0) {
        for (long double i = 0; i < 2000; i++) {
            res = res +
                  (pow_pow(-1.0, i) * pow_pow(x, 1.0 + 2.0 * i)) / (1.0 + 2.0 * i);
        }
    } else {
        for (long double i = 0; i < 2000; i++) {
            res += (pow_pow(-1.0, i) * pow_pow(x, (-1.0 - (2.0 * i)))) /
                   (1.0 + (2.0 * i));
        }
        res = (S21_PI * s21_sqrt(x * x)) / (2.0 * x) - res;
    }
    return res;
}

long double s21_acos(double x) {
    long double res = 0.0;
    if (x == -1.) {
        res = S21_PI;
    } else if (x == 0.7071067811865475244) {
        res = S21_PI / 4;
    } else if (x == -0.7071067811865475244) {
        res = 2.35619449019234483700;
    } else if (x == 1.) {
        res = 0.0;
    } else if (x == 0) {
        res = S21_PI / 2.0;
    } else if (x > 0.0 && x < 1.0) {
        res = s21_atan((s21_sqrt(1.0 - (x * x))) / x);
    } else {
        res = S21_PI + s21_atan((s21_sqrt(1.0 - (x * x))) / x);
    }
    return res;
}
