#ifndef MATH_H
#define MATH_H

#define M_PI_4  0.785398163
#define M_PI_2  1.570796326
#define M_PI    3.141592653
#define M_2_PI  6.283185307

inline float DegToRad(float d)
{
    return d * 0.017453293; // PI / 180
}

inline float RadToDeg(float r)
{
    return r * 57.295779524; // 180 / PI
}

inline int Abs(int v) {
    return __builtin_abs(v);
}

inline float Absf(float v) {
    return __builtin_fabs(v);
}

inline float SquareRootf(float v) {
    return __builtin_sqrt(v);
}

inline double SquareRoot(double v) {
    return __builtin_sqrt(v);
}

inline float Cosf(float a) {
    return __builtin_cos(a);
}

inline double Cos(double a) {
    return __builtin_cos(a);
}

inline float Acosf(float a) {
    return __builtin_acos(a);
}

inline double Acos(double a) {
    return __builtin_acos(a);
}

inline float Sinf(float a) {
    return __builtin_sin(a);
}

inline double Sin(double a) {
    return __builtin_sin(a);
}

inline float Asinf(float a) {
    return __builtin_asin(a);
}

inline double Asin(double a) {
    return __builtin_asin(a);
}

inline float Tanf(float a) {
    return __builtin_tan(a);
}

inline double Tan(double a) {
    return __builtin_tan(a);
}

inline float Atanf(float a) {
    return __builtin_atan(a);
}

inline double Atan(double a) {
    return __builtin_atan(a);
}

inline float Atan2f(float a, float b) {
    return __builtin_atan2(a, b);
}

inline double Atan2(double a, double b) {
    return __builtin_atan2(a, b);
}

#define F_EPS 1e-6

bool float_Equal(float a, float b)
{
    return Absf(Absf(a) - Absf(b)) < F_EPS;
}

#endif // MATH_H
