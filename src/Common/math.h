#ifndef MATH_H
#define MATH_H

#define M_PI_4  0.785398163
#define M_PI_2  1.570796326
#define M_PI    3.141592653
#define M_2_PI  6.283185307

float DegToRad(float d)
{
    return d * 0.017453293; // PI / 180
}

float RadToDeg(float r)
{
    return r * 57.295779524; // 180 / PI
}

int Abs(int v) {
    return __builtin_abs(v);
}

float Absf(float v) {
    return __builtin_fabs(v);
}

float SquareRootf(float v) {
    return __builtin_sqrt(v);
}

double SquareRoot(double v) {
    return __builtin_sqrt(v);
}

float Cosf(float a) {
    return __builtin_cos(a);
}

double Cos(double a) {
    return __builtin_cos(a);
}

float Acosf(float a) {
    return __builtin_acos(a);
}

double Acos(double a) {
    return __builtin_acos(a);
}

float Sinf(float a) {
    return __builtin_sin(a);
}

double Sin(double a) {
    return __builtin_sin(a);
}

float Asinf(float a) {
    return __builtin_asin(a);
}

double Asin(double a) {
    return __builtin_asin(a);
}

float Tanf(float a) {
    return __builtin_tan(a);
}

double Tan(double a) {
    return __builtin_tan(a);
}

float Atanf(float a) {
    return __builtin_atan(a);
}

double Atan(double a) {
    return __builtin_atan(a);
}

float Atan2f(float a, float b) {
    return __builtin_atan2(a, b);
}

double Atan2(double a, double b) {
    return __builtin_atan2(a, b);
}

#define F_EPS 1e-6

bool float_Equal(float a, float b)
{
    return Absf(Absf(a) - Absf(b)) < F_EPS;
}

#endif // MATH_H
