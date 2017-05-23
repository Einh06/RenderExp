#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "common.h"
#include "math.h"

/**********************************
 *             TYPES              *
 **********************************/
typedef struct Vec2f {
    float x, y;
} Vec2f;

typedef struct Vec3f {
    float x, y, z;
} Vec3f;

typedef struct Vec4f {
    float x, y, z, w;
} Vec4f;

/**********************************
 *          2D VECTOR             *
 **********************************/

Vec2f Vec2f_Make(float X, float Y)
{
    return (Vec2f) {.x = X, .y = Y};
}

#define Vec2f_X Vec2f_Make(1, 0)
#define Vec2f_Y Vec2f_Make(0, 1)

Vec2f Vec2f_Add(Vec2f a, Vec2f b) 
{
    return Vec2f_Make( a.x + b.x, a.y + b.y );
}

#define Vec2f_Add3(a, b, c)         Vec2f_Add(Vec2f_Add(a, b), c)
#define Vec2f_Add4(a, b, c, d)      Vec2f_Add(Vec2f_Add3(a, b, c), d)
#define Vec2f_Add5(a, b, c, d, e)   Vec2f_Add(Vec2f_Add4(a, b, c, d), e)

Vec2f Vec2f_Substact(Vec2f a, Vec2f b) 
{
    return Vec2f_Make( a.x - b.x, a.y - b.y );
}

Vec2f Vec2f_MakeInvert(Vec2f v)
{
    return Vec2f_Make( -v.x, -v.y );
}

Vec2f Vec2f_Mult(Vec2f a, float s) 
{
    return Vec2f_Make(a.x * s, a.y * s );
}

float Vec2f_Length(Vec2f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y));
}

Vec2f Vec2f_NormalizeFrom(Vec2f v)
{
    float Ratio = 1 / Vec2f_Length(v);
    Vec2f Result = Vec2f_Mult(v, Ratio);
    return Result;
}

void Vec2f_Normalize(Vec2f* v)
{
    Vec2f Normalized = Vec2f_NormalizeFrom(*v);
    *v = Normalized;
}

float Vec2f_Dot(Vec2f a, Vec2f b)
{
    return (a.x * b.x) + (a.y * b.y);
}

Vec2f Vec2f_Cross(Vec2f a)
{
    Vec2f Result = { -a.y, a.x }; 
    return Result;
}

Vec2f Vec2f_FromVec3f(Vec3f v)
{
    return Vec2f_Make(v.x / v.z, v.y / v.z);
}


/**********************************
 *          3D VECTOR             *
 **********************************/

Vec3f Vec3f_Make(float X, float Y, float Z)
{
    return (Vec3f) {.x = X, .y = Y, .z = Z};
}

#define Vec3f_X Vec3f_Make(1, 0, 0)
#define Vec3f_Y Vec3f_Make(0, 1, 0)
#define Vec3f_Z Vec3f_Make(0, 0, 1)


Vec3f Vec3f_Add(Vec3f a, Vec3f b) 
{
    return Vec3f_Make( a.x + b.x, a.y + b.y, a.z + b.z );
}

#define Vec3f_Add3(a, b, c)         Vec3f_Add(Vec3f_Add(a, b), c)
#define Vec3f_Add4(a, b, c, d)      Vec3f_Add(Vec3f_Add3(a, b, c), d)
#define Vec3f_Add5(a, b, c, d, e)   Vec3f_Add(Vec3f_Add4(a, b, c, d), e)

Vec3f Vec3f_MakeInvert(Vec3f v)
{
    return Vec3f_Make( -v.x, -v.y, v.z);
}

Vec3f Vec3f_Substract(Vec3f a, Vec3f b) 
{
    return Vec3f_Make( a.x - b.x, a.y - b.y, a.z - b.z );
}

Vec3f Vec3f_Mult(Vec3f a, float s) 
{
    return Vec3f_Make( a.x * s, a.y * s, a.z * s );
}

float Vec3f_Length(Vec3f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

Vec3f Vec3f_NormalizeFrom(Vec3f v)
{
    float Ratio = 1.0 / Vec3f_Length(v);
    Vec3f Result = Vec3f_Mult(v, Ratio);
    return Result;
}

void Vec3f_Normalize(Vec3f* v)
{
    Vec3f Normalized = Vec3f_NormalizeFrom(*v);
    *v = Normalized;
}

float Vec3f_Dot(Vec3f a, Vec3f b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Vec3f Vec3f_Cross(Vec3f a, Vec3f b)
{
    Vec3f Result;

    Result.x = (a.y * b.z) - (a.z * b.y);
    Result.y = (a.z * b.x) - (a.x * b.z);
    Result.z = (a.x * b.y) - (a.y * b.x);

    return Result;
}

Vec3f Vec3f_FromVec2f(Vec2f v)
{
    return Vec3f_Make(v.x, v.y, 1);
}

Vec3f Vec3f_FromVec4f(Vec4f v) 
{
    Vec3f Result = {v.x, v.y, v.z}; 
    float Divider = v.w != 0 ? 1.0/v.w : 1.0;
    return Vec3f_Mult(Result, Divider);
}


/**********************************
 *          4D VECTOR             *
 **********************************/

Vec4f Vec4f_Make(float X, float Y, float Z, float W)
{
    return (Vec4f){.x = X, .y = Y, .z = Z, .w = W}; 
}


Vec4f Vec4f_Add(Vec4f a, Vec4f b) 
{
    return Vec4f_Make( a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w );
}

#define Vec4f_Add3(a, b, c)         Vec4f_Add(Vec4f_Add(a, b), c)
#define Vec4f_Add4(a, b, c, d)      Vec4f_Add(Vec4f_Add3(a, b, c), d)
#define Vec4f_Add5(a, b, c, d, e)   Vec4f_Add(Vec4f_Add4(a, b, c, d), e)

Vec4f Vec4f_MakeInvert(Vec4f v)
{
    return Vec4f_Make( -v.x, -v.y, -v.z, -v.w );
}


Vec4f Vec4f_Substract(Vec4f a, Vec4f b) 
{
    return Vec4f_Make( a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w );
}

Vec4f Vec4f_Mult(Vec4f a, float s) 
{
    return Vec4f_Make( a.x * s, a.y * s, a.z * s, a.w * s );
}

float Vec4f_Length(Vec4f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y) + (v.z * v.z) + (v.w * v.w));
}

Vec4f Vec4f_NormalizeFrom(Vec4f v)
{
    float Ratio = 1.0 / Vec4f_Length(v);
    Vec4f Result = Vec4f_Mult(v, Ratio);
    return Result;
}

void Vec4f_Normalize(Vec4f* v)
{
    Vec4f Normalized = Vec4f_NormalizeFrom(*v);
    *v = Normalized;
}

float Vec4f_Dot(Vec4f a, Vec4f b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

Vec4f Vec4f_FromVec3f(Vec3f v, float w)
{
    return Vec4f_Make(v.x, v.y, v.z, w);
}

/**********************************
 *         MATRIX TYPES           *
 **********************************/

typedef struct Mat33 {
    float m[3][3]; //ROW MAJOR
} Mat33;

typedef struct Mat44 {
    float m[4][4]; //RAW MAJOR
} Mat44;

/**********************************
 *        3x3 MATRIX              *
 **********************************/

Mat33 Mat33_MakeEmpty()
{
    Mat33 Result = {.m = {{0, 0, 0},
                          {0, 0, 0},
                          {0, 0, 0}}}; 
    
    return Result;
}

Mat33 Mat33_Make(float m00,float m01,float m02,
                 float m10,float m11,float m12,
                 float m20,float m21,float m22 )
{
    Mat33 Result = {.m = {{ m00, m01, m02,},
                          { m10, m11, m12,},
                          { m20, m21, m22 } }}; 
    return Result;
}

Mat33 Mat33_MakeId()
{
    Mat33 Result = Mat33_Make(1, 0, 0,
                              0, 1, 0,
                              0, 0, 1 ); 
    return Result;
}

Mat33 Mat33_Add(Mat33 a, Mat33 b)
{
    Mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] + b.m[j][i];
        }
    }
    return Result;
}
#define Mat33_Add3(a, b, c)         Mat33_Add(Mat33_Add(a, b), c)
#define Mat33_Add4(a, b, c, d)      Mat33_Add(Mat33_Add3(a, b, c) d)
#define Mat33_Add5(a, b, c, d, e)   Mat33_Add(Mat33_Add4(a, b, c, d), e)

Mat33 Mat33_Substract(Mat33 a, Mat33 b)
{
    Mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] - b.m[j][i];
        }
    }
    return Result;
}

Mat33 Mat33_MultScalar(Mat33 a, float s) 
{
    Mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

Mat33 Mat33_MultMat(Mat33 a, Mat33 b) 
{
    Mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = (a.m[j][0] * b.m[0][i]) + (a.m[j][1] * b.m[1][i]) + (a.m[j][2] * b.m[2][i]);
        }
    }
    return Result;
}

Vec3f Mat33_MultVec(Mat33 a, Vec3f v) 
{
    Vec3f Result;
    Result.x = (a.m[0][0] * v.x) + (a.m[0][1] * v.y) + (a.m[0][2] * v.z);
    Result.y = (a.m[1][0] * v.x) + (a.m[1][1] * v.y) + (a.m[1][2] * v.z);
    Result.z = (a.m[2][0] * v.x) + (a.m[2][1] * v.y) + (a.m[2][2] * v.z);
    return Result;
}

bool Mat33_Equal(Mat33 a, Mat33 b)
{
    float* m1 = (float*)a.m;
    float* m2 = (float*)b.m;
    for (int i = 0; i < 9; ++i) {
        if (!float_Equal(m1[i], m2[i])) {
            return false;
        }
    }
    return true;
}


Mat33 Mat33_Transpose(Mat33 m)
{
    Mat33 Result = Mat33_Make( m.m[0][0], m.m[1][0], m.m[2][0],
                               m.m[0][1], m.m[1][1], m.m[2][1],
                               m.m[0][2], m.m[1][2], m.m[2][2] );
    return Result;
}


Mat33 Mat33_MakeScale(float sx, float sy, float sz)
{
    Mat33 Result = Mat33_Make(sx,  0,  0,
                               0, sy,  0,
                               0,  0, sz ); 
    return Result;
}

Mat33 Mat33_MakeRotation(Vec3f Axis, float Angle)
{

    Assert( float_Equal(Vec3f_Length(Axis), 1.0));

    float x2 = Axis.x * Axis.x;
    float y2 = Axis.y * Axis.y;
    float z2 = Axis.z * Axis.z;

    float S = Sinf(Angle);
    float C = Cosf(Angle);

    Mat33 I = Mat33_MakeId();
    Mat33 J = Mat33_Make(       0, -Axis.z,  Axis.y,
                           Axis.z,       0, -Axis.x,
                          -Axis.y,  Axis.x,       0 );

    Mat33 J2 = Mat33_Make(        -y2 - z2, Axis.x * Axis.y, Axis.z * Axis.x,
                           Axis.x * Axis.y,        -z2 - x2, Axis.y * Axis.z,
                           Axis.z * Axis.x, Axis.y * Axis.z,        -x2 - y2 );
    
    Mat33 Result = Mat33_Add3(Mat33_MultScalar(J, S), Mat33_MultScalar(J2, (1 - C)), I);
    return Result;
}

float Mat33_Trace(Mat33 m)
{
    return m.m[0][0] + m.m[1][1] + m.m[2][2];
}

/**********************************
 *        4x4 MATRIX              *
 **********************************/

Mat44 Mat44_MakeEmpty()
{
    Mat44 Result = {.m = {{0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0} }}; 
    return Result;
}

Mat44 Mat44_Make(float Mat00, float Mat01, float Mat02, float Mat03,
                 float Mat10, float Mat11, float Mat12, float Mat13,
                 float Mat20, float Mat21, float Mat22, float Mat23,
                 float Mat30, float Mat31, float Mat32, float Mat33 )
{
     Mat44 Result = {.m = {{ Mat00,  Mat01,  Mat02,  Mat03},
                           { Mat10,  Mat11,  Mat12,  Mat13},
                           { Mat20,  Mat21,  Mat22,  Mat23},
                           { Mat30,  Mat31,  Mat32,  Mat33} }}; 
    return Result;   
}


Mat44 Mat44_MakeId()
{
    Mat44 Result = Mat44_Make(1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

Mat44 Mat44_FromMat33(Mat33 m)
{
    Mat44 Result = Mat44_Make(m.m[0][0], m.m[0][1], m.m[0][2], 0,
                              m.m[1][0], m.m[1][1], m.m[1][2], 0,
                              m.m[2][0], m.m[2][1], m.m[2][2], 0,
                                      0,         0,         0, 1 ); 
    return Result;
}

Mat44 Mat44_Add(Mat44 a, Mat44 b)
{
    Mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] + b.m[j][i];
        }
    }
    return Result;
}

Mat44 Mat44_Substract(Mat44 a, Mat44 b)
{
    Mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] - b.m[j][i];
        }
    }
    return Result;
}

Mat44 Mat44_MultScalar(Mat44 a, float s) 
{
    Mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

Mat44 Mat44_MultMat(Mat44 a, Mat44 b) 
{
    Mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = (a.m[j][0] * b.m[0][i]) + (a.m[j][1] * b.m[1][i]) + (a.m[j][2] * b.m[2][i]) + (a.m[j][3] * b.m[3][i]);
        }
    }
    return Result;
}

Vec4f Mat44_MultVec(Mat44 a, Vec4f v) 
{
    Vec4f Result;
    Result.x = (a.m[0][0] * v.x) + (a.m[0][1] * v.y) + (a.m[0][2] * v.z) + (a.m[0][3] * v.w);
    Result.y = (a.m[1][0] * v.x) + (a.m[1][1] * v.y) + (a.m[1][2] * v.z) + (a.m[1][3] * v.w);
    Result.z = (a.m[2][0] * v.x) + (a.m[2][1] * v.y) + (a.m[2][2] * v.z) + (a.m[2][3] * v.w);
    Result.w = (a.m[3][0] * v.x) + (a.m[3][1] * v.y) + (a.m[3][2] * v.z) + (a.m[3][3] * v.w);

    return Result;
}

Mat44 Mat44_MakeTranspose(Mat44 m) 
{
    Mat44 Result = Mat44_Make(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
                              m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
                              m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
                              m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3] ); 
    return Result;
}

void Mat44_Transpose(Mat44 *m) {

    Mat44 Result = {.m = {{m->m[0][0], m->m[1][0], m->m[2][0], m->m[3][0]},
                          {m->m[0][1], m->m[1][1], m->m[2][1], m->m[3][1]},
                          {m->m[0][2], m->m[1][2], m->m[2][2], m->m[3][2]},
                          {m->m[0][3], m->m[1][3], m->m[2][3], m->m[3][3]} }}; 
    *m = Result;
}

float Mat44_Trace(Mat44 m)
{
    return m.m[0][0] + m.m[1][1] + m.m[2][2] + m.m[3][3];
}


Mat44 Mat44_MakeScale(float sx, float sy, float sz)
{
    Mat44 Result = Mat44_Make(sx, 0,  0, 0,
                              0, sy,  0, 0,
                              0,  0, sz, 0,
                              0,  0,  0, 1 ); 
    return Result;
}

Mat44 Mat44_MakeTranslate(float dx, float dy, float dz)
{
    Mat44 Result = Mat44_Make(1, 0, 0, dx,
                              0, 1, 0, dy,
                              0, 0, 1, dz,
                              0, 0, 0,  1 ); 
    return Result;
}

Mat44 Mat44_MakeRotationX(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    Mat44 Result = Mat44_Make(1, 0, 0, 0,
                              0, C,-S, 0,
                              0, S, C, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

Mat44 Mat44_MakeRotationY(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    Mat44 Result = Mat44_Make( C, 0, S, 0,
                               0, 1, 0, 0,
                              -S, 0, C, 0,
                               0, 0, 0, 1 ); 
    return Result;
}

Mat44 Mat44_MakeRotationZ(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    Mat44 Result = Mat44_Make(C,-S, 0, 0,
                              S, C, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

Mat44 Mat44_MakeRotation(Vec3f Axis, float Angle)
{
    Mat44 Result = Mat44_FromMat33(Mat33_MakeRotation(Axis, Angle));
    return Result;
}

Mat44 Mat44_LookAt(Vec3f Position, Vec3f LookPosition, Vec3f Up)
{

    Vec3f W = Vec3f_NormalizeFrom(Vec3f_Substract(Position, LookPosition)); //Forward
    Vec3f V = Vec3f_NormalizeFrom(Vec3f_Substract(Up, Vec3f_Mult(W, Vec3f_Dot(Up, W)))); //Up
    Vec3f U = Vec3f_Cross(V, W); //right

    Mat44 M = Mat44_Make(U.x, V.x, W.x, -Position.x,
                         U.y, V.y, W.y, -Position.y,
                         U.z, V.z, W.z, -Position.z,
                           0,   0,   0, 1 ); 
    return M;
}

Mat44 Mat44_Perspective(float Fov, float Near, float Far, float Width, float Height)
{
    float Scale = Near * Tanf(Fov * 0.5 * (M_PI / 180));
    float AspectRatio = Width / Height;
    float r = AspectRatio * Scale, l = -r;
    float t = Scale, b = -t;

    Mat44 Result = Mat44_MakeEmpty();

    Result.m[0][0] = (2 * Near) / (r - l);
    Result.m[1][1] = (2 * Near) / (t - b);
    Result.m[3][2] = -1;

    Result.m[2][2] = - (Far + Near) / (Far - Near);
    Result.m[2][3] = - ( 2 * Far * Near) / (Far - Near);

    return Result;
}



/**********************************
 *            Quaternion          *
 **********************************/

typedef struct Quaternion {
    float r; 
    union {
        struct {
            float x, y, z;
        };
        Vec3f v;
    };
} Quaternion;

Quaternion Quaternion_MakeReal(float R, float X, float Y, float Z)
{
    return (Quaternion){.r = R, .x = X, .y = Y, .z = Z};
}

Quaternion Quaternion_MakeRealVec(float R, Vec3f V)
{
    return Quaternion_MakeReal(R, V.x, V.y, V.z);
}

Mat33 QuaternionToRotation(Quaternion);

#define Quaternion_Identity Quaternion_MakeReal(1, 0, 0, 0)


Quaternion Quaternion_FromAxisAngle(float Angle, Vec3f Axis )
{
    return Quaternion_MakeRealVec(Cosf(Angle / 2), Vec3f_Mult(Axis, Sinf(Angle / 2)));
}

Quaternion Quaternion_Invert(Quaternion q)
{
    Quaternion Result = Quaternion_MakeReal(-q.r, -q.x, -q.y, -q.z);
    return Result;
}

Quaternion Quaternion_Add(Quaternion q, Quaternion r)
{
    return Quaternion_MakeRealVec(q.r + r.r, Vec3f_Add(q.v, r.v));
}

Quaternion Quaternion_Substract(Quaternion q, Quaternion r)
{
    return Quaternion_MakeRealVec(q.r - r.r, Vec3f_Substract(q.v, r.v));
}

Quaternion Quaternion_MultQuaternion(Quaternion q1, Quaternion q2)
{
    float r = (q1.r * q2.r) - Vec3f_Dot(q1.v, q2.v);
    Vec3f v = Vec3f_Add3(Vec3f_Cross(q1.v, q2.v), Vec3f_Mult(q1.v, q2.r), Vec3f_Mult(q2.v, q1.r));
    return Quaternion_MakeRealVec(r,v);
}

Quaternion Quaternion_MultScalar(Quaternion q, float s)
{
    return Quaternion_MakeRealVec(q.r * s, Vec3f_Mult(q.v, s)); 
}


Vec3f Quaternion_MultVec(Quaternion q, Vec3f v)
{
    Mat33 Rotation = QuaternionToRotation(q);
    return Mat33_MultVec(Rotation, v);
}

float Quaternion_Length(Quaternion q)
{
    return SquareRootf((q.r * q.r) + (q.x * q.x) + (q.y * q.y) + (q.z * q.z));
}

Quaternion Quaternion_NormalizeFrom(Quaternion q)
{
    return Quaternion_MultScalar(q, (1 / Quaternion_Length(q)));
}

Quaternion Quaternion_ConjugateFrom(Quaternion q)
{
    return Quaternion_MakeRealVec(q.r, Vec3f_MakeInvert(q.v));
}

Quaternion Quaternion_InverseFrom(Quaternion q)
{
    return Quaternion_NormalizeFrom(Quaternion_ConjugateFrom(q));
}

float Quaternion_Dot(Quaternion q, Quaternion r)
{
   return (q.r * r.r) + (q.x * r.x) + (q.y * r.y) + (q.z * r.z);
}


/**********************************
 *      Non classified Op         *
 **********************************/

typedef struct rotation_axis_angle {
    Vec3f Axis;
    float Angle;
} rotation_axis_angle;

rotation_axis_angle Mat33_AxisAngleFromRotation(Mat33 Rotation)
{
    float Angle = Acosf((Mat33_Trace(Rotation) - 1) / 2);
    if( Angle < F_EPS ) { //near zero
        return (rotation_axis_angle) {.Axis = Vec3f_Make(1, 0, 0), .Angle = Angle};      
    }

    if (float_Equal(Angle, M_PI)) { //near PI
        int Column = 0;
        float Max = 0;
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                float v = Rotation.m[j][i];
                if (v > Max) {
                    Max = v;
                    Column = i;
                }
            }
        }

        return (rotation_axis_angle){.Axis = Vec3f_Make(Rotation.m[0][Column], Rotation.m[1][Column], Rotation.m[2][Column]), .Angle = Angle};
    }

    Mat33 s = Mat33_Substract(Rotation, Mat33_Transpose(Rotation));
    Vec3f v = Vec3f_Make(-s.m[1][2], s.m[0][2], -s.m[1][1]);
    float t = 1.0 / (2 * Sinf(Angle));

    return (rotation_axis_angle) {.Axis = Vec3f_Mult(v, t), .Angle = Angle};
}

Quaternion Slerp(Quaternion Start, Quaternion End, float t)
{
    Assert(float_Equal(Quaternion_Dot(Start, Start), 1.0));
    Assert(float_Equal(Quaternion_Dot(End, End), 1.0));
    Assert(t >= 0.0);
    Assert(t <= 1.0);

    Quaternion q = Quaternion_Substract(End, 
                                        Quaternion_MultScalar(Start, Quaternion_Dot(Start, End))
                                        );
    Quaternion u = Quaternion_NormalizeFrom(q);
    float Angle = Acosf(Quaternion_Dot(Start, End));
    return Quaternion_Add( Quaternion_MultScalar(Start, Cosf(t * Angle)), 
                           Quaternion_MultScalar(u, Sinf(t * Angle))
                           );
}

bool Mat33_IsIdentity(Mat33 m)
{
    return m.m[0][0] == 1 && m.m[0][1] == 0 && m.m[0][2] == 0 &&
           m.m[1][0] == 0 && m.m[1][1] == 1 && m.m[1][2] == 0 &&
           m.m[2][0] == 0 && m.m[2][1] == 0 && m.m[2][2] == 1;
}

Mat33 QuaternionToRotation(Quaternion q)
{
    float r2 = q.r * q.r;
    float x2 = q.x * q.x;
    float y2 = q.y * q.y;
    float z2 = q.z * q.z;

    float rx = q.r * q.x;
    float ry = q.r * q.y;
    float rz = q.r * q.z;
    float xy = q.x * q.y;
    float xz = q.x * q.z;
    float yz = q.y * q.z;

    Mat33 Result = Mat33_Make( r2 + x2 - y2 - z2,     2 * (yz - rx),     2 * (ry + xz) ,
                                   2 * (rz + xy), r2 - x2 + y2 - z2,     2 * (yz - rx) ,
                                   2 * (xz - ry),     2 * (rx + yz), r2 - x2 - y2 + z2  );
    return Result;
}

typedef struct quaternion_pair {
    Quaternion q1, q2;
} quaternion_pair;

quaternion_pair RotationToQuaternion(Mat33 m)
{
    if (Mat33_IsIdentity(m)) {
        Quaternion q = Quaternion_MakeReal(1, 0, 0, 0);
        return (quaternion_pair){.q1 = q, .q2 = Quaternion_Invert(q)};
    }

    rotation_axis_angle AxisAngle = Mat33_AxisAngleFromRotation(m);
    Quaternion q1 = Quaternion_FromAxisAngle(AxisAngle.Angle, AxisAngle.Axis);
    
    return (quaternion_pair){.q1 = q1, .q2 = Quaternion_Invert(q1)};
}

Mat33 RotationInterpolation(Mat33 Start, Mat33 End, float t)
{
    // m1 * m2T == -I
    
    Mat33 mI = Mat33_Make(-1, 0, 0, 0, -1, 0, 0, 0, -1);
    Mat33 Cmp = Mat33_MultMat(Start, Mat33_Transpose(End));
    Assert(!Mat33_Equal(Cmp, mI));
    
    quaternion_pair p1 = RotationToQuaternion(Start);
    quaternion_pair p2 = RotationToQuaternion(End);

    Quaternion q1 = p1.q1;
    Quaternion q2 = p2.q1;
   
    if (Quaternion_Dot(q1, q2) < 0) q2 = p2.q2;
    Quaternion q = Slerp(q1, q2, t);
    return QuaternionToRotation(q);
}

typedef struct euler_angles {
    float Pitch, Yaw, Roll;
} euler_angles;

euler_angles EulerAnglesFromRotation(Mat44 Rotation)
{
    float Pitch, Yaw, Roll;

    Yaw = Asinf(Rotation.m[0][2]);
    if(Abs(Cosf(Yaw)) == 0) {

        Roll = 0;
        Pitch = Atan2f(Rotation.m[2][1], Rotation.m[1][1]);

    } else {

        Roll  = Atan2f(-Rotation.m[0][1], Rotation.m[0][0]);
        Pitch = Atan2f( Rotation.m[1][2], Rotation.m[2][2]);

    }
    euler_angles Result = {Pitch, Yaw, Roll};;
    return Result;
}

#endif // ALGEBRA_H
