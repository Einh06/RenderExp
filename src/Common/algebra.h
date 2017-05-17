#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "common.h"
#include "math.h"

/**********************************
 *          2D VECTOR             *
 **********************************/
struct vec2f {
    float x;
    float y;

    vec2f() : x(0), y(0) {}
    vec2f(float X, float Y) : x(X), y(Y) {}
};

#define vec2f_X vec2f(1, 0);
#define vec2f_Y vec2f(0, 1);

vec2f operator+(vec2f a, vec2f b) 
{
    return { a.x + b.x, a.y + b.y };
}

vec2f operator-(vec2f a, vec2f b) 
{
    return { a.x - b.x, a.y - b.y };
}

vec2f operator-(vec2f v)
{
    return { -v.x, -v.y };
}

vec2f operator*(vec2f a, float s) 
{
    return {a.x * s, a.y * s };
}

vec2f operator*(float s, vec2f a) 
{
    return {a.x * s, a.y * s };
}

float vec2f_Length(vec2f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y));
}

vec2f vec2f_NormalizeFrom(vec2f v)
{
    float Length = vec2f_Length(v);
    vec2f Result = { v.x / Length, v.y/Length};
    return Result;
}

void vec2f_Normalize(vec2f* v)
{
    vec2f Normalized = vec2f_NormalizeFrom(*v);
    *v = Normalized;
}

float vec2f_Dot(vec2f a, vec2f b)
{
    return (a.x * b.x) + (a.y * b.y);
}

vec2f vec2f_Cross(vec2f a)
{
    vec2f Result = { -a.y, a.x }; 
    return Result;
}




/**********************************
 *          3D VECTOR             *
 **********************************/
struct vec3f {
    float x, y, z;
    vec3f() : x(0), y(0), z(0) {}
    vec3f(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
};

#define vec3f_X vec3f(1, 0, 0)
#define vec3f_Y vec3f(0, 1, 0)
#define vec3f_Z vec3f(0, 0, 1)

vec3f operator+(vec3f a, vec3f b) 
{
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

vec3f operator-(vec3f v)
{
    return { -v.x, -v.y, v.z};
}

vec3f operator-(vec3f a, vec3f b) 
{
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

vec3f operator*(vec3f a, float s) 
{
    return { a.x * s, a.y * s, a.z * s };
}

vec3f operator*(float s, vec3f a) 
{
    return { a.x * s, a.y * s, a.z * s };
}

float vec3f_Length(vec3f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

vec3f vec3f_NormalizeFrom(vec3f v)
{
    float Length = vec3f_Length(v);
    vec3f Result = { v.x / Length, v.y/Length, v.z/Length};
    return Result;
}

void vec3f_Normalize(vec3f* v)
{
    vec3f Normalized = vec3f_NormalizeFrom(*v);
    *v = Normalized;
}

float vec3f_Dot(vec3f a, vec3f b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

vec3f vec3f_Cross(vec3f a, vec3f b)
{
    vec3f Result;

    Result.x = (a.y * b.z) - (a.z * b.y);
    Result.y = (a.z * b.x) - (a.x * b.z);
    Result.z = (a.x * b.y) - (a.y * b.x);

    return Result;
}

vec3f vec3f_FromVec2f(vec2f v)
{
    return {v.x, v.y, 1};
}



/**********************************
 *          4D VECTOR             *
 **********************************/
struct vec4f {
    float x;
    float y;
    float z;
    float w;

    vec4f() : x(0), y(0), z(0), w(0) {}
    vec4f(float X, float Y, float Z, float W): x(X), y(Y), z(Z), w(W) {}
};

vec4f operator+(vec4f a, vec4f b) 
{
    return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

vec4f operator-(vec4f v)
{
    return { -v.x, -v.y, -v.z, -v.w };
}

vec4f operator-(vec4f a, vec4f b) 
{
    return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
}

vec4f operator*(vec4f a, float s) 
{
    return { a.x * s, a.y * s, a.z * s, a.w * s };
}

vec4f operator*(float s, vec4f a) 
{
    return { a.x * s, a.y * s, a.z * s, a.w * s };
}

float vec4f_Length(vec4f v)
{
    return SquareRootf((v.x * v.x) + (v.y * v.y) + (v.z * v.z) + (v.w * v.w));
}

vec4f vec4f_NormalizeFrom(vec4f v)
{
    float Ratio = 1.0 / vec4f_Length(v);
    vec4f Result = Ratio * v;
    return Result;
}

void vec4f_Normalize(vec4f* v)
{
    vec4f Normalized = vec4f_NormalizeFrom(*v);
    *v = Normalized;
}

float vec4f_Dot(vec4f a, vec4f b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

vec4f vec4f_FromVec3f(vec3f v, float w)
{
    return {v.x, v.y, v.z, w};
}



vec3f vec3f_FromVec4f(vec4f v) 
{
    vec3f Result = {v.x, v.y, v.z}; 
    float Divider = v.w != 0 ? 1.0/v.w : 1.0;
    return Result * Divider;
}
/**********************************
 *        3x3 MATRIX              *
 **********************************/
struct mat33 {
    float m[3][3]; //ROW MAJOR
};

mat33 operator+(mat33 a, mat33 b)
{
    mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] + b.m[j][i];
        }
    }
    return Result;
}

mat33 operator-(mat33 a, mat33 b)
{
    mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] - b.m[j][i];
        }
    }
    return Result;
}

mat33 operator*(mat33 a, float s) 
{
    mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

mat33 operator*(float s, mat33 a) 
{
    mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

mat33 operator*(mat33 a, mat33 b) 
{
    mat33 Result;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            Result.m[j][i] = (a.m[j][0] * b.m[0][i]) + (a.m[j][1] * b.m[1][i]) + (a.m[j][2] * b.m[2][i]);
        }
    }
    return Result;
}

vec3f operator*(mat33 a, vec3f v) 
{
    vec3f Result;
    Result.x = (a.m[0][0] * v.x) + (a.m[0][1] * v.y) + (a.m[0][2] * v.z);
    Result.y = (a.m[1][0] * v.x) + (a.m[1][1] * v.y) + (a.m[1][2] * v.z);
    Result.z = (a.m[2][0] * v.x) + (a.m[2][1] * v.y) + (a.m[2][2] * v.z);
    return Result;
}

bool mat33_Equal(mat33 a, mat33 b)
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

mat33 mat33_Make()
{
    mat33 Result = {.m = {{0, 0, 0},
                          {0, 0, 0},
                          {0, 0, 0}}}; 
    
    return Result;
}

mat33 mat33_Make(float m00,float m01,float m02,
                 float m10,float m11,float m12,
                 float m20,float m21,float m22 )
{
    mat33 Result = {.m = {{ m00, m01, m02,},
                          { m10, m11, m12,},
                          { m20, m21, m22 } }}; 
    return Result;
}

mat33 mat33_MakeId()
{
    mat33 Result = {.m = {{1, 0, 0},
                          {0, 1, 0},
                          {0, 0, 1}}}; 
    return Result;
}

mat33 mat33_Transpose(mat33 m)
{
    mat33 Result = mat33_Make( m.m[0][0], m.m[1][0], m.m[2][0],
                               m.m[0][1], m.m[1][1], m.m[2][1],
                               m.m[0][2], m.m[1][2], m.m[2][2] );
    return Result;
}


mat33 mat33_MakeScale(float sx, float sy, float sz)
{
    mat33 Result = mat33_Make(sx,  0,  0,
                               0, sy,  0,
                               0,  0, sz ); 
    return Result;
}

mat33 mat33_MakeRotation(vec3f Axis, float Angle)
{

    Assert( float_Equal(vec3f_Length(Axis), 1.0));

    float x2 = Axis.x * Axis.x;
    float y2 = Axis.y * Axis.y;
    float z2 = Axis.z * Axis.z;

    float S = Sinf(Angle);
    float C = Cosf(Angle);

    mat33 I = mat33_MakeId();
    mat33 J = mat33_Make(       0, -Axis.z,  Axis.y,
                           Axis.z,       0, -Axis.x,
                          -Axis.y,  Axis.x,       0 );

    mat33 J2 = mat33_Make(        -y2 - z2, Axis.x * Axis.y, Axis.z * Axis.x,
                           Axis.x * Axis.y,        -z2 - x2, Axis.y * Axis.z,
                           Axis.z * Axis.x, Axis.y * Axis.z,        -x2 - y2 );
    
    mat33 Result = (S * J) + ((1 - C) * J2) + I;
    return Result;
}

float mat33_Trace(mat33 m)
{
    return m.m[0][0] + m.m[1][1] + m.m[2][2];
}

/**********************************
 *        4x4 MATRIX              *
 **********************************/
struct mat44 {
    float m[4][4]; //RAW MAJOR
};

mat44 mat44_Make()
{
    mat44 Result = {.m = {{0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 0, 0} }}; 
    return Result;
}

mat44 mat44_Make(float mat00, float mat01, float mat02, float mat03,
                 float mat10, float mat11, float mat12, float mat13,
                 float mat20, float mat21, float mat22, float mat23,
                 float mat30, float mat31, float mat32, float mat33 )
{
     mat44 Result = {.m = {{ mat00,  mat01,  mat02,  mat03},
                           { mat10,  mat11,  mat12,  mat13},
                           { mat20,  mat21,  mat22,  mat23},
                           { mat30,  mat31,  mat32,  mat33} }}; 
    return Result;   
}


mat44 mat44_MakeId()
{
    mat44 Result = mat44_Make(1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

mat44 mat44_FromMat33(mat33 m)
{
    mat44 Result = mat44_Make(m.m[0][0], m.m[0][1], m.m[0][2], 0,
                              m.m[1][0], m.m[1][1], m.m[1][2], 0,
                              m.m[2][0], m.m[2][1], m.m[2][2], 0,
                                      0,         0,         0, 1 ); 
    return Result;
}

mat44 operator+(mat44 a, mat44 b)
{
    mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] + b.m[j][i];
        }
    }
    return Result;
}

mat44 operator-(mat44 a, mat44 b)
{
    mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] - b.m[j][i];
        }
    }
    return Result;
}

mat44 operator*(mat44 a, float s) 
{
    mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

mat44 operator*(float s, mat44 a) 
{
    mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = a.m[j][i] * s;
        }
    }
    return Result;
}

mat44 operator*(mat44 a, mat44 b) 
{
    mat44 Result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            Result.m[j][i] = (a.m[j][0] * b.m[0][i]) + (a.m[j][1] * b.m[1][i]) + (a.m[j][2] * b.m[2][i]) + (a.m[j][3] * b.m[3][i]);
        }
    }
    return Result;
}

vec4f operator*(mat44 a, vec4f v) 
{
    vec4f Result;
    Result.x = (a.m[0][0] * v.x) + (a.m[0][1] * v.y) + (a.m[0][2] * v.z) + (a.m[0][3] * v.w);
    Result.y = (a.m[1][0] * v.x) + (a.m[1][1] * v.y) + (a.m[1][2] * v.z) + (a.m[1][3] * v.w);
    Result.z = (a.m[2][0] * v.x) + (a.m[2][1] * v.y) + (a.m[2][2] * v.z) + (a.m[2][3] * v.w);
    Result.w = (a.m[3][0] * v.x) + (a.m[3][1] * v.y) + (a.m[3][2] * v.z) + (a.m[3][3] * v.w);

    return Result;
}

mat44 mat44_MakeTranspose(mat44 m) 
{
    mat44 Result = mat44_Make(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
                              m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
                              m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
                              m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3] ); 
    return Result;
}

void mat44_Transpose(mat44 *m) {

    mat44 Result = {.m = {{m->m[0][0], m->m[1][0], m->m[2][0], m->m[3][0]},
                          {m->m[0][1], m->m[1][1], m->m[2][1], m->m[3][1]},
                          {m->m[0][2], m->m[1][2], m->m[2][2], m->m[3][2]},
                          {m->m[0][3], m->m[1][3], m->m[2][3], m->m[3][3]} }}; 
    *m = Result;
}

float mat44_Trace(mat44 m)
{
    return m.m[0][0] + m.m[1][1] + m.m[2][2] + m.m[3][3];
}


mat44 mat44_MakeScale(float sx, float sy, float sz)
{
    mat44 Result = mat44_Make(sx, 0,  0, 0,
                              0, sy,  0, 0,
                              0,  0, sz, 0,
                              0,  0,  0, 1 ); 
    return Result;
}

mat44 mat44_MakeTranslate(float dx, float dy, float dz)
{
    mat44 Result = mat44_Make(1, 0, 0, dx,
                              0, 1, 0, dy,
                              0, 0, 1, dz,
                              0, 0, 0,  1 ); 
    return Result;
}

mat44 mat44_MakeRotationX(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    mat44 Result = mat44_Make(1, 0, 0, 0,
                              0, C,-S, 0,
                              0, S, C, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

mat44 mat44_MakeRotationY(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    mat44 Result = mat44_Make( C, 0, S, 0,
                               0, 1, 0, 0,
                              -S, 0, C, 0,
                               0, 0, 0, 1 ); 
    return Result;
}

mat44 mat44_MakeRotationZ(float Angle) 
{
    float C = Cosf(Angle);
    float S = Sinf(Angle);
    mat44 Result = mat44_Make(C,-S, 0, 0,
                              S, C, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1 ); 
    return Result;
}

mat44 mat44_MakeRotation(vec3f Axis, float Angle)
{
    mat44 Result = mat44_FromMat33(mat33_MakeRotation(Axis, Angle));
    return Result;
}

mat44 mat44_LookAt(vec3f Position, vec3f LookPosition, vec3f Up)
{

    vec3f W = vec3f_NormalizeFrom(Position - LookPosition); //Forward
    vec3f V = vec3f_NormalizeFrom(Up - ((vec3f_Dot(Up, W) * W))); //Left
    vec3f U = vec3f_Cross(V, W); //UP

    mat44 M = mat44_Make(U.x, V.x, W.x, -Position.x,
                         U.y, V.y, W.y, -Position.y,
                         U.z, V.z, W.z, -Position.z,
                           0,   0,   0, 1 ); 
    return M;
}

mat44 mat44_Perspective(float Fov, float Near, float Far, float Width, float Height)
{
    float Scale = Near * Tanf(Fov * 0.5 * (M_PI / 180));
    float AspectRatio = Width / Height;
    float r = AspectRatio * Scale, l = -r;
    float t = Scale, b = -t;

    mat44 Result = mat44_Make();

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
struct quaternion {
    float r; 
    union {
        struct {
            float x, y, z;
        };
        vec3f v;
    };
    quaternion() : r(0), x(0), y(0), z(0) {}
    quaternion(float R, float X, float Y, float Z) : r(R), x(X), y(Y), z(Z) {}
};

#define quaternion_Identity quaternion(1, 0, 0, 0)

quaternion quaternion_MakeFromRealVec(float r, vec3f v)
{
    return quaternion(r, v.x, v.y, v.z);
}

quaternion operator-(quaternion q)
{
    quaternion Result = quaternion(-q.r, -q.x, -q.y, -q.z);
    return Result;
}

quaternion operator+(quaternion q, quaternion r)
{
    return quaternion_MakeFromRealVec(q.r + r.r, q.v + r.v);
}

quaternion operator-(quaternion q, quaternion r)
{
    return quaternion_MakeFromRealVec(q.r - r.r, q.v - r.v);
}

quaternion operator*(quaternion q1, quaternion q2)
{
    float r = (q1.r * q2.r) - vec3f_Dot(q1.v, q2.v);
    vec3f v = vec3f_Cross(q1.v, q2.v) + (q2.r * q1.v) + (q1.r * q2.v);
    return quaternion_MakeFromRealVec(r,v);
}

quaternion operator*(quaternion q, float s)
{
    return quaternion_MakeFromRealVec(q.r * s, q.v * s); 
}

quaternion operator*(float s, quaternion q)
{
    return quaternion_MakeFromRealVec(q.r * s, q.v * s); 
}

float quaternion_Length(quaternion q)
{
    return SquareRootf((q.r * q.r) + (q.x * q.x) + (q.y * q.y) + (q.z * q.z));
}

quaternion quaternion_NormalizeFrom(quaternion q)
{
    return q * (1 / quaternion_Length(q));
}

quaternion quaternion_ConjugateFrom(quaternion q)
{
    return quaternion_MakeFromRealVec(q.r, - q.v);
}

quaternion quaternion_InverseFrom(quaternion q)
{
    return quaternion_NormalizeFrom(quaternion_ConjugateFrom(q));
}

float quaternion_Dot(quaternion q, quaternion r)
{
   return (q.r * r.r) + (q.x * r.x) + (q.y * r.y) + (q.z * r.z);
}


/**********************************
 *      Non classified Op         *
 **********************************/

struct rotation_axis_angle {
    vec3f Axis;
    float Angle;
};

rotation_axis_angle mat33_AxisAngleFromRotation(mat33 Rotation)
{
    float Angle = Acosf((mat33_Trace(Rotation) - 1) / 2);
    if( Angle < F_EPS ) { //near zero
        return {.Axis = vec3f(1, 0, 0), .Angle = Angle};      
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

        return {.Axis = vec3f(Rotation.m[0][Column], Rotation.m[1][Column], Rotation.m[2][Column]), .Angle = Angle};
    }

    mat33 s = Rotation - mat33_Transpose(Rotation);
    float x = -s.m[1][2], y = s.m[0][2], z = -s.m[1][1];
    float t = 2 * Sinf(Angle);
    return {.Axis = vec3f(x/t, y/t, z/t), .Angle = Angle};
}

quaternion Slerp(quaternion Start, quaternion End, float t)
{
    Assert(float_Equal(quaternion_Dot(Start, Start), 1.0));
    Assert(float_Equal(quaternion_Dot(End, End), 1.0));
    Assert(t >= 0.0);
    Assert(t <= 1.0);

    quaternion u = End - (quaternion_Dot(Start, End) * Start);
    u = (1 / quaternion_Length(u)) * u; //Normalize
    float Angle = Acosf(quaternion_Dot(Start, End));
    return (Cosf(t * Angle) * Start) + (Sinf(t * Angle) * u);
}

bool mat33_IsIdentity(mat33 m)
{
    return m.m[0][0] == 1 && m.m[0][1] == 0 && m.m[0][2] == 0 &&
           m.m[1][0] == 0 && m.m[1][1] == 1 && m.m[1][2] == 0 &&
           m.m[2][0] == 0 && m.m[2][1] == 0 && m.m[2][2] == 1;
}

mat33 QuaternionToRotation(quaternion q)
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

    mat33 Result = mat33_Make( r2 + x2 - y2 - z2,     2 * (yz - rx),     2 * (ry + xz) ,
                                   2 * (rz + xy), r2 - x2 + y2 - z2,     2 * (yz - rx) ,
                                   2 * (xz - ry),     2 * (rx + yz), r2 - x2 - y2 + z2  );
    return Result;
}

struct quaternion_pair {
    quaternion q1, q2;
};

quaternion_pair RotationToQuaternion(mat33 m)
{
    if (mat33_IsIdentity(m)) {
        quaternion q = quaternion(1, 0, 0, 0);
        return {.q1 = q, .q2 = -q};
    }

    rotation_axis_angle AxisAngle = mat33_AxisAngleFromRotation(m);
    quaternion q1 = quaternion_MakeFromRealVec(Cosf(AxisAngle.Angle / 2), Sinf(AxisAngle.Angle / 2) * AxisAngle.Axis);
    
    return {.q1 = q1, .q2 = -q1};
}

mat33 RotationInterpolation(mat33 Start, mat33 End, float t)
{
    // m1 * m2T == -I
    Assert(!mat33_Equal(Start * mat33_Transpose(End), -1 * mat33_MakeId()));
    
    quaternion_pair p1 = RotationToQuaternion(Start);
    quaternion_pair p2 = RotationToQuaternion(End);

    quaternion q1 = p1.q1;
    quaternion q2 = p2.q1;
   
    if (quaternion_Dot(q1, q2) < 0) q2 = p2.q2;
    quaternion q = Slerp(q1, q2, t);
    return QuaternionToRotation(q);
}

struct euler_angles {
    float Pitch, Yaw, Roll;
};

euler_angles EulerAnglesFromRotation(mat44 Rotation)
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
    euler_angles Result {Pitch, Yaw, Roll};;
    return Result;
}

#endif // ALGEBRA_H
