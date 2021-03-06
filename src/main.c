#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <unistd.h>    // getcwd
#include <sys/param.h> // MAXPATHLEN
#include <stdio.h>
#include <stdlib.h>

#include "Common/common.h" // some things I like to use
#include "Common/math.h"   // some of the common Math stuff I need
#include "Common/algebra.h" // impo

void MemoryCopy(void *Dest, void *Source, size_t Count)
{
    char* d = Dest;
    char* s = Source;
    for (size_t i = 0; i < Count; ++i) {
        d[i] = s[i]; 
    }
}

typedef struct triangle {
    Vec3f v1;
    Vec3f v2;
    Vec3f v3;
} triangle;

typedef struct triangle_edges_index {
    u32 v1;
    u32 v2;
    u32 v3;
} triangle_edges_index;

typedef struct camera {
    Vec3f Position;
    Quaternion Orientation;
    Vec3f Forward;
    Vec3f Up;
} camera;

// This is a simple mesh generated by a model.
typedef struct mesh3d {

    u32 VerticesBufferSize;
    u32 VerticesCount;
    Vec3f* Vertices;

    u32 IndexesBufferSize;
    u32 IndexesCount;
    u32* Indexes;
} mesh3d;

typedef struct Ray {
    Vec3f Origin;
    Vec3f Direction;
} Ray;

Ray Ray_Make(Vec3f Origin, Vec3f Direction)
{
    return (Ray){.Origin = Origin, .Direction = Direction};
}

float IntersectRayPlane(Ray r, Vec3f v1, Vec3f v2)
{
    //TODO(@Florian): Implement this!
}

float IntersectRaySphere(Ray r, Pos3f SphereCenter, float Radius)
{
    //TODO(@Florian): Implement this!
}

bool IntersectRayTriangle(Ray r, Pos3f p1, Pos3f p2, Pos3f p3)
{
    //TODO(@Florian): Implement this!
}


/**
 * This will make a copy the vertices and triangle buffer.
 */
mesh3d mesh3d_Make(u32 VerticesCount, Vec3f* Vertices, u32 IndexesCount, u32* Indexes)
{
    mesh3d Result;
    Result.VerticesCount = VerticesCount;
    Result.VerticesBufferSize = sizeof(Vec3f) * VerticesCount;
    Result.Vertices = (Vec3f*)malloc(Result.VerticesBufferSize);
    MemoryCopy(Result.Vertices, Vertices, Result.VerticesBufferSize);

    Result.IndexesCount = IndexesCount;
    Result.IndexesBufferSize = sizeof(u32) * IndexesCount;
    Result.Indexes = (u32*)malloc(Result.IndexesBufferSize);
    MemoryCopy(Result.Indexes, Indexes, Result.IndexesBufferSize);

    return Result;
}

u32 StrLength(char *c) {

    u32 Result = 0;
    while (c[Result++]) {}
    return Result - 1;
}

typedef struct string {
    char *Buffer;
    u32   StringSize;
    u32   BufferSize;
} string;

void string_Free(string s) {
    s.BufferSize = 0; 
    s.StringSize = 0; 
    free(s.Buffer);
}

// Assuming 0-termination
string string_MakeFrom(char *Base) {
    string Result;
    Result.StringSize = StrLength(Base);
    Result.BufferSize = Result.StringSize + 1;
    Result.Buffer = (char*)malloc(Result.BufferSize);
    
    for (u32 i = 0; i < Result.StringSize; ++i) {
        Result.Buffer[i] = Base[i];
    }
    Result.Buffer[Result.StringSize] = 0;
    return Result;
}

string string_Append(string A, string B)
{
    string Result;

    Result.StringSize = A.StringSize + B.StringSize;
    Result.BufferSize = Result.StringSize + 1;
    Result.Buffer = (char*)malloc(Result.BufferSize);

    u32 WriteIndex = 0;
    for (u32 i = 0; i < A.StringSize; ++i) {
       Result.Buffer[WriteIndex++] = A.Buffer[i]; 
    }
    for (u32 i = 0; i < B.StringSize; ++i) {
       Result.Buffer[WriteIndex++] = B.Buffer[i]; 
    }
    Result.Buffer[WriteIndex] = 0;

    return Result; 
}

typedef struct file_content {
    size_t BufferSize;
    char* Buffer;
} file_content;

file_content LoadEntireFile(string Location, bool BinaryMode)
{
    file_content Result;

    FILE* File = fopen(Location.Buffer, BinaryMode ? "rb" : "r");
    Assert(File != NULL);

    int Length;
    
    fseek(File, 0, SEEK_END);
    Length = ftell(File); 
    fseek(File, 0, SEEK_SET);

    char* Filebuffer = (char*)malloc(Length);

    int ReadCount = fread(Filebuffer, 1, Length, File);
    Assert(ReadCount == Length);

    Result.BufferSize = Length;
    Result.Buffer = Filebuffer;
    return Result;
}

void FreeFileContent(file_content FileContent)
{
    free(FileContent.Buffer);
}
    
typedef struct ifs_data {

    u32 VerticesCount;
    Vec3f *Vertices;

    u32 FacesCount;
    triangle_edges_index *Faces;
} ifs_data;

ifs_data ParseIFSData(char* Filebuffer) {

#define GET_NEXT_SIZE (*(u32*)(Filebuffer)); Filebuffer += 4

    ifs_data Result;
    //fileheader processing      
    u32 IFSStringSize = GET_NEXT_SIZE;
    Filebuffer += IFSStringSize + sizeof(float) /*skip model version*/;

    u32 ModelNameSize = GET_NEXT_SIZE;
    Filebuffer += ModelNameSize;

    //vertexheader processing
    u32 VERTICESStrSize = GET_NEXT_SIZE; 
    Filebuffer += VERTICESStrSize;

    u32 VerticesCount = GET_NEXT_SIZE;

    //vertex processing
    size_t VerticesBufferSize = (sizeof(Vec3f) * VerticesCount); 
    Vec3f* Vertices = (Vec3f*)malloc(VerticesBufferSize);

    for (u32 i = 0; i < VerticesCount; ++i) {
        Vec3f* VerticesValue = (Vec3f*)(Filebuffer);
        Vertices[i] = *VerticesValue;
        Filebuffer += sizeof(Vec3f);
    }
    //triangle header processing
    u32 TriangleStrSize = GET_NEXT_SIZE; 
    Filebuffer += TriangleStrSize;

    u32 FacesCount = GET_NEXT_SIZE;

    // triangle edges processing 
    size_t FacesBufferSize = (sizeof(triangle_edges_index) * FacesCount);
    triangle_edges_index* Faces = (triangle_edges_index*)malloc( FacesBufferSize);

    for (u32 i = 0; i < FacesCount; ++i) {

        triangle_edges_index *IndexValue = (triangle_edges_index*)(Filebuffer);
        Faces[i] = *IndexValue;

        Filebuffer += sizeof(triangle_edges_index);
    }

    Result.VerticesCount = VerticesCount;
    Result.Vertices = Vertices;

    Result.FacesCount = FacesCount;
    Result.Faces = Faces;

    return Result;
#undef GET_NEXT_SIZE
}


struct {
    mesh3d ModelMesh;
} Global;


file_content LoadShaderSource(string Location)
{
    return LoadEntireFile(Location, false);
}

/********************
 *  MAIN SHIT
 ********************/
void KeyCallback(GLFWwindow* Window, int Key, int Scancode __unused, int Action, int Mode __unused)
{
    if(Key == GLFW_KEY_ESCAPE && Action == GLFW_PRESS) {
        glfwSetWindowShouldClose(Window, GL_TRUE);
    }
}


void Initialisation(string ModelFilepath)
{
    //Get file buffer here and send it to the parser
    file_content ifsContent = LoadEntireFile(ModelFilepath, true);
    
    ifs_data IFSData = ParseIFSData(ifsContent.Buffer);

    mesh3d NewMesh = mesh3d_Make(IFSData.VerticesCount, IFSData.Vertices, IFSData.FacesCount * 3, (u32*)IFSData.Faces);
    Global.ModelMesh = NewMesh;

    FreeFileContent(ifsContent);
}

Mat44 ViewMatrixFromCamera(camera Camera) 
{
    Vec3f Forward = Vec3f_NormalizeFrom(Camera.Orientation.v);
    return Mat44_LookAt(Camera.Position, Vec3f_Add(Camera.Position, Forward), Vec3f_Y);
}

int main(int Argc __unused, char** Argv __unused)
{
    char CurrentWorkingDirectoryPath[MAXPATHLEN];
    getcwd(CurrentWorkingDirectoryPath, MAXPATHLEN);

    string Cwd = string_MakeFrom(CurrentWorkingDirectoryPath);
    string TeapotPath = string_MakeFrom("/data/teapot.ifs");
    string AbsoluteTeapotPath = string_Append(Cwd, TeapotPath);
     
    Initialisation(AbsoluteTeapotPath);

    camera Camera = (camera) {.Position = Vec3f_Make(0, 0, -1), .Orientation = Quaternion_MakeReal(0, 0, 0, 1), .Forward = Vec3f_Make(0, 0, 1), .Up = Vec3f_Y};

    string_Free(TeapotPath);
    string_Free(AbsoluteTeapotPath);

    // I want to do some flat shading. My current understanding of OpenGL is that there are
    // no ways to get a normal per triangle, all normals are stored per vertices. This 
    // prevents me of using the vertices/indexes (called elements) based drawing.
    // The approch is to then create a expended vertex buffer and use the DrawArrays method, 
    // and calculate the assign the triangle normal as a vertex normal. Since the buffer 
    // contains duplicated vertices, a vertex might have multiple normal assigned

    // Generate the expended vertices buffer and its associated normal buffer
    u32 VerticesCount = Global.ModelMesh.IndexesCount;
    Vec3f* Vertices = (Vec3f*)malloc(VerticesCount * sizeof(Vec3f));

    u32 NormalsCount = VerticesCount;
    Vec3f* Normals = (Vec3f*)malloc(NormalsCount * sizeof(Vec3f));

    for (u32 i = 0; i < Global.ModelMesh.IndexesCount; i+=3) {
        
        u32 i1 = Global.ModelMesh.Indexes[i + 0];
        u32 i2 = Global.ModelMesh.Indexes[i + 1];
        u32 i3 = Global.ModelMesh.Indexes[i + 2];

        Vec3f v1 = Global.ModelMesh.Vertices[i1]; 
        Vec3f v2 = Global.ModelMesh.Vertices[i2]; 
        Vec3f v3 = Global.ModelMesh.Vertices[i3]; 

        Vec3f e1 = Vec3f_Substract(v2, v1);
        Vec3f e2 = Vec3f_Substract(v3, v2);

        Vertices[i + 0] = v1;
        Vertices[i + 1] = v2;
        Vertices[i + 2] = v3;

        Vec3f Normal = Vec3f_NormalizeFrom(Vec3f_Cross(e2, e1));
        Normals[i + 0] = Normal;
        Normals[i + 1] = Normal;
        Normals[i + 2] = Normal;
    }
    //End of manually computing vertices and their normals


    // Init drawing stuff
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); //macOS
    
    GLFWwindow* Window = glfwCreateWindow(800, 600, "Experiment", NULL, NULL);    
    if(Window == NULL) {
        printf("Cannot open window"); 
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(Window);

    glfwSetKeyCallback(Window, KeyCallback);
    glfwSetInputMode(Window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        
        printf("Cannot init glew"); 
        return -1; 
    }

    int Width, Height;
    glfwGetFramebufferSize(Window, &Width, &Height);
    glViewport(0, 0, Width, Height);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    // end of drawing initlisation


    // Creating shader program
    string VertexShaderSourcePath = string_MakeFrom("/data/shaders/default.vert");
    string FragmentShaderSourcePath = string_MakeFrom("/data/shaders/default.frag");

    string VSAbsolutePath = string_Append(Cwd, VertexShaderSourcePath);
    string FSAbsolutePath = string_Append(Cwd, FragmentShaderSourcePath);

    file_content  VertexShaderSource = LoadShaderSource(VSAbsolutePath);
    file_content  FragmentShaderSource = LoadShaderSource(FSAbsolutePath);

    string_Free(VertexShaderSourcePath);
    string_Free(FragmentShaderSourcePath);
    string_Free(VSAbsolutePath);
    string_Free(FSAbsolutePath);


    GLuint VertexShader, FragmentShader;

    VertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(VertexShader, 1, &VertexShaderSource.Buffer, NULL);
    glCompileShader(VertexShader);

    GLint Success;
    GLchar Info[512];

    glGetShaderiv(VertexShader, GL_COMPILE_STATUS, &Success);
    if (!Success) {
        glGetShaderInfoLog(VertexShader, 512, NULL, Info);
        printf("Could not compile vertex shader: %s",Info);
        exit(1);
    }

    FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(FragmentShader, 1, &FragmentShaderSource.Buffer, NULL);
    glCompileShader(FragmentShader);

    glGetShaderiv(FragmentShader, GL_COMPILE_STATUS, &Success);
    if (!Success) {
        glGetShaderInfoLog(FragmentShader, 512, NULL, Info);
        printf("Could not compile fragment shader: %s",Info);
        exit(1);
    }

    GLuint ShaderProgram;

    ShaderProgram = glCreateProgram();
    glAttachShader(ShaderProgram, VertexShader);
    glAttachShader(ShaderProgram, FragmentShader);
    glLinkProgram(ShaderProgram);

    glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &Success);
    if(!Success) {
        glGetProgramInfoLog(ShaderProgram, 512, NULL, Info);
        printf("Could not link program: %s", Info);
        exit(1);
    }

    glUseProgram(ShaderProgram);
    glDeleteShader(VertexShader);
    glDeleteShader(FragmentShader);
    FreeFileContent(VertexShaderSource);
    FreeFileContent(FragmentShaderSource);
    //Done creating shader
    
    //Generating buffers for vertices
    GLuint VAO, VBO, NBO;
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &NBO);

    glGenVertexArrays(1, &VAO);

    glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, VerticesCount * sizeof(Vec3f), (void*)Vertices, GL_STATIC_DRAW); // Computer not happy when you put this in the loop
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

        glBindBuffer(GL_ARRAY_BUFFER, NBO);
        glBufferData(GL_ARRAY_BUFFER, NormalsCount * sizeof(Vec3f), (void*)Normals, GL_STATIC_DRAW);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
        
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    //Done generating vertices

    // Some stuff that doesn't need to be recomputed every frame
    // Not that it Matters anyway, lets pretend we are efficient here
    Mat33 StartRotationInterpolation = Mat33_MakeRotation(Vec3f_Y, M_PI_4);
    Mat33 EndRotationInterpolation   = Mat33_MakeRotation(Vec3f_Y, M_PI_2 + M_PI_4);
    Mat44 TranslationMatrix = Mat44_MakeTranslate(0, 0, 3);

    Mat44 Projection = Mat44_Perspective(35, 0.1, 100, Width, Height);

    double LastTimeStamp = glfwGetTime();
    double TimeForFrame = 0.01666;


    double CursorX, CursorY, PrevCursorX, PrevCursorY;
    glfwGetCursorPos(Window, &CursorX, &CursorY); 
    PrevCursorX = CursorX;
    PrevCursorY = CursorY;

    while(!glfwWindowShouldClose(Window)) {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
        glfwPollEvents();

        double CurrentTime = glfwGetTime();
        double DeltaTime = CurrentTime - LastTimeStamp;


        glfwGetCursorPos(Window, &CursorX, &CursorY);
        float DiffX = CursorX - PrevCursorX;
        float DiffY = CursorY - PrevCursorY;
        if ((Absf(DiffX) > F_EPS) ||
            (Absf(DiffY) > F_EPS) )
        {
            const float Sensitivity = 0.05;

            Mat33 MatRotationY = Mat33_MakeRotationY(-DiffX * Sensitivity * DeltaTime);
            Mat33 MatRotationX = Mat33_MakeRotationX(-DiffY * Sensitivity * DeltaTime);
            Mat33 MatRotation = Mat33_MultMat(MatRotationY, MatRotationX);

            quaternion_pair Pair = RotationToQuaternion(MatRotation);

            Quaternion Rotation = Pair.q1;
            Camera.Orientation = Quaternion_NormalizeFrom(Quaternion_MultQuaternion(Rotation, Camera.Orientation));
        }  
        PrevCursorX = CursorX;
        PrevCursorY = CursorY;

        Vec3f MoveDirection = Vec3f_Make(0, 0, 0);

        if (glfwGetKey(Window, GLFW_KEY_W) == GLFW_PRESS) {
            MoveDirection.z += 1;
        }
        if (glfwGetKey(Window, GLFW_KEY_S) == GLFW_PRESS) {
            MoveDirection.z -= 1;
        }
        if (glfwGetKey(Window, GLFW_KEY_A) == GLFW_PRESS) {
            MoveDirection.x -=1;
        }
        if (glfwGetKey(Window, GLFW_KEY_D) == GLFW_PRESS) {
            MoveDirection.x +=1;
        }
        if (MoveDirection.x != 0 || MoveDirection.z != 0) {
            Vec3f Forward = Vec3f_NormalizeFrom(Quaternion_MultVec(Camera.Orientation, MoveDirection));
             
            const float Speed = 10.0;
            Camera.Position = Vec3f_Add(Camera.Position, Vec3f_Mult(Forward, Speed * DeltaTime));
        }

        float t = (Cosf(glfwGetTime()) + 1.0) / 2.0;

        Mat33 Rotation = RotationInterpolation(StartRotationInterpolation, EndRotationInterpolation, t);
        Mat44 RotationMatrix = Mat44_FromMat33(Rotation);

        Mat44 Model = Mat44_MultMat(TranslationMatrix, RotationMatrix);
        Mat44 View = ViewMatrixFromCamera(Camera);//*/Mat44_LookAt(Vec3f_Make(0, 0, -1), Vec3f_Make(0, 0, 2.5), Vec3f_Make(0, 1, 0));
        Mat44 MVP = Mat44_MultMat3(Projection, View, Model);

        //my Matrices are row-major, opengl wants column-major, therefore we need to transpose
        GLint ModelLocation = glGetUniformLocation(ShaderProgram, "Model");
        glUniformMatrix4fv(ModelLocation, 1, GL_TRUE, (float*)Model.m);
        
        GLint ViewLocation = glGetUniformLocation(ShaderProgram, "View");
        glUniformMatrix4fv(ViewLocation, 1, GL_TRUE, (float*)View.m); 

        GLint MVPLocation = glGetUniformLocation(ShaderProgram, "MVP");
        glUniformMatrix4fv(MVPLocation, 1, GL_TRUE, (float*)MVP.m);

        // Draw some shit
        glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, VerticesCount);
        glBindVertexArray(0);

        glfwSwapBuffers(Window);

        LastTimeStamp = CurrentTime;

        // Me showing off i can manually manage a 60fps so that i don't use 80% because i'm drawing like a mad man
        // Obviously, it barelly works, sometime it just hangs
        while ((glfwGetTime() - CurrentTime) < TimeForFrame) { usleep(200); }
        
    }

    glfwTerminate();
    return 0;
}
