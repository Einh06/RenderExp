#version 330 core 

layout (location = 0) in vec3 position; 
layout (location = 1) in vec3 normal;

uniform mat4 Model;
uniform mat4 View;

uniform mat4 MVP; 

out vec3 Normal_cameraspace;
out vec3 LightDirection_cameraspace;

void main() { 

    gl_Position = MVP * vec4(position, 1.0); 

    Normal_cameraspace = (View * Model * vec4(normal, 1.0)).xyz;
    LightDirection_cameraspace = (vec4(0, 0, 1, 1)).xyz;
}
