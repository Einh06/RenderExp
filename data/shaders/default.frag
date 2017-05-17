#version 330 core

out vec4 color;

in vec3 Normal_cameraspace;
in vec3 LightDirection_cameraspace;

void main() { 

    vec3 l = normalize(-LightDirection_cameraspace);
    vec3 n = normalize(Normal_cameraspace);

    float cosTheta = abs(dot(l, n));
    color = vec4(1.0, 1.0, 1.0, 1.0) * cosTheta;
}
