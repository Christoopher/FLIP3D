#version 400

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;

in vec4 vertex;
inout vec3 velocity; //might not work with inout
 
void main(void)
{
	gl_Position = projectionMatrix*modelViewMatrix*vertex;
}