#version 400

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;

in vec3 velocity;
in vec3 vertex;

out vec3 color;
 
void main(void)
{	
	color = velocity;
	gl_Position = projectionMatrix*modelViewMatrix*vec4(vertex,1.0);
}