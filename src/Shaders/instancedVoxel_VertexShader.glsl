#version 400

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;

in vec3 vertex;
in vec3 position;
in float isFluid;

out float isFLuid0;
 
void main(void)
{	
	isFLuid0 = isFluid;
	gl_Position = projectionMatrix*modelViewMatrix*vec4(vertex + position,1.0);
}