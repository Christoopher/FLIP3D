#version 400

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform int dimz;
uniform float h;

in vec3 velocity;
in vec3 vertex;

out vec3 color;
 
void main(void)
{	
	color = velocity;
	vec3 pos = vertex;
	pos.z = h*dimz - pos.z;
	pos /= h;
	//pos *= h/h;
	gl_Position = projectionMatrix*modelViewMatrix*vec4(pos,1.0);
}