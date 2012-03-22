#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in vec3 color;
uniform float edge;

void main(void)
{
	//Do velocity based shading here
	fragColor = mix(vec4(0.25f,0.25f,1.0f,1.0f),vec4(1.0f),smoothstep(0.0f,edge,length(color)));
}