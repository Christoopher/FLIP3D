#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in vec3 color;

void main(void)
{
	//Do velocity based shading here
	//Code
	vec3 poo = mix(vec3(0.25f,0.25f,1.0f),vec3(1.0f),smoothstep(0.0,5.0,length(color)));
	fragColor = vec4(poo,1.0);
}