#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in vec3 color;

void main(void)
{
	//Do velocity based shading here
	//Code
	vec3 poo = color;
	fragColor = vec4(poo,1.0);
}