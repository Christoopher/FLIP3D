#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in vec3 velocity;

void main(void)
{
	//Do velocity based shading here
	//Code

	fragColor = vec4(vec3(smoothstep(0.0f,1.0f,sqrt(length(velocity))),1.0f);
}
