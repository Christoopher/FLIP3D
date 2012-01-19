#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in float isFLuid0;

void main(void)
{
	fragColor.rgb = mix(vec3(0.4f,0.4f,1.0),vec3(1), isFLuid0);
	fragColor.a = 1.0;
}