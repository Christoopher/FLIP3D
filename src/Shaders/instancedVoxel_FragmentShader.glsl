#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in float isFLuid0;

void main(void)
{
	fragColor.rgb = vec3(isFLuid0);
	fragColor.a = isFLuid0;
}