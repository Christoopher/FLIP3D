#version 400

layout(location = 0, index = 0) out vec4 fragColor;

in float isFLuid0;

void main(void)
{
	if(isFLuid0 < 0.5f) // AIR
	{
		fragColor.rgb = vec3(1.0f);
		fragColor.a = 0.0f;
	}
	else if(isFLuid0 < 1.5f && isFLuid0 > 0.5f) // FLUID == 1.0
	{
		fragColor.rgb = vec3(0.25f,0.25f,1.0f);
		fragColor.a = 1.0f;
	}
	else //SOLID
	{
		fragColor.a = 1.0f;
		fragColor.rgb = vec3(0.5f,0.25f,0.25f);
	}

}