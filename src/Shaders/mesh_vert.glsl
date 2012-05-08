#version 120


uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform mat3 normalMatrix;

varying vec3 v_N;
varying vec3 vLightDir;
 
void main(void)
{	
	v_N = normalize(normalMatrix*gl_Normal);	
	// Get vertex position in eye coordinates
	vec4 vPosition4 = modelViewMatrix * gl_Vertex;
	vec3 vPosition3 = vPosition4.xyz / vPosition4.w;
	vec3 lightpos = gl_LightSource[0].position.xyz;
	vLightDir = normalize(lightpos - vPosition3);

	vec3 pos = gl_Vertex.xyz;
	//pos.z = 16.0 - pos.z;
	gl_Position = projectionMatrix*modelViewMatrix*vec4(pos,1.0);
}
