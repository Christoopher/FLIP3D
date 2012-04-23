
uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
 
void main(void)
{	
	vec3 pos = gl_Vertex.xyz;
	//HÅRDKODAT 32
	pos.z = 32.0-pos.z;
	gl_Position = projectionMatrix*modelViewMatrix*vec4(pos,1.0);
}
