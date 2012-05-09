#version 120

varying vec3 v_N;
varying vec3 vLightDir;

void main( void )
{
	vec3 L = vLightDir;
	vec3 N = normalize(v_N);
	vec3 R = normalize(2.0*N*dot(N,L)-L);
	vec3 color = vec3(0.2,0.2,0.5);
	float amb = 0.5;
	vec3 ambient = color * amb;
	vec3 diffuse = color * (1.0-amb) * max(dot(L,N), 0.0);
	vec3 specular = vec3(1.0) * pow(max(dot(R,N), 0.0), 40.0);
	vec4 finalColor = vec4(ambient + diffuse + specular, 1.0);

	gl_FragColor = finalColor;

}

