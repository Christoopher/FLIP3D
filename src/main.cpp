#include <iostream>
#include "OpenGLViewer.h"

bool running = true;
const int N = 100000;

void moveParticles(float * verticies, float dt)
{
	int stride = 3;
	
	float* pos = verticies;

	for (int i = 0; i<N; ++i)
	{

		*pos++ = *pos + dt* sin(2.0*((float(rand()) / RAND_MAX) - 0.5)); //x
		*pos++ = *pos + dt* sin(2.0*((float(rand()) / RAND_MAX) - 0.5)); //y
		*pos++ = *pos + dt* sin(2.0*((float(rand()) / RAND_MAX) - 0.5)); //z		
	}
}

int main(void)
{
	float * verticies = new float[3*N];
	float * velocities = new float[3*N];

	for(int i =0; i<3*N; i++)
	{
		verticies[i] = 2*2.0*((float(rand()) / RAND_MAX) - 0.5);
		velocities[i] = float(rand()) / RAND_MAX;
	}

	initViewer(600, 600);
	initParticles(verticies, velocities, sizeof(float)*N*3, N);

	
	while(running)
	{

		drawAndUpdate();

		moveParticles(verticies, 0.01);
		updateParticleLocation(verticies, sizeof(float)*N*3);

	}

	delete [] verticies;
	delete [] velocities;


	return 0;
}