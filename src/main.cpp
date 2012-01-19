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

void initVoxels(float * voxelPositions, float * voxelFlags, int k)
{
	float * posItr = voxelPositions;
	float * flagItr = voxelFlags;
	for (int x = 0; x < k; ++x)
	{
		for (int y = 0; y < k; ++y)
		{
			for (int z = 0; z < k; ++z)
			{
				*posItr++ = x;
				*posItr++ = y;
				*posItr++ = z;

				*flagItr++ = float(rand()) / RAND_MAX < 0.2 ? 1.0f : 0.0f;
			}
		}
	}
}

int main(void)
{
	float * verticies = new float[3*N];
	float * velocities = new float[3*N];

	float * voxelPositions = new float[3*10*10*10];
	float * voxelFlags = new float[10*10*10]; 
	initVoxels(voxelPositions,voxelFlags,10);

	for(int i =0; i<3*N; i++)
	{
		verticies[i] = 2*2.0*((float(rand()) / RAND_MAX) - 0.5);
		velocities[i] = float(rand()) / RAND_MAX;
	}

	OpenGl_initViewer(600, 600);
	//OpenGl_initParticles(verticies, velocities, sizeof(float)*N*3, N);
	OpenGl_initWireframeCube(voxelPositions,voxelFlags,10);
	
	while(running)
	{

		OpenGl_drawAndUpdate(running);

		//moveParticles(verticies, 0.1);
		//OpenGl_updateParticleLocation(verticies, sizeof(float)*N*3);

	}

	delete [] verticies;
	delete [] velocities;
	delete [] voxelFlags;
	delete [] voxelPositions;

	return 0;
}