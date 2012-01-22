
#include <iostream>
#include "OpenGLViewer.h"
#include "Vector3.h"
#include "Particles.h"
#include <stdio.h>

#include <xmmintrin.h>

#include <time.h>
#include "hr_time.h"
#include <limits>

bool running = true;
const int Nparticles = 1000000;
const int Nvoxels = 10;

void initVoxels(float * voxelPositions, float * voxelFlags, int k)
{
  float * posItr = voxelPositions;
  float * flagItr = voxelFlags;

  for (int x = 0; x < k; ++x) {
	for (int y = 0; y < k; ++y) {
	  for (int z = 0; z < k; ++z) {
	*posItr++ = x;
	*posItr++ = y;
	*posItr++ = z;	  
	*flagItr++ = float(rand()) / RAND_MAX < 0.2 ? 1.0f : 0.0f;
	  }
	}
  }

}

void TestSSE() 
{
	vec3f v1(0.0f);
	vec3f v2(1.0f);

	CStopWatch time;
	INT64 iterations = 0, additions = 100000000;
	time.startTimer();

	while(iterations < additions)
	{
		v1+=v2;
		iterations++;
	}
	time.stopTimer();
	std::cout << time.getElapsedTime() << std::endl;
	std::cout << v1 << std::endl;
	std::cout << v2 << std::endl;
	std::cin.get();
}

int main(void)
{
   float * verticies = new float[3*Nparticles];
   float * velocities = new float[3*Nparticles];

  Particles particles(Nparticles);
  
  float * voxelPositions = new float[3*Nvoxels*Nvoxels*Nvoxels];
  float * voxelFlags = new float[Nvoxels*Nvoxels*Nvoxels]; 
  initVoxels(voxelPositions,voxelFlags,Nvoxels);


  for(int i =0; i<Nparticles; i++){
	  vec3f newpos(5*2.0*((float(rand()) / RAND_MAX) - 0.5),5*2.0*((float(rand()) / RAND_MAX) - 0.5),5*2.0*((float(rand()) / RAND_MAX) - 0.5));
	  vec3f newvel(2.0*((float(rand()) / RAND_MAX) - 0.5), 2.0*((float(rand()) / RAND_MAX) - 0.5), 2.0*((float(rand()) / RAND_MAX) - 0.5));
	  
	  add_particle(particles, newpos, newvel);
  }

  OpenGl_initViewer(600, 600);
  get_position_larray(particles, verticies);
  get_velocity_larray(particles, velocities);
  OpenGl_initParticles(verticies, velocities, sizeof(float)*Nparticles*3, Nparticles);
  OpenGl_initWireframeCube(voxelPositions,voxelFlags,Nvoxels);
  
  int iterations = 0;
  while(running) {
	OpenGl_drawAndUpdate(running);
	  
	if(iterations > 200) {
	  initVoxels(voxelPositions,voxelFlags,Nvoxels);
	  OpenGl_updateVoxels(voxelPositions, voxelFlags, Nvoxels);
	  iterations = 0;
	}

	advect_particles(particles,0.1);
	get_position_larray(particles, verticies);
	OpenGl_updateParticleLocation(verticies, sizeof(float)*Nparticles*3);
	++iterations;
  }
  
  delete [] verticies;
  delete [] velocities;
  delete [] voxelFlags;
  delete [] voxelPositions;
  
  return 0;
}
