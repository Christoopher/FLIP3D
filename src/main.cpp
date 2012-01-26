
#include <iostream>
#include "OpenGLViewer.h"
#include "Vector3.h"
#include "SSEVector3.h"
#include "Particles.h"
#include <stdio.h>
#include "Grid.h"
#include "Array3D.h"

#include <xmmintrin.h>

#include <time.h>
#include "hr_time.h"
#include <limits>

bool running = true;
const int Nparticles = 4;
const int Nvoxels = 10;
const float h = 1.0f/Nvoxels;

void initVoxels(float * voxelPositions, Array3f & voxelFlags, int k)
{
	float * posItr = voxelPositions;

	for (int x = 0; x < k; ++x) {
		for (int y = 0; y < k; ++y) {
			for (int z = 0; z < k; ++z) {
				*posItr++ = x;
				*posItr++ = y;
				*posItr++ = z;	  
				voxelFlags(x,y,z) = 1.0;
			}
		}
	}

}

void TestSSE() 
{
	vec3f v1(1.0,2.0,1.0);
	SSE::vec3f v2(1.0,2.0,1.0);

	CStopWatch time;
	INT64 iterations = 0, additions = 10000000;
	float result;
	time.startTimer();
	while(iterations < additions)
	{
		result = mag2(v1);
		iterations++;
	}
	time.stopTimer();
	std::cout << time.getElapsedTime() << std::endl;
	std::cout << result  << std::endl;

	iterations = 0;
	time.startTimer();
	while(iterations < additions)
	{
		result = mag2(v2);
		iterations++;
	}
	time.stopTimer();
	std::cout << time.getElapsedTime() << std::endl;
	std::cout << result  << std::endl;
	std::cin.get();
}

void update_voxel_flags(Grid & grid, Array3f & flags)
{
	for(int i = 0; i < grid.marker.size; ++i)
	{
		flags.data[i] = grid.marker.data[i] == FLUIDCELL ? 1.0f : 0.0f;
	}
}

//----------------------------------------------------------------------------//
// Iterates fluid solver on timestep
//----------------------------------------------------------------------------//
void Solve(float timestep, Grid & grid, Particles & particles)
{
	int iterations = 0;
	std::cout << std::scientific;
	for(float elapsed = 0; elapsed < timestep;)
	{
		//Fetch max timestep
		float dt = grid.CFL();
		//std::cout << dt << "\n";
		if( dt > timestep - elapsed)
			dt = timestep - elapsed;

		elapsed += dt;

		transfer_to_grid(particles,grid);
		grid.save_velocities();
		grid.add_gravity(dt);
		grid.apply_boundary_conditions();
		grid.get_velocity_update();
		update_from_grid(particles,grid);
		advect_particles(particles,dt);
		iterations++;
	}

	//std::cout << "iterations: " << iterations << " timestep:" << timestep <<  "\n";
}





int main(void)
{

	float * verticies = new float[3*Nparticles];
	float * velocities = new float[3*Nparticles];


	Array3f voxelFlags(Nvoxels,Nvoxels,Nvoxels);
	float voxelPositions[3*Nvoxels*Nvoxels*Nvoxels];
	initVoxels(voxelPositions,voxelFlags,Nvoxels);


	Grid grid(Nvoxels,Nvoxels,Nvoxels,h,9.82f);
	Particles particles(Nparticles,grid);
	
		
// 	for(int i =0; i<Nparticles; i++){
// 		vec3f newpos(5*2.0*((float(rand()) / RAND_MAX) - 0.5),5*2.0*((float(rand()) / RAND_MAX) - 0.5),5*2.0*((float(rand()) / RAND_MAX) - 0.5));
// 		vec3f newvel(2.0*((float(rand()) / RAND_MAX) - 0.5), 2.0*((float(rand()) / RAND_MAX) - 0.5), 2.0*((float(rand()) / RAND_MAX) - 0.5));
// 
// 		add_particle(particles, newpos, newvel);
// 	}

	add_particle(particles,vec3f(0.5f*h,h*9.5f,0.5f*h), vec3f(0.0f));
	add_particle(particles,vec3f(0.5f*h,h*9.5f,h*9.5f), vec3f(0.0f));
	add_particle(particles,vec3f(h*9.5f,h*9.5f,0.5f*h), vec3f(0.0f));
	add_particle(particles,vec3f(h*9.5f,h*9.5f,h*9.5f), vec3f(0.0f));

	OpenGl_initViewer(600, 600);
	get_position_larray(particles, verticies);
	get_velocity_larray(particles, velocities);
	OpenGl_initParticles(verticies, velocities, sizeof(float)*3*Nparticles, Nparticles);
	OpenGl_initWireframeCube(voxelPositions,voxelFlags.data,Nvoxels);


	float dt = 10e-5;
	while(running) {
		OpenGl_drawAndUpdate(running);

		Solve(dt,grid, particles);

		get_position_larray(particles, verticies);
		OpenGl_updateParticleLocation(verticies, sizeof(float)*3*Nparticles);
		update_voxel_flags(grid,voxelFlags);
		OpenGl_updateVoxels(voxelPositions, voxelFlags.data, Nvoxels);

	}

	delete [] verticies;
	delete [] velocities;
	delete [] voxelPositions;

	return 0;
}



