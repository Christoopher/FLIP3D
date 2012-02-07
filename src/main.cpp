
#include <limits>
#include <stdio.h>
#include <iostream>

#include "OpenGLViewer.h"
#include "Vector3.h"
#include "Particles.h"
#include "Grid.h"
#include "Array3D.h"
#include <time.h>
#include "hr_time.h"
#include "Sparse_Matrix.h"
#include "Fluid_Solver.h"
#include "Unconditioned_CG_Solver.h"


bool running = true;
const int perCell = 8;
const int box = 20;
const int Nparticles = box*box*box*perCell;
const int Nvoxels = 32;

void initVoxels(float * voxelPositions, Array3f & voxelFlags, int k)
{
	float * posItr = voxelPositions;

	for (int x = 0; x < k; ++x) {
		for (int y = 0; y < k; ++y) {
			for (int z = 0; z < k; ++z) {
				*posItr++ = (float)x;
				*posItr++ = (float)y;
				*posItr++ = (float)z;	  
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

		flags.data[i] = grid.marker.data[i];

		//if(grid.marker.data[i] == FLUIDCELL)
		//	flags.data[i] = FLUIDCELL;
		//else
		//	flags.data[i] = 0;
	}
}

/*void TestSparse() 
{
	int dim = 4;
	Sparse_Matrix mtx(dim,dim,dim);
	VectorN vec(dim*dim*dim);
	VectorN res(dim*dim*dim);

	size_t size = sizeof(mtx);

	for(int k=1; k<dim-1; ++k)
	{
		for(int j = 1; j < dim-1; ++j)
		{
			for (int i = 1; i < dim-1; ++i)
			{
				mtx(i,j,k,0) += 1;
				mtx(i,j,k,1) += 1;
				mtx(i,j,k,2) += 1;
				mtx(i,j,k,3) += 1;
			}
		}
	}

	int row = 0;
	for(int k=0; k<dim; ++k)
	{
		for(int j = 0; j < dim; ++j)
		{
			for (int i = 0; i < dim; ++i)
			{
				std::cout << row << ": " << mtx(i,j,k,0) << ", " << mtx(i,j,k,1) << ", " << mtx(i,j,k,2) << ", " << mtx(i,j,k,2)  << ", " << mtx(i,j,k,3) << "\n"; 
				row++;
			}
		}
	}


	for(int i = 0; i<dim*dim*dim; ++i)
	{
		vec.data[i] = 1;
	}

	mtx_mult_vectorN(mtx,vec,res,marker);

	row = 0;
	for(int i = 0; i<dim*dim*dim; ++i)
	{
		std::cout << row << ": " << res.data[i] << "\n";
		row++;
	}
}*/





int main(void)
{

	float * verticies = new float[3*Nparticles];
	float * velocities = new float[3*Nparticles];

	//Array3f voxelFlags(Nvoxels,Nvoxels,Nvoxels);
	//float * voxelPositions  = new float[3*Nvoxels*Nvoxels*Nvoxels];
	//initVoxels(voxelPositions,voxelFlags,Nvoxels);

	Fluid_Solver fluid_solver(Nvoxels,Nvoxels,Nvoxels,1.0/100,9.82f,1000,Nparticles);
	

	float h = fluid_solver.grid.h;
// 	add_particle(fluid_solver.particles,vec3f(0.5f*h,h*9.5f,0.5f*h), vec3f(0.0f));	
// 	add_particle(fluid_solver.particles,vec3f(0.5f*h,h*9.5f,h*9.5f), vec3f(0.0f));
// 	add_particle(fluid_solver.particles,vec3f(h*9.5f,h*9.5f,0.5f*h), vec3f(0.0f));
// 	add_particle(fluid_solver.particles,vec3f(h*9.5f,h*9.5f,h*9.5f), vec3f(0.0f));

	fluid_solver.init_box();

	OpenGl_initViewer(600, 600);
	get_position_larray(fluid_solver.particles, verticies);
	get_velocity_larray(fluid_solver.particles, velocities);
	OpenGl_initParticles(verticies, velocities, sizeof(float)*3*Nparticles, Nparticles);
	//OpenGl_initWireframeCube(voxelPositions,voxelFlags.data,Nvoxels);


	while(running) {
		if(reset)
			fluid_solver.reset();
		reset = false;
		
		OpenGl_drawAndUpdate(running);

		if(step || play)
			fluid_solver.step_frame();

		get_position_larray(fluid_solver.particles, verticies);
		OpenGl_updateParticleLocation(verticies, sizeof(float)*3*Nparticles);
		get_velocity_larray(fluid_solver.particles,velocities);
		OpenGl_updateParticleVelocity(velocities,sizeof(float)*3*Nparticles);
		
// 		if(showgrid)
// 		{
// 			fluid_solver.grid.classify_voxel();
// 			update_voxel_flags(fluid_solver.grid,voxelFlags);
// 			OpenGl_updateVoxels(voxelPositions, voxelFlags.data, Nvoxels);
// 		}

	}

	delete [] verticies;
	delete [] velocities;
	//delete [] voxelPositions;

	return 0;
}



