
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
const int box = 10;
const int Nparticles = 70000*perCell;
const int dimx = 64, dimy = 64, dimz = 32;

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

/*void TestSSE() 
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
}*/

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

int main(void)
{
// 	int Nvoxels = dimx*dimy*dimz;
// 	Array3f voxelFlags(dimx,dimy,dimz);
// 	float * voxelPositions  = new float[3*Nvoxels*Nvoxels*Nvoxels];
// 	initVoxels(voxelPositions,voxelFlags,Nvoxels);

	Fluid_Solver fluid_solver(dimx,dimy,dimz,1.0f/32.0f,1.0f/30.0f,9.82f,1000,Nparticles);
	fluid_solver.init_box();

	OpenGl_initViewer(600, 600,fluid_solver.grid);
	OpenGl_initParticles(fluid_solver.particles.pos, fluid_solver.particles.vel, sizeof(vec3f)*fluid_solver.particles.currnp, fluid_solver.particles.currnp);
	//OpenGl_initWireframeCube(voxelPositions,voxelFlags.data,Nvoxels);
	
	while(running) {
		if(reset)
			fluid_solver.reset();
		reset = false;
		
		OpenGl_drawAndUpdate(running);

		if(step || play)
		{
			fluid_solver.step_frame();
		//	write_paricle_pos_binary(fluid_solver.particles);
		}

		OpenGl_updateParticleLocation(fluid_solver.particles.pos, sizeof(vec3f)*fluid_solver.particles.currnp);
		OpenGl_updateParticleVelocity(fluid_solver.particles.vel,sizeof(vec3f)*fluid_solver.particles.currnp);
		
// 		if(showgrid)
// 		{
// 			fluid_solver.grid.classify_voxel();
// 			update_voxel_flags(fluid_solver.grid,voxelFlags);
// 			OpenGl_updateVoxels(voxelPositions, voxelFlags.data, Nvoxels);
// 		}

	}
	
	return 0;
}



