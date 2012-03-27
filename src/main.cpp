
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

#include "ObjLoader.h"


const int perCell = 8;
const int box = 10;
const int Nparticles = 30000*perCell;
const int dimx = 64, dimy = 32, dimz = 32;

void initVoxels(float * voxelPositions, int dx, int dy, int dz)
{
	//OBS Ordningen på loopen är viktig här
	float * posItr = voxelPositions;
	for (int z = 0; z < dz; ++z) {
		for (int y = 0; y < dy; ++y) {
			for (int x = 0; x < dx; ++x) {
				*posItr++ = (float)x;
				*posItr++ = (float)y;
				*posItr++ = dz - (float)z;	  
			}
		}
	}
}

void update_voxel_flags(Grid & grid, Array3f & flags)
{

	for(int k = 0; k < grid.Nz; ++k)
		for(int j = 0; j < grid.Ny; ++j)
			for(int i = 0; i < grid.Nx; ++i)
			{
				flags(i,j,k) = (float)grid.marker(i,j,k);
			}
}

CStopWatch stopwatch;

int numframes;
double avgtime = 0;
int main(void)
{
	Fluid_Solver fluid_solver(dimx,dimy,dimz,1.0f,1.0f/30.0f,9.82f,1.0f,Nparticles);
	fluid_solver.init_box();

	OpenGl_initViewer(600, 600,fluid_solver.grid);
	OpenGl_initParticles(fluid_solver.particles.pos, fluid_solver.particles.vel, sizeof(vec3f)*fluid_solver.particles.currnp, fluid_solver.particles.currnp);
	
	
	int Nvoxels = dimx*dimy*dimz;
	Array3f voxelFlags(dimx,dimy,dimz);	
	float * voxelPositions  = new float[3*dimx*dimy*dimz];
	
	initVoxels(voxelPositions,dimx,dimy,dimz);

	OpenGl_initWireframeCube(voxelPositions,voxelFlags.data,Nvoxels);
	update_voxel_flags(fluid_solver.grid,voxelFlags);
	
	while(running) {
		if(reset)
			fluid_solver.reset();
		reset = false;

		if(showgrid)
		{
			testMesh.mesh_to_grid(fluid_solver.grid);
			update_voxel_flags(fluid_solver.grid,voxelFlags);
			OpenGl_updateVoxels(voxelPositions, voxelFlags, Nvoxels);
		}
		
		OpenGl_drawAndUpdate(running);

		if(step || play)
		{
//			stopwatch.startTimer();
			fluid_solver.step_frame();
// 			stopwatch.stopTimer();
// 			std::cout << std::scientific;
// 			std::cout << stopwatch.getElapsedTime() << "\n";
// 			if(numframes >= 6)
// 				avgtime += stopwatch.getElapsedTime();
// 			numframes++;
			
		//	write_paricle_pos_binary(fluid_solver.particles);
		}

		

		OpenGl_updateParticleLocation(fluid_solver.particles.pos, sizeof(vec3f)*fluid_solver.particles.currnp);
		OpenGl_updateParticleVelocity(fluid_solver.particles.vel,sizeof(vec3f)*fluid_solver.particles.currnp);
		
		

	}
	

	TerminateViewer();
// 	std::cout << std::scientific;
// 	std::cout << "Avg. time" << avgtime/(numframes-6) << "\n";
// 	std::cout << "Press any key to quit...\n";
// 	std::cin.get();

	

	

	return 0;
}



