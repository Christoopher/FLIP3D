
#include <limits>
#include <stdio.h>
#include <iostream>
#include <vector>

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

#include "MarchingCubes/MarchingCubes.h"



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



void initLevelset(int dim, float r)
{
	srand(int(time(NULL)));
	//levelset.resize(dim*dim*dim);
	levelset = new mp4Vector[dim*dim*dim];
	float center = dim/2.0f + 0.5f;
	for(int k = 0; k < dim; ++k)
		for(int j = 0; j < dim; ++j)
			for(int i = 0; i < dim; ++i)
			{
				float dx = (center - i);
				float dy = (center - j);
				float dz = (center - k);
				float dist = sqrtf( dx*dx + dy*dy + dz*dz )/dim;
				mp4Vector vert(10*float(i-center)/dim, 10*float(j-center)/dim,10*float(k-center)/dim,0);
				vert.val = r-dist;
				int off = i + dim*(j + dim *k);
				levelset[off] = vert; //Negative on the inside	
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

	initLevelset(128,0.3);

	levelset = new mp4Vector[dimx*dimy*dimz];
	for(int k = 0; k < dimz; ++k)
		for(int j = 0; j < dimy; ++j)
			for(int i = 0; i < dimx; ++i)
			{
				levelset[(i + dimx*(j + dimy *k))] = mp4Vector(i,j,k,1000.0);
			}

	//tri = MarchingCubes(128,128,128,0,levelset,LinearInterp,numOfTriangles);
	
	
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
			for(int k = 0; k < dimz; ++k)
				for(int j = 0; j < dimy; ++j)
					for(int i = 0; i < dimx; ++i)
					{
						levelset[(i + dimx*(j + dimy *k))] = mp4Vector(i,j,k,1000.0);
					}
			fluid_solver.step_frame();
// 			stopwatch.stopTimer();
// 			std::cout << std::scientific;
// 			std::cout << stopwatch.getElapsedTime() << "\n";
// 			if(numframes >= 6)
// 				avgtime += stopwatch.getElapsedTime();
// 			numframes++;
			
		//	write_paricle_pos_binary(fluid_solver.particles);
			delete [] tri;		
			tri = MarchingCubes(dimx,dimy,dimz,0,levelset,LinearInterp,numOfTriangles);
		}

		

		OpenGl_updateParticleLocation(fluid_solver.particles.pos, sizeof(vec3f)*fluid_solver.particles.currnp);
		OpenGl_updateParticleVelocity(fluid_solver.particles.vel,sizeof(vec3f)*fluid_solver.particles.currnp);
		
		

	}
	

	TerminateViewer();
// 	std::cout << std::scientific;
// 	std::cout << "Avg. time" << avgtime/(numframes-6) << "\n";
// 	std::cout << "Press any key to quit...\n";
// 	std::cin.get();

	

	
	delete [] levelset;
	return 0;
}



