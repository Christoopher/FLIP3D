
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
const int dimx = 32, dimy = 32, dimz = 32;
const float gridh = 0.1;

void initVoxels(float * voxelPositions, int dx, int dy, int dz)
{
	//OBS Ordningen på loopen är viktig här
	float * posItr = voxelPositions;
	for (int z = 0; z < dz; ++z) {
		for (int y = 0; y < dy; ++y) {
			for (int x = 0; x < dx; ++x) {
				*posItr++ = (float)x;
				*posItr++ = (float)y;
				*posItr++ = (float)z;	  
			}
		}
	}
}

void update_voxel_flags(Grid & grid, Array3f & flags)
{

	for(int k = 1; k < grid.Nz-1; ++k)
		for(int j = 1; j < grid.Ny-1; ++j)
			for(int i = 1; i < grid.Nx-1; ++i)
			{
				flags(i,j,k) = (float)grid.marker(i,j,k);
			}
}


/*
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
*/

CStopWatch stopwatch;

int numframes;
double avgtime = 0;




void
runFluidSim()
{
	Fluid_Solver fluid_solver(dimx,dimy,dimz,0.1,1.0f/30.0f,9.82f,1.0f,Nparticles);
	fluid_solver.init_box();

	OpenGl_initViewer(600, 600, dimx, dimy, dimz, gridh);
	OpenGl_initParticles(&fluid_solver.particles.pos[0], &fluid_solver.particles.vel[0], sizeof(vec3f)*fluid_solver.particles.currnp, fluid_solver.particles.currnp);	
	
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
			
			//write_paricle_pos_binary(fluid_solver.particles);

			//fluid_solver.createSurface();
			//openGl_setMesh(fluid_solver.tri,fluid_solver.nrofTriangles);
		}

		OpenGl_updateParticles(fluid_solver.particles);
	}
	

	TerminateViewer();
// 	std::cout << std::scientific;
// 	std::cout << "Avg. time" << avgtime/(numframes-6) << "\n";
// 	std::cout << "Press any key to quit...\n";
//	std::cin.get();
}

void
runSurfaceReconstruction(int frame)
{
	Particles particles;
	TRIANGLE * tri;
	int nrofTriangles;
	
	OpenGl_initViewer(600, 600, dimx, dimy, dimz, gridh);

	
	while(running) {

		
		/*if(reset)
			fluid_solver.reset();
		reset = false;*/

		OpenGl_drawAndUpdate(running);

		if(step || play)
		{
			
			tri = new TRIANGLE[1];
			read_paricle_pos_binary(particles, frame);

			mesh(particles,dimx, dimy, dimz, h, 1, nrofTriangles, tri);
			openGl_setMesh(tri, nrofTriangles);
		}
	}
	

	TerminateViewer();

}



int main(void)
{
	//runFluidSim();
	
	runSurfaceReconstruction(1); //Read specific frame
	
	return 0;
}



