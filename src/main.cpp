
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


//CStopWatch stopwatch;

int numframes;
double avgtime = 24;

void
runFluidSim()
{
	Fluid_Solver fluid_solver(dimx,dimy,dimz,gridh,1.0f/30.0f,9.82f,1.0f,Nparticles);
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
			
			write_paricle_pos_binary(fluid_solver.particles);

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


struct IFL
{
	struct Tri
	{
		int i,j,k;
		vec3f norm;
	};

	struct Normal
	{
		Normal() : n(0) {}
		vec3f norm;
		int n;
	};

	std::vector<Tri> tris;
	std::vector<Normal> normals;
	std::vector<vec3f> vecs;

	int insertPoint(vec3f & v)
	{
		int i;
		for(i = 0; i < vecs.size(); ++i)
		{
			if(v == vecs[i])
				return i;
		}
		vecs.push_back(v);
		return i;
	}

	void init(TRIANGLE * t, int n)
	{
		for(int i = 0; i < n; ++i)
		{
			vec3f v1(t[i].p[0].x,t[i].p[0].y,t[i].p[0].z);
			vec3f v2(t[i].p[1].x,t[i].p[1].y,t[i].p[1].z);
			vec3f v3(t[i].p[2].x,t[i].p[2].y,t[i].p[2].z);



			Tri tri;
			tri.i = insertPoint(v1);
			tri.j = insertPoint(v2);
			tri.k = insertPoint(v3);
			tri.norm = vec3f(t[i].norm.x,t[i].norm.y, t[i].norm.z);
			tris.push_back(tri);

		}

		std::cout << "Inserted " << vecs.size() << " unique points from " << n*3 << "original points\n";
		calcNormals();
	}

	void
	calcNormals()
	{
		normals.resize(vecs.size());
		for(int i = 0; i < tris.size(); ++i)
		{
			normals[tris[i].i].norm += vecs[tris[i].i];
			++normals[tris[i].i].n;
			
			normals[tris[i].j].norm += vecs[tris[i].j];
			++normals[tris[i].j].n;

			normals[tris[i].k].norm += vecs[tris[i].k];
			++normals[tris[i].k].n;
		}

		for(int i = 0; i < normals.size(); ++i)
		{
			normals[i].norm /= normals[i].n;
			normals[i].norm /= mag(normals[i].norm);
		}

	}
};

void
writeObj(TRIANGLE * tri, int & ntriangles, int frameNr)
{
	IFL ifl;
	ifl.init(tri,ntriangles);


	std::ofstream file;
	std::stringstream ss;
	ss << "objs/mesh_" << frameNr << ".obj";
	std::string filename;
	ss >> filename;
	file.open (filename);

	//Print points
	for(int i = 0; i < ifl.vecs.size(); ++i)
	{
		file << "v " << ifl.vecs[i][0] << ' ' << ifl.vecs[i][1] << ' ' << ifl.vecs[i][2] << "\n";
	}
	
	file << "\n";

	//Print normals
	/*for(int i = 0; i < ifl.vecs.size(); ++i)
	{
		file << "vn " << ifl.normals[i].norm[0] << ' ' << ifl.normals[i].norm[1] << ' ' << ifl.normals[i].norm[2] << "\n";
		
	}
	
	file << "\n";
	
	//print vecs
	for(int i = 0; i < ifl.tris.size(); ++i)
	{
		file << "f ";
		file << ifl.tris[i].i+1 << "\\\\" << ifl.tris[i].i+1 << ' ';
		file << ifl.tris[i].j+1 << "\\\\" << ifl.tris[i].j+1 << ' ';
		file << ifl.tris[i].k+1 << "\\\\" << ifl.tris[i].k+1 << ' ';
		file << "\n";
	}*/
	for(int i = 0; i < ifl.tris.size(); ++i)
	{
		file << "f ";
		file << ifl.tris[i].i+1 << ' ';
		file << ifl.tris[i].j+1 << ' ';
		file << ifl.tris[i].k+1;
		file << "\n";
	}
	file.close();
}

void
runSurfaceReconstruction(int frame)
{
	Particles particles;
	TRIANGLE * tri;
	int nrofTriangles;
	
	OpenGl_initViewer(600, 600, dimx, dimy, dimz, gridh);

	tri = new TRIANGLE[1];
	read_paricle_pos_binary(particles, frame);


	OpenGl_initParticles(&particles.pos[0], &particles.pos[0], sizeof(vec3f)*particles.currnp, particles.currnp);	
	
	while(running) {

		
		if(reset)
		{
			std::cout << "setting changed: press \'s\' to genereate surface \n";
			calcmesh = true;
		}
		reset = false;

		OpenGl_drawAndUpdate(running);

		if(step || play)
		{
			mesh(particles,dimx, dimy, dimz, gridh, 1, nrofTriangles, tri);
			openGl_setMesh(tri, nrofTriangles);
		}
	}
	
	writeObj(tri,nrofTriangles, frame);
	TerminateViewer();

}

void
runManySurfaceReconstructions(int frame_begin, int frame_end)
{
	Particles particles;
	TRIANGLE * tri;
	int nrofTriangles;

	tri = new TRIANGLE[1];

	for(int currFrame = frame_begin; currFrame <= frame_end; ++currFrame)
	{
		read_paricle_pos_binary(particles, currFrame);

		//if(reset)
		//{
			std::cout << "setting changed: press \'s\' to genereate surface \n";
			calcmesh = true;
		//}
		reset = false;

		//if(step || play)
		//{
			mesh(particles,dimx, dimy, dimz, gridh, 2, nrofTriangles, tri);
		//}
		writeObj(tri,nrofTriangles, currFrame);
		std::cout << "Wrote frame: " << currFrame << "\n\n\n";
	}
	
}

int main(void)
{
	//runFluidSim();
	
	runSurfaceReconstruction(35); //Read specific frame
	//runManySurfaceReconstructions(10, 40);
	
	return 0;
}



