#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <iostream>
#include <sstream>
#include <fstream>


#include "Vector3.h"
#include "Grid.h"
#include "Array3D.h"
#include <omp.h>

#include "OpenGLViewer.h"

struct Particles
{
	//maximum nr of particles and the number of particles in use
	int maxnp, currnp; 
	vec3f *vel, *pos;
	Array3f weightsumx,weightsumy,weightsumz;

	Particles() {}

	Particles(int maxParticles,Grid & grid)
	{
		init(maxParticles,grid);
	}

	void init(int maxParticles, Grid & grid)
	{
		maxnp = maxParticles;
		currnp = 0;
		pos = new vec3f[maxnp];
		vel = new vec3f[maxnp];
		weightsumx.init(grid.u.nx, grid.u.ny, grid.u.nz);
		weightsumy.init(grid.v.nx, grid.v.ny, grid.v.nz);
		weightsumz.init(grid.w.nx, grid.w.ny, grid.w.nz);
	}

	~Particles() 
	{
		delete [] vel;
		delete [] pos;
	}

	void clear()
	{
		std::memset(vel, 0, maxnp*sizeof(vec3f));
		std::memset(pos, 0, maxnp*sizeof(vec3f));
		currnp = 0;
	}


};

void move_particles_in_grid(Particles & particles, Grid & grid, float dt)
{
	vec3f vel;
	int ui, i, vj, j, wk, k;
	float ufx, fx, vfy, fy, wfz, fz;

	float xmax = (grid.Nx-1.001)*grid.h, xmin = 1.001*grid.h;
	float ymax = (grid.Ny-1.001)*grid.h, ymin = 1.001*grid.h;
	float zmax = (grid.Nz-1.001)*grid.h, zmin = 1.001*grid.h;

#pragma omp parallel for private(ufx, fx, vfy, fy, wfz, fz,ui, i, vj, j, wk, k,vel) shared(xmax,ymax,zmax)
	for(int p = 0; p < particles.currnp; p++)
	{
		//Trilerp from grid
		grid.bary_x(particles.pos[p][0], ui, ufx);
		grid.bary_x_centre(particles.pos[p][0], i, fx);

		grid.bary_y(particles.pos[p][1], vj, vfy);
		grid.bary_y_centre(particles.pos[p][1], j, fy);

		grid.bary_z(particles.pos[p][2], wk, wfz);
		grid.bary_z_centre(particles.pos[p][2], k, fz);

		vel = vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz));			

		//Move particle one step with forward euler
		particles.pos[p] += dt*vel;

		clamp(particles.pos[p][0], xmin,xmax);
		clamp(particles.pos[p][1], ymin,ymax);
		clamp(particles.pos[p][2], zmin,zmax);

		//particles.pos[p] = newpos;
		
// 		if(grid.marker(floor(newpos[0]), floor(newpos[1]), floor(newpos[2])) == SOLIDCELL)
// 		{
// 			//Particles has enetered a solid voxel. clamp back to surface.
// 
// 			//Direction
// 			//vec3f dir = newpos - particles.pos[p];
// 			//float length = mag(dir);
// 			newpos = particles.pos[p];// + 0.05*dir;
// 
// 		}
// 		//else
// 		particles.pos[p] = newpos;
			
		
		//Clamp pos to be inside of the solid walls
		//clamp(particles.pos[p][0], xmin,xmax);
		//clamp(particles.pos[p][1], ymin,ymax);
		//clamp(particles.pos[p][2], zmin,zmax);


	}
}

// void write_to_file(const char * filename);

void update_from_grid(Particles & particles, Grid & grid)
{
	int i, ui, j, vj, k, wk;
	float fx, ufx, fy, vfy, fz, wfz;
#pragma omp parallel for private(fx, ufx, fy, vfy, fz, wfz, i, ui, j, vj, k, wk)
	for (int p = 0; p < particles.currnp; ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0], ui, ufx);
		grid.bary_x_centre(particles.pos[p][0], i, fx);

		grid.bary_y(particles.pos[p][1], vj, vfy);
		grid.bary_y_centre(particles.pos[p][1], j, fy);

		grid.bary_z(particles.pos[p][2], wk, wfz);
		grid.bary_z_centre(particles.pos[p][2], k, fz);

		//Pic
		//particles.vel[p] = vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz)); 
		//Flip
		//particles.vel[p] += vec3f(grid.du.trilerp(ui,j,k,ufx,fy,fz), grid.dv.trilerp(i,vj,k,fx,vfy,fz), grid.dw.trilerp(i,j,wk,fx,fy,wfz));
		
		//PIC/FLIP
		float alpha = 0.05f;
		particles.vel[p] =  alpha*vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz))
			+ (1.0f - alpha)*(particles.vel[p] + vec3f(grid.du.trilerp(ui,j,k,ufx,fy,fz), grid.dv.trilerp(i,vj,k,fx,vfy,fz), grid.dw.trilerp(i,j,wk,fx,fy,wfz)));

	}
}

void accumulate(Array3f &macvel,Array3f &sum, float & pvel, int i, int j, int k, float fx, float fy, float fz)
{
	float weight = 0;

	weight = (1-fx)*(1-fy)*(1-fz);
	macvel(i,j,k) += weight*pvel;
	sum(i,j,k) += weight;

	weight = fx*(1-fy)*(1-fz);
	macvel(i+1,j,k) += weight*pvel;
	sum(i+1,j,k) += weight;

	weight = (1-fx)*fy*(1-fz);
	macvel(i,j+1,k) += weight*pvel;
	sum(i,j+1,k) += weight;

	weight = fx*fy*(1-fz);
	macvel(i+1,j+1,k) += weight*pvel;
	sum(i+1,j+1,k) += weight;

	weight = (1-fx)*(1-fy)*fz;
	macvel(i,j,k+1) += weight*pvel;
	sum(i,j,k+1) += weight;

	weight = fx*(1-fy)*fz;
	macvel(i+1,j,k+1) += weight*pvel;
	sum(i+1,j,k+1) += weight;

	weight = (1-fx)*fy*fz;
	macvel(i,j+1,k+1) += weight*pvel;
	sum(i,j+1,k+1) += weight;

	weight = fx*fy*fz;
	macvel(i+1,j+1,k+1) += weight*pvel;
	sum(i+1,j+1,k+1) += weight;


}

void transfer_to_grid(Particles & particles, Grid & grid)
{
	int ui, vj, wk,i,j,k;
	float fx, ufx, fy, vfy, fz, wfz;
	int tmpi, tmpj, tmpk;

	particles.weightsumx.zero();
	particles.weightsumy.zero();
	particles.weightsumz.zero();

#pragma omp parallel for private(ui, vj, wk,i,j,k,fx, ufx, fy, vfy, fz, wfz,tmpi, tmpj, tmpk)
	for (int p = 0; p < particles.currnp; ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0],ui,ufx); tmpi = ui;
		grid.bary_y_centre(particles.pos[p][1],j,fy);
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.u,particles.weightsumx,particles.vel[p][0],ui,j,k,ufx,fy,fz);


		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_y(particles.pos[p][1],vj,vfy); tmpj = vj;
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.v,particles.weightsumy,particles.vel[p][1],i,vj,k,fx,vfy,fz);


		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_y_centre(particles.pos[p][1],j,fy);
		grid.bary_z(particles.pos[p][2],wk,wfz); tmpk = wk;
		accumulate(grid.w,particles.weightsumz,particles.vel[p][2],i,j,wk,fx,fy,wfz);

		grid.marker(tmpi,tmpj,tmpk) = FLUIDCELL;
		mp4Vector val(tmpi,tmpj,tmpk,-1.0);
		levelset[(i + grid.Nx*(j + grid.Ny*k))] = val;

	}


	//Scale u velocities with weightsumx

	#pragma omp parallel for
	for (int i = 0; i < grid.u.size; i++)
	{
		if(grid.u.data[i] != 0)
			grid.u.data[i] /= particles.weightsumx.data[i];
	}

	//Scale v velocities with weightsumy
	#pragma omp parallel for
	for (int i = 0; i < grid.v.size; i++)
	{
		if(grid.v.data[i] != 0)
			grid.v.data[i] /= particles.weightsumy.data[i];
	}

	//Scale w velocities with weightsumz
	#pragma omp parallel for
	for (int i = 0; i < grid.w.size; i++)
	{
		if(grid.w.data[i] != 0)
			grid.w.data[i] /= particles.weightsumz.data[i];
	}

}

//----------------------------------------------------------------------------//
// Adds a particle to the particles struct
//----------------------------------------------------------------------------//
inline void add_particle(Particles & particles, vec3f & pos, vec3f  & vel)
{
	particles.pos[particles.currnp] = pos;
	particles.vel[particles.currnp] = vel;
	++particles.currnp;
}

void write_paricle_pos_binary(Particles & particles)
{
	static int frame = 0;
	std::stringstream ss;
	ss << "frame_" << frame << ".dat";
	std::string str = ss.str();
	FILE *file = fopen(str.c_str(),"wb");

	for (int i = 0; i < particles.currnp; ++i)
	{
		fwrite((const void *) (particles.pos[i].v),sizeof(float),3,file);
	}

	fclose(file);

	frame++;
}

#endif