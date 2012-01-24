#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include "Vector3.h"
#include "Grid.h"
#include "Array3D.h"



struct Particles
{
	//maximum nr of particles and the number of particles in use
	int maxnp, currnp; 
	vec3f *vel, *pos;
	Array3f weightsumx,weightsumy,weightsumz;

	Particles(int maxParticles, Grid & grid)
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

	
};



// void write_to_file(const char * filename);

void update_from_grid(Particles & particles, Grid & grid)
{
	int p;
	int i, ui, j, vj, k, wk;
	float fx, ufx, fy, vfy, fz, wfz;
	for (int p = 0; p < particles.currnp; ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0], ui, ufx);
		grid.bary_x_centre(particles.pos[p][0], i, fx);

		grid.bary_y(particles.pos[p][1], vj, vfy);
		grid.bary_y_centre(particles.pos[p][1], j, fy);

		grid.bary_z(particles.pos[p][2], wk, wfz);
		grid.bary_z_centre(particles.pos[p][2], k, fz);

		//Flip
		particles.vel[p] += vec3f(grid.du.trilerp(ui,j,k,ufx,fy,fz), grid.dv.trilerp(i,vj,k,fx,vfy,fz), grid.dw.trilerp(i,j,wk,fx,fy,wfz));
		
		//Pic
		//particles.vel[p] = vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz)); 
	}
}

void accumulate(Array3f &accum,Array3f &sum, float vel, int i, int j, int k, float fx, float fy, float fz)
{
	float weight;

	weight = (1-fx)*(1-fy)*(1-fz);
	accum(i,j,k) += weight*vel;
	sum(i,j,k) += weight;

	weight = fx*(1-fy)*(1-fz);
	accum(i+1,j,k) += weight*vel;
	sum(i+1,j,k) += weight;

	weight = (1-fx)*fy*(1-fz);
	accum(i,j+1,k) += weight*vel;
	sum(i,j+1,k) += weight;

	weight = fx*fy*(1-fz);
	accum(i+1,j+1,k) += weight*vel;
	sum(i+1,j+1,k) += weight;

	weight = (1-fx)*(1-fy)*fz;
	accum(i,j,k+1) += weight*vel;
	sum(i,j,k+1) += weight;

	weight = fx*(1-fy)*fz;
	accum(i+1,j,k+1) += weight*vel;
	sum(i+1,j,k+1) += weight;

	weight = (1-fx)*fy*fz;
	accum(i,j+1,k+1) += weight*vel;
	sum(i,j+1,k+1) += weight;

	weight = fx*fy*fz;
	accum(i+1,j+1,k+1) += weight*vel;
	sum(i+1,j+1,k+1) += weight;


}

void transfer_to_grid(Particles & particles, Grid & grid)
{
	int vi, vj, vk,i,j,k;
	float fx, fy, fz;
	grid.marker.zero();
	for (int p = 0; p < particles.currnp; ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0],vi,fx);
		i = vi;
		grid.bary_y_centre(particles.pos[p][1],vj,fy);
		grid.bary_z_centre(particles.pos[p][2],vk,fz);
		accumulate(grid.u,particles.weightsumx,particles.vel[p][0],vi,vj,vk,fx,fy,fz);
		

		grid.bary_x_centre(particles.pos[p][0],vi,fx);
		grid.bary_y(particles.pos[p][1],vj,fy);
		j = vj;
		grid.bary_z_centre(particles.pos[p][2],vk,fz);
		accumulate(grid.v,particles.weightsumy,particles.vel[p][1],vi,vj,vk,fx,fy,fz);

		
		grid.bary_x_centre(particles.pos[p][0],vi,fx);
		grid.bary_y_centre(particles.pos[p][1],vj,fy);
		grid.bary_z(particles.pos[p][2],vk,fz);
		k = vk;
		accumulate(grid.w,particles.weightsumz,particles.vel[p][2],vi,vj,vk,fx,fy,fz);

		grid.marker(i,j,k) = FLUIDCELL;
	}
	
	//Scale u velocities with weightsumx
	for (int i = 0; i < grid.u.size; i++)
	{
		if(grid.u.data[i] > 0)
			grid.u.data[i] /= particles.weightsumx.data[i];
	}

	//Scale v velocities with weightsumy
	for (int i = 0; i < grid.v.size; i++)
	{
		if(grid.v.data[i] > 0)
			grid.v.data[i] /= particles.weightsumy.data[i];
	}

	//Scale w velocities with weightsumz
	for (int i = 0; i < grid.w.size; i++)
	{
		if(grid.w.data[i] > 0)
			grid.w.data[i] /= particles.weightsumz.data[i];
	}

	
}

//----------------------------------------------------------------------------//
// Advect particles with stored velocity
//----------------------------------------------------------------------------//
void advect_particles(Particles & particles, float dt)
{
	for (int i = 0; i < particles.currnp; ++i)
	{
		particles.pos[i] += particles.vel[i]*dt;
	}

// 	vec3f * velItr = particles.vel;
// 	for (vec3f * ppos = particles.pos; ppos < ppos + particles.currnp; ppos++ )
// 	{
// 		*ppos += (*velItr++)*dt;
// 	}
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

//----------------------------------------------------------------------------//
// Converts particle velocities into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_velocity_larray(Particles & particles, float velArray[])
{
	float * itr = velArray;

	for (int i = 0; i<particles.currnp; i++)
	{
		*itr++ = particles.vel[i].v[0];
		*itr++ = particles.vel[i].v[1];
		*itr++ = particles.vel[i].v[2];

	}
	
// 	for (vec3f * pvel = particles.vel; pvel < pvel + particles.currnp; pvel++ )
// 	{
// 		*itr++ = (*pvel)[0];
// 		*itr++ = (*pvel)[1];
// 		*itr++ = (*pvel)[2];		
// 	}
}

//----------------------------------------------------------------------------//
// Converts particle positions into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_position_larray(Particles & particles, float posArray[])
{
	float * itr = posArray;
	int i = 0;
// 	for (vec3f * ppos = particles.pos; ppos < ppos + particles.currnp; ++ppos )
// 	{
// 		*itr++ = (*ppos)[0];
// 		*itr++ = (*ppos)[1];
// 		*itr++ = (*ppos)[2];
// 		i++;
// 	}
	for (int i = 0; i<particles.currnp; i++)
	{
		*itr++ = particles.pos[i].v[0]*10;
		*itr++ = particles.pos[i].v[1]*10;
		*itr++ = particles.pos[i].v[2]*10;

	}
}



#endif