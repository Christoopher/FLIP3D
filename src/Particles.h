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

// void write_to_file(const char * filename);

void update_from_grid(Particles & particles, Grid & grid)
{
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

		//Pic
		//particles.vel[p] = vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz)); 
		//Flip
		//particles.vel[p] += vec3f(grid.du.trilerp(ui,j,k,ufx,fy,fz), grid.dv.trilerp(i,vj,k,fx,vfy,fz), grid.dw.trilerp(i,j,wk,fx,fy,wfz));
		
		//PIC/FLIP
		float alpha = 0.2f;
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
	
	particles.weightsumx.zero();
	particles.weightsumy.zero();
	particles.weightsumz.zero();

	for (int p = 0; p < particles.currnp; ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0],ui,ufx);
		grid.bary_y_centre(particles.pos[p][1],j,fy);
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.u,particles.weightsumx,particles.vel[p][0],ui,j,k,ufx,fy,fz);
		

		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_y(particles.pos[p][1],vj,vfy);
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.v,particles.weightsumy,particles.vel[p][1],i,vj,k,fx,vfy,fz);

		
		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_y_centre(particles.pos[p][1],j,fy);
		grid.bary_z(particles.pos[p][2],wk,wfz);
		accumulate(grid.w,particles.weightsumz,particles.vel[p][2],i,j,wk,fx,fy,wfz);

	}
	
	//Scale u velocities with weightsumx
	for (int i = 0; i < grid.u.size; i++)
	{
		if(grid.u.data[i] != 0)
			grid.u.data[i] /= particles.weightsumx.data[i];
	}

	//Scale v velocities with weightsumy
	for (int i = 0; i < grid.v.size; i++)
	{
		if(grid.v.data[i] != 0)
			grid.v.data[i] /= particles.weightsumy.data[i];
	}

	//Scale w velocities with weightsumz
	for (int i = 0; i < grid.w.size; i++)
	{
		if(grid.w.data[i] != 0)
			grid.w.data[i] /= particles.weightsumz.data[i];
	}

	for (int p = 0; p < particles.currnp; ++p)
	{
		grid.bary_x(particles.pos[p][0],i,fx);
		grid.bary_y(particles.pos[p][1],j,fy);
		grid.bary_z(particles.pos[p][2],k,fz);
		grid.marker(i,j,k) = FLUIDCELL;
	}
	
}

//----------------------------------------------------------------------------//
// Advect particles with stored velocity
//----------------------------------------------------------------------------//
void advect_particles(Particles & particles, float dt)
{
	vec3f * velItr = particles.vel;
	for (vec3f * ppos = particles.pos; ppos < particles.pos + particles.currnp; ++ppos )
	{
		*ppos += (*velItr++)*dt; // position_n+1 = position_n + velocity * dt
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

//----------------------------------------------------------------------------//
// Converts particle velocities into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_velocity_larray(Particles & particles, float velArray[])
{
	float * itr = velArray;
	for (vec3f * pvel = particles.vel; pvel < particles.vel + particles.currnp; ++pvel )
	{
		*itr++ = (*pvel)[0];
		*itr++ = (*pvel)[1];
		*itr++ = (*pvel)[2];		
	}
}

//----------------------------------------------------------------------------//
// Converts particle positions into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_position_larray(Particles & particles, float posArray[])
{
	float * itr = posArray;
	for (vec3f * ppos = particles.pos; ppos < particles.pos + particles.currnp; ++ppos )
	{
		*itr++ = (*ppos)[0]*11;
		*itr++ = (*ppos)[1]*11;
		*itr++ = (*ppos)[2]*11;
	}
}



#endif