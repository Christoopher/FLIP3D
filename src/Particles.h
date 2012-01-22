#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include "Vector3.h"
#include "Grid.h"

struct Particles
{
	//maximum nr of particles and the number of particles in use
	int maxnp, currnp; 
	vec3f *vel, *pos;

	Particles(int maxParticles)
	{
		maxnp = maxParticles;
		currnp = 0;
		pos = new vec3f[maxnp];
		vel = new vec3f[maxnp];
	}

	~Particles() 
	{
		delete [] vel;
		delete [] pos;
	}
};

// void write_to_file(const char * filename);
// void transfer_to_grid(Particles & particles, Grid & grid);
// void update_from_grid(Particles & particles, Grid & grid);

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
		*itr++ = particles.pos[i].v[0];
		*itr++ = particles.pos[i].v[1];
		*itr++ = particles.pos[i].v[2];

	}
}



#endif