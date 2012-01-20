#include "Particles.h"

//----------------------------------------------------------------------------//
// Converts particle positions into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_position_larray(Particles & particles, float posArray[])
{
	float * itr = posArray;
	for (vec3f * ppos = particles.pos; ppos < ppos + particles.currnp; ppos++ )
	{
		*itr++ = (*ppos)[0];
		*itr++ = (*ppos)[1];
		*itr++ = (*ppos)[2];
	}
}

//----------------------------------------------------------------------------//
// Converts particle velocities into a linear array to use with e.g OpenGL
//----------------------------------------------------------------------------//
void get_velocity_larray(Particles & particles, float velArray[])
{
	float * itr = velArray;
	for (vec3f * pvel = particles.vel; pvel < pvel + particles.currnp; pvel++ )
	{
		*itr++ = (*pvel)[0];
		*itr++ = (*pvel)[1];
		*itr++ = (*pvel)[2];
	}
}

//----------------------------------------------------------------------------//
// Adds a particle to the particles struct
//----------------------------------------------------------------------------//
inline void add_particle(Particles & particles, vec3f pos, vec3f vel)
{
	particles.pos[particles.currnp] = pos;
	particles.vel[particles.currnp] = vel;
	++particles.currnp;
}

//----------------------------------------------------------------------------//
// Advect particles with stored velocity
//----------------------------------------------------------------------------//
void advect_particles(Particles & particles, float dt)
{
	vec3f * velItr = particles.vel;
	for (vec3f * ppos = particles.pos; ppos < ppos + particles.currnp; ppos++ )
	{
		*ppos += (*velItr++)*dt;
	}
}





