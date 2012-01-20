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
		pos = new vec3f[maxnp];
		vel = new vec3f[maxnp];
	}

	~Particles() 
	{
		delete [] vel;
		delete [] pos;
	}
};

void advect_particles(Particles & particles, float dt);
void write_to_file(const char * filename);
void transfer_to_grid(Particles & particles, Grid & grid);
void update_from_grid(Particles & particles, Grid & grid);
void add_particle(Particles & particles, vec3f pos, vec3f vel);
void get_position_larray(Particles & particles, float posArray[]);
void get_velocity_larray(Particles & particles, float velArray[]);





#endif