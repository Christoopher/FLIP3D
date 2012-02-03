#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "Particles.h"
#include "Grid.h"


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Unconditioned_CG_Solver.h"

struct Fluid_Solver
{
	Particles particles;
	Grid grid;
	
	int dimx,dimy,dimz;
	float timestep;

	Fluid_Solver(int dimx, int dimy, int dimz, float timestep, float gravity,float rho, int max_particles) 
		: dimx(dimx), dimy(dimy), dimz(dimz), timestep(timestep)
	{
		grid.init(dimx,dimy,dimz,1.0f/dimx,gravity, rho);
		particles.init(max_particles,grid);
	}

	void reset();

	void step_frame();
	void step(float dt);

	void init_box(int dim = 1, int particles_per_cell = 1);

};

void Fluid_Solver::reset()
{
	particles.clear();
	init_box();

}

void Fluid_Solver::init_box(int dim, int particles_per_cell)
{
	const static int pcel = particles_per_cell;
	const static int d = dim;
	int halfz = dimz / 2;
	int halfx = dimx / 2;
	float hdim = d/2.0f;
	

	srand ( time(NULL) );
	float r1,r2,r3;
	for(float k = halfz-hdim; k < halfz + hdim; ++k)
		for(float j = dimz-d; j < dimy; ++j)
			for(float i = halfx-hdim; i < halfx+hdim; ++i)
			{
				for (int p = 0; p < pcel; ++p)
				{
					r1 = 0.99f*(float(rand()) / RAND_MAX - 0.5);
					r2 = 0.99f*(float(rand()) / RAND_MAX - 0.5);
					r3 = 0.99f*(float(rand()) / RAND_MAX - 0.5);
					
					add_particle(particles,vec3f(i+r1,j+r2,k+r3)*grid.h,vec3f(0.0f));
				}

			}
}

void Fluid_Solver::step_frame()
{
	for(float elapsed = 0; elapsed < timestep;)
	{

		float dt = grid.CFL();
		if( dt > timestep - elapsed)
			dt = timestep - elapsed;

		elapsed += dt;

		step(dt);
	}
}

void Fluid_Solver::step(float dt)
{

	transfer_to_grid(particles,grid);

	grid.save_velocities();
	grid.add_gravity(dt);
	grid.apply_boundary_conditions();
	//grid.classify_voxel();
	grid.form_poisson(dt);
	grid.calc_divergence();
	grid.solve_pressure(100,10e-6);
	grid.project(dt);
	grid.apply_boundary_conditions();
	grid.get_velocity_update();
	update_from_grid(particles,grid);

	for (int i = 0; i < 5; i++)
		advect_particles(particles,0.2*dt);
}












#endif