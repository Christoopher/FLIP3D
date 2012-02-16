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

	Fluid_Solver(int dimx, int dimy, int dimz, float h, float timestep, float gravity,float rho, int max_particles) 
		: dimx(dimx), dimy(dimy), dimz(dimz), timestep(timestep)
	{
		grid.init(dimx,dimy,dimz,h,gravity, rho);
		particles.init(max_particles,grid);
	}

	void reset();

	void step_frame();
	void step(float dt);

	void init_box();

};

void Fluid_Solver::reset()
{
	particles.clear();
	init_box();

}

void Fluid_Solver::init_box()
{
	srand ( time(NULL) );
	float r1,r2,r3;
	float subh = grid.h/2.0f;
	vec3f pos(0);
	for(float k = 1; k < 31; ++k)
		for(float j = 1; j < 31; ++j)
			for(float i = 1; i < 10; ++i)
			{
				for (int kk = -1; kk < 1; ++kk)
					for(int jj = -1; jj < 1; ++jj)
						for (int ii = -1; ii < 1; ++ii)
						{
							r1 = float(rand()) / RAND_MAX - 0.5; //[-0.5, 0.5]
							r2 = float(rand()) / RAND_MAX - 0.5;
							r3 = float(rand()) / RAND_MAX - 0.5;
							pos[0] = (i+0.5f)*grid.h + (ii + 0.5f + 0.95*r1)*subh;
							pos[1] = (j+0.5f)*grid.h + (jj + 0.5f + 0.95*r2)*subh;
							pos[2] = (k+0.5f)*grid.h + (kk + 0.5f + 0.95*r3)*subh;
							add_particle(particles,pos,vec3f(0.0f));
						}

			}

// 			for(float k = 1; k < 31; ++k)
// 				for(float j = 1; j < 31; ++j)
// 					for(float i = 59; i < 63; ++i)
// 					{
// 						for (int kk = -1; kk < 1; ++kk)
// 							for(int jj = -1; jj < 1; ++jj)
// 								for (int ii = -1; ii < 1; ++ii)
// 								{
// 									r1 = float(rand()) / RAND_MAX - 0.5; //[-0.5, 0.5]
// 									r2 = float(rand()) / RAND_MAX - 0.5;
// 									r3 = float(rand()) / RAND_MAX - 0.5;
// 									pos[0] = (i+0.5f)*grid.h + (ii + 0.5f + 0.95*r1)*subh;
// 									pos[1] = (j+0.5f)*grid.h + (jj + 0.5f + 0.95*r2)*subh;
// 									pos[2] = (k+0.5f)*grid.h + (kk + 0.5f + 0.95*r3)*subh;
// 									add_particle(particles,pos,vec3f(0.0f));
// 								}
// 
// 					}


// 			for(float k = 20; k < 30; ++k)
// 				for(float j = 10; j < 30; ++j)
// 					for(float i = 20; i < 30; ++i)
// 					{
// 						for (int kk = -1; kk < 1; ++kk)
// 							for(int jj = -1; jj < 1; ++jj)
// 								for (int ii = -1; ii < 1; ++ii)
// 								{
// 									r1 = float(rand()) / RAND_MAX - 0.5; //[-0.5, 0.5]
// 									r2 = float(rand()) / RAND_MAX - 0.5;
// 									r3 = float(rand()) / RAND_MAX - 0.5;
// 									pos[0] = (i+0.5f)*grid.h + (ii + 0.5f + 0.95*r1)*subh;
// 									pos[1] = (j+0.5f)*grid.h + (jj + 0.5f + 0.95*r2)*subh;
// 									pos[2] = (k+0.5f)*grid.h + (kk + 0.5f + 0.95*r3)*subh;
// 									add_particle(particles,pos,vec3f(0.0f));
// 								}
// 
// 					}
}

void Fluid_Solver::step_frame()
{
	static int frame = 0;
	for(float elapsed = 0; elapsed < timestep;)
	{

		float dt = grid.CFL();
		if( dt > timestep - elapsed)
			dt = timestep - elapsed;

		elapsed += dt;

		step(dt);
	}
	frame++;
}

void Fluid_Solver::step(float dt)
{

	for (int i = 0; i < 5; i++)
		move_particles_in_grid(particles,grid,0.2*dt);

	grid.zero();
	transfer_to_grid(particles,grid);
	grid.save_velocities();
	grid.add_gravity(dt);
	grid.apply_boundary_conditions();
	grid.classify_voxel();
	grid.form_poisson(dt); //Funkar med hög sannolikhet
	grid.calc_divergence();
	grid.solve_pressure(100,1e-6);
	grid.project(dt);
	//grid.extend_velocity();
	grid.apply_boundary_conditions();
	grid.get_velocity_update();
	update_from_grid(particles,grid);
	

	
}












#endif