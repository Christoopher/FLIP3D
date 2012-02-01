#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "Particles.h"
#include "Grid.h"
#include "Sparse_Matrix.h"

struct Fluid_Solver
{
	Particles particles;
	Grid grid;
	Sparse_Matrix poisson;
	VectorN r;

	int dimx,dimy,dimz;
	float timestep;

	Fluid_Solver(int dimx, int dimy, int dimz, float timestep, float gravity, int max_particles) 
		: dimx(dimx), dimy(dimy), dimz(dimz), timestep(timestep)
	{
		grid.init(dimx,dimy,dimz,1.0/dimx,gravity);
		particles.init(max_particles,grid);
		poisson.init(dimx,dimy,dimz);
		r.init(dimx*dimy*dimz);
	}

	void step_n(int n);
	void step_once();

	void form_poisson();

	//void reset(); Later
};

void Fluid_Solver::form_poisson()
{

}

void Fluid_Solver::step_once()
{
	int iterations = 0;
	for(float elapsed = 0; elapsed < timestep;)
	{
		transfer_to_grid(particles,grid);
		//Fetch max timestep	
		float dt = grid.CFL();

		if( dt > timestep - elapsed)
			dt = timestep - elapsed;

		elapsed += dt;

		grid.save_velocities();
		grid.add_gravity(dt);

		/*
		form_poisson()
		CG.solve(max_iterations, error);
		project();

		*/


		grid.apply_boundary_conditions();
		grid.get_velocity_update();
		update_from_grid(particles,grid);
		advect_particles(particles,dt);
		iterations++;
	}

}

void Fluid_Solver::step_n(int n)
{
	for (int i = 0; i < n; i++)
	{
		step_once();
		//write_to_disk()
	}		
}










#endif