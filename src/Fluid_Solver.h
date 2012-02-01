#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "Particles.h"
#include "Grid.h"
#include "Sparse_Matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct Fluid_Solver
{
	Particles particles;
	Grid grid;
	Sparse_Matrix poisson;
	VectorN rhs;

	int dimx,dimy,dimz;
	float timestep, rho;

	Fluid_Solver(int dimx, int dimy, int dimz, float timestep, float gravity,float rho, int max_particles) 
		: dimx(dimx), dimy(dimy), dimz(dimz), rho(rho), timestep(timestep)
	{
		grid.init(dimx,dimy,dimz,1.0f/dimx,gravity);
		particles.init(max_particles,grid);
		poisson.init(dimx,dimy,dimz);
		rhs.init(dimx*dimy*dimz);
	}

	void step_n(int n);
	void step_once();

	void form_poisson(float dt);
	void calc_divergence();
	void project(float dt);

	void init_box(int dim, int particles_per_cell);

};

void Fluid_Solver::init_box(int dim, int particles_per_cell)
{
	float dist = 0;
	srand ( time(NULL) );
	float r1,r2,r3;
	for(int k = 3; k < 8; ++k)
		for(int j = 3; j < 8; ++j)
			for(int i = 3; i < 8; ++i)
			{
				for (int p = 0; p < particles_per_cell; ++p)
				{
					r1 = 0.49*(float(rand()) / RAND_MAX - 0.5);
					r2 = 0.49*(float(rand()) / RAND_MAX - 0.5);
					r3 = 0.49*(float(rand()) / RAND_MAX - 0.5);
					add_particle(particles,vec3f(i+r1,j+r2,k+r3)*grid.h,vec3f(0.0f));
				}

			}
}

void Fluid_Solver::project(float dt)
{
	float scale = grid.overh * dt / rho;
	int offset;
	float val;
	for(int k = 0; k < dimz; ++k)
		for(int j = 0; j < dimy; ++j)
			for(int i = 0; i < dimx; ++i)
			{
				if(grid.marker(i,j,k) == FLUIDCELL)
				{
					offset = i + dimx*(j + dimy*k); //Offset into rhs.data array
					val = scale * float(rhs.data[offset]);
					
					grid.u(i,j,k) -= val;
					grid.u(i+1,j,k) += val;

					grid.v(i,j,k) -= val;
					grid.v(i,j+1,k) += val;

					grid.w(i,j,k) -= val;
					grid.w(i,j+1,k) += val;

				}
				else if(grid.marker(i,j,k) == SOLIDCELL)
				{
					grid.u(i,j,k) = 0;
					grid.u(i+1,j,k) = 0;

					grid.v(i,j,k) = 0;
					grid.v(i,j+1,k) = 0;

					grid.w(i,j,k) = 0;
					grid.w(i,j,k+1) = 0;
				}
			}
}

void Fluid_Solver::calc_divergence()
{
	rhs.zero();
	int offset;
	float scale = grid.overh;
	for(int k = 0; k < dimz; ++k)
		for(int j = 0; j < dimy; ++j)
			for(int i = 0; i < dimx; ++i)
			{
				if(grid.marker(i,j,k) == FLUIDCELL)
				{
					offset = i + dimx*(j + dimy*k); //Offset into rhs.data array
					
					if(grid.marker(i+1,j,k) == SOLIDCELL)		//If cell(i+1,j,k) remove u(i+1/2,j,k)
						rhs.data[offset] -= grid.u(i-1,j,k);
					else if(grid.marker(i-1,j,k) == SOLIDCELL)	//If cell(i-1,j,k) remove u(i-1/2,j,k)
						rhs.data[offset] += grid.u(i+1,j,k);
					else
						rhs.data[offset] += grid.u(i+1,j,k) - grid.u(i-1,j,k);

					if(grid.marker(i,j+1,k) == SOLIDCELL)		//If cell(i,j+1,k) remove u(i,j+1/2,k)
						rhs.data[offset] -= grid.v(i,j-1,k);
					else if(grid.marker(i-1,j,k) == SOLIDCELL)	//If cell(i,j-1,k) remove u(i,j-1/2,k)
						rhs.data[offset] += grid.v(i,j+1,k);
					else
						rhs.data[offset] += grid.v(i,j+1,k) - grid.v(i,j-1,k);

					if(grid.marker(i,j,k+1) == SOLIDCELL)		//If cell(i,j,k+1) remove u(i,j,k+1/2)
						rhs.data[offset] -= grid.w(i,j,k-1);
					else if(grid.marker(i-1,j,k) == SOLIDCELL)	//If cell(i,j,k-1) remove u(i,j,k-1/2)
						rhs.data[offset] += grid.w(i,j,k+1);
					else
						rhs.data[offset] += grid.w(i,j,k+1) - grid.w(i,j,k-1);
					
					rhs.data[offset] *= -scale;
				}
			}
}

void Fluid_Solver::form_poisson(float dt)
{
	float scale = grid.overh*grid.overh * timestep / rho;
	for(int k = 1; k < dimz-1; ++k)
		for(int j = 1; j < dimy-1; ++j)
			for(int i = 1; i < dimx-1; ++i)
				if(grid.marker(i,j,k) == FLUIDCELL)
				{
					if(grid.marker(i-1,j,k) != SOLIDCELL)		//Cell(i-1,j,k) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(grid.marker(i+1,j,k) != SOLIDCELL)		//Cell(i+1,j,k) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(grid.marker(i+1,j,k) == FLUIDCELL)	//Cell(i+1,j,k) Is fluid
							poisson(i,j,k,1) -= scale; 
					}

					if(grid.marker(i,j-1,k) != SOLIDCELL)		//Cell(i,j-1,k) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(grid.marker(i,j+1,k) != SOLIDCELL)		//Cell(i,j+1,k) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(grid.marker(i,j+1,k) == FLUIDCELL)	//Cell(i,j+1,k) Is fluid
							poisson(i,j,k,2) -= scale; 
					}

					if(grid.marker(i,j,k-1) != SOLIDCELL)		//Cell(i,j,k-1) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(grid.marker(i,j,k+1) != SOLIDCELL)		//Cell(i,j,k+1) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(grid.marker(i,j+1,k) == FLUIDCELL)	//Cell(i,j,k+1) Is fluid
							poisson(i,j,k,3) -= scale; 
					}					
				} //End if CELL(i,j,k) == FLUIDCELL
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