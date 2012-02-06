#ifndef _GRID_H_
#define _GRID_H_

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

#include "Array3D.h"
#include "Sparse_Matrix.h"
#include "Unconditioned_CG_Solver.h"

struct Grid
{
	int Nx, Ny,Nz;
	float h,overh, gravity, rho;

	Array3f u,v,w,du,dv,dw; //Staggered u,v,w velocities
	Array3c marker; //Voxel classification
	Sparse_Matrix poisson; // The matrix for pressure stage
	VectorN rhs; //Right hand side of the poisson equation
	VectorN x; //Right hand side of the poisson equation
	Uncondioned_CG_Solver cg;

	Grid() {}
	
	Grid(int Nx_, int Ny_, int Nz_, float h_, float gravity_, float rho_) : Nx(Nx_), Ny(Ny_), Nz(Nz_), h(h_), overh(1.0f/h), gravity(gravity_), rho(rho_)
	{
		init(Nx_,Ny_,Nz_,h_,gravity_,rho_);
	}

	void init(int Nx_, int Ny_, int Nz_, float h_, float gravity_, float rho_)
	{
		Nx = Nx_; Ny = Ny_; Nz = Nz_; h = h_; gravity = gravity_;
		overh = 1.0f/h;
		rho = rho_;
		u.init(Nx_+1,Ny_,Nz_);
		v.init(Nx_,Ny_+1,Nz_);
		w.init(Nx_,Ny_,Nz_+1);
		du.init(Nx_+1,Ny_,Nz_);
		dv.init(Nx_,Ny_+1,Nz_);
		dw.init(Nx_,Ny_,Nz_+1);
		marker.init(Nx_,Ny_,Nz_);
		poisson.init(Nx_,Ny_,Nz_);
		rhs.init(Nx_,Ny_,Nz_);
		x.init(Nx_,Ny_,Nz_);
		cg.init(Nx_,Ny_,Nz_);
	}

	void zero()
	{
		u.zero();
		v.zero();
		w.zero();
		du.zero();
		dv.zero();
		dw.zero();
		poisson.zero();
		rhs.zero();
		x.zero();
		marker.zero();
	}

	void bary_x(float x, int &i, float &fx)
	{
		float sx=x*overh;
		i=(int)sx;
		fx=sx-floor(sx);
	}

	void bary_x_centre(float x, int &i, float &fx)
	{
		float sx=x*overh-0.5f;
		i=(int)sx;
		if(i<0){ 
			i=0; fx=0.0; 
		}
		else if(i>Nx-2) { 	//Why?
			i=Nx-2; fx=1.0; 
		}
		else{ 
			fx=sx-floor(sx); 
		}
	}

	void bary_y(float y, int &j, float &fy)
	{
		float sy=y*overh;
		j=(int)sy;
		fy=sy-floor(sy);
	}

	void bary_y_centre(float y, int &j, float &fy)
	{
		float sy=y*overh-0.5f;
		j=(int)sy;
		if(j<0){ j=0; fy=0.0; }
		else if(j>Ny-2){ j=Ny-2; fy=1.0; }
		else{ fy=sy-floor(sy); }
	}

	void bary_z(float z, int &k, float &fz)
	{
		float sz=z*overh;
		k=(int)sz;
		fz=sz-floor(sz);
	}

	void bary_z_centre(float z, int &k, float &fz)
	{
		float sz=z*overh-0.5f;
		k=(int)sz;
		if(k<0){ k=0; fz=0.0; }
		else if(k>Nz-2){ k=Nz-2; fz=1.0; }
		else{ fz=sz-floor(sz); }
	}

	void save_velocities()
	{
		u.copy_to(du);
		v.copy_to(dv);
		w.copy_to(dw);
	}

	void get_velocity_update()
	{
		du.zero();
		dv.zero();
		dw.zero();
		//du,dv,dw holds the saved velocites
		//u, v, w hols the new velocities
		//Thus the change in velocity in e.g u is: u - du
// 		float *pold,*pnew;
// 		for(pnew = u.data, pold = du.data; pnew < u.data + u.size; ++pnew, ++pold)
// 			*pold = *pnew - *pold;
// 		for(pnew = v.data, pold = dv.data; pnew < v.data + v.size; ++pnew, ++pold)
// 			*pold = *pnew - *pold;
// 		for(pnew = w.data, pold = dw.data; pnew < w.data + w.size; ++pnew, ++pold)
// 			*pold = *pnew - *pold;

		int i;
		for(i=0; i<u.size; ++i)
			du.data[i]=u.data[i]-du.data[i];
		for(i=0; i<v.size; ++i)
			dv.data[i]=v.data[i]-dv.data[i];
		for(i=0; i<w.size; ++i)
			dw.data[i]=w.data[i]-dw.data[i];
	}

	void add_gravity(float dt)
	{
		float gdt = gravity*dt;
		for (float * itr = v.data; itr < v.data + v.size; ++itr)
			*itr -= gdt;
	}

	void classify_voxel()
	{
// 		for(int k = 0; k < Nz; ++k)
// 			for(int j = 0; j < Ny; ++j)
// 				for(int i = 0; i < Nx; ++i)
// 				{
// 					if( i == 0 || i == Nx - 1 || j == 0 || k == 0 || k == Nz - 1)
// 						marker(i,j,k) = SOLIDCELL;
// 					//if(j == 0)
// 				}

		for(int k = 0; k < Nz; ++k)
			for(int j = 0; j < Ny; ++j)
				marker(0,j,k) = marker(Nx-1,j,k) = SOLIDCELL; //Left, right cells

		for(int k = 0; k < Nz; ++k)
			for(int i = 0; i < Nx; ++i)
				marker(i,0,k) = SOLIDCELL; //floor cells

		for(int j = 0; j < Ny; ++j)
			for(int i = 0; i < Nx; ++i)
				marker(i,j,0) = marker(i,j,Nz-1) = SOLIDCELL; //Front back wall cells
	}
	
	void apply_boundary_conditions() //As of know we set zero on the boundary "floor"
	{
		//Project velocities

		for(int k = 0; k<u.nz; ++k)
			for(int j = 0; j < u.ny; ++j)
				u(0,j,k) = u(1,j,k) = u(u.nx-1,j,k) = u(u.nx-2,j,k) = 0; //Left, right wall u component

		for(int k = 0; k<v.nz; ++k)
			for(int i = 0; i < v.nx; ++i)
				v(i,0,k) = v(i,1,k) = 0; //floor v component

		for(int j = 0; j < w.ny; ++j)
			for(int i = 0; i < w.nx; ++i)
				w(i,j,0) = w(i,j,1) = w(i,j,w.nz-1) = w(i,j,w.nz-2) = 0; //Front back wall w component

		
	}
	
	float CFL()
	{
		float maxu = u.infnorm();
		float maxv = v.infnorm();
		float maxw = w.infnorm();
		float maxvel = max(h*gravity,sqr(maxu) + sqr(maxv) + sqr(maxw));
		if(maxvel < 10e-16f)
			maxvel = 10e-16f;
		return h/sqrtf(maxvel);
	}

	void form_poisson(float dt);
	void calc_divergence();
	void project(float dt);
	void solve_pressure(int maxiterations, double tolerance);

};


void Grid::solve_pressure(int maxiterations, double tolerance)
{
	//poisson.write_file_to_matlab("test.txt");

	cg.solve(poisson,rhs,maxiterations,tolerance,x,marker);
}

//----------------------------------------------------------------------------//
// Subtracts the pressure gradient from the velocities making the velocity field 
// divergence free
//----------------------------------------------------------------------------//
void Grid::project(float dt)
{
	float scale = overh * dt / rho;
	int offset;
	float val;
	for(int k = 0; k < Nz; ++k)
		for(int j = 0; j < Ny; ++j)
			for(int i = 0; i < Nx; ++i)
			{
				if(marker(i,j,k) == FLUIDCELL)
				{
					offset = i + Nx*(j + Ny*k); //Offset into rhs.data array
					val = scale * float(x.data[offset]);

					u(i,j,k) -= val;
					u(i+1,j,k) += val;

					v(i,j,k) -= val;
					v(i,j+1,k) += val;

					w(i,j,k) -= val;
					w(i,j,k+1) += val;

				}
				else if(marker(i,j,k) == SOLIDCELL)
				{
					u(i,j,k) = 0;
					u(i+1,j,k) = 0;

					v(i,j,k) = 0;
					v(i,j+1,k) = 0;

					w(i,j,k) = 0;
					w(i,j,k+1) = 0;
				}
			}
}


//----------------------------------------------------------------------------//
// Calculates the divergence in the velocity field
// and fills the the b of the poisson equation
//----------------------------------------------------------------------------//
void Grid::calc_divergence()
{
	rhs.zero();
	int offset;
	double scale = overh;
	for(int k = 0; k < Nx; ++k)
		for(int j = 0; j < Ny; ++j)
			for(int i = 0; i < Nz; ++i)
			{
				if(marker(i,j,k) == FLUIDCELL)
				{
					offset = i + Nx*(j + Ny*k); //Offset into rhs.data array

					if(marker(i+1,j,k) == SOLIDCELL)		//If cell(i+1,j,k) remove u(i+1/2,j,k)
						rhs.data[offset] -= u(i,j,k);
					else if(marker(i-1,j,k) == SOLIDCELL)	//If cell(i-1,j,k) remove u(i-1/2,j,k)
						rhs.data[offset] += u(i+1,j,k);
					else
						rhs.data[offset] += u(i+1,j,k) - u(i,j,k);

					if(marker(i,j+1,k) == SOLIDCELL)		//If cell(i,j+1,k) remove u(i,j+1/2,k)
						rhs.data[offset] -= v(i,j,k);
					else if(marker(i,j-1,k) == SOLIDCELL)		//If cell(i,j-1,k) remove u(i,j-1/2,k)
						rhs.data[offset] += v(i,j+1,k);
					else
						rhs.data[offset] += v(i,j+1,k) - v(i,j,k);

					if(marker(i,j,k+1) == SOLIDCELL)		//If cell(i,j,k+1) remove u(i,j,k+1/2)
						rhs.data[offset] -= w(i,j,k);
					else if(marker(i,j,k-1) == SOLIDCELL)	//If cell(i,j,k-1) remove u(i,j,k-1/2)
						rhs.data[offset] += w(i,j,k+1);
					else
						rhs.data[offset] += w(i,j,k+1) - w(i,j,k);

					rhs.data[offset] *= -scale;
				}
			}
}

//----------------------------------------------------------------------------//
// Sets up the coefficients in the matrix of the poisson equation
//----------------------------------------------------------------------------//
void Grid::form_poisson(float dt)
{
	poisson.zero();
	double scale = overh * overh * dt / rho; // dt / (rho * dx^2) = (1/dx^2) * dt / rho
	for(int k = 1; k < Nz-1; ++k)
		for(int j = 1; j < Ny-1; ++j)
			for(int i = 1; i < Nx-1; ++i)
			{
				if(marker(i,j,k) == FLUIDCELL)
				{
					if(marker(i-1,j,k) != SOLIDCELL)		//Cell(i-1,j,k) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(marker(i+1,j,k) != SOLIDCELL)		//Cell(i+1,j,k) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(marker(i+1,j,k) == FLUIDCELL)	//Cell(i+1,j,k) Is fluid
							poisson(i,j,k,1) -= scale; 
					}

					if(marker(i,j-1,k) != SOLIDCELL)		//Cell(i,j-1,k) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(marker(i,j+1,k) != SOLIDCELL)		//Cell(i,j+1,k) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(marker(i,j+1,k) == FLUIDCELL)	//Cell(i,j+1,k) Is fluid
							poisson(i,j,k,2) -= scale; 
					}

					if(marker(i,j,k-1) != SOLIDCELL)		//Cell(i,j,k-1) Is air or fluid
						poisson(i,j,k,0) += scale; 
					if(marker(i,j,k+1) != SOLIDCELL)		//Cell(i,j,k+1) Is air or fluid
					{
						poisson(i,j,k,0) += scale; 
						if(marker(i,j,k+1) == FLUIDCELL)	//Cell(i,j,k+1) Is fluid
							poisson(i,j,k,3) -= scale; 
					}					
				} //End if CELL(i,j,k) == FLUIDCELL
			}
}





#endif