#ifndef _GRID_H_
#define _GRID_H_

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

#include "Array3D.h"

struct Grid
{
	int Nx, Ny,Nz;
	float h,overh, gravity;

	Array3f u,v,w,du,dv,dw; //Staggered u,v,w velocities
	Array3d p; //Pressure
	Array3c marker; //Voxel classification

	Grid() {}
	
	Grid(int Nx_, int Ny_, int Nz_, float h_, float gravity_) : Nx(Nx_), Ny(Ny_), Nz(Nz_), h(h_), overh(1.0f/h), gravity(gravity_)
	{
		init(Nx_,Ny_,Nz_,h_,gravity_);
	}

	void init(int Nx_, int Ny_, int Nz_, float h_, float gravity_)
	{
		Nx = Nx_; Ny = Ny_; Nz = Nz_; h = h_; gravity = gravity_;
		overh = 1.0f/h;
		u.init(Nx_+1,Ny_,Nz_);
		v.init(Nx_,Ny_+1,Nz_);
		w.init(Nx_,Ny_,Nz_+1);
		du.init(Nx_+1,Ny_,Nz_);
		dv.init(Nx_,Ny_+1,Nz_);
		dw.init(Nx_,Ny_,Nz_+1);
		marker.init(Nx_,Ny_,Nz_);
		p.init(Nx_,Ny_,Nz_);
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
		else if(i>p.nx-2) { 	//Why?
			i=p.nx-2; fx=1.0; 
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
		else if(j>p.ny-2){ j=p.ny-2; fy=1.0; }
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
		else if(k>p.nz-2){ k=p.nz-2; fz=1.0; }
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
		
		//du,dv,dw holds the saved velocites
		//u, v, w hols the new velocities
		//Thus the change in velocity in e.g u is: 
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
		{
			*itr -= gdt;
		}
	}

	void apply_boundary_conditions() //As of know we set zero on the boundary "floor"
	{
		for(int k=0; k<v.nz; ++k)
			for(int i=0; i < v.nx; ++i)
				{
					v(i,0,k) = 0;
				}
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

};





#endif