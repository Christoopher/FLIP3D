#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "Vector3.h"
#include "Grid.h"
#include "Array3D.h"
#include <omp.h>

//#include "OpenGLViewer.h"

struct Particles
{
	//maximum nr of particles and the number of particles in use
	int maxnp, currnp; 
	//vec3f *vel, *pos;

	std::vector<vec3f> vel, pos;
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
		//pos = new vec3f[maxnp];
		//vel = new vec3f[maxnp];
		//pos.resize(maxnp);
		//vel.resize(maxnp);
		weightsumx.init(grid.u.nx, grid.u.ny, grid.u.nz);
		weightsumy.init(grid.v.nx, grid.v.ny, grid.v.nz);
		weightsumz.init(grid.w.nx, grid.w.ny, grid.w.nz);
	}

	~Particles() 
	{
		//delete [] vel;
		//delete [] pos;
	}

	void clear()
	{
		//std::memset(vel, 0, maxnp*sizeof(vec3f));
		//std::memset(pos, 0, maxnp*sizeof(vec3f));
		pos.clear();
		vel.clear();
		currnp = 0;
	}

	void
	remove(int i )
	{
		//vel.erase(vel.begin()+i);
		//pos.erase(pos.begin()+i);
		std::swap(vel[i],vel.back());
		std::swap(pos[i],pos.back());
		vel.pop_back();
		pos.pop_back();
		currnp = vel.size();
		//std::cout << "REMOVED PARTICLE " << i << "\n";
	}

};

void move_particles_in_grid(Particles & particles, Grid & grid, float dt)
{
	vec3f vel;
	int ui, i, vj, j, wk, k;
	float ufx, fx, vfy, fy, wfz, fz;

	float xmax = (grid.Nx-1.001)*grid.h, xmin = 1.001*grid.h;
	float ymax = (grid.Ny-1.001)*grid.h, ymin = 1.001*grid.h;
	float zmax = (grid.Nz-1.001)*grid.h, zmin = 1.001*grid.h;

#if NDEBUG
	std::vector< std::vector<int> > removeIndices;
	removeIndices.resize(omp_get_max_threads());
#endif

#pragma omp parallel for private(ufx, fx, vfy, fy, wfz, fz,ui, i, vj, j, wk, k,vel) shared(xmax,ymax,zmax)
	for(int p = 0; p < particles.pos.size(); p++)
	{
		//Trilerp from grid
		grid.bary_x(particles.pos[p][0], ui, ufx);
		grid.bary_x_centre(particles.pos[p][0], i, fx);

		grid.bary_y(particles.pos[p][1], vj, vfy);
		grid.bary_y_centre(particles.pos[p][1], j, fy);

		grid.bary_z(particles.pos[p][2], wk, wfz);
		grid.bary_z_centre(particles.pos[p][2], k, fz);

		vel = vec3f(grid.u.trilerp(ui,j,k,ufx,fy,fz), grid.v.trilerp(i,vj,k,fx,vfy,fz), grid.w.trilerp(i,j,wk,fx,fy,wfz));			


		//Move particle one step with forward euler
		if(grid.marker(ui,vj,wk) == SOLIDCELL)
		{
			//std::cout << "particle in solid voxel: [" << ui << "," << vj << "," << wk << "]\n";
			continue;
		}

		vec3f newpos = particles.pos[p] + dt*vel;

		//Clamp to border
		//clamp(newpos[0], xmin,xmax);
		//clamp(newpos[1], ymin,ymax);
		//clamp(newpos[2], zmin,zmax);

		grid.bary_x(newpos[0], ui, ufx);		
		grid.bary_y(newpos[1], vj, vfy);		
		grid.bary_z(newpos[2], wk, wfz);

		if(ui < 0 || vj < 0 || wk < 0)
			continue;

		//Push particle out
		float scale = 1.0;
		bool moveX = true, moveY = true, moveZ = true;
		vec3f movevec(0.0);
		if(grid.marker(ui,vj,wk) == SOLIDCELL)
		{
			// X - AXIS
			if(grid.marker(ui+1,vj,wk) != SOLIDCELL) //Push left
			{
				//newpos[0] += grid.h*scale;
				movevec[0] += 1.0;
			}
			else if(grid.marker(ui-1,vj,wk) != SOLIDCELL) //Push right
			{
				//newpos[0] -= grid.h*scale;
				movevec[0] -= 1.0;
			}
			else
			{				
				moveX = false;
			}

			// Y - AXIS 
			
			if(grid.marker(ui,vj+1,wk) != SOLIDCELL) //Push up
			{
				//newpos[1] += grid.h*scale;
				movevec[1] += 1.0;
			}
			else if(grid.marker(ui,vj-1,wk) != SOLIDCELL) //Push down
			{
				movevec[1] -= 1.0;
				//newpos[1] -= grid.h*scale;
			}
			else
			{
				moveY = false;
			}

			// Z - AXIS 
			if(grid.marker(ui,vj,wk+1) != SOLIDCELL) //Push backwards
			{
				//newpos[2] += grid.h*scale;
				movevec[2] += 1.0;
			}
			else if(grid.marker(ui,vj,wk-1) != SOLIDCELL) //Push forward
			{
				//newpos[2] -= grid.h*scale;
				movevec[2] += 1.0;
			}
			else
			{
				//std::cout << "else in Z_AXIS\n";
				//newpos = particles.pos[p];
				moveZ = false;
			}
				
		}

		//if surrounded by solid
		if(!moveZ && !moveY && !moveX)
		{
			//std::cout << "could not move particle\n";
			newpos = particles.pos[p];
		}
		else
		{
			newpos += movevec*grid.h*scale;
		}

	
		particles.pos[p] = newpos;
	}


}

//void write_to_file(const char * filename);

void update_from_grid(Particles & particles, Grid & grid)
{
	int i, ui, j, vj, k, wk;
	float fx, ufx, fy, vfy, fz, wfz;
#pragma omp parallel for private(fx, ufx, fy, vfy, fz, wfz, i, ui, j, vj, k, wk)
	for (int p = 0; p < particles.pos.size(); ++p) //Loop over all particles
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
		float alpha = 0.05f;
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
	int tmpi, tmpj, tmpk;

	particles.weightsumx.zero();
	particles.weightsumy.zero();
	particles.weightsumz.zero();

#if NDEBUG
	std::vector< std::vector<int> > removeIndices;
	removeIndices.resize(omp_get_max_threads());
#else
	std::vector< int > removeIndices;
#endif

#pragma omp parallel for private(ui, vj, wk,i,j,k,fx, ufx, fy, vfy, fz, wfz,tmpi, tmpj, tmpk)
	for (int p = 0; p < particles.pos.size(); ++p) //Loop over all particles
	{
		grid.bary_x(particles.pos[p][0],ui,ufx); tmpi = ui;
		grid.bary_y(particles.pos[p][1],vj,vfy); tmpj = vj;
		grid.bary_z(particles.pos[p][2],wk,wfz); tmpk = wk;

		if(grid.marker(ui,vj,wk) == SOLIDCELL)
		{
			//std::cout << "solid particle in transf2grid.. removign...\n";
#if NDEBUG
			removeIndices[omp_get_thread_num()].push_back(p);
#else
			removeIndices.push_back(p);
#endif
			continue;
		}
		else
			grid.marker(ui,vj,wk) = FLUIDCELL;
		

		grid.bary_y_centre(particles.pos[p][1],j,fy);
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.u,particles.weightsumx,particles.vel[p][0],ui,j,k,ufx,fy,fz);


		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_z_centre(particles.pos[p][2],k,fz);
		accumulate(grid.v,particles.weightsumy,particles.vel[p][1],i,vj,k,fx,vfy,fz);


		grid.bary_x_centre(particles.pos[p][0],i,fx);
		grid.bary_y_centre(particles.pos[p][1],j,fy);
		accumulate(grid.w,particles.weightsumz,particles.vel[p][2],i,j,wk,fx,fy,wfz);

		

		//mp4Vector val(tmpi,tmpj,tmpk,-1.0);
		//levelset[(i + grid.Nx*(j + grid.Ny*k))] = val;

	}



#if NDEBUG
	for(int i = 0; i < removeIndices.size(); ++i)
	{
		for(int j = 0; j < removeIndices[i].size(); ++j)
		{
			particles.remove(removeIndices[i].at(j));
		}
	}
#else
	for(int j = 0; j < removeIndices.size(); ++j)
	{
		particles.remove(removeIndices[j]);
	}
#endif



	//Scale u velocities with weightsumx

	#pragma omp parallel for
	for (int i = 0; i < grid.u.size; i++)
	{
		if(grid.u.data[i] != 0)
			grid.u.data[i] /= particles.weightsumx.data[i];
	}

	//Scale v velocities with weightsumy
	#pragma omp parallel for
	for (int i = 0; i < grid.v.size; i++)
	{
		if(grid.v.data[i] != 0)
			grid.v.data[i] /= particles.weightsumy.data[i];
	}

	//Scale w velocities with weightsumz
	#pragma omp parallel for
	for (int i = 0; i < grid.w.size; i++)
	{
		if(grid.w.data[i] != 0)
			grid.w.data[i] /= particles.weightsumz.data[i];
	}

}

//----------------------------------------------------------------------------//
// Adds a particle to the particles struct
//----------------------------------------------------------------------------//
inline void add_particle(Particles & particles, vec3f & pos, vec3f  & vel)
{
	//particles.pos[particles.currnp] = pos;
	//particles.vel[particles.currnp] = vel;
	particles.pos.push_back(pos);
	particles.vel.push_back(vel);
	++particles.currnp;
}

void write_paricle_pos_binary(Particles & particles)
{
	static int frame = 0;
	std::stringstream ss;
	ss << "frame_" << frame << ".dat";
	std::string str = ss.str();
	FILE *file = fopen(str.c_str(),"wb");

	int count = particles.pos.size();
	fwrite((const void *) (&count),sizeof(int),1,file);		
	fwrite((const void *) (&particles.pos[0]),sizeof(vec3f),count,file);

	fclose(file);

	frame++;
}

void read_paricle_pos_binary(Particles & particles, int frame)
{
	int count;
	std::string str;
	std::stringstream ss;
	ss << "frame_" << frame << ".dat";
	ss >> str;
	FILE *file = fopen(str.c_str(),"rb");
	if(file == NULL)
	{
		std::cout << "Hej";
		system("pause");
	}
		
	fread((void *) (&count),sizeof(int),1,file);
	
	particles.maxnp = count;
	particles.currnp = count; 
	particles.pos.clear();
	particles.pos.resize(count);
	fread((void *) (&particles.pos[0]),sizeof(vec3f),particles.pos.size(),file);

	fclose(file);
}

#endif