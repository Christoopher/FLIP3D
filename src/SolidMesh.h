#ifndef SOLID_MESH_H
#define SOLID_MESH_H

#include "ObjLoader.h"
#include <vector>
#include <string>
#include "Vector3.h"

struct SolidMesh
{
	SolidMesh() {}
	SolidMesh(std::string filename, vec3f pos, float h) : position(pos) , velocity(0.0)
	{
		mesh = LoadObjMesh(filename);		
	}

	void init(std::string filename, vec3f pos, vec3f vel, float h)
	{
		mesh = LoadObjMesh(filename);
		scale = h;
		position = pos;
		velocity = vel;
		seed_random_points(h);
	}

	void mesh_to_grid(Grid & grid)
	{
		std::vector<ObjMeshFace>::iterator itr;
			
		int i,j,k;
		for (itr = mesh.faces.begin(); itr != mesh.faces.end(); itr++)
		{
			vec3f pos = itr->vertices->pos;
			pos += position;

					

			i = floor(pos[0]);
			j = floor(pos[1]);
			k = floor(pos[2]);

			//check x
			if(i < 0 || i > grid.Nx-1 )
				continue;
			//check y
			//check x
			if(j < 0 || j > grid.Ny-1 )
				continue;
			//check z
			//check x
			if(k < 0 || k > grid.Nz-1 )
				continue;

			grid.marker(i,j,k) = SOLIDCELL;
			grid.u(i,j,k) = grid.u(i+1,j,k) = velocity[0];
			grid.v(i,j,k) = grid.v(i,j+1,k) = velocity[1];
			grid.w(i,j,k) = grid.w(i,j,k+1) = velocity[2];
		}

		int shrinks = 10;
		float shrinkage = 1.0;
		for(int itr = 0; itr < shrinks; ++itr)
		{
			std::vector<vec3f>::iterator itr2;
			for (itr2 = randomPoints.begin(); itr2 != randomPoints.end(); itr2++)
			{
				vec3f pos(*itr2*shrinkage);
				pos += position;
	
				i = floor(pos[0]);
				j = floor(pos[1]);
				k = floor(pos[2]);

				//check x
				if(i < 0 || i > grid.Nx )
					continue;
				//check y
				//check x
				if(j < 0 || j > grid.Ny )
					continue;
				//check z
				//check x
				if(k < 0 || k > grid.Nz )
					continue;

		
				grid.marker(i,j,k) = SOLIDCELL;
				grid.u(i,j,k) = grid.u(i+1,j,k) = velocity[0];
				grid.v(i,j,k) = grid.v(i,j+1,k) = velocity[1];
				grid.w(i,j,k) = grid.w(i,j,k+1) = velocity[2];
			}
			shrinkage *= (1.0 - scale);
		}
	}

	void seed_random_points(float h)
	{
		std::vector<ObjMeshFace>::iterator itr;

		int i,j,k;
		for (itr = mesh.faces.begin(); itr != mesh.faces.end(); itr++)
		{
			vec3f p1 = itr->vertices[0].pos;
			vec3f p2 = itr->vertices[1].pos;
			vec3f p3 = itr->vertices[2].pos;
			
			float d1 = dist(p1,p2);
			float d2 = dist(p1,p3);
			float d3 = dist(p2,p3);
			
			//If the dist between all pairs of vertives are less than gridsize than 
			//we dont seed points
			//if(d1 < h || d2 < h || d3 < h)
			//	continue;

			vec3f v1 = p2-p1;
			vec3f v2 = p3-p1;
			float area = mag(cross(v1,v2))*0.5;

			int points = 100;//2.0 * area / h;

			double a,b,c;
			float step = 1.0/sqrtf(points);
			for (float aa = 0; aa < 1.0; aa+=step)
			{
				for (float bb = 0; bb < 1.0; bb+=step)
				{
					if(aa+bb > 1.0)
					{
						a = 1-aa;
						b = 1-bb;
					}
					else
					{
						a = aa;
						b = bb;
					}
					c = 1 - a - b;

					vec3f p = a*p1 + b*p2 + c*p3;

					randomPoints.push_back(p);
				}
			}

		}
	}

	void draw()
	{
		glBegin(GL_POINTS);
		std::vector<ObjMeshFace>::iterator itr;

		int i,j,k;
		for (itr = mesh.faces.begin(); itr != mesh.faces.end(); itr++)
		{
			vec3f p1 = ((*itr).vertices[0].pos + position)*scale;
			vec3f p2 = ((*itr).vertices[1].pos + position)*scale;
			vec3f p3 = ((*itr).vertices[2].pos + position)*scale;
			glVertex3f(p1[0],p1[1],p1[2]);
			glVertex3f(p2[0],p2[1],p2[2]);
			glVertex3f(p3[0],p3[1],p3[2]);
			//glVertex3f((*itr).vertices[0].pos[0], (*itr).vertices[0].pos[1], itr->vertices[0].pos[2]);
			//glVertex3f((*itr).vertices[1].pos[0], (*itr).vertices[1].pos[1], itr->vertices[1].pos[2]);
			//glVertex3f((*itr).vertices[2].pos[0], (*itr).vertices[2].pos[1], itr->vertices[2].pos[2]);
		}

		std::vector<vec3f>::iterator itr2;
		for (itr2 = randomPoints.begin(); itr2 != randomPoints.end(); itr2++)
		{
			vec3f p1 = (*itr2 + position)*scale;
			glVertex3f(p1[0],p1[1],p1[2]);
		}
		glEnd();
		glFlush();
	}

	void
	move(const float &dt)
	{
		position += dt*velocity;
	}

	ObjMesh mesh;
	vec3f position;
	vec3f velocity;
	float scale;
	std::vector<vec3f> randomPoints;

	
};

#endif