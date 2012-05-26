
#ifndef ANI_MESHER_H
#define ANI_MESHER_H

#include "MarchingCubes/MarchingCubes.h"
#include "Particles.h"

#include <cmath>
#include <iostream>
#include <armadillo>
#include "Array3D.h"

#include <time.h>
#include "hr_time.h"

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

typedef arma::mat armaMat;
typedef arma::vec armaVec;


bool calcmesh = true;
armaMat * Gs;
vec3f * smoothPos;
float * density;
mp4Vector * phi;

float isovalue = 34;
float ks = 1400;
float kr = 4.0;	
float H = 0.05;
float kspecial = 0.9;
float amp = 1.0;

int currentParticles;
int ldx, ldy, ldz;


void 
weightedPos(Particles & p, int i,float r, vec3f & xiW)
{
	float weightsum = 0.0;
	float weighti;
	float dr;
	vec3f sumweightpos(0.0);

	//calc the weighted sum
	for (int j = 0; j < p.currnp; ++j) //alla andra partiklar
	{
		dr = dist(p.pos[i], p.pos[j])/r;
		weighti = max(1.0 - dr,0.0);            ///CHANGED THIS FROM max(1.0 - dr*dr*dr,0.0); 
		weightsum += weighti;
		sumweightpos += weighti*p.pos[j];
	}
	xiW = sumweightpos/weightsum;
}

void 
smoothParticles(const Particles & p,float *density, vec3f * smoothpos, float lam,float h)
{
	float weightsum;
	float weighti;
	float dr;
	float r = 2*h;
	vec3f weightpos;

	for (int i = 0; i < p.currnp; ++i) //i �r den partiklen vi kollar p�
	{
		//Reset the sum and pos
		weightsum = 0;
		weightpos[0] = weightpos[1] = weightpos[2] = 0.0;
		//calc the weighted sum
		for (int j = 0; j < p.currnp; ++j) //alla andra partiklar
		{
			dr = dist(p.pos[i], p.pos[j])/r;
			weighti = max(1.0 - dr*dr*dr,0.0);
			weightsum += weighti;
			weightpos += weighti*p.pos[j];
		}
		smoothpos[i] = (1.0 - lam)*p.pos[i] + lam*weightpos/weightsum;
		density[i] = weightsum;
	}
}

inline void 
vecvecMult(vec3f v,float w, armaMat & res)
{
	res(0,0) = w*v[0]*v[0];
	res(0,1) = w*v[0]*v[1];
	res(0,2) = w*v[0]*v[2]; 

	res(1,0) = w*v[1]*v[0];
	res(1,1) = w*v[1]*v[1]; 
	res(1,2) = w*v[1]*v[2];

	res(2,0) = w*v[2]*v[0];
	res(2,1) = w*v[2]*v[1];
	res(2,2) = w*v[2]*v[2];
}

void 
calcCovMtx(Particles & p, int i,float r, armaMat & Ci, vec3f & xiW)
{
	
	float weightsum = 0.0;
	float weighti = 0.0;
	float dr = 0.0;
	armaMat matSum; matSum.zeros(3,3);
	armaMat tmpMat; tmpMat.zeros(3,3);
	
	for (int j = 0; j < p.currnp; ++j) //alla andra partiklar
	{
		dr = dist(p.pos[i], p.pos[j])/r;
		weighti = max(1.0 - dr*dr*dr,0.0);

		if(weighti <= 10e-16)
			continue;

		weightsum += weighti;
		vecvecMult(p.pos[j]-xiW,weighti,tmpMat);
		matSum += tmpMat;
	}
	Ci = matSum / weightsum;
	
	/*for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			Ci(i,j) = matSum(i,j) / weightsum;*/
					
}

void 
calcGi(const armaMat & Ci, armaMat & Gi)
{
	armaVec eigenv; 
	armaMat R, RT;
	
	

	svd(R, eigenv, RT, Ci);
	
	//float ks = 1.0/arma::norm(Ci,"inf");//1.0/eigenv(0); // gives norm(ks*Ci) approx 1
	//Check if distribution is uniform
	//if(eigenv(2) <= eigenv(0)*kspecial) 
	//{		
		//Check if to much stretch
		if(eigenv(0) >= kr*eigenv(2))
		{
			/*float alpha = ( kr*eigenv(2)) / (eigenv(0) + kr*eigenv(2));
			float gamma = eigenv(1) / eigenv(0);
			eigenv(1) *= gamma;
			eigenv(0) *= alpha;
			eigenv(2) *= (1.0f-alpha);*/
			eigenv(0) = max(eigenv(0),eigenv(0)/kr);
			eigenv(1) = max(eigenv(1),eigenv(0)/kr);
			eigenv(2) = max(eigenv(2),eigenv(0)/kr);
			//ks = 1.0/eigenv(0);
		}

		//eigenv(0) = eigenv(0);
		//eigenv(1) = eigenv(1);
		//eigenv(2) = eigenv(2);
		
	//}
	//else
	//	eigenv(0) = eigenv(1) = eigenv(2) = kn;

	//inverse of eigenmatrix
	eigenv(0) = 1.0/(eigenv(0)*ks);
	eigenv(1) = 1.0/(eigenv(1)*ks);
	eigenv(2) = 1.0/(eigenv(2)*ks);

	Gi = (1.0/(2.0*H))*R*diagmat(eigenv)*trans(RT);
	

	//matrixMultiplyToGetGi(eigenv, R, RT, Gi);
}
	
/*
 * This method constructs G_i matrices for each particle i
 * By first calculating weighted covariance matrices C_i
*/
void 
calcAniMatrices(Particles & p, armaMat * Gs)
{
	float r = 2*H;
#pragma omp parallel for
	for (int i = 0; i < p.currnp; ++i) //i �r den partikeln vi kollar p�
	{		
		vec3f xiW;
		armaMat Ci; 
		Ci.zeros(3,3);
		//First calc xiW
		weightedPos(p,i,r,xiW);
		//Calc Ci
		calcCovMtx(p,i,r,Ci,xiW);
		
		//Do SVD and check eigenvalues
		Gs[i].zeros(3,3); //Must be allocated
		
		calcGi(Ci,Gs[i]);
	}
}

void 
getPhi(Particles & p, armaMat * G,float *density ,vec3f * smoothPos, mp4Vector & point)
{
	point.val = 0.0;
	float sum = 0.0;
//#pragma omp parallel for reduction( +: sum)
	for(int i = 0; i < p.currnp; ++i)
	{
		float Gnorm, Grnorm;
		vec3f Gr, r;
		Gnorm = -100.0;
		r[0] = point.x - smoothPos[i][0]; 
		r[1] = point.y - smoothPos[i][1];
		r[2] = point.z - smoothPos[i][2];

		Gr[0] = G[i](0,0)*r[0] + G[i](0,1)*r[1] + G[i](0,2)*r[2];
		Gr[1] = G[i](1,0)*r[0] + G[i](1,1)*r[1] + G[i](1,2)*r[2];
		Gr[2] = G[i](2,0)*r[0] + G[i](2,1)*r[1] + G[i](2,2)*r[2];

		Grnorm = mag(Gr);	
		//std::cout << "rnorm: " << mag(r) << "\n";
		//std::cout << "Grnorm: " << Grnorm << "\n\n";

		//Get norm (i.e. max) of G
		Gnorm = norm(G[i],"inf");

		float gg = Grnorm*Grnorm;
		float P = max((1.0-gg)*(1.0-gg)*(1.0-gg),0.0);
		sum += (1.0/density[i])*Gnorm * P * amp;
	}

	point.val = sum;
}
void
createPhi(Particles & p, vec3f * smoothpos,float * density, armaMat * Gs, const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{

	Array3c marker;
	marker.init(Nx,Ny,Nz);

	std::cout << "Build marker grid\n";
	#pragma omp parallel for
	for(int n = 0; n < p.currnp; ++n)
	{
		int i = floor(p.pos[n][0]/h);
		int j = floor(p.pos[n][1]/h);
		int k = floor(p.pos[n][2]/h);

		marker(i,j,k) = FLUIDCELL;
	}

		CStopWatch stopwatch;
	stopwatch.startTimer();

	std::cout << "Buildlevelset\n";
#pragma omp parallel for
	for(int k = 1; k < Nz-1; ++k)
	{
		for(int j = 1; j < Ny-1; ++j)
		{
			for(int i = 1; i < Nx-1; ++i)
			{
				
				/*if(	marker(i,j,k) == FLUIDCELL || 
					marker(i-1,j,k) == FLUIDCELL || marker(i+1,j,k) == FLUIDCELL || 
					marker(i,j-1,k) == FLUIDCELL || marker(i,j+1,k) == FLUIDCELL || 
					marker(i,j,k-1) == FLUIDCELL || marker(i,j,k+1) == FLUIDCELL ||
					marker(i+1,j+1,k+1) == FLUIDCELL || marker(i-1,j-1,k-1) == FLUIDCELL
				)		
				{*/
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = (i)*h;
					phi[phioffset].y = (j)*h;
					phi[phioffset].z = (k)*h;
					getPhi(p,Gs,density,smoothpos,phi[phioffset]);
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				/*}
				else
				{
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = (i)*h;
					phi[phioffset].y = (j)*h;
					phi[phioffset].z = (k)*h;
					phi[phioffset].val = -10000;
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				}
				*/
			}	
		}	
	}
	stopwatch.stopTimer();
	std::cout << std::scientific;
	std::cout << stopwatch.getElapsedTime() << "\n";
}

void 
getSimplePhi2(Particles & p,armaMat * G, mp4Vector & point)
{
	vec3f x, point_;
	point_[0] = point.x; point_[1] = point.y; point_[2] = point.z; 

	float s, ss, wi, R = 1.0, weightsum = 0.0f;
	point.val = 0.0;

	for(int i = 0; i < p.currnp; ++i)
	{
		float Gnorm, Grnorm;
		vec3f Gr, r;
		r = point_ - p.pos[i];

		Gr[0] = G[i](0,0)*r[0] + G[i](0,1)*r[1] + G[i](0,2)*r[2];
		Gr[1] = G[i](1,0)*r[0] + G[i](1,1)*r[1] + G[i](1,2)*r[2];
		Gr[2] = G[i](2,0)*r[0] + G[i](2,1)*r[1] + G[i](2,2)*r[2];
		
		Gnorm = norm(G[i],"inf");
		Grnorm = mag(Gr);
		std::cout << "r: " << mag(r) << "\n";
		std::cout << "GrNorm: " << Grnorm << "\n";
		s = Grnorm*Grnorm/R;
		ss = s*s;
		wi = max(0.0, (1.0 - ss)*(1.0 - ss)*(1.0 - ss));
		x += point_*wi; //Averaged point
		weightsum += wi;
	}
	float r = 0.1;
	if(weightsum < 10e-7)
	{
		point.val = 1000.0;
	}
	else
	{
		x /= weightsum;
		point.val = mag(point_ - x) - r;
	}
	std::cout << "isovalue: " << point.val << "\n";
}
void 
createSimplePhi2(Particles & p, armaMat * G,const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{
	Array3c marker;
	marker.init(Nx,Ny,Nz);
	vec3f phiPos; 
	std::cout << "Build marker grid\n";
#pragma omp parallel for
	for(int n = 0; n < p.currnp; ++n)
	{
		int i = floor(p.pos[n][0]/h);
		int j = floor(p.pos[n][1]/h);
		int k = floor(p.pos[n][2]/h);

		marker(i,j,k) = FLUIDCELL;
	}

	std::cout << "Buildlevelset\n";
#pragma omp parallel for
	for(int k = 1; k < Nz-1; ++k)
	{
		for(int j = 1; j < Ny-1; ++j)
		{
			for(int i = 1; i < Nx-1; ++i)
			{
				if(	marker(i,j,k) == FLUIDCELL || 
					marker(i-1,j,k) == FLUIDCELL || marker(i+1,j,k) == FLUIDCELL || 
					marker(i,j-1,k) == FLUIDCELL || marker(i,j+1,k) == FLUIDCELL || 
					marker(i,j,k-1) == FLUIDCELL || marker(i,j,k+1) == FLUIDCELL ||
					marker(i+1,j+1,k+1) == FLUIDCELL || marker(i-1,j-1,k-1) == FLUIDCELL)		
				{
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = phiPos[0] = (i+0.5)*h;
					phi[phioffset].y = phiPos[1] = (j+0.5)*h;
					phi[phioffset].z = phiPos[2] = (k+0.5)*h;
					getSimplePhi2(p,Gs,phi[phioffset]);
				}
				else
				{
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = (i+0.5)*h;
					phi[phioffset].y = (j+0.5)*h;
					phi[phioffset].z = (k+0.5)*h;
					phi[phioffset].val = 0.0;
				}
			}
		}
	}
}


void 
getSimplePhi(Particles & p, mp4Vector & point)
{
	vec3f x, point_;
	point_[0] = point.x; point_[1] = point.y; point_[2] = point.z; 

	float magni, s, ss, wi, r = 0.0, wi_den, R = 0.5, ri = 0.1;
	point.val = 0.0;

	for(int i = 0; i < p.currnp; ++i)
	{
		s = mag(point_ - p.pos[i])/R;
		ss = s*s;
		wi += max(0.0, (1.0 - ss)*(1.0 - ss)*(1.0 - ss));
	}

	r = wi*ri;
	x = wi*point_;
	point.val = mag(point_ - x) - r;
}
void 
createSimplePhi(Particles & p, const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{
	Array3c marker;
	marker.init(Nx,Ny,Nz);
	vec3f phiPos; 
	std::cout << "Build marker grid\n";
#pragma omp parallel for
	for(int n = 0; n < p.currnp; ++n)
	{
		int i = floor(p.pos[n][0]/h);
		int j = floor(p.pos[n][1]/h);
		int k = floor(p.pos[n][2]/h);

		marker(i,j,k) = FLUIDCELL;
	}

	std::cout << "Buildlevelset\n";
#pragma omp parallel for
	for(int k = 1; k < Nz-1; ++k)
	{
		for(int j = 1; j < Ny-1; ++j)
		{
			for(int i = 1; i < Nx-1; ++i)
			{
				if(	marker(i,j,k) == FLUIDCELL || 
					marker(i-1,j,k) == FLUIDCELL || marker(i+1,j,k) == FLUIDCELL || 
					marker(i,j-1,k) == FLUIDCELL || marker(i,j+1,k) == FLUIDCELL || 
					marker(i,j,k-1) == FLUIDCELL || marker(i,j,k+1) == FLUIDCELL ||
					marker(i+1,j+1,k+1) == FLUIDCELL || marker(i-1,j-1,k-1) == FLUIDCELL)		
				{
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = phiPos[0] = (i+0.5)*h;
					phi[phioffset].y = phiPos[1] = (j+0.5)*h;
					phi[phioffset].z = phiPos[2] = (k+0.5)*h;
					getSimplePhi(p,phi[phioffset]);
				}
				else
				{
					int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
					phi[phioffset].x = (i+0.5)*h;
					phi[phioffset].y = (j+0.5)*h;
					phi[phioffset].z = (k+0.5)*h;
					phi[phioffset].val = 0.0;
				}
			}
		}
	}
}


void 
initMeshData(int nParts, int Nx, int Ny, int Nz, int res)
{

	if(ldx != Nx || ldy != Ny || ldz != Nz) //if changed size
	{
		std::cout << "allocating new space for phi\n";
		if(phi != NULL)
		{
			delete [] phi;
			phi = new mp4Vector[(Nx-2)*(Ny-2)*(Nz-2)];
		}
		else //first time
			phi = new mp4Vector[(Nx-2)*(Ny-2)*(Nz-2)];				
	}


	if(currentParticles != nParts)
	{
		std::cout << "allocating new space for smoothpos and Gs\n";
		if(smoothPos != NULL)
		{
			delete [] smoothPos;
			smoothPos = new vec3f[nParts];

			delete density;
			density = new float[nParts];
			
		}
		else
		{
			smoothPos = new vec3f[nParts];
			density = new float[nParts];
		}

		if(Gs != NULL)
		{
			delete [] Gs;
			Gs = new armaMat[nParts];
		}
		else
			Gs = new armaMat[nParts];
	}


	ldx = Nx;
	ldy = Ny;
	ldz = Nz;
	currentParticles = nParts;
}

void 
mesh(Particles & p, const int Nx_, const int Ny_, const int Nz_, const float h_, int res, int & numOfTriangles, TRIANGLE *& tri)
{
	int Nx = res*Nx_, Ny = res*Ny_, Nz = res*Nz_;
	initMeshData(p.currnp, Nx, Ny, Nz, res);
	float h = h_/(float)res;

	if(calcmesh)
	{
		std::cout << "Started meshing\n";			
		
		//std::cout << "Smooth particle positions\n";
		smoothParticles(p,density, smoothPos,0.9, h*0.5);

		
		std::cout << "Calculate Anisostropic matrices G\n";
		calcAniMatrices(p,Gs);
		
		std::cout << "Create the levelset\n";
		createPhi(p,smoothPos,density,Gs,Nx,Ny,Nz,h,phi);
		//createSimplePhi2(p,Gs,Nx,Ny,Nz,h,phi);
	}

	std::cout << "Run marching cubes\n";
	delete [] tri;
	tri = MarchingCubes(Nx-2, Ny-2, Nz-2, isovalue, phi, LinearInterp, numOfTriangles);
	std::cout << "Marching cubes creates: " << numOfTriangles << " triangles\n";

	calcmesh = false;
}





#endif