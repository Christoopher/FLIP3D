
#ifndef ANI_MESHER_H
#define ANI_MESHER_H

#include "MarchingCubes/MarchingCubes.h"
#include "Particles.h"

#include <cmath>
#include <iostream>
#include <armadillo>


typedef arma::mat armaMat;
typedef arma::vec armaVec;

float H = 0.5;


void weightedPos(Particles & p, int i,float r, vec3f & xiW)
{
	float weightsum = 0.0;
	float weighti;
	float dr;
	vec3f sumweightpos(0.0);

	//calc the weighted sum
	for (int j = 0; j < p.currnp; ++j) //alla andra partiklar
	{
		dr = dist(p.pos[i], p.pos[j])/r;
		weighti = max(1.0 - dr*dr*dr,0.0);
		weightsum += weighti;
		sumweightpos += weighti*p.pos[j];
	}
	xiW = sumweightpos/weightsum;
}

void smoothParticles(const Particles & p, vec3f * smoothpos, float lam,float h)
{
	float weightsum;
	float weighti;
	float dr;
	float r = h;
	vec3f weightpos;

	for (int i = 0; i < p.currnp; ++i) //i är den partiklen vi kollar på
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
	}
}

inline void vecvecMult(vec3f v,float w, armaMat & res)
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

void calcCovMtx(Particles & p, int i,float r, armaMat & Ci, vec3f & xiW)
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

/*void matrixMultiplyToGetGi(const armaVec & eigenv, const armaMat & R, const armaMat & RT, armaMat & Gi)
{
	Gi(0,0) = R(0,0)*RT(0,0)*eigenv(0)   +   R(1,0)*RT(0,1)*eigenv(1)   +   R(2,0)*RT(0,2)*eigenv(2); //OK
	Gi(0,1) = R(0,1)*RT(0,0)*eigenv(0)   +   R(1,1)*RT(0,1)*eigenv(1)   +   R(2,1)*RT(0,2)*eigenv(2); //OK
	Gi(0,2) = R(0,2)*RT(0,0)*eigenv(0)   +   R(1,2)*RT(0,1)*eigenv(1)   +   R(2,2)*RT(0,2)*eigenv(2); //OK
	Gi(1,0) = R(0,0)*RT(1,0)*eigenv(0)   +   R(1,0)*RT(1,1)*eigenv(1)   +   R(2,0)*RT(1,2)*eigenv(2); //OK
	Gi(1,1) = R(0,1)*RT(1,0)*eigenv(0)   +   R(1,1)*RT(1,1)*eigenv(1)   +   R(2,1)*RT(1,2)*eigenv(2); //OK
	Gi(1,2) = R(0,2)*RT(1,0)*eigenv(0)   +   R(1,2)*RT(1,1)*eigenv(1)   +   R(2,2)*RT(1,2)*eigenv(2); //OK
	Gi(2,0) = R(0,0)*RT(2,0)*eigenv(0)   +   R(1,0)*RT(2,1)*eigenv(1)   +   R(2,0)*RT(2,2)*eigenv(2); //OK
	Gi(2,1) = R(0,1)*RT(2,0)*eigenv(0)   +   R(1,1)*RT(2,1)*eigenv(1)   +   R(2,1)*RT(2,2)*eigenv(2); //OK
	Gi(2,2) = R(0,2)*RT(2,0)*eigenv(0)   +   R(1,2)*RT(2,1)*eigenv(1)   +   R(2,2)*RT(2,2)*eigenv(2); //OK
}*/

void calcGi(const armaMat & Ci, armaMat & Gi)
{
	armaVec eigenv; 
	armaMat R, RT;
	float kr = 4.0;	
	float kspecial = 0.8;	//(1-kn) -> Hur mycket i procent får egenvärdena skilja sig i axlarna.
	float kn = 1.0;
	
	float Cnorm = norm(Ci,"inf");
	
	float ks = 1.0/Cnorm; // gives norm(ks*Ci) approx 1

	svd(R, eigenv, RT, Ci);
	
	bool adjust = false;
	if(eigenv(2) / eigenv(0) < kn)
	{		
		if(eigenv(0) >= kr*eigenv(2))
		{
			eigenv(0) = max(eigenv(0), eigenv(0)/ks);
			eigenv(1) = max(eigenv(1), eigenv(0)/ks);
			eigenv(2) = max(eigenv(2), eigenv(0)/ks);
		}		

		eigenv(0) = eigenv(0) *ks*H;
		eigenv(1) = eigenv(1) *ks*H;
		eigenv(2) = eigenv(2) *ks*H;
	}
	else
		eigenv(0) = eigenv(1) = eigenv(2) = kn;

	//inverse of eigenmatrix
	eigenv(0) = 1.0/eigenv(0);
	eigenv(1) = 1.0/eigenv(1);
	eigenv(2) = 1.0/eigenv(2);

	Gi = (1/H)*R*diagmat(eigenv)*trans(RT);
	//matrixMultiplyToGetGi(eigenv, R, RT, Gi);
}
	
/*
 * This method constructs G_i matrices for each particle i
 * By first calculating weighted covariance matrices C_i
*/
void calcAniMatrices(Particles & p, armaMat * Gs)
{
	vec3f xiW;
	float weightsum;
	float weighti;
	armaMat Ci; 
	Ci.zeros(3,3);

	float r = 2.0*H;
	for (int i = 0; i < p.currnp; ++i) //i är den partikeln vi kollar på
	{		
		//First calc xiW
		weightedPos(p,i,r,xiW);
		//Calc Ci
		calcCovMtx(p,i,r,Ci,xiW);
		
		//Do SVD and check eigenvalues
		Gs[i].zeros(3,3); //Must be allocated
		
		calcGi(Ci,Gs[i]);

	}
}

inline void getPhi(Particles & p, armaMat * G, vec3f * smoothPos, mp4Vector & point)
{
	float Gnorm, Grnorm;
	vec3f Gr, r;
	point.val = 0.0;
	for(int i = 0; i < p.currnp; ++i)
	{
		Gnorm = -100.0;
		r[0] = point.x - smoothPos[i][0]; 
		r[1] = point.y - smoothPos[i][1];
		r[2] = point.z - smoothPos[i][2];

		Gr[0] = G[i](0,0)*r[0] + G[i](0,1)*r[1] + G[i](0,2)*r[2];
		Gr[1] = G[i](1,0)*r[0] + G[i](1,1)*r[1] + G[i](1,2)*r[2];
		Gr[2] = G[i](2,0)*r[0] + G[i](2,1)*r[1] + G[i](2,2)*r[2];

		Grnorm = mag(Gr);	

		//Get norm (i.e. max) of G
		Gnorm = norm(G[i],"inf");
		/*for(int n = 0; n < 3; ++n)
		{
			for(int m = 0; m<3; ++m)
			{
			    if(G[i](n,m) > Gnorm) //WE assume that G has only positive entries
					Gnorm = G[i](n,m);
			}  
		}*/

		float gg = Grnorm*Grnorm;
		float P = max((1.0-gg)*(1.0-gg)*(1.0-gg),0.0);
		point.val += Gnorm * P;
	}
#if (DEBUG)
	if(point.val > 0)
		std::cout << "[" << point.x << "," << point.y << "," << point.z << "] : " << point.val << "\n";
#endif
}

/*inline void getSimplePhi(Particles & p, mp4Vector & point)
{
	float XminusXnorm, r0, wi;
	vec3f r;
	point.val = 0.0;
	for(int i = 0; i < p.currnp; ++i)
	{
		Grnorm = mag(Gr);	

		//Get norm (i.e. max) of G
		Gnorm = max(max(G[i]));


		float gg = Grnorm*Grnorm;
		float P = max((1.0-gg)*(1.0-gg)*(1.0-gg),0.0);
		point.val += Gnorm * P;
	}
}*/

void createPhi(Particles & p, vec3f * smoothpos, armaMat * Gs, const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{
	for(int k = 1; k < Nz-1; ++k)
		for(int j = 1; j < Ny-1; ++j)
			for(int i = 1; i < Nx-1; ++i)
			{
				int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
				phi[phioffset].x = i*h;
				phi[phioffset].y = j*h;
				phi[phioffset].z = k*h;
				getPhi(p,Gs,smoothpos,phi[phioffset]);
			}
}

/*void createSimplePhi(Particles & p, const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{
	for(int k = 1; k < Nz-1; ++k)
		for(int j = 1; j < Ny-1; ++j)
			for(int i = 1; i < Nx-1; ++i)
			{
				int phioffset = i-1 + (Nx-2)*(j-1 + (Ny-2)*(k-1));
				phi[phioffset].x = i*h;
				phi[phioffset].y = j*h;
				phi[phioffset].z = k*h;
				getSimplePhi(p,phi[phioffset]);
			}
}*/

void 
mesh(Particles & p, const int Nx, const int Ny, const int Nz, const float h, int res, int & numOfTriangles, TRIANGLE *& tri)
{
	std::cout << "Started meshing\n";
	vec3f * smoothPos = new vec3f[p.currnp];
	std::cout << "Smooth particle positions\n";
	smoothParticles(p, smoothPos,0.0, h);
	
	armaMat * Gs = new armaMat[p.currnp];
	std::cout << "Calculate Anisostropic matrices G\n";
	calcAniMatrices(p,Gs);

	mp4Vector * phi;
	phi = new mp4Vector[(Nx-2)*(Ny-2)*(Nz-2)];
	std::cout << "Create the levelset\n";
	createPhi(p,smoothPos,Gs,Nx,Ny,Nz,h,phi);
	//createSimplePhi(p,Nx,Ny,Nz,h,phi);

	delete [] Gs;
	delete [] smoothPos;

	std::cout << "Run marching cubes\n";
	delete [] tri;
	tri = MarchingCubes(Nx-2, Ny-2, Nz-2, 2, phi, LinearInterp, numOfTriangles);
	std::cout << "Marching cubes creates: " << numOfTriangles << " triangles\n";
	delete [] phi;
}





#endif