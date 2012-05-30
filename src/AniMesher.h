
#ifndef ANI_MESHER_H
#define ANI_MESHER_H

#include "MarchingCubes/MarchingCubes.h"
#include "Particles.h"

#include <cmath>
#include <iostream>
#include <armadillo>
#include "Array3D.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <limits>

#include <time.h>
#include "hr_time.h"

#define OUTSIDE 0
#define FLUIDCELL 1
#define INSIDE 2

typedef arma::mat armaMat;
typedef arma::vec armaVec;


bool calcmesh = true;
armaMat * Gs;
vec3f * smoothPos;
float * density;
mp4Vector * phi;


float isovalue = 0.06; //if simple
//float isovalue = 3.5; //anisotropic
float ks = 3500;//1400;
float kr = 4.0;	
float H = 0.05; //0.07;
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
calcWeight(const Particles & p,const int & i,const float & r, vec3f & weightpos, float & weightsum) 
{
	weightpos[0]=weightpos[1]=weightpos[2] = 0;
	weightsum = 0;
	float weighti;
	float dr;
	for (int j = 0; j < p.currnp; ++j) //alla andra partiklar
	{
		dr = dist(p.pos[i], p.pos[j])/r;
		weighti = max(1.0 - dr*dr*dr,0.0);
		weightsum += weighti;
		weightpos += weighti*p.pos[j];
	}
}

void 
smoothParticles(const Particles & p,float *density, vec3f * smoothpos, float lam,float h)
{
	CStopWatch stopwatch;
	stopwatch.startTimer();
	float r = 2*h;
	vec3f weightpos(0);
	float weightsum;	
#pragma omp parallel for schedule(dynamic) private(weightpos,weightsum)
	for (int i = 0; i < p.currnp; ++i) //i är den partiklen vi kollar på
	{
		calcWeight(p,i,r, weightpos, weightsum); 
		smoothpos[i] = (1.0 - lam)*p.pos[i] + lam*weightpos/weightsum;
		density[i] = weightsum;
	}
	stopwatch.stopTimer();
	std::cout << std::scientific;
	std::cout << stopwatch.getElapsedTime() << "\n";
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
	if(eigenv(0) < 1e-7 || eigenv(1) < 1e-7 || eigenv(2) < 1e-7)
	{
		Gi.eye(3,3);
		Gi*= (1.0/(2.0*H));
	}
	else
	{
		eigenv(0) = 1.0/(eigenv(0)*ks);	
		eigenv(1) = 1.0/(eigenv(1)*ks);
		eigenv(2) = 1.0/(eigenv(2)*ks);
		Gi = (1.0/(2.0*H))*R*diagmat(eigenv)*trans(RT);
	}


	//matrixMultiplyToGetGi(eigenv, R, RT, Gi);
}

/*
* This method constructs G_i matrices for each particle i
* By first calculating weighted covariance matrices C_i
*/
void 
	calcAniMatrices(Particles & p, armaMat * Gs)
{
	CStopWatch stopwatch;
	stopwatch.startTimer();
	float r = 2*H;
#pragma omp parallel for
	for (int i = 0; i < p.currnp; ++i) //i är den partikeln vi kollar på
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
	stopwatch.stopTimer();
	std::cout << std::scientific;
	std::cout << stopwatch.getElapsedTime() << "\n";
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
		sum += (1.0/density[i])*Gnorm * P;

	}

	point.val = sum;
}
void
createPhi(Particles & p, Array3c & marker,vec3f * smoothpos,float * density, armaMat * Gs, const int Nx, const int Ny, const int Nz, const float h,const float res, mp4Vector * phi)
{

	CStopWatch stopwatch;
	stopwatch.startTimer();

	float h2 = h/res;
	std::cout << "Buildlevelset\n";
	int Nzres = Nz*res;
	int Nyres = Ny*res;
	int Nxres = Nx*res;
#pragma omp parallel for
	for(int k = 0; k < Nzres; ++k)
	{
		for(int j = 0; j < Nyres; ++j)
		{
			for(int i = 0; i < Nxres; ++i)
			{
				int phioffset = i + Nxres*(j + Nyres*k);
				/*int i_ = (i * h2) / h + 0.001;
				int j_ = (j * h2) / h + 0.001;
				int k_ = (k * h2) / h + 0.001;*/
				if(	marker(i,j,k) == FLUIDCELL)		
				{
					phi[phioffset].x = (i)*h2;
					phi[phioffset].y = (j)*h2;
					phi[phioffset].z = (k)*h2;
					getPhi(p,Gs,density,smoothpos,phi[phioffset]);
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
					
				}
			
				else if(marker(i,j,k) == OUTSIDE)
				{
					phi[phioffset].x = (i)*h2;
					phi[phioffset].y = (j)*h2;
					phi[phioffset].z = (k)*h2;
					phi[phioffset].val = 0;
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				}
				else//INSIDE
				{
					phi[phioffset].x = (i)*h2;
					phi[phioffset].y = (j)*h2;
					phi[phioffset].z = (k)*h2;
					phi[phioffset].val = -std::numeric_limits<float>::max();
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				}
				//std::cout << "phi: " << phi[phioffset].val << "\n";
				
			}	
			
		}	
	}
	stopwatch.stopTimer();
	std::cout << std::scientific;
	std::cout << stopwatch.getElapsedTime() << "\n";
}

void 
	getSimplePhi(Particles & p,const float h ,mp4Vector & point)
{
	vec3f x(0.0), point_;
	point_[0] = point.x; point_[1] = point.y; point_[2] = point.z; 

	float s, ss, weightsum=0.0,wi, r = h*0.5;
	float R = 2*r;
	point.val = 0.0;

	for(int i = 0; i < p.currnp; ++i)
	{
		s = mag2(point_ - p.pos[i])/R;
		ss = s;
		wi = max(0.0, (1.0 - ss)*(1.0 - ss)*(1.0 - ss));
		weightsum += wi;
		x += p.pos[i]*wi;
	}
	

	if(weightsum > 1e-7)
		x = x/weightsum;

	point.val = mag(point_ - x) - r;

}
void 
	createSimplePhi(Particles & p,Array3c & marker, const int Nx, const int Ny, const int Nz, const float h, mp4Vector * phi)
{
	CStopWatch stopwatch;
	stopwatch.startTimer();

	std::cout << "Buildlevelset\n";
#pragma omp parallel for
	for(int k = 0; k < Nz; ++k)
	{
		for(int j = 0; j < Ny; ++j)
		{
			for(int i = 0; i < Nx; ++i)
			{
				int phioffset = i + Nx*(j + Ny*k);
				/*int i_ = (i * h2) / h + 0.001;
				int j_ = (j * h2) / h + 0.001;
				int k_ = (k * h2) / h + 0.001;*/
				if(	marker(i,j,k) == FLUIDCELL)		
				{
					phi[phioffset].x = (i)*h;
					phi[phioffset].y = (j)*h;
					phi[phioffset].z = (k)*h;
					getSimplePhi(p,h,phi[phioffset]);
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";

				}			
				else if(marker(i,j,k) == OUTSIDE)
				{
					phi[phioffset].x = (i)*h;
					phi[phioffset].y = (j)*h;
					phi[phioffset].z = (k)*h;
					phi[phioffset].val = std::numeric_limits<float>::max();
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				}
				else//INSIDE
				{
					phi[phioffset].x = (i)*h;
					phi[phioffset].y = (j)*h;
					phi[phioffset].z = (k)*h;
					phi[phioffset].val = -std::numeric_limits<float>::max();
					//std::cout << "isovalue: " << phi[phioffset].val << "\n";
				}
			}
		}
	}

	stopwatch.stopTimer();
	std::cout << std::scientific;
	std::cout << stopwatch.getElapsedTime() << "\n";
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
			phi = new mp4Vector[Nx*res*Ny*res*Nz*res];
		}
		else //first time
			phi = new mp4Vector[Nx*res*Ny*res*Nz*res];				
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
	//Ska inte allokera plats för Gs etc innan vi tagit bort partiklar!!!!!!!!!!!
	int Nx = Nx_+1, Ny = Ny_+1, Nz = Nz_+1;
		
	float h = h_;
	if(calcmesh)
	{
		Array3c marker1;
		Array3c marker;
		marker.init(Nx,Ny,Nz);
		marker1.init(Nx,Ny,Nz);

		std::cout << "Build marker grid\n";
#pragma omp parallel for
		for(int n = 0; n < p.currnp; ++n)
		{
			int i = floor(p.pos[n][0]/h + 0.5);
			int j = floor(p.pos[n][1]/h + 0.5);
			int k = floor(p.pos[n][2]/h + 0.5);
			marker1(i,j,k) = INSIDE;
		}

		//Set inside cell to INSIDE
#pragma omp parallel for
		for(int k = 1; k < Nz-1; ++k)
			for(int j = 1; j < Ny-1; ++j)
				for(int i = 1; i < Nx-1; ++i)
				{
					//IF BORDER CELL
					if(marker1(i,j,k) == INSIDE)
					{
					
						if( marker1(i-1,j,k) == OUTSIDE || marker1(i+1,j,k) == OUTSIDE || 
							marker1(i,j-1,k) == OUTSIDE || marker1(i,j+1,k) == OUTSIDE || 
							marker1(i,j,k-1) == OUTSIDE || marker1(i,j,k+1) == OUTSIDE)	
						{
							marker(i,j,k) = FLUIDCELL; //SPARAS
						}
						else
						{
							marker(i,j,k) = marker1(i,j,k);
						}
						
					}
					else
						marker(i,j,k) = marker1(i,j,k);

				}

#pragma omp parallel for
		for(int k = 1; k < Nz-1; ++k)
			for(int j = 1; j < Ny-1; ++j)
				for(int i = 1; i < Nx-1; ++i)
				{
					if(marker(i,j,k) == INSIDE)
					{
						if( marker(i-1,j,k) == FLUIDCELL || marker(i+1,j,k) == FLUIDCELL || 
							marker(i,j-1,k) == FLUIDCELL || marker(i,j+1,k) == FLUIDCELL || 
							marker(i,j,k-1) == FLUIDCELL || marker(i,j,k+1) == FLUIDCELL)	
						{
							marker1(i,j,k) = FLUIDCELL; //SPARAS
						}
						else
							marker1(i,j,k) = marker(i,j,k);
					}
					else
						marker1(i,j,k) = marker(i,j,k);
				}
#pragma omp parallel for
		for(int k = 1; k < Nz-1; ++k)
			for(int j = 1; j < Ny-1; ++j)
				for(int i = 1; i < Nx-1; ++i)
				{
					marker(i,j,k) = marker1(i,j,k);
				}





		// REMOVE INSIDE PARTICLES
		std::cout << "currnp: " << p.currnp << "\n";

#if NDEBUG
		std::vector< std::vector<int> > indices;
		indices.resize(omp_get_max_threads());
#else
		std::vector< int > indices;
#endif

#pragma omp parallel for
		for (int n = 0; n < p.pos.size(); ++n) //Loop over all particles
		{
			int i = floor(p.pos[n][0]/h + 0.5);
			int j = floor(p.pos[n][1]/h + 0.5);
			int k = floor(p.pos[n][2]/h + 0.5);

#if NDEBUG
			if(marker(i,j,k) == FLUIDCELL) //If we want to save
				indices[omp_get_thread_num()].push_back(n);
#else
			if(marker(i,j,k) == FLUIDCELL)
				indices.push_back(n);
#endif

		}

std::vector<vec3f> pos;
#if NDEBUG
		int sum = 0;
		for(int i = 0; i < indices.size(); ++i)
		{
			//sum += removeIndices.size();
			for(int j = 0; j < indices[i].size(); ++j)
			{
				++sum;
				//p.remove(removeIndices[i][j]);
				pos.push_back( p.pos[indices[i][j] ]);
			}
		}
		
#else
		for(int j = 0; j < indices.size(); ++j)
		{
			pos.push_back( p.pos[indices[j] ]);
		}
#endif

		p.pos.clear();
		copy(pos.begin(),pos.end(),std::back_inserter(p.pos));
		p.currnp = p.pos.size();
		std::cout << "currnp: " << p.pos.size() << "\n";

		Array3c * from = &marker1;
		Array3c * to = &marker;
		for(int x = 0; x<2; ++x)
		{
			std::swap(from,to);

#pragma omp parallel for
			for(int k = 1; k < Nz-1; ++k)
			{
				for(int j = 1; j < Ny-1; ++j)
				{
					for(int i = 1; i < Nx-1; ++i)
					{
						if((*from)(i,j,k) == OUTSIDE) //ARICELL will be expanced to FLUIDCELL if connected to border
						{
							if( (*from)(i-1,j,k) == FLUIDCELL || (*from)(i+1,j,k) == FLUIDCELL || 
								(*from)(i,j-1,k) == FLUIDCELL || (*from)(i,j+1,k) == FLUIDCELL || 
								(*from)(i,j,k-1) == FLUIDCELL || (*from)(i,j,k+1) == FLUIDCELL)		
							{
								(*to)(i,j,k) = FLUIDCELL;
							}
							else
								(*to)(i,j,k) = (*from)(i,j,k);
						}
						else
							(*to)(i,j,k) = (*from)(i,j,k);

					}	
				}	
			}


		}


				int fluids=0, solids = 0, outsides= 0;

		for(int k = 1; k < Nz-1; ++k)
			for(int j = 1; j < Ny-1; ++j)
				for(int i = 1; i < Nx-1; ++i)
					if(marker(i,j,k) == INSIDE)
						++solids;
					else if(marker(i,j,k) == FLUIDCELL)
						++fluids;
					else if(marker(i,j,k) == OUTSIDE)
						++outsides;
		

		std::cout << "TOTALEN: " << solids+fluids << "\n";
		std::cout << "inside: " << solids << "\n";
		std::cout << "fluids: " << fluids << "\n";
		std::cout << "outsides: " << outsides << "\n";


		initMeshData(p.currnp, Nx, Ny, Nz, res);

		std::cout << "Started meshing\n";			

		//std::cout << "Smooth particle positions\n";
		//smoothParticles(p,density, smoothPos,0.9, h_*0.5);


		//std::cout << "Calculate Anisostropic matrices G\n";
		//calcAniMatrices(p,Gs);

		std::cout << "Create the levelset\n";
		//createPhi(p,marker,smoothPos,density,Gs,Nx,Ny,Nz,h,res,phi); isovalue = 30;
		createSimplePhi(p,marker,Nx,Ny,Nz,h,phi); isovalue = 0.06;
		
	}

	std::cout << "Run marching cubes\n";
	delete [] tri;
	tri = MarchingCubes(Nx, Ny, Nz, isovalue, phi, LinearInterp, numOfTriangles);
	std::cout << "Marching cubes creates: " << numOfTriangles << " triangles\n";

	calcmesh = false;
}





#endif