#ifndef UNCONDITIONED_CG_SOLVER_H
#define UNCONDITIONED_CG_SOLVER_H

#include "Sparse_Matrix.h"
#include "Array3D.h"
#include <cmath>

struct Uncondioned_CG_Solver
{
	VectorN d; //Search vector
	VectorN r;
	VectorN Adj;
	double beta, alpha;
		
	Uncondioned_CG_Solver() {}

	Uncondioned_CG_Solver(int dimx, int dimy, int dimz)
	{
		init(dimx,dimy,dimz);
	}

	void init(int dimx, int dimy, int dimz)
	{
		d.init(dimx, dimy, dimz);
		r.init(dimx, dimy, dimz);
		Adj.init(dimx, dimy, dimz);
	}

	void solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & x, Array3c & marker);

};

void Uncondioned_CG_Solver::solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & pressure, Array3c &marker)
{
	//write_system_to_matlab(A,b,marker); //Debug
	b.copy_to(r);
	if(r.infnorm() == 0.0) // we should check sqrt(rnorm) but we don't want to so we check rnorm < tol^2
		return;

	tol = tol*r.infnorm();

	double rnorm = vectorN_norm2(r); //r(0) = b - Ax(0) = b
	if(rnorm == 0.0)
		return;

	Adj.zero();
	b.copy_to(d); // d(0) = r(0) = b

	int i = 0;
	double rnextnorm = 0;

	while(true)
	{
		//Calc alpha(i): alpha(i) = dot(r,r) / dot( d(i), A*d(i) );
		mtx_mult_vectorN(A,d,Adj,marker);
		//write_system_to_matlab(A,b,marker); //Debug
		//d.write_file_to_matlab("d.txt",marker);
		//tmp == A*d(i)
		double dnorm  = vectorN_dot(d,Adj);

		if(dnorm == 0 || dnorm != dnorm)
			int apa = 0;

		alpha = rnorm/dnorm;
		
		if(alpha != alpha)
			int apa = 0;
		//Calc new position x(i+1): x(i+1) = x(i) + alpha(i)*d(i);
		vectorN_add_scale(pressure, d, alpha); //void

		//Calc new residual r(i+1) = r(i) - alpha(i)*A*d(i);
		vectorN_sub_scale(r, Adj, alpha); 
		


		i++; //We have now moved one step
		if(r.infnorm() <= tol || i == maxiterations)
		{
			std::cout << std::scientific;
			std::cout << "CG: " << i << " iterations, " << "norm_squared = " << rnextnorm <<  "\n";
			pressure.write_file_to_matlab("x.txt",marker);
			break;
		}

		rnextnorm = vectorN_norm2(r); //r = r(i+1)

		// *** Calculate new search direction; d(j+1) ***

		//Calc beta(i+1) = dot(r(i+1),r(i+1)) / dot(r(i),r(i));
		beta = rnextnorm / rnorm;

		//Calc new search vector d(i+1) = r(i+1) + beta(i+1)*d(i);
		vectorN_scale_add(d, r, beta);
		rnorm = rnextnorm;
	}
}









#endif