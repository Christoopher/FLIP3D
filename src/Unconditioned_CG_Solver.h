#ifndef UNCONDITIONED_CG_SOLVER_H
#define UNCONDITIONED_CG_SOLVER_H

#include "Sparse_Matrix.h"
#include <cmath>

struct Uncondioned_CG_Solver
{
	VectorN d; //Search vector
	VectorN r, rprev;
	VectorN tmp;
	double beta, alpha;

	Uncondioned_CG_Solver() {}

	Uncondioned_CG_Solver(int size)
	{
		init(size);
	}

	void init(int size)
	{
		d.init(size);
		r.init(size);
		rprev.init(size);
		tmp.init(size);
	}

	void solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & x);

};

void Uncondioned_CG_Solver::solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & x)
{
	//Calc d(0): d(0) = r(0) = b - A*x(0); 
	x.zero();
	//mtx_mult_vectorN(A,x,tmp);
	//vectorN_VectorN_sub(b,tmp,d);
	b.copy_to(r);
	b.copy_to(d);
	int i = 0;
	double rnorm;
	double rnextnorm;
	while(true)
	{
		//Calc alpha(i): alpha(i) = dot(r,r) / dot( d(i), A*d(i) );
		mtx_mult_vectorN(A,d,tmp);

		for (int i = 0; i < tmp.size; ++i)
		{
			std::cout << tmp.data[i] << "\n";
		}

		rnorm = vectorN_norm2(r);
		double dnorm = vectorN_norm2(d);
		if(rnorm > 10e-6)
			alpha = rnorm/vectorN_dot(d,tmp);
		else
			break;

		//Calc new position x(i+1): x(i+1) = x(i) + alpha(i)*d(i);
		vectorN_add_scale(x, d, alpha); //void

		//Calc new residual r(i+1) = r(i) - alpha(i)*A*d(i);
		// tmp = A*d(i)
		r.copy_to(rprev);
		vectorN_sub_scale(r, tmp, alpha);
		i++;
		rnextnorm = vectorN_norm2(r);
		if(std::sqrt(rnextnorm) < tol || i == maxiterations) 
			break;

		//Calc beta(i+1) = dot(r(i+1),r(i+1)) / dot(r(i),r(i));
		beta = rnextnorm / rnorm;

		//Calc new search vector d(i+1) = r(i+1) + beta(i+1)*d(i);
		vectorN_scale_add(d, r, beta);
	}
}









#endif