#ifndef UNCONDITIONED_CG_SOLVER_H
#define UNCONDITIONED_CG_SOLVER_H

#include "Sparse_Matrix.h"
#include "Array3D.h"
#include <cmath>

struct Uncondioned_CG_Solver
{
	VectorN d; //Search vector
	VectorN z;
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
		z.init(dimx, dimy, dimz);
		Adj.init(dimx, dimy, dimz);
	}

	void apply_precond(const Sparse_Matrix & A, const Sparse_Matrix & precond,const VectorN &r, VectorN &z,const Array3c & marker);

	void solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & x, Array3c & marker);
	void solve_precond(const Sparse_Matrix & A,const VectorN & b,const Sparse_Matrix & precond,int maxiterations, double tol, VectorN & pressure, Array3c &marker);

};

void Uncondioned_CG_Solver::solve(const Sparse_Matrix & A,const VectorN & b,int maxiterations, double tol, VectorN & pressure, Array3c &marker)
{
	b.copy_to(r);
	double rinfnorm = r.infnorm();
	if(rinfnorm == 0.0)
		return;

	tol = tol*rinfnorm;

	double rnorm = vectorN_norm2(r);
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
		
		alpha = rnorm/vectorN_dot(d,Adj);
		
		//Calc new position x(i+1): x(i+1) = x(i) + alpha(i)*d(i);
		vectorN_add_scale(pressure, d, alpha); //void

		//Calc new residual r(i+1) = r(i) - alpha(i)*A*d(i);
		vectorN_sub_scale(r, Adj, alpha); 
		


		i++; //We have now moved one step
		if(r.infnorm() <= tol || i == maxiterations)
		{
			std::cout << std::scientific;
			std::cout << "CG: " << i << " iterations, " << "norm_squared = " << rnextnorm <<  "\n";
			//pressure.write_file_to_matlab("x.txt",marker);
			return;
		}

		// *** Calculate new search direction; d(j+1) ***

		//Calc beta(i+1) = dot(r(i+1),r(i+1)) / dot(r(i),r(i));
		//the norm of the new residual
		rnextnorm = vectorN_norm2(r); //r = r(i+1)
		beta = rnextnorm / rnorm;

		//Calc new search vector d(i+1) = r(i+1) + beta(i+1)*d(i);
		vectorN_scale_add(d, r, beta);
		rnorm = rnextnorm;
	}
}

void Uncondioned_CG_Solver::apply_precond(const Sparse_Matrix & A, const Sparse_Matrix & precond,const VectorN &r, VectorN &z,const Array3c & marker)
{
	//Solve Lq = r
	double t = 0;
	Adj.zero();
	for(int k = 1; k < A.dimz-1; ++k)
		for(int j = 1; j < A.dimy-1; ++j)
			for(int i = 1; i < A.dimx-1; ++i)
			{
				if(marker(i,j,k) == FLUIDCELL)
				{
					t = r(i,j,k) - A(i-1,j,k,1)*precond(i-1,j,k,0) * Adj(i-1,j,k)
								 - A(i,j-1,k,2)*precond(i,j-1,k,0) * Adj(i,j-1,k)
								 - A(i,j,k-1,3)*precond(i,j,k-1,0) * Adj(i,j,k-1);
					Adj(i,j,k) = t * precond(i,j,k,0);
				}
			}

	//Solve Lt z = q	

	z.zero();
	for(int k = A.dimz-2; k > 0; --k)
		for(int j = A.dimy-2; j > 0; --j)
			for(int i = A.dimx-2; i > 0; --i)
			{
				if(marker(i,j,k) == FLUIDCELL)
				{
					t = Adj(i,j,k) - A(i,j,k,1)*precond(i,j,k,0) * z(i+1,j,k)
								   - A(i,j,k,2)*precond(i,j,k,0) * z(i,j+1,k)
								   - A(i,j,k,3)*precond(i,j,k,0) * z(i,j,k+1);
					z(i,j,k) = t * precond(i,j,k,0);						
						
				}
			}
}


void Uncondioned_CG_Solver::solve_precond(const Sparse_Matrix & A,const VectorN & b,const Sparse_Matrix & precond,int maxiterations, double tol, VectorN & pressure, Array3c &marker)
{
	b.copy_to(r);
	double rinfnorm = r.infnorm();
	if(rinfnorm == 0.0)
		return;

	tol = tol*rinfnorm;
	pressure.zero();


	//z(0) = precond * r0
	apply_precond(A,precond,r,z,marker);
	z.copy_to(d); // d(0) = r(0) = b

	double rznorm = vectorN_dot(z,r);
	if(rznorm == 0.0)
		return;


	int i = 0;
	double rznextnorm = 0;

	while(true)
	{
		mtx_mult_vectorN(A,d,z,marker);
		alpha = rznorm/vectorN_dot(d,z);

		//Calc new position x(i+1): x(i+1) = x(i) + alpha(i)*d(i);
		vectorN_add_scale(pressure, d, alpha); //void

		//Calc new residual r(i+1) = r(i) - alpha(i)*A*d(i);
		vectorN_sub_scale(r, z, alpha); 

		i++; //We have now moved one step
		if(r.infnorm() <= tol || i == maxiterations)
		{
			//std::cout << std::scientific;
			//std::cout << "CG: " << i << " iterations, " << "norm_squared = " << rznextnorm <<  "\n";
			//pressure.write_file_to_matlab("x.txt",marker);
			return;
		}

		// *** Calculate new search direction; d(j+1) ***
		apply_precond(A,precond,r,z,marker);
		//Calc beta(i+1) = dot(r(i+1),r(i+1)) / dot(r(i),r(i));
		//the norm of the new residual
		rznextnorm = vectorN_dot(r,z); //r = r(i+1)
		beta = rznextnorm / rznorm;

		//Calc new search vector d(i+1) = r(i+1) + beta(i+1)*d(i);
		vectorN_scale_add(d, z, beta);
		rznorm = rznextnorm;
	}
}









#endif