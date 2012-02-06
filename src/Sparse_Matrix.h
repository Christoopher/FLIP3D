#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "Grid.h"

struct VectorN
{
	VectorN() : size(0) {}

	VectorN(int dimx_, int dimy_, int dimz_) 
	{
		init(dimx_, dimy_, dimz_);
	}

	void init(int dimx_, int dimy_, int dimz_)
	{
		size = dimx_*dimy_*dimz_;
		dimx = dimx_;
		dimy = dimy_;
		dimz = dimz_;
		data = new double[size];
		zero();
	}

	void copy_to(VectorN & vec) const
	{
		std::memcpy(vec.data, data, size*sizeof(double)); 
	}


	~VectorN()
	{
		delete [] data;
	}

	double & operator()(int i, int j, int k)
	{
		return data[i + dimx*(j + dimy*k)];
	}

	const double & operator()(int i, int j, int k) const
	{
		return data[i + dimx*(j + dimy*k)];
	}

	void zero()
	{
		std::memset(data, 0, size*sizeof(double));
	}

	void write_file_to_matlab(std::string filename, Array3c & marker);

	double norm();


	int size,dimx,dimy,dimz;
	double * data;
};

struct Sparse_Matrix
{
	Sparse_Matrix() : data(NULL)
	{

	}

	Sparse_Matrix(int dimx, int dimy, int dimz) : size(4*dimx*dimy*dimz), dimx(dimx), dimy(dimy), dimz(dimz), stride_y(4*dimx), stride_z(4*dimx*dimy)
	{
		init(dimx,dimy,dimz);
	}

	void init(int dimx_, int dimy_, int dimz_)
	{
		dimx = dimz_; dimy = dimy_; dimz = dimz_;
		size = 4*dimx*dimy*dimz;
		stride_y = 4*dimx;
		stride_z = stride_y*dimy;
		data = new double[size];
		zero();
	}

	~Sparse_Matrix()
	{
		delete [] data;
	}

	double & operator()(int i, int j, int k, int offset)
	{
		return data[i*4 + j*stride_y + k*stride_z + offset];
	}

	const double & operator()(int i, int j, int k, int offset) const
	{
		return data[i*4 + j*stride_y + k*stride_z + offset];
	}


	void zero()
	{
		std::memset(data, 0, size*sizeof(double));
	}

	void write_file_to_matlab(std::string filename, VectorN & vec);


	int size, dimx, dimy, dimz;
	int stride_y, stride_z;
	double * data;
};

void VectorN::write_file_to_matlab(std::string filename, Array3c & marker)
{
	std::ofstream file;
	file.open(filename);
	for(int k = 0; k < dimz; ++k)
		for(int j = 0; j < dimy; ++j)
			for(int i = 0; i < dimx; ++i)
				if(marker(i,j,k) == FLUIDCELL)
					file << (*this)(i,j,k) << ",";
	
	file.close();
}

void write_system_to_matlab(const Sparse_Matrix & A,const VectorN & vec, Array3c & marker)
{

	std::ofstream diag, plusi,plusj,plusk,b;
	diag.open("diag.txt");
	plusi.open("plusi.txt");
	plusj.open("plusj.txt");
	plusk.open("plusk.txt");
	b.open("b.txt");

	int dx = 0,dy = 0,dz = 0;

	int itr = 1;
	for(int k = 0; k < marker.nz; ++k)
	{
		for(int j = 0; j < marker.ny; ++j)
		{
			for(int i = 0; i < marker.nx; ++i)
			{
				if(marker(i,j,k) != FLUIDCELL)
					continue;

				dx++;
				dy++;
				dz++;

				diag << A(i,j,k,0) << ",";
				plusi << A(i,j,k,1) << ",";
				plusj << A(i,j,k,2) << ",";
				plusk << A(i,j,k,3) << ",";
				b << vec(i,j,k) << ",";
// 				}
// 				else
// 				{
// 					diag << A(i,j,k,0);
// 					plusi << A(i,j,k,1);
// 					plusj << A(i,j,k,2);
// 					plusk << A(i,j,k,3);
// 					b << vec(i,j,k);
// 				}
				++itr;
			}
		}
	}
	diag.close();
	plusi.close();
	plusj.close();
	plusk.close();
	b.close();
}


void mtx_mult_vectorN(const Sparse_Matrix & mtx, const VectorN & vec, VectorN & res, Array3c &marker)
{

	//All boundary cells are SOLIDS
	//Thus the boundary cell rows in the Poisson matrix are ZERO: No need to operate on them
	res.zero();
	int offset;
	for(int k=1; k<mtx.dimz-1; ++k)
		for(int j = 1; j < mtx.dimy-1; ++j)
			for (int i = 1; i < mtx.dimx-1; ++i)
			{
				if(marker(i,j,k) != FLUIDCELL)
					continue;

				offset = i-1 + mtx.dimx*(j + mtx.dimy*k); //What row in matrix and/or vector
				res.data[offset] += mtx(i-1,j,k,1)*vec.data[offset]; //i-1,j,k

				offset += 1;
				res.data[offset] += mtx(i,j,k,0)*vec.data[offset]; //i,j,k

				offset += 1;
				res.data[offset] += mtx(i,j,k,1)*vec.data[offset]; //i+1,j,k

				offset = i + mtx.dimx*(j-1 + mtx.dimy*k);
				res.data[offset] += mtx(i,j-1,k,2)*vec.data[offset]; //i,j-1,k

				offset += 2*mtx.dimx;
				res.data[offset] = mtx(i,j,k,2)*vec.data[offset]; //i,j+1,k

				offset = i + mtx.dimx*(j + mtx.dimy*(k-1));
				res.data[offset] += mtx(i,j,k-1,3)*vec.data[offset]; //i,j,k-1

				offset += 2*mtx.dimx*mtx.dimy;
				res.data[offset] += mtx(i,j,k,3)*vec.data[offset]; //i,j,k+1

			}

}

void vectorN_VectorN_add(const VectorN & lhs, const VectorN & rhs, VectorN & res)
{

}

void vectorN_add(VectorN & lhs, const VectorN & rhs)
{
	double * rhsItr = rhs.data;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		*lhsitr = *lhsitr + *rhsItr++;
	}
}

void vectorN_add_scale(VectorN & lhs, const VectorN & rhs, double scale)
{
	double * rhsItr = rhs.data;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		*lhsitr = *lhsitr + scale*(*rhsItr++);
	}
}

void vectorN_scale_add(VectorN & lhs, const VectorN & rhs, double scale)
{
	double * rhsItr = rhs.data;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		*lhsitr = (*lhsitr) * scale + (*rhsItr++);
	}
}

void vectorN_VectorN_sub(const VectorN & lhs, const VectorN & rhs, VectorN & res)
{
	double * resItr = res.data;
	double * rhsItr = rhs.data;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		*resItr++ = *lhsitr - *rhsItr++;
	}
}

void vectorN_sub_scale(VectorN & lhs, const VectorN & rhs, double scale)
{
	double * rhsItr = rhs.data;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		*lhsitr = *lhsitr - scale*(*rhsItr++);
	}
}

double vectorN_dot(const VectorN & lhs, const VectorN & rhs)
{
	double * rhsitr = rhs.data;
	double sum = 0;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		sum += (*lhsitr) * (*rhsitr++);
	}
	return sum;
}

double vectorN_norm2(const VectorN & lhs)
{
	double sum = 0;
	for(double *lhsitr = lhs.data; lhsitr < lhs.data + lhs.size; ++lhsitr)
	{
		sum += (*lhsitr) * (*lhsitr);
	}
	return sum;
}


#endif