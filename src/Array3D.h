#ifndef ARRAY_3D_H
#define ARRAY_3D_H

#include <cstdio>
#include <cmath>

template<class T>
struct Array3{
	int nx, ny,nz;
	int size;
	T *data;

	Array3()
		:nx(0), ny(0), nz(0), size(0), data(0)
	{}

	Array3(int nx_, int ny_, int nz_) : nx(0), ny(0), nz(0), size(0), data(0) 
	{ 
		init(nx_,ny_,nz_); 
	}

	void init(int nx_, int ny_, int nz_)
	{
		delete_memory();
		nx = nx_; ny = ny_; nz = nz_;
		size=nx*ny*nz;
		data=new T[size];
		zero();
	}

	~Array3()
	{
		delete_memory(); 
	}

	void delete_memory()
	{
		delete[] data; data=0;
		nx=ny=size=0;
	}

	const T &operator() (int i, int j, int k) const
	{ 
		return data[ i+nx*(j + k*ny) ]; 
	}

	T &operator() (int i, int j, int k)
	{ 
		return data[ i+nx*(j + k*ny) ]; 
	}

	T trilerp(int i, int j, int k, T fx, T fy, T fz)
	{ 
		T fval = (1-fx)*((1-fy)*(*this)(i,j,k)+fy*(*this)(i,j+1,k))+fx*((1-fy)*(*this)(i+1,j,k)+fy*(*this)(i+1,j+1,k)); 
		T bval = (1-fx)*((1-fy)*(*this)(i,j,k+1)+fy*(*this)(i,j+1,k+1))+fx*((1-fy)*(*this)(i+1,j,k+1)+fy*(*this)(i+1,j+1,k+1)); 
		return (1-fz)*fval + fz*bval;		
		
	}

	void copy_to(Array3 &a) const
	{ 
		std::memcpy(a.data, data, size*sizeof(T)); 
	}

	T infnorm() const
	{ 
		T r=0;
		for(int i=0; i<size; ++i)
			if(!(std::fabs(data[i])<=r)) 
				r=std::fabs(data[i]);
		return r;
	}

	void zero()
	{ 
		std::memset(data, 0, size*sizeof(T)); 
	}

	double dot(const Array3 &a) const
	{
		double r=0;
		for(int i=0; i<size; ++i)
			r+=data[i]*a.data[i];
		return r;
	}

	void increment(double scale, const Array3 &a)
	{ 
		for(int i=0; i<size; ++i) 
			data[i]+=scale*a.data[i]; 
	}

	void scale_and_increment(double scale, const Array3 &a)
	{ 
		for(int i=0; i<size; ++i) 
			data[i]=scale*data[i]+a.data[i]; 
	}
};

typedef Array3<float> Array3f;
typedef Array3<double> Array3d;
typedef Array3<char> Array3c;




#endif