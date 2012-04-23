/////////////////////////////////////////////////////////////////////////////////////////////
//	FileName:	MarchingCubes.h
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//	website	:	www.angelfire.com/linux/myp
//	date	:	July 2002
//	
//	Description:	Marching Cubes Algorithm
/////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include "mpVector.h"
#include "MCTable.h"

typedef struct {
	mpVector p[3];
	mpVector norm;
} TRIANGLE;

//does Linear Interpolation between points p1 and p2 (they already contain their computed values)
mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value);

////////////////////////////////////////////////////////////////////////////////////////
//POINTERS TO FUNCTIONS
//pointer to function which computes if point p is outside the surface
typedef bool (*OUTSIDE)(mpVector);

//pointer to function which determines the point of intersection of the edge with 
//the isosurface between points p1 and p2
//any other information is passed in the void array mcInfo
typedef mpVector (*INTERSECTION)(mp4Vector, mp4Vector, float);

//pointer to function which computes the value at point p
typedef float (*FORMULA)(mpVector);

///// the MARCHING CUBES algorithm itself /////

//	1A).
//takes number of cells (ncellsX, ncellsY, ncellsZ) to subdivide on each axis
// minValue is the third argument for INTERSECTION function
// array of length (ncellsX+1)(ncellsY+1)(ncellsZ+1) of mp4Vector points containing coordinates and values
// function of type mpVector (mp4Vector p1, mp4Vector p2) intersection, which determines the 
//  point of intersection of the surface and the edge between points p1 and p2
//returns pointer to triangle array and the number of triangles in numTriangles
//note: array of points is first taken on x axis, then y and then z. So for example, if u iterate through it in a
//       for loop, have indexes i, j, k for x, y, z respectively, then to get the point you will have to make the
//		 following index: i*(ncellsY+1)*(ncellsZ+1) + j*(ncellsZ+1) + k .
//		Also, the array starts at the minimum on all axes.
//TODO: another algorithm which takes array of JUST values. Coordinates then start at farthest, lower-right corner.
TRIANGLE* MarchingCubes(int ncellsX, int ncellsY, int ncellsZ, float minValue, mp4Vector * points,  
									INTERSECTION intersection, int &numTriangles);
//  1B).
//same as above only does linear interpolation so no INTERSECTION function is needed
TRIANGLE* MarchingCubesLinear(int ncellsX, int ncellsY, int ncellsZ, float minValue, 
									mp4Vector * points, int &numTriangles);


//	2A).
//takes dimensions (minx,maxx,miny,...) and the number of cells (ncellsX,...) to subdivide on each axis
// minValue is the third argument for INTERSECTION function
// function of type float (mpVector p) formula, which computes value of p at its coordinates
// function of type mpVector (mp4Vector p1, mp4Vector p2) intersection, which determines the 
//  point of intersection of the surface and the edge between points p1 and p2
// saves number of triangles in numTriangles and the pointer to them is returned
// (note: mins' and maxs' are included in the algorithm)
TRIANGLE* MarchingCubes(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
									int ncellsX, int ncellsY, int ncellsZ, float minValue, 
									FORMULA formula, INTERSECTION intersection, int &numTriangles);
//	2B).
//same as above only does linear interpolation to determine intersection of edge and surface
// INTERSECTION function is no more needed
TRIANGLE* MarchingCubesLinear(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
									int ncellsX, int ncellsY, int ncellsZ, float minValue, 
									FORMULA formula, int &numTriangles);
#endif