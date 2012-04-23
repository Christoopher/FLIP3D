/////////////////////////////////////////////////////////////////////////////////////////////
//	FileName:	MarchingCubes.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//	website	:	www.angelfire.com/linux/myp
//	date	:	July 2002
//	
//	Description:	Marching Cubes Algorithm
/////////////////////////////////////////////////////////////////////////////////////////////

#include "MarchingCubes.h"

mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value)
{
	mpVector p;
	if(p1.val != p2.val && value != p1.val)
		p = (mpVector)p1 + ((mpVector)p2 - (mpVector)p1)/(p2.val - p1.val)*(value - p1.val);
	else 
		p = (mpVector)p1;
	return p;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//	MARCHING CUBES	//

int dimx = 0, dimy = 0;

inline int offset(int x, int y, int z)
{
	return (x + dimx*(y + dimy*z));
}

//  VERSION  1A).  //
TRIANGLE* MarchingCubes(int ncellsX, int ncellsY, int ncellsZ, float minValue, mp4Vector * points,  
										INTERSECTION intersection, int &numTriangles)
{
	dimx = ncellsX;
	dimy = ncellsY;
	TRIANGLE * triangles = new TRIANGLE[3*ncellsX*ncellsY*ncellsZ];//this should be enough space, if not change 4 to 5
	numTriangles = int(0);

	//go through all the points
	for(int k=0; k < ncellsZ-1; k++)	//z axis
		for(int j=0; j < ncellsY-1; j++)		//y axis
			for(int i=0; i < ncellsX-1; i++)			//x axis
			{
				//initialize vertices
				mp4Vector verts[8];
  /*(step 3)*/  verts[0] = points[offset(i,j,k)];
				verts[1] = points[offset(i+1,j,k)];
				verts[2] = points[offset(i+1,j,k+1)];
				verts[3] = points[offset(i,j,k+1)];
				verts[4] = points[offset(i,j+1,k)];
				verts[5] = points[offset(i+1,j+1,k)];
				verts[6] = points[offset(i+1,j+1,k+1)];
				verts[7] = points[offset(i,j+1,k+1)];
				
				//get the index
				int cubeIndex = int(0);
				for(int n=0; n < 8; n++)
				{
   /*(step 4)*/		if(verts[n].val <= minValue) 
						cubeIndex |= (1 << n);
				}

				//check if its completely inside or outside
   /*(step 5)*/ if(!edgeTable[cubeIndex]) 
					continue;
			
				//get intersection vertices on edges and save into the array
				mpVector intVerts[12];
   /*(step 6)*/ if(edgeTable[cubeIndex] & 1) intVerts[0] = intersection(verts[0], verts[1], minValue);
				if(edgeTable[cubeIndex] & 2) intVerts[1] = intersection(verts[1], verts[2], minValue);
				if(edgeTable[cubeIndex] & 4) intVerts[2] = intersection(verts[2], verts[3], minValue);
				if(edgeTable[cubeIndex] & 8) intVerts[3] = intersection(verts[3], verts[0], minValue);
				if(edgeTable[cubeIndex] & 16) intVerts[4] = intersection(verts[4], verts[5], minValue);
				if(edgeTable[cubeIndex] & 32) intVerts[5] = intersection(verts[5], verts[6], minValue);
				if(edgeTable[cubeIndex] & 64) intVerts[6] = intersection(verts[6], verts[7], minValue);
				if(edgeTable[cubeIndex] & 128) intVerts[7] = intersection(verts[7], verts[4], minValue);
				if(edgeTable[cubeIndex] & 256) intVerts[8] = intersection(verts[0], verts[4], minValue);
				if(edgeTable[cubeIndex] & 512) intVerts[9] = intersection(verts[1], verts[5], minValue);
				if(edgeTable[cubeIndex] & 1024) intVerts[10] = intersection(verts[2], verts[6], minValue);
				if(edgeTable[cubeIndex] & 2048) intVerts[11] = intersection(verts[3], verts[7], minValue);

				//now build the triangles using triTable
				for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
   /*(step 7)*/ 	triangles[numTriangles].p[0] = intVerts[triTable[cubeIndex][n+2]];
					triangles[numTriangles].p[1] = intVerts[triTable[cubeIndex][n+1]];
					triangles[numTriangles].p[2] = intVerts[triTable[cubeIndex][n]];
   /*(step 8)*/ 	triangles[numTriangles].norm = ((triangles[numTriangles].p[1] - 
						triangles[numTriangles].p[0]).Cross(triangles[numTriangles].p[2] - 
						triangles[numTriangles].p[0])).Normalize();
					numTriangles++;
				}
			
			}	//END OF FOR LOOP
		
		//free all the wasted space
		TRIANGLE * retTriangles = new TRIANGLE[numTriangles];
		for(int i=0; i < numTriangles; i++)
			retTriangles[i] = triangles[i];
		delete [] triangles;
	
	return retTriangles;
}


//	VERSION  1B).  //
TRIANGLE* MarchingCubesLinear(int ncellsX, int ncellsY, int ncellsZ, float minValue, 
									mp4Vector * points, int &numTriangles)
{
	return MarchingCubes(ncellsX, ncellsY, ncellsZ, minValue, points, LinearInterp, numTriangles);
}


//	VERSION  2A).  //
TRIANGLE* MarchingCubes(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
							int ncellsX, int ncellsY, int ncellsZ, float minValue, 
							FORMULA formula, INTERSECTION intersection, int &numTriangles)
{
	//space is already defined and subdivided, staring with step 3
	//first initialize the points
	mp4Vector * mcDataPoints = new mp4Vector[(ncellsX+1)*(ncellsY+1)*(ncellsZ+1)];
	mpVector stepSize((mcMaxX-mcMinX)/ncellsX, (mcMaxY-mcMinY)/ncellsY, (mcMaxZ-mcMinZ)/ncellsZ);
	
	int YtimesZ = (ncellsY+1)*(ncellsZ+1);	//for extra speed
	for(int i=0; i < ncellsX+1; i++) {
		int ni = i*YtimesZ;						//for speed
		float vertX = mcMinX + i*stepSize.x;
		for(int j=0; j < ncellsY+1; j++) {
			int nj = j*(ncellsZ+1);				//for speed
			float vertY = mcMinY + j*stepSize.y;
			for(int k=0; k < ncellsZ+1; k++) {
				mp4Vector vert(vertX, vertY, mcMinZ + k*stepSize.z, 0);
				vert.val = formula((mpVector)vert);
   /*(step 3)*/ mcDataPoints[ni + nj + k] = vert;
			}
		}
	}
	//then run Marching Cubes (version 1A) on the data
	return MarchingCubes(ncellsX, ncellsY, ncellsZ, minValue, mcDataPoints, intersection, numTriangles);
}

//	VERSION  2B).  //
TRIANGLE* MarchingCubesLinear(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ, 
								int ncellsX, int ncellsY, int ncellsZ, float minValue, 
								FORMULA formula, int &numTriangles)
{
	return MarchingCubes(mcMinX, mcMaxX, mcMinY, mcMaxY, mcMinZ, mcMaxZ, ncellsX, ncellsY, ncellsZ, minValue,
		formula, LinearInterp, numTriangles);
}