/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	FileName:	MCExample.h and MCExample.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu	or  mikepolyakov@hotmail.com
//	Website	:	www.angelfire.com/linux/myp
//	Date	:	7/16/2002
//	
//
//	This demonstrates using the Marching Cubes Algorithm (MarchingCubes.h and MarchingCubes.cpp).
//	needed files: mpVector.cpp and mpVector.h (visit my website to get these files)
//	(Note: Uses wxWindows (www.wxwindows.org))
//
//	Controls:	Up/Down		- walks forwards and backwards
//				Left/Right	- strafes to the side
//				Del/End		- strafes up and down
//				Space		- runs Marching Cubes algorithm with minValue
//				Ctrl		- decrements minValue and runs Marching Cubes
//				Insert		- increments minValue and runs Marching Cubes
//				F12			- toggles the fullscreen mode
//				q/w			- decrements/increments number of subdivisions on X-axis
//				a/s			- decrements/increments number of subdivisions on Y-axis
//				z/x			- decrements/increments number of subdivisions on Z-axis
//				e/r			- decreses/increases number of subdivisions on all axis simultaneously
//
//	
//	email me with any suggestions or bugs...
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TEST_H_
#define TEST_H_

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/glcanvas.h>
#include "MarchingCubes.h"

class GLCanvas;

class MyApp : public wxApp
{
public:
    virtual bool OnInit();

	
};

class GLCanvas : public wxGLCanvas
{
public:
	GLCanvas(wxFrame* parent);
	void RunMarchingCubesTest();
	void InitData();

	void OnPaint(wxPaintEvent& event);
	void OnEraseBackground(wxEraseEvent& event) { /* does nothing to avoid flashing */ }
	void OnSize(wxSizeEvent& event);
	void OnKey(wxKeyEvent& event);
	void OnMouseMotion(wxMouseEvent& event);

private:
	//camera movement variables
	float stepZ, stepX, stepY;		//increment step for moving camera	
	float transZ, transX, transY;	//translates camera by this amount on each axis
	float angleX, angleY, angleZ;	//angles of rotation
	float saX, saY, saZ;			//increments on angles
	
	//minimal value and number of cells on each axis: passed into Marching Cubes
	float minValue;
	int nX, nY, nZ;
	//data points passed to Marching Cubes
	mp4Vector * mcPoints;
	//data returned by Marching Cubes
	TRIANGLE * Triangles;
	int numOfTriangles;
	
	wxFrame * parent;

DECLARE_EVENT_TABLE()	
};

//boundary values for Marching Cubes
#define MINX -4.0
#define MAXX 5.0
#define MINY -8.0
#define MAXY 8.0
#define MINZ -8.0
#define MAXZ 8.0

//function which computes a value at point p
float Potential(mpVector p);

//draws coordinates of length s on all axis
void drawCoordinates(float s);

#endif
