/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	FileName:	MCExample.h and MCExample.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu	or  mikepolyakov@hotmail.com
//	Website	:	www.angelfire.com/linux/myp	
//	Date	:	7/16/2002
//
//	This demonstrates using the Marching Cubes Algorithm. 
//	Refer to MCExample.h for more information
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "MCExample.h"

// MyApp //////////////////////////////////////////////////////////////////////////////////////
//starts the main application
IMPLEMENT_APP(MyApp)
//open a frame and glcanvas
bool MyApp::OnInit()
{
	wxFrame * frame = new wxFrame((wxWindow*)NULL, -1, "Marching Cubes Example", wxPoint(0, 60), wxSize(600, 600));
	GLCanvas * canvas = new GLCanvas(frame);
	frame->Show(TRUE);
	return TRUE;
}
// END MyApp //////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////

// Glcanvas ///////////////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(GLCanvas, wxGLCanvas)
	EVT_PAINT(GLCanvas::OnPaint)
	EVT_SIZE(GLCanvas::OnSize)
	EVT_ERASE_BACKGROUND(GLCanvas::OnEraseBackground)
	EVT_KEY_DOWN(GLCanvas::OnKey)
	EVT_MOTION(GLCanvas::OnMouseMotion)
END_EVENT_TABLE()

//constructor
GLCanvas::GLCanvas(wxFrame* parent) : wxGLCanvas((wxWindow*)parent, -1, wxPoint(0, 0), parent->GetSize()),
	stepZ(0.2), stepX(0.2), stepY(0.2), transZ(-8), transX(0), transY(0), angleX(0), angleY(0), angleZ(0),
	saX(5), saY(5), saZ(5), minValue(1.8), nX(20), nY(20), nZ(20)
{ 
	//init openGL
	SetCurrent(); 
	glClearColor(0, 0, 0, 1); 
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	float light_color[] = {1, 1, 1, 1};
	float light_position[] = {5, 2, 7, 1};
	glLightfv(GL_LIGHT0, GL_COLOR, light_color);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	//init material properties
	float ambient_color[] = {0, 1, 1, 1};
	float diff_color[] = {1, 1, 1, 1};
	float spec_color[] = {1, 1, 1, 1};
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diff_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_color);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0);
	glColor4fv(ambient_color);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
	this->parent = parent;
	
	//initialize data to be passed to Marching Cubes
	mcPoints = new mp4Vector[(nX+1)*(nY+1)*(nZ+1)];
	mpVector stepSize((MAXX-MINX)/nX, (MAXY-MINY)/nY, (MAXZ-MINZ)/nZ);
	for(int i=0; i < nX+1; i++)
		for(int j=0; j < nY+1; j++)
			for(int k=0; k < nZ+1; k++) {
				mp4Vector vert(MINX+i*stepSize.x, MINY+j*stepSize.y, MINZ+k*stepSize.z, 0);
				vert.val = Potential((mpVector)vert);
				mcPoints[i*(nY+1)*(nZ+1) + j*(nZ+1) + k] = vert;
			}
	//runs Marching Cubes
	Triangles = MarchingCubes(nX, nY, nZ, minValue, mcPoints, LinearInterp, numOfTriangles);
}

//re-initializes the data to be used by Marching Cubes
void GLCanvas::InitData()
{
	delete [] mcPoints;	//first free the previous allocated memory
	mcPoints = new mp4Vector[(nX+1)*(nY+1)*(nZ+1)];
	mpVector stepSize((MAXX-MINX)/nX, (MAXY-MINY)/nY, (MAXZ-MINZ)/nZ);
	for(int i=0; i < nX+1; i++)
		for(int j=0; j < nY+1; j++)
			for(int k=0; k < nZ+1; k++) {
				mp4Vector vert(MINX+i*stepSize.x, MINY+j*stepSize.y, MINZ+k*stepSize.z, 0);
				vert.val = Potential((mpVector)vert);
				mcPoints[i*(nY+1)*(nZ+1) + j*(nZ+1) + k] = vert;
	}
}

void GLCanvas::RunMarchingCubesTest()
{
	delete [] Triangles;	//first free the previous allocated memory
	Triangles = MarchingCubes(nX, nY, nZ, minValue, mcPoints, LinearInterp, numOfTriangles);
}

//moves the camera around
//increments and decrements the minimum value, runs Marching Cubes and switches to full screen mode
void GLCanvas::OnKey(wxKeyEvent& event)
{
	switch(event.GetKeyCode()) {
		case WXK_UP: transZ += stepZ; Refresh(FALSE); break;
		case WXK_DOWN:transZ -= stepZ; Refresh(FALSE); break;
		case WXK_INSERT: minValue += 0.05; RunMarchingCubesTest(); Refresh(FALSE); break;
		case WXK_CONTROL: minValue -= 0.05; RunMarchingCubesTest(); Refresh(FALSE); break;
		case WXK_SPACE: RunMarchingCubesTest(); Refresh(FALSE); break;
		case WXK_LEFT: transX -= stepX; Refresh(FALSE); break;
		case WXK_RIGHT: transX += stepX; Refresh(FALSE); break;
		case WXK_DELETE: transY += stepY; Refresh(FALSE); break;
		case WXK_END: transY -= stepY; Refresh(FALSE); break;
		case WXK_F12: parent->IsFullScreen() ? parent->ShowFullScreen(FALSE) : parent->ShowFullScreen(TRUE); 
					  Refresh(FALSE); break;
		case 'Q': nX > 1 ? nX--:nX; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'W': nX++; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'A': nY > 1 ? nY--:nY; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'S': nY++; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'Z': nZ > 1 ? nZ--:nZ; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'X': nZ++; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		case 'E': if(nX > 3 && nY > 3 && nZ > 3) {	//make sure non are less or equal to zero
						nX -= 2; nY -= 2; nZ -= 2; 
						InitData(); 
						RunMarchingCubesTest(); 
						Refresh(FALSE);} 
					break;
		case 'R': nX += 2; nY += 2; nZ += 2; InitData(); RunMarchingCubesTest(); Refresh(FALSE); break;
		default: break;
	}
}

//controls the mouse (camera rotation)
void GLCanvas::OnMouseMotion(wxMouseEvent& event)
{
	static int x = event.m_x, y = event.m_y;
	if (event.Dragging())
	{
		if(event.m_shiftDown) {
			angleZ += event.m_x < x ? -saZ : (event.m_x == x ? 0 : saZ);
			angleX += event.m_y < y ? -saX : (event.m_y == y ? 0 : saX);
		} else {
			angleY += event.m_x < x ? -saY : (event.m_x == x ? 0 : saY);
			angleX += event.m_y < y ? -saX : (event.m_y == y ? 0 : saX);
		}
		x = event.m_x; 
		y = event.m_y;
		Refresh(FALSE);
	}
}

//called by resize events
void GLCanvas::OnSize(wxSizeEvent& event)
{
	int w, h;
    GetClientSize(&w, &h);
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	float ratio = (float)h/(float)w;
	glFrustum(-1,1, -ratio, ratio, 1, 100);
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

//does rotation and translation before drawing the triangles returned by Marching Cubes
void GLCanvas::OnPaint(wxPaintEvent& event)
{
	wxPaintDC dc(this);
	//clears the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glPushMatrix();
		glTranslatef(transX, transY, transZ);
		glRotatef(angleX, 1, 0, 0);
		glRotatef(angleY, 0, 1, 0);
		glRotatef(angleZ, 0, 0, 1);
		//draws triangles
		glBegin(GL_TRIANGLES);
			for(int i=0; i < numOfTriangles; i++){
				glNormal3f(Triangles[i].norm.x, Triangles[i].norm.y, Triangles[i].norm.z);
				for(int j=0; j < 3; j++)	
					glVertex3f(Triangles[i].p[j].x,Triangles[i].p[j].y,Triangles[i].p[j].z);
			}
		glEnd();
		//draws coordinates
		drawCoordinates(6);
	glPopMatrix();
	
	//tells OpenGL to execute the above commands
	glFlush();
	SwapBuffers();
}
// END GLCanvas //////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


float Potential(mpVector p)
{
	mpVector dp1 = mpVector( 0.0, -2.0,  0.0)-p;
	mpVector dp2 = mpVector( 0.0,  2.0,  0.0)-p;
	mpVector dp3 = mpVector( 2.0,  2.0,  0.0)-p;
	mpVector dp4 = mpVector( 0.0,  0.0,  4.0)-p;
	mpVector dp5 = mpVector(-0.5,  3.1, -1.0)-p;
	mpVector dp6 = mpVector( 0.0,  0.0, -4.0)-p;
	return 1/dp1.Magnitude() + 1/dp2.Magnitude() + 1/dp3.Magnitude() + 1/dp4.Magnitude() + 1/dp5.Magnitude() + 
		1/dp6.Magnitude();
}

void drawCoordinates(float s)
{
	glBegin(GL_LINES);
		glVertex3f(0.0, 0.0, 0.0);		//x coordinate
		glVertex3f(s, 0.0, 0.0);

		glVertex3f(0.0, 0.0, 0.0);		//y coordinate
		glVertex3f(0.0, s, 0.0);
		
		glVertex3f(0.0, 0.0, 0.0);		//z coordinate
		glVertex3f(0.0, 0.0, s);
	glEnd();
}