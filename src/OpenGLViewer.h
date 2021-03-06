#ifndef _OPENGL_VIEWER_
#define _OPENGL_VIEWER_


#include "config.h"
#include "Particles.h"

#include "MarchingCubes/MarchingCubes.h"
#include "LoadShaderUtility.h"

//C/C++ headers
#include <stdlib.h>
#include <stdio.h>

//GlEW and GLFW
#include "glew.h"
#include "wglew.h"
#include "glfw.h"

//GLTools
#include "math3d.h"
#include "GLFrame.h"
#include "GLMatrixStack.h"
#include "GLFrustum.h"
#include "GLGeometryTransform.h"
#include "Vector3.h"
#include "Grid.h"
#include "SolidMesh.h"



#include "AniMesher.h"



TRIANGLE * meshTriangles;
int numOfTriangles;


void openGl_setMesh(TRIANGLE *& tri, int n)
{
	meshTriangles = tri;
	numOfTriangles = n;
}


//----------------------------------------------------------------------------//
// Variables declaration
//----------------------------------------------------------------------------//

//Window properties and GUI-related
int winw, winh;
int lastmousex, lastmousey, lastwheelpos;
float posDx = 0.0f, posDy = 0.0f, zoom = 0.0f, rotDx = 0.0f, rotDy = 0.0f;
int nrOfParticles, nrOfVoxels;
double t0 = 0;
int frames = 0;
int Nx, Ny, Nz;
float h;
bool running = true;
GLfloat edge = 3.0f;
bool step = false, reset = false, showgrid = false, play = false;
bool up_is_down;
bool down_is_down;
SolidMesh testMesh, testMesh2;

bool isokeyleft = false, isokeyright = false;

//----------------------------------------------------------------------------//
// Shaders
//----------------------------------------------------------------------------//

//Shader for drawing particles
GLhandleARB particle_SP;
GLhandleARB particle_VS;
GLhandleARB particle_FS;
GLuint particle_ModelViewMtxLocation;
GLuint particle_ProjectionMtxLocation;
GLuint particle_VertexLocation;
GLuint particle_VelocityLocation;
GLuint particle_dimxLocation;
GLuint particle_dimyLocation;
GLuint particle_dimzLocation;
GLuint particle_hLocation;
GLuint particle_edgeLocation;

//Shader for drawing solid particles
GLhandleARB solid_SP;
GLhandleARB solid_VS;
GLhandleARB solid_FS;
GLuint solid_ModelViewMtxLocation;
GLuint solid_ProjectionMtxLocation;

//Shader for drawing surface mesh
GLhandleARB mesh_SP;
GLhandleARB mesh_VS;
GLhandleARB mesh_FS;
GLuint mesh_ModelViewMtxLocation;
GLuint mesh_ProjectionMtxLocation;
GLuint mesh_NormalMtxLocation;


//Shader for drawing instanced wireframe cubes
//to display grid
GLhandleARB instancedVoxel_SP;
GLhandleARB instancedVoxel_VS;
GLhandleARB instancedVoxel_FS;
GLuint instancedVoxel_ModelViewMtxLocation;
GLuint instancedVoxel_ProjectionMtxLocation;
GLuint instancedVoxel_VertexLocation;
GLuint instancedVoxel_PositionLocation;
GLuint instancedVoxel_FlagLocation;

//VAO's and VBO's
GLuint unitWFCube_VAO;
GLuint unitWFCube_vert_VBO;
GLuint unitWFCube_flags_and_postions_VBO;
GLuint unitWFCube_ind_VBO;
GLuint particle_VAO;
GLuint particle_VBO;

//OpenGLTools
GLFrame cameraFrame;
GLMatrixStack modelViewMatrix;
GLMatrixStack projectionMatrix;
GLFrustum viewFrustum;
GLGeometryTransform transformPipeline;

//----------------------------------------------------------------------------//
// Unit wireframe cube: vertices indices
//----------------------------------------------------------------------------//
GLfloat unitWFCubeVertices[] = 
{ 
	1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,
	1.0f,0.0f,0.0f,  1.0f,1.0f,0.0f,	0.0f,1.0f,0.0f,  0.0f,0.0f,0.0f
};

GLuint unitWFCubeindices[] = { 0,1, 1,2, 2,3, 0,3, 3,4, 4,7, 2,7, 4,5, 5,6, 7,6, 0,5, 1,6};


//----------------------------------------------------------------------------//
// Print opengl error to console
//----------------------------------------------------------------------------//
void getOpenGLError() 
{
	GLenum errCode; 
	const GLubyte *errString;

	if ((errCode = glGetError()) != GL_NO_ERROR) {
		errString = gluErrorString(errCode);
		fprintf (stderr, "OpenGL Error: %s\n", errString);
	}
}


float getFPS(void)
{
	double t = glfwGetTime(); 
	return (float)((double)frames / (t - t0));

}
//----------------------------------------------------------------------------//
// Display fps, winh, winw, zoom and title
//----------------------------------------------------------------------------//
void showFPS(int winw, int winh, float zoom) {

	static char titlestr[200];
	double t, fps;
	// Get current time
	t = glfwGetTime();  // Number of seconds since glfwInit()
	// If one second has passed, or if this is the very first frame
	if( (t - t0) > 1.0 || frames == 0 )
	{
		fps = (double)frames / (t - t0);
		sprintf(titlestr, "FLIP3D, %dx%d pixels, %.1fx zoom, %.1f FPS -> %.1f Mpixels/s",
			winw, winh, zoom, fps, winw*winh*fps*1e-6);
		glfwSetWindowTitle(titlestr);
		t0 = t;
		frames = 0;
	}
	frames ++;
}


//----------------------------------------------------------------------------//
// Cleans up the OpenGL viewer: Destroy window, destroys shaders, fbos, vbos, etc...
//----------------------------------------------------------------------------//
void TerminateViewer()
{
	//Detach and destroy shader's
	glDetachShader(particle_SP,particle_FS);
	glDetachShader(particle_SP,particle_VS);
	glDeleteShader(particle_FS);
	glDeleteShader(particle_VS);
	glDeleteProgram(particle_SP);

	glDetachShader(instancedVoxel_SP,instancedVoxel_FS);
	glDetachShader(instancedVoxel_SP,instancedVoxel_VS);
	glDeleteShader(instancedVoxel_FS);
	glDeleteShader(instancedVoxel_VS);
	glDeleteProgram(instancedVoxel_SP);

	//Delete VBO's and VAO's
	glDeleteVertexArrays(1,&particle_VAO);
	glDeleteBuffers(1,&particle_VBO);

	//Terminate the window
	glfwTerminate();
}

//----------------------------------------------------------------------------//
// Parses and compiles shaders
//----------------------------------------------------------------------------//
void initShaders()
{
	particle_SP = glCreateProgramObjectARB();
	particle_VS = compileShader("particle_VertexShader.glsl", VERTEX_SHADER);
	particle_FS = compileShader("particle_FragmentShader.glsl",FRAGMENT_SHADER);
	glAttachObjectARB(particle_SP,particle_VS);
	glAttachObjectARB(particle_SP,particle_FS);
	link_shaders(particle_SP);
	particle_ModelViewMtxLocation = glGetUniformLocationARB(particle_SP,"modelViewMatrix");
	particle_ProjectionMtxLocation = glGetUniformLocationARB(particle_SP,"projectionMatrix");
	particle_VertexLocation = glGetAttribLocation(particle_SP,"vertex");
	particle_VelocityLocation = glGetAttribLocation(particle_SP,"velocity");
	particle_dimzLocation = glGetUniformLocationARB(particle_SP, "dimz");
	particle_hLocation = glGetUniformLocationARB(particle_SP, "h");
	particle_edgeLocation = glGetUniformLocationARB(particle_SP,"edge");
	getOpenGLError();

	solid_SP = glCreateProgramObjectARB();
	solid_VS = compileShader("solidParticle_vert.glsl", VERTEX_SHADER);
	solid_FS = compileShader("solidParticle_frag.glsl",FRAGMENT_SHADER);
	glAttachObjectARB(solid_SP,solid_VS);
	glAttachObjectARB(solid_SP,solid_FS);
	link_shaders(solid_SP);
	solid_ModelViewMtxLocation = glGetUniformLocationARB(solid_SP,"modelViewMatrix");
	solid_ProjectionMtxLocation = glGetUniformLocationARB(solid_SP,"projectionMatrix");
	getOpenGLError();

	mesh_SP = glCreateProgramObjectARB();
	mesh_VS = compileShader("mesh_vert.glsl", VERTEX_SHADER);
	mesh_FS = compileShader("mesh_frag.glsl",FRAGMENT_SHADER);
	glAttachObjectARB(mesh_SP,mesh_VS);
	glAttachObjectARB(mesh_SP,mesh_FS);
	link_shaders(mesh_SP);
	mesh_ModelViewMtxLocation = glGetUniformLocationARB(mesh_SP,"modelViewMatrix");
	mesh_ProjectionMtxLocation = glGetUniformLocationARB(mesh_SP,"projectionMatrix");
	mesh_NormalMtxLocation = glGetUniformLocationARB(mesh_SP,"normalMatrix");
	getOpenGLError();

	instancedVoxel_SP = glCreateProgramObjectARB();
	instancedVoxel_VS = compileShader("instancedVoxel_VertexShader.glsl", VERTEX_SHADER);
	instancedVoxel_FS = compileShader("instancedVoxel_FragmentShader.glsl",FRAGMENT_SHADER);
	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_VS);
	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_FS);
	link_shaders(instancedVoxel_SP);
	instancedVoxel_ModelViewMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"modelViewMatrix");
	instancedVoxel_ProjectionMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"projectionMatrix");
	instancedVoxel_VertexLocation = glGetAttribLocationARB(instancedVoxel_SP,"vertex");
	instancedVoxel_PositionLocation = glGetAttribLocationARB(instancedVoxel_SP,"position");
	instancedVoxel_FlagLocation = glGetAttribLocationARB(instancedVoxel_SP,"isFluid");
	getOpenGLError();
}

//----------------------------------------------------------------------------//
// Initialize vao, fbo for particles
//----------------------------------------------------------------------------//
void OpenGl_initParticles(void * vertices, void * velocities, size_t size, int nrofparticles_)
{
	nrOfParticles = nrofparticles_;
	glGenVertexArrays(1,&particle_VAO);
	glGenBuffers(1,&particle_VBO);
	glBindVertexArray(particle_VAO);
	glBindBuffer(GL_ARRAY_BUFFER ,particle_VBO);
	glBufferData(GL_ARRAY_BUFFER,2*size,NULL,GL_STREAM_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER,0, size,vertices);
	glBufferSubData(GL_ARRAY_BUFFER,size,size,velocities);
	glVertexAttribPointer(particle_VertexLocation, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(particle_VelocityLocation, 3, GL_FLOAT, GL_FALSE,0,(GLvoid *) size);
	glEnableVertexAttribArray(particle_VertexLocation);
	glEnableVertexAttribArray(particle_VelocityLocation);	
	glBindVertexArray(0);	
}

//----------------------------------------------------------------------------//
// Initialize vao, fbo for wireframe cube
//----------------------------------------------------------------------------//
void OpenGl_initWireframeCube(void * positions, void * flags, int nrOfVoxels_)
{
	nrOfVoxels = nrOfVoxels_;
	size_t positionSize = 3*sizeof(float)*nrOfVoxels;
	size_t flagSize = sizeof(float)*nrOfVoxels;


	glGenVertexArrays(1,&unitWFCube_VAO);
	glGenBuffers(1,&unitWFCube_vert_VBO);
	glGenBuffers(1,&unitWFCube_ind_VBO);
	glGenBuffers(1,&unitWFCube_flags_and_postions_VBO);
	getOpenGLError();

	//Vertices
	glBindVertexArray(unitWFCube_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, unitWFCube_vert_VBO);
	glBufferData(GL_ARRAY_BUFFER,sizeof(unitWFCubeVertices),unitWFCubeVertices,GL_STATIC_DRAW);

	glVertexAttribPointer(instancedVoxel_VertexLocation, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(instancedVoxel_VertexLocation);
	getOpenGLError();

	//Instanced color and position
	glBindBuffer(GL_ARRAY_BUFFER ,unitWFCube_flags_and_postions_VBO);
	//3 floats for position , 1 for the flag
	glBufferData(GL_ARRAY_BUFFER,positionSize + flagSize,0,GL_STREAM_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER,0,positionSize,positions); //Positions
	glBufferSubData(GL_ARRAY_BUFFER,positionSize,flagSize,flags); //Flags
	getOpenGLError();

	glVertexAttribPointer(instancedVoxel_PositionLocation, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(instancedVoxel_FlagLocation, 1, GL_FLOAT, GL_FALSE,0, (GLvoid *)positionSize );
	glEnableVertexAttribArray(instancedVoxel_PositionLocation);
	glEnableVertexAttribArray(instancedVoxel_FlagLocation);
	glVertexAttribDivisor(instancedVoxel_PositionLocation,1);
	glVertexAttribDivisor(instancedVoxel_FlagLocation,1);
	getOpenGLError();

	//Indicies
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, unitWFCube_ind_VBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(unitWFCubeindices), unitWFCubeindices,GL_STATIC_DRAW);
	getOpenGLError();

	glBindVertexArray(0);
}

////////////////////////
// Update particles
////////////////////////
void OpenGl_updateParticles(Particles & p)
{
	nrOfParticles = p.pos.size();
	int size = p.pos.size()*sizeof(vec3f);
	glBindVertexArray(particle_VAO);
	glBindBuffer(GL_ARRAY_BUFFER ,particle_VBO);
	glBufferData(GL_ARRAY_BUFFER,2*size,NULL,GL_STREAM_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER,0, size,&p.pos[0]);
	glBufferSubData(GL_ARRAY_BUFFER,size,size,&p.vel[0]);
	glVertexAttribPointer(particle_VertexLocation, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(particle_VelocityLocation, 3, GL_FLOAT, GL_FALSE,0,(GLvoid *) size);
	glEnableVertexAttribArray(particle_VertexLocation);
	glEnableVertexAttribArray(particle_VelocityLocation);	
	glBindVertexArray(0);
}

//----------------------------------------------------------------------------//
// Update voxels locations and flags
//----------------------------------------------------------------------------//
void OpenGl_updateVoxels(void * positions,Array3f & data, int nrOfVoxels_)
{
	nrOfVoxels = nrOfVoxels_;
	size_t PosSize = 3*sizeof(float)*nrOfVoxels;
	size_t FlagSize = sizeof(float)*nrOfVoxels;
	glBindBuffer(GL_ARRAY_BUFFER ,unitWFCube_flags_and_postions_VBO);
	glBufferSubData(GL_ARRAY_BUFFER,0, PosSize,positions);
	glBufferSubData(GL_ARRAY_BUFFER,PosSize, FlagSize,data.data);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glVertexAttribPointer(instancedVoxel_FlagLocation, 1, GL_FLOAT, GL_FALSE,0, (GLvoid *)PosSize );
	getOpenGLError();
}

//----------------------------------------------------------------------------//
// Responsible for updating OpenGL due to screen resizing
//----------------------------------------------------------------------------//
void Resize()
{
	glfwGetWindowSize(&winw, &winh);
	glViewport(0, 0, winw, winh);

	// Calculate the projection matrix and bind it to shader
	viewFrustum.SetPerspective(45.0f, (GLfloat) winw / (GLfloat) winh, 1.0f, 10000000.0f);	

	projectionMatrix.LoadMatrix(viewFrustum.GetProjectionMatrix());

	transformPipeline.SetMatrixStacks(modelViewMatrix,projectionMatrix);
}

//----------------------------------------------------------------------------//
// Draws particles to screen
//----------------------------------------------------------------------------//
void DrawParticles() 
{
	glUseProgram(particle_SP);
	glBindVertexArray(particle_VAO);
	glUniformMatrix4fv(particle_ModelViewMtxLocation, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
	glUniformMatrix4fv(particle_ProjectionMtxLocation, 1, GL_FALSE, transformPipeline.GetProjectionMatrix());
	glUniform1iv(particle_dimzLocation,1, &Nz);
	glUniform1fv(particle_hLocation,1, &h);
	if(up_is_down)
		edge += 0.025;
	if(down_is_down)
		edge -= 0.025;
	glUniform1fv(particle_edgeLocation,1,&edge);
	glPointSize(4.0);
	glDrawArrays(GL_POINTS,0,nrOfParticles);

	glBindVertexArray(0);
}


//----------------------------------------------------------------------------//
// Draw Voxels to screen
//----------------------------------------------------------------------------//
void DrawVoxels() 
{
	glUseProgram(instancedVoxel_SP);
	glBindVertexArray(unitWFCube_VAO);
	glUniformMatrix4fv(instancedVoxel_ModelViewMtxLocation, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
	glUniformMatrix4fv(instancedVoxel_ProjectionMtxLocation, 1, GL_FALSE, transformPipeline.GetProjectionMatrix());
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDrawElementsInstanced(GL_LINES, 2*12, GL_UNSIGNED_INT, 0, nrOfVoxels);
	glBindVertexArray(0);
	glDisable(GL_BLEND);
}

void DrawSolidVoxels() 
{
	glUseProgram(solid_SP);
	glUniformMatrix4fv(solid_ModelViewMtxLocation, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
	glUniformMatrix4fv(solid_ProjectionMtxLocation, 1, GL_FALSE, transformPipeline.GetProjectionMatrix());
	glPointSize(0.1);

	testMesh.draw();
	//testMesh2.draw();
}

void DrawMesh()
{
	glUseProgram(mesh_SP);
	
	glUniformMatrix3fv(mesh_NormalMtxLocation, 1, GL_FALSE, transformPipeline.GetNormalMatrix());
	glUniformMatrix4fv(mesh_ModelViewMtxLocation, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
	glUniformMatrix4fv(mesh_ProjectionMtxLocation, 1, GL_FALSE, transformPipeline.GetProjectionMatrix());

	//KOM IH�G ATT S�TTA lightpos
	float lightpos0[4]={10.0f, 0.0f, 0.0f, 1.0f};  
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos0);




	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	
	glBegin(GL_TRIANGLES);
	for(int i=0; i < numOfTriangles; i++){
		for(int j=0; j < 3; j++)	
		{
			glNormal3f(meshTriangles[i].norm.x, meshTriangles[i].norm.y, meshTriangles[i].norm.z);
			glVertex3f(meshTriangles[i].p[j].x,meshTriangles[i].p[j].y,meshTriangles[i].p[j].z);
		}
	}
	glEnd();
}

//----------------------------------------------------------------------------//
// Draws all content
//----------------------------------------------------------------------------//

void OpenGl_drawAndUpdate(bool &running)
{
	if(isokeyleft)
	{
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
		{
			if(glfwGetKey(GLFW_KEY_LCTRL) == GLFW_PRESS)
				isovalue -= 0.01;
			else
				--isovalue;
		}
		else
			isovalue -= 0.1;
		std::cout << "iso: " << isovalue << "\n";
	}
	if(isokeyright)
	{
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
		{
			if(glfwGetKey(GLFW_KEY_LCTRL) == GLFW_PRESS)
				isovalue += 0.01;
			else
				++isovalue;

		}
		else
			isovalue += 0.1;
		std::cout << "iso: " << isovalue << "\n";
	}



	//glEnable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	
	running = !glfwGetKey( GLFW_KEY_ESC ) &&
		glfwGetWindowParam( GLFW_OPENED );

	if(!running)
		return;

	showFPS(winw, winh, zoom);
	Resize(); //Update viewport and projection matrix

	//Clear the buffer color and depth
	glClearColor(0.5f,0.5f,0.5f,1.0f);
	//glClearColor(0.0f,0.0f,0.0f,1.0f);
	//glClearColor(1.0f,1.0f,1.0f,1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Disables vsync
	//wglSwapIntervalEXT(0);

	//Camera Matrix
	M3DMatrix44f mCamera;
	cameraFrame.GetCameraMatrix(mCamera);		
	modelViewMatrix.PushMatrix(mCamera);
	modelViewMatrix.Translate(posDx,posDy,zoom);
	modelViewMatrix.Rotate(-rotDx,1.0f,0.0f,0.0f);
	modelViewMatrix.Rotate(-rotDy,0.0f,1.0f,0.0f);
	modelViewMatrix.Translate(-Nx*0.5f*h,-Ny*0.5f*h,-Nz*0.5f*h);

	
#ifdef SOLIDS
	DrawSolidVoxels();
#endif

	DrawParticles();
	DrawMesh();
	

	if(showgrid)
		DrawVoxels();


	modelViewMatrix.PopMatrix();

	glfwSwapBuffers();

}



//----------------------------------------------------------------------------//
// GLFW Keyboard callback
//----------------------------------------------------------------------------//
void GLFWCALL KeyboardFunc( int key, int action )
{
	if(key == 'P')
	{
		if(action == GLFW_PRESS)
			play = !play;
	}

	if(key == 'G')
	{
		if(action == GLFW_PRESS)
			showgrid = !showgrid;

	}

	if(key == 'R' && action == GLFW_PRESS)
		reset = true;

	if(key == 'S' && action == GLFW_PRESS)
		step = true;
	if(key == 'S' && action == GLFW_RELEASE)
		step = false;

	if(key == 'C' && action == GLFW_PRESS) //Cam reset
	{
		rotDy = 0;
		rotDx = 0;
		posDx = 0;
		posDy = 0;
		zoom = 0;

	}

	if(key == GLFW_KEY_ESC)
	{
		running = false;
	}

	if (key == GLFW_KEY_UP && action == GLFW_PRESS)
		up_is_down = true;
	if(key == GLFW_KEY_UP && action == GLFW_RELEASE)
		up_is_down = false;

	if(key == GLFW_KEY_DOWN && action == GLFW_PRESS)
		down_is_down = true;
	if(key == GLFW_KEY_DOWN  && action == GLFW_RELEASE)
		down_is_down = false;


	// CHANGE ISOVALUE IN ANIMESHER
	if(key == GLFW_KEY_RIGHT && action == GLFW_PRESS)
	{
			isokeyright = true;
	}
	if(key == GLFW_KEY_RIGHT && action == GLFW_RELEASE)
	{
			isokeyright = false;
	}

	
	if(key == GLFW_KEY_LEFT && action == GLFW_PRESS)
	{
			isokeyleft = true;
	}
	if(key == GLFW_KEY_LEFT && action == GLFW_RELEASE)
	{
			isokeyleft = false;
	}


	///////////////////////////////////

	/*
	kr = 4.0;	
	kspecial = 0.8;	//(1-kn) -> Hur mycket i procent f�r egenv�rdena skilja sig i axlarna.
	kn = 1.0;
	H = 0.5;
	*/

	if(key == 'B' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			kr += 0.5;
		else
			kr -= 0.5;

		std::cout << "kr = " << kr << "\n";
	}

	if(key == 'M' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			ks += 100;
		else
			ks -= 100;

		std::cout << "ks = " << ks << "\n";
	}

	if(key == 'V' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			H += 0.01;
		else
			H -= 0.01;

		std::cout << "H = " << H << "\n";
	}

	/*if(key == 'N' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			kn += 0.1;
		else
			kn -= 0.1;

		std::cout << "kn = " << H << "\n";
	}

	if(key == 'N' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			kn += 0.1;
		else
			kn -= 0.1;

		std::cout << "kn = " << kn << "\n";
	}*/

	if(key == 'K' && action == GLFW_PRESS)
	{	
		if(glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS)
			kspecial += 0.1;
		else
			kspecial -= 0.1;

		//std::cout << "kspecial = " << kspecial << "\n";
	}


	if(key == 'L' && action == GLFW_PRESS)
	{	

		std::cout << "H = " << H << ": key \'V\'\n";
		//std::cout << "kn = " << kn << ": key \'N\'\n";
		std::cout << "kr = " << kr << ": key \'B\'\n";
		std::cout << "ks = " << ks  << ": key \'M\'\n";
		//std::cout << "kspecial = " << kspecial <<  ": key \'K\'\n";
	}



	

}


//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MouseButtonFunc( int button, int action )
{

}

//----------------------------------------------------------------------------//
// GLFW MouseWheel callback
//----------------------------------------------------------------------------//
void GLFWCALL MouseWheelFunc( int pos )
{
	zoom += (pos - lastwheelpos) *1.2;

	lastwheelpos = pos;
}

//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MousePosFunc( int x, int y )
{
	if(glfwGetKey(GLFW_KEY_LCTRL) == GLFW_PRESS && glfwGetMouseButton(GLFW_MOUSE_BUTTON_1) == GLFW_PRESS)
	{
		rotDy += (lastmousex - x) * 0.5;
		rotDx += (lastmousey - y) * 0.5;
	}
	else if(glfwGetMouseButton(GLFW_MOUSE_BUTTON_1) == GLFW_PRESS)
	{
		posDx -= (lastmousex - x) * 0.05;
		posDy += (lastmousey - y) * 0.05;
	}
	lastmousex = x;
	lastmousey = y;
}


//----------------------------------------------------------------------------//
// Creates and sets up a window
//----------------------------------------------------------------------------//
void OpenGl_initViewer(int width_, int height_, const int dimx_, const int dimy_, const int dimz_, const float gridh)
{
	winw = width_;
	winh = height_;
	Nx = dimx_;
	Ny = dimy_;
	Nz = dimz_;
	h = gridh;

	// Initialize GLFW
	if( !glfwInit() )
	{
		exit( EXIT_FAILURE );
	}

	// Open an OpenGL window using glfw
	if( !glfwOpenWindow( winw,winh, 8,8,8,8,32,0, GLFW_WINDOW ) )
	{
		glfwTerminate();
		exit( EXIT_FAILURE );
	}

	//Init glew!
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	//Setup callbacks
	glfwSetKeyCallback(KeyboardFunc);
	glfwSetMouseButtonCallback(MouseButtonFunc);
	glfwSetMousePosCallback(MousePosFunc);
	glfwSetMouseWheelCallback(MouseWheelFunc);

	glfwGetMousePos(&lastmousex, &lastmousey);
	lastwheelpos = glfwGetMouseWheel();

	initShaders();

	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_POINT_SMOOTH);
	//glCullFace(GL_FRONT);
	

	//Move the camera back 5 units
	cameraFrame.SetOrigin(0.0f,0.0f,10.0);

#ifdef SOLIDS
	//testMesh.init("cube10.obj", vec3f(35,-1,33),vec3f(0.0.0f,20.0f,0.0f) , h);
	testMesh.init("cube10.obj", vec3f(15,-10,15),vec3f(0.0f,1.0f,0.0f) , h);
	//testMesh2.init("cube10.obj", vec3f(0,-1,9), vec3f(0.0f) , h);
#endif

}













#endif