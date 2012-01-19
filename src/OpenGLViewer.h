#ifndef _OPENGL_VIEWER_
#define _OPENGL_VIEWER_

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

//----------------------------------------------------------------------------//
// Variables declaration
//----------------------------------------------------------------------------//

//Window properties and GUI-related
int winw, winh;
int lastmousex, lastmousey, lastwheelpos;
float posDx = 0.0f, posDy = 0.0f, zoom = 0.0f, rotDx = 0.0f, rotDy = 0.0f;
int nrofparticles;

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

//Shader for drawing instanced wireframe cubes
//to display grid
GLhandleARB instancedVoxel_SP;
GLhandleARB instancedVoxel_VS;
GLhandleARB instancedVoxel_FS;
GLuint instancedVoxel_ModelViewMtxLocation;
GLuint instancedVoxel_ProjectionMtxLocation;
GLuint instancedVoxel_VertexLocation;
GLuint instancedVoxel_PositionLocation;
GLuint instancedVoxel_ColorLocation;

//VAO's and VBO's
GLuint unitWFCube_VAO;
GLuint unitWFCube_VBO;
GLuint unitWFCube_ind_VBO;
GLuint particle_VAO;
GLuint particle_VBO;

//OpenGLTools
GLFrame cameraFrame;
GLMatrixStack modelViewMatrix;
GLMatrixStack projectionMatrix;
GLFrustum viewFrustum;
GLGeometryTransform transformPipeline;



void getOpenGLError() 
{
	GLenum errCode; 
	const GLubyte *errString;

	if ((errCode = glGetError()) != GL_NO_ERROR) {
		errString = gluErrorString(errCode);
		fprintf (stderr, "OpenGL Error: %s\n", errString);
	}
}

void showFPS(int winw, int winh, float zoom) {

	static double t0 = 0.0;
	static int frames = 0;
	static char titlestr[200];
	double t, fps;
	// Get current time
	t = glfwGetTime();  // Number of seconds since glfwInit()
	// If one second has passed, or if this is the very first frame
	if( (t - t0) > 1.0 || frames == 0 )
	{
		fps = (double)frames / (t - t0);
		sprintf(titlestr, "%dx%d pixels, %.1fx zoom, %.1f FPS -> %.1f Mpixels/s",
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

// 	glDetachShader(instancedVoxel_SP,instancedVoxel_FS);
// 	glDetachShader(instancedVoxel_SP,instancedVoxel_VS);
// 	glDeleteShader(instancedVoxel_FS);
// 	glDeleteShader(instancedVoxel_VS);
// 	glDeleteProgram(instancedVoxel_SP);

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
	getOpenGLError();

// 	instancedVoxel_SP = glCreateProgramObjectARB();
// 	instancedVoxel_VS = compileShader("instancedVoxel_VertexShader.glsl", VERTEX_SHADER);
// 	instancedVoxel_FS = compileShader("instancedVoxel_FragmentShader.glsl",FRAGMENT_SHADER);
// 	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_VS);
// 	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_FS);
// 	link_shaders(instancedVoxel_SP);
// 	instancedVoxel_ModelViewMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"ModelViewMatrix");
// 	instancedVoxel_ProjectionMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"ProjectionMatrix");
// 	instancedVoxel_VertexLocation = glGetAttribLocationARB(instancedVoxel_SP,"Vertex");
// 	instancedVoxel_PositionLocation = glGetAttribLocationARB(instancedVoxel_SP,"Position");
// 	instancedVoxel_ColorLocation = glGetAttribLocationARB(instancedVoxel_SP,"Color");

}

//----------------------------------------------------------------------------//
// Initialize vao, fbo for particles
//----------------------------------------------------------------------------//
void initParticles(float * vertices, float * velocities, size_t size, int nrofparticles_)
{
	nrofparticles = nrofparticles_;
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
// Update particle locations
//----------------------------------------------------------------------------//
void updateParticleLocation(void * vertices, size_t size)
{
	glBindBuffer(GL_ARRAY_BUFFER ,particle_VBO);
	glBufferSubData(GL_ARRAY_BUFFER,0, size,vertices);
	glBindBuffer(GL_ARRAY_BUFFER,0);
}

//----------------------------------------------------------------------------//
// Responsible for updating OpenGL due to screen resizing
//----------------------------------------------------------------------------//
void Resize()
{
	glfwGetWindowSize(&winw, &winh);
	glViewport(0, 0, winw, winh);

	// Calculate the projection matrix and bind it to shader
	viewFrustum.SetPerspective(45.0f, (GLfloat) winw / (GLfloat) winh, 1.0f, 100.0f);	

	projectionMatrix.LoadMatrix(viewFrustum.GetProjectionMatrix());

	transformPipeline.SetMatrixStacks(modelViewMatrix,projectionMatrix);
}

//----------------------------------------------------------------------------//
// Draws all content
//----------------------------------------------------------------------------//
void drawAndUpdate()
{
	showFPS(winw, winh, zoom);
	Resize();
	glClearColor(0.0,0.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//Disable vsync
	wglSwapIntervalEXT(0);

	//Camera
	M3DMatrix44f mCamera;
	//Cameraframe is located in (0,0,5)

	cameraFrame.GetCameraMatrix(mCamera);
	
	modelViewMatrix.PushMatrix(mCamera);
	modelViewMatrix.Translate(0.0f,0.0f,zoom);
	modelViewMatrix.Rotate(-rotDx,1.0f,0.0f,0.0f);
	modelViewMatrix.Rotate(-rotDy,0.0f,1.0f,0.0f);
	

	glUseProgram(particle_SP);
	glBindVertexArray(particle_VAO);
	glUniformMatrix4fv(particle_ModelViewMtxLocation, 1, GL_FALSE, transformPipeline.GetModelViewMatrix());
	glUniformMatrix4fv(particle_ProjectionMtxLocation, 1, GL_FALSE, transformPipeline.GetProjectionMatrix());
	getOpenGLError();

	glPointSize(1.0);
	glDrawArrays(GL_POINTS,0,nrofparticles);
	getOpenGLError();
	
	modelViewMatrix.PopMatrix();
	glBindVertexArray(0);
	glfwSwapBuffers();

}



//----------------------------------------------------------------------------//
// GLFW Keyboard callback
//----------------------------------------------------------------------------//
void GLFWCALL KeyboardFunc( int key, int action )
{

}


//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MouseButtonFunc( int button, int action )
{

}


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
	if(glfwGetKey(GLFW_KEY_LCTRL) == GLFW_PRESS && glfwGetMouseButton(GLFW_MOUSE_BUTTON_1) == GLFW_PRESS){
		rotDy += (lastmousex - x) * 0.5;
		rotDx += (lastmousey - y) * 0.5;
		
	}
	lastmousex = x;
	lastmousey = y;
}


//----------------------------------------------------------------------------//
// Creates and sets up a window
//----------------------------------------------------------------------------//
void initViewer(int width_, int height_)
{
	winw = width_;
	winh = height_;

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

	//Move the camera back 5 units
	cameraFrame.SetOrigin(0.0f,0.0f,10.0f);
}













#endif