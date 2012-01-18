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
int lastmousex, lastmousey;
float posDx = 0.0f, posDy = 0.0f, zoom = 1.0f, rotDx = 0.0f, rotDy = 0.0f;


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

	glfwGetMousePos(&lastmousex, &lastmousey);

	initShaders();
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
	particle_ModelViewMtxLocation = glGetUniformLocationARB(particle_SP,"ModelViewMatrix");
	particle_ProjectionMtxLocation = glGetUniformLocationARB(particle_SP,"ProjectionMatrix");
	particle_VertexLocation = glGetAttribLocationARB(particle_SP,"Vertex");
	particle_VelocityLocation = glGetAttribLocationARB(particle_SP,"Velocity");

	instancedVoxel_SP = glCreateProgramObjectARB();
	instancedVoxel_VS = compileShader("instancedVoxel_VertexShader.glsl", VERTEX_SHADER);
	instancedVoxel_FS = compileShader("instancedVoxel_FragmentShader.glsl",FRAGMENT_SHADER);
	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_VS);
	glAttachObjectARB(instancedVoxel_SP,instancedVoxel_FS);
	link_shaders(instancedVoxel_SP);
	instancedVoxel_ModelViewMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"ModelViewMatrix");
	instancedVoxel_ProjectionMtxLocation = glGetUniformLocationARB(instancedVoxel_SP,"ProjectionMatrix");
	instancedVoxel_VertexLocation = glGetAttribLocationARB(instancedVoxel_SP,"Vertex");
	instancedVoxel_PositionLocation = glGetAttribLocationARB(instancedVoxel_SP,"Position");
	instancedVoxel_ColorLocation = glGetAttribLocationARB(instancedVoxel_SP,"Color");

}

//----------------------------------------------------------------------------//
// Initialize vao, fbo for particles
//----------------------------------------------------------------------------//
void initParticles(void * vertices, size_t vertSize, void * velocities , size_t velocitySize)
{
	glGenVertexArrays(1,&particle_VAO);
	glGenBuffers(1,&particle_VBO);
	glBindVertexArray(particle_VAO);
	glBindBuffer(GL_ARRAY_BUFFER,particle_VBO);
	glBufferData(GL_ARRAY_BUFFER,vertSize + velocitySize,vertices,GL_STREAM_DRAW); //Will be updated every frame
	glBufferSubData(GL_ARRAY_BUFFER,0,vertSize,vertices);
	glBufferSubData(GL_ARRAY_BUFFER,vertSize,velocitySize,velocities);
	glVertexAttribPointer(particle_VertexLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(particle_VelocityLocation, 3, GL_FLOAT, GL_FALSE, (void *) vertSize);
	glEnableVertexAttribArray(particle_VertexLocation);
	glEnableVertexAttribArray(particle_VelocityLocation);
	glBindVertexArray(0);



}


//----------------------------------------------------------------------------//
// Update particle locations
//----------------------------------------------------------------------------//
void updateParticleLocation(void * locations, size_t size)
{

}

//----------------------------------------------------------------------------//
// Draws all content
//----------------------------------------------------------------------------//
void drawAndUpdate()
{

}

//----------------------------------------------------------------------------//
// Responsible for updating OpenGL due to screen resizing
//----------------------------------------------------------------------------//
void Resize()
{
	glViewport(0, 0, winw, winh);

	// Calculate the projection matrix and bind it to shader
	viewFrustum.SetPerspective(45.0f, (GLfloat) winw / (GLfloat) winh, 1.0f, 100.0f);	

	projectionMatrix.LoadMatrix(viewFrustum.GetProjectionMatrix());

	transformPipeline.SetMatrixStacks(modelViewMatrix,projectionMatrix);
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

//----------------------------------------------------------------------------//
// GLFW MouseButton callback
//----------------------------------------------------------------------------//
void GLFWCALL MousePosFunc( int x, int y )
{
	posDx = (lastmousex - x) / zoom;
	posDy = (lastmousey - y) / zoom;
	rotDx = (lastmousex - x) * 0.5;
	posDy = (lastmousey - y) * 0.5;
	lastmousex = x;
	lastmousey = y;
}














#endif