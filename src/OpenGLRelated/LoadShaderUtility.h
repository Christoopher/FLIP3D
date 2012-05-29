#ifndef	_LOAD_SHADER_UTILITY_
#define _LOAD_SHADER_UTILITY_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <sstream>
#include <string>

#include "glew.h"
#include "wglew.h"
#include "glfw.h"


enum ShaderType
{
	VERTEX_SHADER = 0,
	FRAGMENT_SHADER,
	GEOMETRY_SHADER,
};

//----------------------------------------------------------------------------//
// printError - Signal an error.
// Authour: Stefan Gustavson
//----------------------------------------------------------------------------//
void printError(const char *errtype, const char *errmsg) {
  fprintf(stderr, "%s: %s\n", errtype, errmsg);
}


/*
 * filelength - Determine the number of bytes in a file.
 * This is a lazy hack to avoid calling stat(), but it works.
 */
int filelength(const char *filename) {
// 	FILE *ifp;
// 	int length = 0;
// 
// 	ifp = fopen(filename, "r");
// 	fseek(ifp, 0, SEEK_END);
// 	length = (int)ftell(ifp);
// 	fclose(ifp);
// 	return length;

	struct stat results;
	if (stat(filename, &results) == 0)
		return results.st_size;
	else
	{
		return NULL;
	}
}

//----------------------------------------------------------------------------//
//readShaderFile - read shader source from a file to a string.
//Author: Stefan Gustavson
//----------------------------------------------------------------------------//
char* readShaderFile(const char *filename) 
{
	std::stringstream ss;
	ss << "../src/Shaders/" << filename;
	std::string fullPath;
	ss >> fullPath;

	FILE *file = fopen( fullPath.c_str(), "r");
	if(file == NULL)
	{
		printError("I/O error", "Cannot open shader file!");
		return 0;
	}
	int bytesinfile = filelength(fullPath.c_str());
	
	if(bytesinfile == NULL)
	{
		printError("I/O error", "Cannot get shader file size!");
		return 0;
	}

	unsigned char *buffer = (unsigned char*)malloc(bytesinfile+1);
	int bytesread = fread( buffer, 1, bytesinfile, file);
	buffer[bytesread] = 0; // Terminate the string with 0
	fclose(file);

	return (char *)buffer;
}


//Got this from http://www.lighthouse3d.com/opengl/glsl/index.php?oglinfo
// it prints out shader info (debugging!)
void printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;
	glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);
	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("printShaderInfoLog: %s\n",infoLog);
		free(infoLog);
	}else{
		printf("Shader Info Log: OK\n");
	}
}

//Got this from http://www.lighthouse3d.com/opengl/glsl/index.php?oglinfo
// it prints out shader info (debugging!)
void printProgramInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;
	glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);
	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("printProgramInfoLog: %s\n",infoLog);
		free(infoLog);
	}else{
		printf("Program Info Log: OK\n");
	}
}

//----------------------------------------------------------------------------//
// Links shaders, paramlist
//----------------------------------------------------------------------------//
GLhandleARB create_and_link_shaders(int count, ...)
{
	GLhandleARB programObj = glCreateProgramObjectARB();
	GLint shadersLinked;
	char str[4096]; // For error messages from the GLSL compiler and linker

	va_list ap;
	va_start(ap, count); //Requires the last fixed parameter (to get the address)
	for(int j = 0; j < count; j++)
		glAttachObjectARB(programObj,va_arg(ap, GLhandleARB));
	va_end(ap);

	// Link the program object and print out the info log
	glLinkProgramARB( programObj );
	glGetObjectParameterivARB( programObj, GL_OBJECT_LINK_STATUS_ARB, &shadersLinked );

	//Check if linkage was OK!
	if( shadersLinked == GL_FALSE )
	{
		glGetInfoLogARB( programObj, sizeof(str), NULL, str );
		printError("Program object linking error", str);
	}

	//printProgramInfoLog(programObj);

	return programObj;
}

//----------------------------------------------------------------------------//
// Links shaders, paramlist
//----------------------------------------------------------------------------//
GLhandleARB link_shaders(GLhandleARB & programObj)
{
	GLint shadersLinked;
	char str[4096]; // For error messages from the GLSL compiler and linker

	// Link the program object and print out the info log
	glLinkProgramARB( programObj );
	glGetObjectParameterivARB( programObj, GL_OBJECT_LINK_STATUS_ARB, &shadersLinked );

	//Check if linkage was OK!
	if( shadersLinked == GL_FALSE )
	{
		glGetInfoLogARB( programObj, sizeof(str), NULL, str );
		printError("Program object linking error", str);
	}

	return programObj;
}

//----------------------------------------------------------------------------//
// Compiles a vertex shader and return the ID
//----------------------------------------------------------------------------//
GLhandleARB compileShader(const char *shaderfilename, const ShaderType Type)
{
	
	GLhandleARB shader;
	const char * shaderStrings[1];
	GLint compiled;
	char str[4096]; // For error messages from the GLSL compiler and linker
	char *shaderAssembly;

	//Create the shader object
	if(Type ==  VERTEX_SHADER)
		shader = glCreateShaderObjectARB( GL_VERTEX_SHADER_ARB );
	else if(Type == FRAGMENT_SHADER)
		shader = glCreateShaderObjectARB( GL_FRAGMENT_SHADER_ARB );
	else if(Type == GEOMETRY_SHADER)
		shader = glCreateShaderObjectARB( GL_GEOMETRY_SHADER_ARB );

	//Build the shader
	shaderAssembly = readShaderFile( shaderfilename );
	shaderStrings[0] = shaderAssembly;
	glShaderSourceARB( shader, 1, shaderStrings, NULL );
	glCompileShaderARB(shader);
	free((void *)shaderAssembly);

	//Check if compilation succeeded!
	glGetObjectParameterivARB( shader, GL_OBJECT_COMPILE_STATUS_ARB, 
		&compiled );
	if(compiled == GL_FALSE)
	{
		glGetInfoLogARB( shader, sizeof(str), NULL, str );
		if(Type ==  VERTEX_SHADER)
			printError("Vertex shader compile error", str);
		else if(Type == FRAGMENT_SHADER)
			printError("Fragment shader compile error", str);
		else if(Type == GEOMETRY_SHADER)
			printError("Geometry shader compile error", str);
		
	}

	//printShaderInfoLog(shader);

	return shader;
	
}


#endif