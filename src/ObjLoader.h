#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include "Vector3.h"

#define TOKEN_VERTEX_POS "v"
#define TOKEN_VERTEX_NOR "vn"
#define TOKEN_VERTEX_TEX "vt"
#define TOKEN_FACE "f"

struct ObjMeshVertex{
	vec3f pos;
	vec3f texcoord;
	vec3f normal;
};

/* This is a triangle, that we can render */
struct ObjMeshFace{
	ObjMeshVertex vertices[3];
};

/* This contains a list of triangles */
struct ObjMesh{
	std::vector<ObjMeshFace> faces;
};



/* Internal structure */
struct _ObjMeshFaceIndex{
	int pos_index[3];
};

/* Call this function to load a model, only loads triangulated meshes */
ObjMesh LoadObjMesh(std::string filename){
	ObjMesh myMesh;

	std::vector<vec3f>           positions;
	std::vector<_ObjMeshFaceIndex>  faces;
	/**
	 * Load file, parse it
	 * Lines beginning with:
	 * '#'  are comments can be ignored
	 * 'v'  are vertices positions (3 floats that can be positive or negative)
	 * 'vt' are vertices texcoords (2 floats that can be positive or negative)
	 * 'vn' are vertices normals   (3 floats that can be positive or negative)
	 * 'f'  are faces, 3 values that contain 3 values which are separated by / and <space>
	 */

	std::ifstream filestream;
	filestream.open(filename.c_str());

	if (filestream.is_open() != true)
	{
		std::cout << "could not fnid file\n";
		exit(0);
	}

	std::string line_stream;	// No longer depending on char arrays thanks to: Dale Weiler
	while(std::getline(filestream, line_stream)){
		std::stringstream str_stream(line_stream);
		std::string type_str;
		str_stream >> type_str;
		if(type_str == TOKEN_VERTEX_POS){
			vec3f pos;
			str_stream >> pos[0] >> pos[1] >> pos[2];
			positions.push_back(pos);
		}
		else if(type_str == TOKEN_FACE){
			_ObjMeshFaceIndex face_index;
			char interupt;
			for(int i = 0; i < 3; ++i){
				str_stream >> face_index.pos_index[i];
			}
			faces.push_back(face_index);
		}
	}
	// Explicit closing of the file
	filestream.close();

	for(size_t i = 0; i < faces.size(); ++i){
		ObjMeshFace face;
		for(size_t j = 0; j < 3; ++j){
			face.vertices[j].pos = positions[faces[i].pos_index[j] - 1];
		}
		myMesh.faces.push_back(face);
	}

	return myMesh;
}



#endif
