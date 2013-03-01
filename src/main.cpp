//============================================================================
//
//
//============================================================================
#include <iostream>

#include "ObjIO.h"
#include "Mesh.h"
#include "MeshEx.h"
#include "RetileSurface.h"






int main(int argc, char ** argv)
{
	// load mesh
	dk::Mesh *m = dk::io::importFromObjFile( "path/to/file.obj" );

	// generate winged edge datastructure
	MeshEx *mx = new MeshEx(m);

	delete m;

	// perform retiling
	int numVertices = 100; // manually specify the number of vertices (== level of detail)
	RetileSurface retiler;
	RetileSurface::ConstrainedVertexSet cvs; // this is a list of vertices which are fixed in their position (and dont follow repulsion)
	retiler.retile( mx, numVertices, cvs );

	// mx is now retiled and should have numVertices vertices etc.
	m = mx->getMesh();
	dk::io::exportToObjFile( m, "path/to/outputfile.obj" );

	delete m;
	delete mx;

	return 0;
}
