// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

A very straightforward trianglemesh class

----------------------------------------------------------------------*/
#include "Mesh.h"
#include <map>
#include <algorithm>

namespace dk
{
	//
	// constructor
	//
	Mesh::Mesh()
	{
	}

	//
	// destructor
	//
	Mesh::~Mesh()
	{
		for( std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it )
			delete *it;
		for( std::vector<Triangle*>::iterator it = triangles.begin(); it != triangles.end(); ++it )
			delete *it;

		vertices.clear();
		triangles.clear();
	}

	//
	// constructor
	//
	Mesh::Mesh( const std::vector<math::Vec3f> &_vertices, const std::vector<int> &indicees )
	{
		for( std::vector<math::Vec3f>::const_iterator it=_vertices.begin(); it != _vertices.end(); ++it )
		{
			// create a vertex and
			// add it to the list of vertices
			vertices.push_back( new Vertex( *it ) );
		}

		std::vector<Vertex *> temp;

		temp.resize( vertices.size() );

		// resolve the indicees to vertex references
		int currentIndex = 0;
		for( std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it, ++currentIndex )
			temp[currentIndex] = *it;


		// add triangles
		for( std::vector<int>::const_iterator it=indicees.begin(); it != indicees.end(); it += 3 )
		{
			int i1 = *( it+0 );
			int i2 = *( it+1 );
			int i3 = *( it+2 );

			// createadd triangle
			triangles.push_back( new Triangle( temp[i1], temp[i2], temp[i3] ) );
		}
	}


	//
	//
	//
	void Mesh::computeVertexIndicees( void )
	{
		// resolve the vertex references to indicees
		int currentIndex = 0;
		for( std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it, ++currentIndex )
			(*it)->index = currentIndex;
	}


	//
	//
	//
	void Mesh::computeNormals( void )
	{
		// reset vertex normals
		for( std::vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it )
			(*it)->normal = math::Vec3f();

		// compute face normals
		for( std::vector<Triangle *>::iterator it = triangles.begin(); it != triangles.end(); ++it )
		{
			//compute face normal
			math::Vec3f faceNormal = math::normalize( math::crossProduct( (*it)->v1->position - (*it)->v0->position, (*it)->v2->position - (*it)->v0->position ) );

			// add face normal to all vertices which are touched by the triangle
			(*it)->v0->normal += faceNormal;
			(*it)->v1->normal += faceNormal;
			(*it)->v2->normal += faceNormal;
		}

		// normalize vertex normals
		for( std::vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it )
			(*it)->normal.normalize();
		
	}


	void Mesh::assertManifold()
	{
		std::map<Vertex*, std::vector<Triangle *> > triRing;

		// for each triangle
		for( std::vector<Triangle *>::iterator tit = triangles.begin(); tit != triangles.end(); ++tit )
		{
			Triangle *t = *tit;

			for( size_t i=0; i<3; ++i )
				triRing[ t->v[i] ].push_back( t );
		}

		for( std::map<Vertex*, std::vector<Triangle *> >::iterator it = triRing.begin(); it != triRing.end(); ++it )
		{
			Vertex *v = it->first;

			if( it->second.size() < 3 )
			{
				printf( "degenerated vertex detected (referenced by less then 3 triangles) -> removed\n" );

				// remove all triangles referencing v
				for( std::vector<Triangle *>::iterator tit = it->second.begin(); tit != it->second.end(); ++tit )
				{
					Triangle *t = *tit;

					// remove triangle from triangleRings from all reference vertices
					for( size_t i =0;i<3; ++i )
						if( t->v[i] != v )
						{
							triRing[t->v[i]].erase( std::find( triRing[t->v[i]].begin(), triRing[t->v[i]].end(), t ) );
							// TODO: check if new degenerate vertices are created
						}

					// remove the triangle
					removeTriangle( t );
				}

				// remove v itsself
				removeVertex( v );
			}
		}


	}

	void Mesh::removeVertex( Vertex *v )
	{
		std::vector<Vertex *>::iterator it = std::find( vertices.begin(), vertices.end(), v );

		if( it != vertices.end() )
			vertices.erase( it );

		delete v;
	}

	void Mesh::removeTriangle( Triangle *t )
	{
		std::vector<Triangle *>::iterator it = std::find( triangles.begin(), triangles.end(), t );

		if( it != triangles.end() )
			triangles.erase( it );

		delete t;
	}
}