// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This class performs the algorithm which is proposed in the paper
"Retiling Polygonal Surfaces" from G.Turk. The Algorithm uses
the CGAL for computational geometry stuff.

----------------------------------------------------------------------*/
#include "RetileSurface.h"

#include <iostream>
#include <fstream>

// CGAL ---------------------------------------------------------------
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>


typedef CGAL::Filtered_kernel< CGAL::Simple_cartesian<double> > K;


#include <CGAL/Constrained_Delaunay_triangulation_2.h>
typedef CGAL::Triangulation_vertex_base_2<K>                                                               Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K, CGAL::Triangulation_face_base_with_info_2<int, K> > Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                                       TDS;
typedef CGAL::Exact_predicates_tag                                                                       Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>                                Triangulation;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Face_iterator                 Face_iterator;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Face_handle                     Face_handle;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Vertex_handle                 Vertex_handle;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Point                                 Point;



Point get2d_2( math::Vec3f o, math::Vec3f e1,  math::Vec3f e2, math::Vec3f point )
{
	return Point( math::dotProduct( e1, point - o ), math::dotProduct( e2, point - o ) );
}

math::Vec3f get3d_2( math::Vec3f o, math::Vec3f e1,  math::Vec3f e2, Point point )
{
	return o + e1*point.x() + e2*point.y();
}

//
// constructor
//
RetileSurface::RetileSurface()
{
	// initialize members
	m_radius        = 0.0f;
	m_damp         = 0.15f; // damping factor of the repulsion forces
	m_iterationCount = 100;

	mesh               = 0;
}

//
// destructor
//
RetileSurface::~RetileSurface()
{
	clearVertices();
}

//
//
//
void RetileSurface::clearVertices()
{
	for(std::vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
		delete *it;
	vertices.clear();
}

//
// generates a specified number of vertices which will be randomly
// distributed (not evenly) over the mesh
//
void RetileSurface::generateVertices( size_t number )
{
	std::map< float, MeshEx::Triangle* > areas; // used for randomly distributing vertices over the mesh

	float area_sum = 0.0f;
	for( size_t i=0; i<mesh->m_triangles.size(); ++i )
	{
		float area = math::area( mesh->m_triangles[i]->v0->position, mesh->m_triangles[i]->v1->position, mesh->m_triangles[i]->v2->position );
		area_sum += area;
		areas[area_sum] = mesh->m_triangles[i];
	}

	m_radius = 2.0f*sqrt( area_sum / float(number)  );

	// distribute vertex candidates ------------------------------------------------------------
	// scatter randomly distributed vertices across the surface of the orginigal mesh
	for( size_t i=0; i<number; ++i )
	{
		// randomly select triangle
		float s = math::g_randomNumber() * area_sum;
		MeshEx::Triangle *selected = 0;
		float area_i_minus_1 = 0.0f;
		float area_i = 0.0f;

	

		for( std::map<float, MeshEx::Triangle *>::iterator it = areas.begin(); it != areas.end(); area_i_minus_1 = it->first, ++it )
			if( s < it->first )
			{
				selected = it->second;
				area_i = it->first;
				break;
			}
		// TODO: optimize searching
		//std::map<float, MeshEx::Triangle *>::iterator it = std::lower_bound( areas.begin(), areas.end(), s );

		if( !selected )
			printf( "error: no triangle randomly selected!\n" );

		// TODO: make sure these points are not too close to an vertex or an edge (?)

		// create random barycentric coords
		float a = math::mapValueTo0_1( area_i_minus_1, area_i, s );
		float b = (1.0f - a) * math::g_randomNumber();
		float c = 1.0f - a - b;

		float x = a * selected->v0->position.x + b * selected->v1->position.x + c * selected->v2->position.x;
		float y = a * selected->v0->position.y + b * selected->v1->position.y + c * selected->v2->position.y;
		float z = a * selected->v0->position.z + b * selected->v1->position.z + c * selected->v2->position.z;

		math::Vec3f vertexPosition( x, y, z );

		// make sure that the created point lies within the polygon and not outside or on the boundary
		math::Vec3f                     planeNormal; // normal of the plane on which we move the vertex
		math::Vec3f                           base1; // base vectors which build the 2d coordinate system within |R³
		math::Vec3f                           base2;
		float                         planeDistance; // distance of the plane on which we move the vertex to the origin

		planeNormal = selected->normal;
		planeDistance = -math::dotProduct( vertexPosition, planeNormal );
		base1 = math::normalize( selected->u );
		base2 = math::crossProduct( planeNormal, base1 );


		Point t_v1 = get2d_2( vertexPosition, base1, base2, selected->v[0]->position );
		Point t_v2 = get2d_2( vertexPosition, base1, base2, selected->v[1]->position );
		Point t_v3 = get2d_2( vertexPosition, base1, base2, selected->v[2]->position );

		Point vert = get2d_2( vertexPosition, base1, base2, vertexPosition );

		std::vector<Point> ps;

		ps.push_back( t_v1 );
		ps.push_back( t_v2 );
		ps.push_back( t_v3 );


		switch( CGAL::bounded_side_2( ps.begin(), ps.end(), vert, K()) )
		{
		case CGAL::ON_BOUNDED_SIDE :
		  break;
		case CGAL::ON_BOUNDARY:
		case CGAL::ON_UNBOUNDED_SIDE:
			--i;
			continue;
		}

		vertices.push_back( new Vertex( x, y, z ) );
		// store the triangle this vertex belongs to
		vertices.back()->t = selected;

		vertices.back()->path.push_back( vertices.back()->position );
		vertices.back()->index = i;
	}
}

//
//
//
void RetileSurface::moveVertexOnMesh( Vertex *vertex, const math::Vec3f &direction )
{
	math::Vec3f  targetPosition = vertex->position + direction; // the current target position
	MeshEx::Edge                                 *lastEdge = 0; // edge which was last crossed


	int count = 0; // tells how many times the point was rotated around an edge
	do
	{
		MeshEx::Edge       *edge = 0;
		math::Vec3f     intersection;
		bool     pointOnLine = false; // tells whether the point lies exactly on a edge

		math::Vec3f      planeNormal; // normal of the plane on which we move the vertex
		math::Vec3f       planePoint; // a point on the plane
		math::Vec3f            base1; // base vectors which build the 2d coordinate system within |R³
		math::Vec3f            base2;
		float          planeDistance; // distance of the plane on which we move the vertex to the origin

		planeNormal = vertex->t->normal;

		planePoint = vertex->t->v0->position;
		planeDistance = -math::dotProduct( planePoint, planeNormal );
		base1 = vertex->t->u;
		base2 = vertex->t->v_v;


		if(0)
		//if( count == 0 )
		//if( vertex->index == 257  )
		{

			Point t_v1 = get2d_2( planePoint, base1, base2, vertex->t->v[0]->position );
			Point t_v2 = get2d_2( planePoint, base1, base2, vertex->t->v[1]->position );
			Point t_v3 = get2d_2( planePoint, base1, base2, vertex->t->v[2]->position );

			Point vert = get2d_2( planePoint, base1, base2, vertex->position );

			std::vector<Point> ps;

			ps.push_back( t_v1 );
			ps.push_back( t_v2 );
			ps.push_back( t_v3 );


			switch( CGAL::bounded_side_2( ps.begin(), ps.end(), vert, K()) )
			{
			case CGAL::ON_BOUNDED_SIDE :
				//printf( "count %i: is inside the polygon.\n", count);
			  break;
			case CGAL::ON_BOUNDARY:
			  //printf( "count %i: is on the polygon boundary. (%i)\n", count ,vertex->index);
			  break;
			case CGAL::ON_UNBOUNDED_SIDE:
			  //printf( "count %i: is outside the polygon. (%i)\n", count ,vertex->index);
			  break;
			}
		}

		// intersect the vector from currentPosition to targetPosition with the edges
		// to see whether the points moves out of the polygon
		for( size_t j=0; j<3; ++j )
		{
			// if the vertex has already travelled over an edge, than we wont
			// test against that edge
			if( vertex->t->e[j] != lastEdge )
			{
				// project edge-vertices, current vertex and target position into the 2d-plane of the triangle
				Point e_v1 = get2d_2( planePoint, base1, base2, vertex->t->e[j]->v1->position );
				Point e_v2 = get2d_2( planePoint, base1, base2, vertex->t->e[j]->v2->position );
				Point vert = get2d_2( planePoint, base1, base2, vertex->position );
				Point targ = get2d_2( planePoint, base1, base2, targetPosition );

				// test if point lies on that edge
				if( CGAL::orientation( e_v1, e_v2, vert ) == CGAL::COLLINEAR )
				{
					printf( "point lies on the line! (%i)  count: %i\n", vertex->index, count );
					/*
					pointOnLine = true;
					edge = lastEdge = vertex->t->e[j];
					intersection = vertex->position;
					vertex->t->e[j]->tag = true;
					
					break;
					*/
					// this is a bit nasty, because the point lies on a vertex of the triangle, touching
					// multiple edges

					// currently we handle this by replacing the vertex to the center of the triangle
					vertex->position = vertex->t->center;
					//points.push_back( std::make_pair( vertex->position, math::Color::Yellow() ) );
					return;
				}else
				{
					// do intersection test with cgal
					CGAL::Segment_2<K> line_1( e_v1, e_v2 );
					CGAL::Segment_2<K> line_2( vert, targ );
					if( CGAL::do_intersect( line_1, line_2 ) )
					{
						Point tempIntersection;
						if( CGAL::assign(tempIntersection, CGAL::intersection(line_1, line_2)) )
						{
							// point
							intersection = get3d_2( planePoint, base1, base2, tempIntersection );
						}else
						{
							// segment!
							printf( "error : intersection during vertex movement is a line\n" );
						}

						edge = lastEdge = vertex->t->e[j];
						break;
					}
					// project the target position onto the plane
					//targetPosition = get3d_2( vertex->position, base1, base2, targ ); // doesnt help
				}
			}
		}

		// if there was no intersection...
		if( !edge )
		{
			// ...then we are done - the targetposition doesnt move out of the triangle
			vertex->position = targetPosition;
			return;
		}

		// get the adjacent Triangle
		MeshEx::Triangle *adjTri = edge->getOtherTriangle( vertex->t );

		// if the edge is a boundary edge (which is the case when there is no
		// triangle on the other side of that edge) then move the point onto
		// the boundary edge and quit
		if( !adjTri )
		{
			vertex->position = intersection;
			return;
		}


		// the point leaves the triangle over an edge to another triangle
		// -> rotate the point around the edge which is crossed to keep the
		// point planar to the adjacent triangle
		// -> this procedure is similar to the one in the mapPointToPlane function
		
		// axis, the point will be rotated about
		math::Vec3f axis = math::normalize(edge->v2->position - edge->v1->position);

		float dp = math::dotProduct( vertex->t->normal, adjTri->normal );


		if( dp < -1.0f )
			dp = -1.0f;
		if( dp > 1.0f )
			dp = 1.0f;

		// rotation angle around edge
		float dihedralAngle = acosf( dp );

		float check = math::dotProduct( math::crossProduct(axis, vertex->t->normal) , math::normalize(vertex->position - edge->v1->position) );
		float check2 = math::dotProduct( vertex->t->normal , math::normalize(adjTri->getOtherVertex(edge)->position - edge->v1->position) );


		if( check < 0.0f )
			dihedralAngle = -dihedralAngle;

		if( check2 > 0.0f )
			dihedralAngle = -dihedralAngle;




		math::Matrix44f trans = math::Matrix44f::TranslationMatrix( -edge->v1->position );
		math::Matrix44f rot = math::Matrix44f::RotationMatrix( axis, -dihedralAngle );
		math::Matrix44f itrans = math::Matrix44f::TranslationMatrix( edge->v1->position );

		// move the current position to the intersection point
		vertex->position = intersection;
		// rotate the target position onto the plane of the adjacent triangle
		targetPosition = math::transform( targetPosition, trans*rot*itrans );
		// set the adjacent triangle to the current triangle
		vertex->t = adjTri;

		// check if the target position falls directly on an edge of the new polygon

		// check if the vertex is still within the polygon


		// we track the path of the vertex for debugging purposes
		//vertex->path.push_back( vertex->position );
	}while( ++count < 100 );

	//vertex->position = targetPosition;
}


//
// maps the given point into the plane defined by t
//
// optional an edge may be specified around which the point q should be rotated
// if such an edge is not given, the edge is guesstimated
//
math::Vec3f RetileSurface::mapPointToPlane( const math::Vec3f &q, MeshEx::Triangle *t, MeshEx::Edge *hingeEdge )
{
	MeshEx::Edge *edge = hingeEdge;

	// if no edge is defined
	if( !edge )
	{
		// guesstimate the edge of the direction the point q lies in respect to the triangle t
		// we do this by checking q with each plane defined by the edge-normals
		float   distances[3];

		// for each vertex-edge
		for( size_t i=0; i<3; ++i )
		{
			MeshEx::Vertex *v1 = t->v[i];
			MeshEx::Vertex *v2 = t->v[(i+1)%3];

			math::Vec3f planeNormal = math::normalize( math::crossProduct( t->normal, v2->position - v1->position ) );
			float planeDistance = -math::dotProduct( v1->position, planeNormal );
		
			distances[i] = math::distancePointPlane( q, planeNormal, planeDistance );

			// if the distance is below 0, then the point lies outside the triangle on that side of the edge
		}

		// we take the one with the smallest value
		if( (distances[0] <= distances[1])&&(distances[0] <= distances[2]) )
			edge = t->getEdge( t->v0, t->v1 );
		else
			if( (distances[1] <= distances[0])&&(distances[1] <= distances[2]) )
				edge = t->getEdge( t->v1, t->v2 );
			else
				if( (distances[2] <= distances[0])&&(distances[2] <= distances[1]) )
					edge = t->getEdge( t->v2, t->v0 );
	}

	// axis, the point will be rotated about
	math::Vec3f axis = math::normalize(edge->v2->position - edge->v1->position);

	// normal of the adjacent triangle
	math::Vec3f adjNormal;

	// if there is an adjacent triangle, we can get the normal from that
	MeshEx::Triangle *adjTri = edge->getOtherTriangle( t );
	if( adjTri )
		adjNormal = adjTri->normal;
	else
		// else we have to construct the normal
		adjNormal = math::crossProduct( axis, math::normalize( q - edge->v1->position ) );

	// rotation angle around edge
	float dihedralAngle = acosf( math::dotProduct( t->normal, adjNormal ) );

	float check = math::dotProduct( math::crossProduct(axis, t->normal) , math::normalize( t->center - edge->v1->position) );
	float check2 = math::dotProduct( t->normal , math::normalize( q - edge->v1->position) );

	if( check < 0.0f )
		dihedralAngle = -dihedralAngle;

	if( check2 > 0.0f )
		dihedralAngle = -dihedralAngle;


	math::Matrix44f trans = math::Matrix44f::TranslationMatrix( -edge->v1->position );
	math::Matrix44f rot = math::Matrix44f::RotationMatrix( axis, dihedralAngle );
	math::Matrix44f itrans = math::Matrix44f::TranslationMatrix( edge->v1->position );

	return math::transform( q, trans*rot*itrans );
}


//
// Retriangulates each of the new vertices into the triangle on which they rest
//
void RetileSurface::doMutualTesselation()
{
	printf( "mutual retriangulation...\n");

	// 1. do-retriangulation including the vertex-candidates

	// build hash
	std::map<MeshEx::Triangle *, std::vector<math::Vec3f> > tris;
	std::map<MeshEx::Triangle *, std::vector<Vertex *> > tris2;

	for( size_t i=0; i<vertices.size(); ++i )
	{
		Vertex *v = vertices[i];
		tris[ v->t ].push_back( v->position );
		tris2[ v->t ].push_back( v );
	}

	// retriangulate triangles
	for( std::map<MeshEx::Triangle *, std::vector<math::Vec3f> >::iterator it = tris.begin(); it != tris.end(); ++it )
		mesh->reTriangulate( it->first, it->second );
}


//
// Moves the new vertex over the meshsurface depending on repulsion forces excerted from nearby
// neighbours.
//
void RetileSurface::performIteration( float radius, float damp )
{
	// for k iterations
	//  for each point P on surface
	//    determine nearby points to P
	//    map these nearby points onto the plane of the triangle P lies on
	//    compute repulsion forces of theses points exerted on P
	//  for each point of P on surface
	//    compute new position of P on surface

	// for each point on surface - compute superposition of all neighbouring points
	for( size_t i=0; i<vertices.size(); ++i)
	{
		Vertex *p = vertices[i];

		// reset force
		p->force = math::Vec3f( 0.0f, 0.0f, 0.0f );

		// determine nearby points
		std::vector<Vertex *>   neighbours;
		for( size_t j=0; j<vertices.size(); ++j)
			if( i != j )
				if( math::distance( p->position, vertices[j]->position ) < radius )
					neighbours.push_back( vertices[j] );


		std::vector<math::Vec3f> projectedNeighbours;

		// for each nearby point - compute repulsing force
		for( size_t j=0; j < neighbours.size(); ++j )
		{
			Vertex *q = neighbours[j];
			math::Vec3f    qProjected;
			MeshEx::Edge      *sharedEdge = 0;

			// if q lies on the same triangle as p
			if( p->t == q->t )
			{
				// point q lies on the same triangle as p
				qProjected = q->position;
			}else
			// if the triangle of q is adjacent to the triangle of p
			if( sharedEdge = p->t->getSharedEdge( q->t ) )
			{
				// rotate the vertex position around the edge between
				// those 2 triangles
				qProjected = mapPointToPlane( q->position, p->t, sharedEdge );
			}else
			{
				// rotate q around the edge which points into the direction of q->position
				qProjected = mapPointToPlane( q->position, p->t );

				// project q onto the plane of the triangle in which p lies
				qProjected = math::projectPointOnPlane( p->t->normal, -math::dotProduct(p->t->normal, p->t->center), qProjected );
			}

			// compute surface distance approximation to the point
			float d = math::distance( p->position, qProjected );

			// if the projected point q is still within range
			if( d<radius )
			{
				// compute repulsing force for the current neighbour
				// depending on approximated surface distance and the
				// vector from projected q to p
				p->force += (radius-d)*math::normalize( p->position - qProjected );

				projectedNeighbours.push_back( qProjected );
			}
		}

		// check if all neighbours lie in the plane of the triangle of p
		for( size_t n = 0; n<projectedNeighbours.size(); ++n )
		{
			float dis = math::distancePointPlane( projectedNeighbours[n], p->t->normal, -math::dotProduct( p->t->normal, p->t->center ) );

			if( fabs(dis) > 0.0001f )
			{
				printf( "error neighbour points dont lie on the plane of t\n" );
				p->tag = true;
			}
		}
	}

	// for each point on surface - move point depending on the forces excerted from neighbours
	for( size_t i=0; i<vertices.size(); ++i)
	{
		Vertex *p = vertices[i];

		// we dont move the vertex if it is not fixed
		if( p->fixed )
			continue;

		// move the vertex on the mesh for the specified distance
		moveVertexOnMesh( p, damp*p->force );

		// we project the vertex back to the plane of the triangle to counter instablities
		p->position = math::projectPointOnPlane( p->t->normal, -math::dotProduct( p->t->center, p->t->normal ), p->position );
	}
}


//
// performs the algorithm on the given mesh with a set of contrained vertices
//
void RetileSurface::retile( MeshEx *_mesh, size_t vertexCount, ConstrainedVertexSet &cvs )
{
	mesh = _mesh;

	// clear the list of vertices
	clearVertices();

	// store current mesh vertices
	std::vector<MeshEx::Vertex *>  oldVertices;
	for( size_t i=0; i<mesh->m_vertices.size(); ++i )
		oldVertices.push_back( mesh->m_vertices[i] );

	// create random vertices all over the mesh -----------------------------------
	printf( "generating vertices....\n" );
	generateVertices( vertexCount - cvs.size() );

	// add constrained vertex
	for( ConstrainedVertexSet::iterator it = cvs.begin(); it != cvs.end(); ++it )
	{
		vertices.push_back( new Vertex( it->first.x, it->first.y, it->first.z ) );
		// store the triangle this vertex belongs to
		vertices.back()->t = it->second;

		// fix this vertex so that it wont be moved
		vertices.back()->fixed = true;
	}

	// perform relaxation ---------------------------------------------------------
	printf( "performing relaxation...\n" );
	for( size_t i=0; i<m_iterationCount; ++i )
		performIteration( m_radius, m_damp );

	// mutual tesselation --------------------------------------------------------
	printf( "performing mutual tesselation...\n" );
	doMutualTesselation();


	// 2. remove old vertices
	std::vector<MeshEx::Vertex *> retainedVertices; // vertices which could not be removed

	printf( "retriangulating...\n" );
	for( size_t i=0; i<oldVertices.size(); ++i )
	{
		bool removed = mesh->removeVertexAndReTriangulateNeighbourhood( oldVertices[i] );

		// if vertex could not be removed
		if( !removed )
			// try again later...
			retainedVertices.push_back( oldVertices[i] );
	}

	{
		int count = 0; // this variable counts the number of tries

		// as long as there are any retained vertices and
		// we have not hit the maximum trial count yet:
		while( !retainedVertices.empty() && (count < 5) )
		{
			// get the first entry
			std::vector<MeshEx::Vertex *>::iterator it = retainedVertices.begin();

			// as long as we have not reached the list of all retained vertices...
			while( it != retainedVertices.end() )
			{
				// try to remove the current vertex
				if( mesh->removeVertexAndReTriangulateNeighbourhood( *it ) )
					// if it was successfull, then remove the vertex
					// from the list of retained vertices
					retainedVertices.erase( it ); // it now points to the next element within the vector
				else
					// vertex could not be removed -> go ahead
					++it;
			};

			++count;
		};
	}

	printf( "done (%i vertices retained from removal)\n", retainedVertices.size() );


	// postprocess the mesh
	for( std::vector<MeshEx::Triangle *>::iterator it = mesh->m_triangles.begin(); it != mesh->m_triangles.end(); ++it )
		(*it)->computeNormal();
	for( std::vector<MeshEx::Vertex *>::iterator it = mesh->m_vertices.begin(); it != mesh->m_vertices.end(); ++it )
		(*it)->computeNormal();
}




//
// this is a utility routine for debug and visualization purposes
// writes the marker Particles out to a binary file
//
void RetileSurface::writeVertices( const std::vector<RetileSurface::Vertex *> &vertices, const char *filenameFormat, ... )
{
	char filename[2048];

	// assemble filename from given format and arguments
	va_list argumentList;
	va_start( argumentList, filenameFormat );
	vsprintf_s( filename, 2048, filenameFormat, argumentList);
	va_end(argumentList);

	// createopen the file
	std::ofstream file;
	file.open( filename, std::ios::out | std::ios::binary );


	// file creation not successfull?
	if( !file )
		// quit
		return;

	// write into the file

	// number of particles
	int vertexNum = (int)vertices.size();
	file.write( (char *) &vertexNum, sizeof(int) );

	// now write the position of each particle
	for( unsigned int i=0; i<(unsigned int)vertexNum; ++i )
	{
		float x = vertices[i]->position.x;
		float y = vertices[i]->position.y;
		float z = vertices[i]->position.z;

		file.write( (char *) &x, sizeof(float) );
		file.write( (char *) &y, sizeof(float) );
		file.write( (char *) &z, sizeof(float) );
	}

	// done
	file.close();
}