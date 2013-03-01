/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#include "MeshEx.h"
#include <algorithm>



#include <list>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polygon_2.h>

//typedef CGAL::Cartesian<double> K;
//struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Filtered_kernel< CGAL::Simple_cartesian<double> > K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K, std::list<Point> > Polygon_2;
//typedef CGAL::Polygon_2<K, std::list<Point> >::Vertex_iterator VertexIterator;
//typedef CGAL::Polygon_2<K, std::list<Point> >::Edge_const_iterator EdgeIterator; 

size_t g_iterations = 0;

//
// constructor which takes a mesh as input
//
MeshEx::MeshEx( Mesh *mesh )
{
	stop = false;

	mesh->computeVertexIndicees();
	mesh->computeNormals();

	// convert the mesh into the local problem domain
	// create a node for each vertex
	for( std::vector<Mesh::Vertex *>::iterator it=mesh->vertices.begin(); it != mesh->vertices.end(); ++it )
	{
		Mesh::Vertex *vert = *it;
		Vertex *n = createVertex(vert->position);
		n->normal = vert->normal;
	}

	/*
	// create a element for each triangle
	for( std::vector<Mesh::Triangle *>::iterator it=mesh->triangles.begin(); it != mesh->triangles.end(); ++it )
	{
		Mesh::Triangle *tri = *it;
		Triangle *e = createTriangle( tri->v[0]->index,tri->v[1]->index, tri->v[2]->index );
	}
	*/

	std::map< Vertex *, std::vector<Edge *> > edgeHash; // stores all incoming and outgoing edges for each vertex

	//int tempCount = 0;
	for( std::vector<Mesh::Triangle *>::iterator it=mesh->triangles.begin(); it != mesh->triangles.end(); ++it )
	{
		Mesh::Triangle *tri = *it;

		//if( tempCount++>1000)
		//	break;

		// get vertices of the triangle
		Vertex *v0 = m_vertices[tri->v[0]->index];
		Vertex *v1 = m_vertices[tri->v[1]->index];
		Vertex *v2 = m_vertices[tri->v[2]->index];

		// edges which will build the mesh
		Edge *e0, *e1, *e2;

		e0 = e1 = e2 = 0;

		// look for edge v0->v1 and v2->v0
		std::vector<Edge *> &v0_edges = edgeHash[v0];

		for( std::vector<Edge *>::iterator eit = v0_edges.begin(); eit != v0_edges.end(); ++eit )
		{
			Edge *e = *eit;
			if( e->getOtherVertex( v0 ) == v1 )
			{
				e0 = e;
			}else
			if( e->getOtherVertex( v0 ) == v2 )
			{
				e2 = e;
			}
		}

		// edge not found?
		if( !e0 )
		{
			// create it
			e0 = createEdge( v0, v1 );
			// register edge with the hash
			edgeHash[v0].push_back( e0 );
			edgeHash[v1].push_back( e0 );
		}

		// edge not found?
		if( !e2 )
		{
			// create it
			e2 = createEdge( v2, v0 );
			// register edge with the hash
			edgeHash[v2].push_back( e2 );
			edgeHash[v0].push_back( e2 );
		}

		// look for edge v1->v2
		std::vector<Edge *> &v1_edges = edgeHash[v1];

		for( std::vector<Edge *>::iterator eit = v1_edges.begin(); eit != v1_edges.end(); ++eit )
		{
			Edge *e = *eit;
			if( e->getOtherVertex( v1 ) == v2 )
			{
				e1 = e;
				break;
			}
		}

		// edge not found?
		if( !e1 )
		{
			// create it
			e1 = createEdge( v1, v2 );
			// register edge with the hash
			edgeHash[v1].push_back( e1 );
			edgeHash[v2].push_back( e1 );
		}

		Triangle *e = createTriangle( v0, v1, v2, e0, e1, e2 );

		if( stop )
		{
			//e->tag = true;
			//return;
		}
	}


	std::vector<MeshEx::Vertex *> tempVerts( m_vertices.begin(), m_vertices.end() );
	
	for( std::vector<MeshEx::Vertex *>::iterator it=tempVerts.begin(); it != tempVerts.end(); ++it )
	{
		Vertex *v = *it;
		if( v->triangleRing.size() < 3 )
		{
			printf( "error degenerated vertex detected\n" );
			removeVertex( v );
		}
	}

	//// remove all triangles with boundary edges
	//for( std::vector<MeshEx::Triangle *>::iterator it=tempTris.begin(); it != tempTris.end(); ++it )
	//{
	//	Triangle *t = *it;

	//	if( t->tag2 )
	//	{
	//		int debugHit =0;
	//		debugHit++;
	//	}

	//	bool remove = false;

	//	for( size_t i=0; i<3; ++i )
	//		if( t->e[i]->isBoundaryEdge() )
	//		{
	//			remove = true;
	//			break;
	//		}

	//	if( remove )
	//	{
	//		printf( "removing degenerate triangle\n" );
	//		removeTriangle( t );
	//	}

	//}

	detectAndFillHoles();
}


//
// destructor
//
MeshEx::~MeshEx()
{
	for( std::vector<Vertex *>::iterator it = m_vertices.begin(); it != m_vertices.end(); ++it )
		delete *it;
	for( std::vector<Edge *>::iterator it = m_edges.begin(); it != m_edges.end(); ++it )
		delete *it;
	for( std::vector<Triangle *>::iterator it = m_triangles.begin(); it != m_triangles.end(); ++it )
		delete *it;

	m_vertices.clear();
	m_edges.clear();
	m_triangles.clear();
}


//
// creates a mesh which represents the MeshEx
//
Mesh *MeshEx::getMesh()
{
	std::map<Vertex *, int>        indicees; // the index into the vertexvector for each vertex
	std::vector<math::Vec3f> vertexPosition; // position for each vertex
	std::vector<int>       triangleIndicees; // 3 indicees into the vertexvector for each triangle


	vertexPosition.reserve( m_vertices.size() );
	triangleIndicees.reserve( m_triangles.size()*3 );

	int count = 0;
	for( std::vector<Vertex *>::iterator it = m_vertices.begin(); it != m_vertices.end(); ++it )
	{
		Vertex *v = *it;
		vertexPosition.push_back( v->position );
		indicees[ v ] = count++;
	}

	for( std::vector<Triangle *>::iterator it = m_triangles.begin(); it != m_triangles.end(); ++it )
	{
		Triangle *t = *it;
		triangleIndicees.push_back( indicees[t->v0] );
		triangleIndicees.push_back( indicees[t->v1] );
		triangleIndicees.push_back( indicees[t->v2] );
	}

	return new Mesh( vertexPosition, triangleIndicees );
}

//
//
//
void MeshEx::clear()
{
	for( std::vector<Triangle *>::iterator it = m_triangles.begin(); it != m_triangles.end(); ++it )
		delete *it;
	m_triangles.clear();
	for( std::vector<Edge *>::iterator it = m_edges.begin(); it != m_edges.end(); ++it )
		delete *it;
	m_edges.clear();
	for( std::vector<Vertex *>::iterator it = m_vertices.begin(); it != m_vertices.end(); ++it )
		delete *it;
	m_vertices.clear();
}

//
// splits the given edge into 2 and returns the splitnode which was created to seperate the edge
//
MeshEx::Vertex *MeshEx::splitEdge( Edge *edge, const math::Vec3f &splitPosition )
{
	// if the edge had an element on the left then it must be replaced with 2 new m_triangles
	// same counts for the right side
	Triangle   *left = 0; // indicates whether the edge had a reference to a left element
	Triangle    leftCopy; // used to hold the values of the parent element after it has been removed by calling removeEdge
	Vertex *l1[3], *l2[3]; // m_vertices for the new sub-m_triangles which appear from dividing the left element

	Triangle  *right = 0; // indicates whether the edge had a reference to a right element
	Triangle   rightCopy; // used to hold the values of the parent element after it has been removed by calling removeEdge
	Vertex *r1[3], *r2[3]; // m_vertices for the new sub-m_triangles which appear from dividing the right element

	// add the new node
	Vertex *sn = this->createVertex( splitPosition );

	// is there a element on the left side of the edge?
	if( edge->left )
	{
		left = edge->left;
		leftCopy = *edge->left;

		// edgeindex of the edge which has to be split
		int ei;
		// which edge has to be split?
		for( ei=0; ei<3; ++ei )
			if( left->e[ei] == edge )
				break;
		l1[0] = left->v[ei];
		l1[1] = sn;
		l1[2] = left->v[(ei+2)%3];
		l2[0] = sn;
		l2[1] = left->v[(ei+1)%3];
		l2[2] = left->v[(ei+2)%3];
	}

	// is there a element on the right side of the edge?
	if( edge->right )
	{
		right = edge->right;
		rightCopy = *edge->right;

		// edgeindex of the edge which has to be split
		int ei;
		// which edge has to be split?
		for( ei=0; ei<3; ++ei )
			if( right->e[ei] == edge )
				break;
		r1[0] = right->v[ei];
		r1[1] = sn;
		r1[2] = right->v[(ei+2)%3];
		r2[0] = sn;
		r2[1] = right->v[(ei+1)%3];
		r2[2] = right->v[(ei+2)%3];
	}

	// remove edge from list (which will remove the 2 wing m_triangles also )
	removeEdge( edge );

	// add the new m_triangles (beware that left and right point to invalid locations since these m_triangles
	// have been removed by removeEdge)
	if( left )
	{
		Triangle *el1 = createTriangle( l1[0], l1[1], l1[2] );
		Triangle *el2 = createTriangle( l2[0], l2[1], l2[2] );
	}
	if( right )
	{
		Triangle *er1 = createTriangle( r1[0], r1[1], r1[2] );
		Triangle *er2 = createTriangle( r2[0], r2[1], r2[2] );
	}

	return sn;
}

//
// looks for an edge between the given m_vertices - returns 0 if it could not be found
//
MeshEx::Edge *MeshEx::findEdge( Vertex *v1, Vertex *v2 )
{
	for( std::vector<Edge *>::iterator eit = m_edges.begin(); eit != m_edges.end(); ++eit )
	{
		Edge *e = *eit;
		if( ((e->v1 == v1)&&(e->v2 == v2))||((e->v1 == v2)&&(e->v2 == v1)) )
			return e;
	}

	return 0;
}

//
// creates a new edge from 2 given m_vertices
//
MeshEx::Edge *MeshEx::createEdge( Vertex *v1, Vertex *v2 )
{
	Edge *edge = new Edge( v1, v2 );
	// store edge
	m_edges.push_back( edge );
	return edge;
}

//
// removes edge from the list of m_edges
//
void MeshEx::removeEdge( Edge *edge )
{
	for( std::vector<Edge *>::iterator eit = m_edges.begin(); eit != m_edges.end(); ++eit )
		if( edge == *eit )
		{
			m_edges.erase( eit );

			// remove the m_triangles which neighboured the edge
			if( edge->left )
				removeTriangle( edge->left );
			if( edge->right )
				removeTriangle( edge->right );

			delete edge;			

			return;
		}
}

//
// creates a Vertex and adds it to the node list with the given world position
//
MeshEx::Vertex *MeshEx::createVertex( const math::Vec3f &position )
{
	Vertex *n = new Vertex();
	n->position = position;
	m_vertices.push_back( n );
	return n;
}

//
// removes given node from the node list
//
void MeshEx::removeVertex( Vertex *vertex )
{
	std::vector<Triangle *> triangleRing_copy( vertex->triangleRing.begin(), vertex->triangleRing.end() );
	// remove the triangles which contain v
	for( std::vector<Triangle *>::iterator it=triangleRing_copy.begin(); it != triangleRing_copy.end(); ++it )
		// remove Triangle and dont remove isolated elements
		removeTriangle( *it, false, true );

	// check
	if( !vertex->isDesolate() )
		// strange
		printf( "error: vertex still referenced although all triangles sharing the vertex had been removed\n" );
		
	m_vertices.erase( std::remove( m_vertices.begin(), m_vertices.end(), vertex ) );

	delete vertex;
}


//
// creates a Triangle and adds it to the element list with the given node indices
//
MeshEx::Triangle *MeshEx::createTriangle( const size_t &index0, const size_t &index1, const size_t &index2 )
{
	return createTriangle( m_vertices[index0], m_vertices[index1], m_vertices[index2] );
}

//
// creates a Triangle and adds it to the triangle list from given vertices
//
MeshEx::Triangle *MeshEx::createTriangle( Vertex *v0, Vertex *v1, Vertex *v2 )
{
	Triangle *e = new Triangle();

	e->v0 = v0;
	e->v1 = v1;
	e->v2 = v2;

	// register the element with its ith node
	e->v0->registerTriangle( e );
	e->v1->registerTriangle( e );
	e->v2->registerTriangle( e );


	// create winged edge information
	Edge *edge = 0;

	edge = findEdge( e->v0, e->v1 );
	// if edge between n0 and n1 does not exist
	if( !edge )
		// create edge
		edge = createEdge( e->v0, e->v1 );
	edge->registerTriangle( e );
	e->e1 = edge;

	edge = findEdge( e->v1, e->v2 );
	// if edge between n1 and n2 does not exist
	if( !edge )
		// create edge
		edge = createEdge( e->v1, e->v2 );
	edge->registerTriangle( e );
	e->e2 = edge;


	edge = findEdge( e->v2, e->v0 );
	// if edge between n2 and n0 does not exist
	if( !edge )
		// create edge
		edge = createEdge( e->v2, e->v0 );
	edge->registerTriangle( e );
	e->e3 = edge;

	e->computeProperties();


	m_triangles.push_back( e );
	return e;
}

//
// creates a Triangle and adds it to the triangle list from given edges
//
MeshEx::Triangle     *MeshEx::createTriangle( Vertex *v0, Vertex *v1, Vertex *v2, Edge *e0, Edge *e1, Edge *e2 )
{
	Triangle *t = new Triangle();

	// find the 3 vertices from edges


	t->v0 = v0;
	t->v1 = v1;
	t->v2 = v2;

	t->e1 = e0;
	t->e2 = e1;
	t->e3 = e2;

	// register the element with its ith node
	t->v0->registerTriangle( t );
	t->v1->registerTriangle( t );
	t->v2->registerTriangle( t );

	for( int i=0; i<3; ++i )
	{
		if( !t->e[i]->isBoundaryEdge() )
		{
			stop = true;
			//t->e[i]->tag = true;
		}
		t->e[i]->registerTriangle( t );
	}

	t->computeProperties();

	m_triangles.push_back( t );
	return t;
}

MeshEx::Triangle *MeshEx::createTriangle( MeshEx::Vertex *v0, MeshEx::Vertex *v1, MeshEx::Vertex *v2, std::vector<Edge *> &edges )
{
	Triangle *t = new Triangle();

	// find the 3 vertices from edges
	t->v0 = v0;
	t->v1 = v1;
	t->v2 = v2;

	Edge *e0, *e1, *e2;

	e0 = e1 = e2 = 0;

	for( std::vector<Edge *>::iterator it = edges.begin(); it != edges.end(); ++it )
	{
		Edge *e = *it;
		if( e->contains( v0, v1 ) )
			e0 = e;
		else
		if( e->contains( v1, v2 ) )
			e1 = e;
		else
		if( e->contains( v2, v0 ) )
			e2 = e;
	}

	t->e1 = e0;
	t->e2 = e1;
	t->e3 = e2;

	// register the element with its ith node
	t->v0->registerTriangle( t );
	t->v1->registerTriangle( t );
	t->v2->registerTriangle( t );

	for( int i=0; i<3; ++i )
		t->e[i]->registerTriangle( t );

	t->computeProperties();

	m_triangles.push_back( t );
	return t;
}


//
// removes element from the list of m_triangles
//
void MeshEx::removeTriangle( Triangle *element, bool removeVertices, bool removeEdges )
{
	// remove m_m_triangles from the list
	for( std::vector<Triangle *>::iterator eit = m_triangles.begin(); eit != m_triangles.end(); ++eit )
		if( element == *eit )
		{
			m_triangles.erase( eit );

			for( size_t i=0; i<3; ++i )
			{
				// remove references to the element from neighbouring m_m_edges
				element->e[i]->unregisterTriangle( element );
				// remove the element from the elementring of neighbouring m_m_vertices
				element->v[i]->unRegisterTriangle( element );
			}

			// TODO: take removeIsolated into account

			// if isolated edges have to be removed
			if( removeEdges )
				for( size_t i=0; i<3; ++i )
					if( element->e[i]->isDesolate() )
						removeEdge( element->e[i] );

			// if isolated vertices have to be removed
			if( removeVertices )
				for( size_t i=0; i<3; ++i )
					if( element->v[i]->isDesolate() )
						removeVertex( element->v[i] );

			delete element;

			return;
		}
}

//
// retriangulates the given triangle taking the
// given points into acount (which have to be coplanar with the triangle)
//
void MeshEx::reTriangulate( Triangle *t, const std::vector<math::Vec3f> &pointsWithin )
{
	struct helper
	{
		helper( Vertex *_v1, Vertex *_v2 ) : v1(_v1), v2(_v2) { sqDistance = (v2->position - v1->position).getSquaredLength(); }
		float sqDistance;
		Vertex *v1;
		Vertex *v2;

		bool operator<( const helper &rhs )
		{
			return sqDistance < rhs.sqDistance;
		}
	};





	// prepare -------------------------------------------------------------
	// store the m_m_vertices and m_m_edges of the original triangle
	Vertex                    *v[3];
	Edge                      *e[3];
	math::Vec3f              normal;

	normal = t->normal;

	for( size_t i=0; i<3; ++i )
	{
		v[i] = t->v[i];
		e[i] = t->e[i];
	}

	// remove the triangle and keep isolated vertices and edges
	removeTriangle( t, false, false );

	// if there is only one point within the triangle (which will be the case to 99.99 percent)
	if( pointsWithin.size() == 1 )
	{
		// then we triangulate manually for speed reasons
		Vertex *vn = createVertex( pointsWithin[0] );

		Edge *e0 = createEdge( v[0], vn );
		Edge *e1 = createEdge( v[1], vn );
		Edge *e2 = createEdge( v[2], vn );

		createTriangle( v[0], v[1], vn, e[0], e1, e0 );
		createTriangle( v[1], v[2], vn, e[1], e2, e1 );
		createTriangle( v[2], v[0], vn, e[2], e0, e2 );

		// done
		return;
	}

	// -----------------------------------------------------------------------------------------------------------
	// there are more than 1 points within the triangle so we have to perform                                    -
	// greedy triangulation                                                                                      -
	// -----------------------------------------------------------------------------------------------------------

	// create vertices for each point which has to be taken into account while triangulating
	std::vector<Vertex *>    vertices;
	std::vector<Edge *>         edges;

	// add the candidate vertices to the local vertex list
	for( std::vector<math::Vec3f>::const_iterator it=pointsWithin.begin(); it != pointsWithin.end(); ++it )
	{
		vertices.push_back( createVertex( *it ) );
		// the vertices normal is identical to the normal of the original polygon
		vertices.back()->normal = normal;
	}

	// precompute sorted distances from each point to each other point (without the vertices of the original
	// trianlge which form the boundary)

	std::vector<helper> sqDistances;  // sorted distances of each pair of points

	for( size_t i=0; i<vertices.size()-1; ++i )
		for( size_t j=i+1; j<vertices.size(); ++j )
			sqDistances.push_back( helper( vertices[i], vertices[j] ) );

	// add the distance of each point to each of the boundary vertices comming from the
	// original triangle, we do this after computing the distance of each vertex to each other
	// vertex so that we dont check-create the boundary edges from the original triangle
	// which already exist
	for( size_t i=0; i<3; ++i )
		for( size_t j=0; j<vertices.size(); ++j )
			sqDistances.push_back( helper( v[i], vertices[j] ) );


	// sort
	std::sort( sqDistances.begin(), sqDistances.end() );

	// triangulate the area of the original triangle -----------------------------------------

	// perform greedy triangluation
	for( std::vector<helper>::iterator it=sqDistances.begin(); it != sqDistances.end(); ++it )
	{
		helper *h = &(*it);

		// test for intersection with all existing edges
		bool intersection = false;
		for( std::vector<Edge *>::iterator eit = edges.begin(); eit != edges.end(); ++eit)
		{
			Edge *e = (*eit);

			// intersection?
			if( !(e->contains(h->v1) || e->contains(h->v2)) && e->intersectsLine( h->v1->position, h->v2->position ) )
			{
				intersection = true;
				break;
			}
		}

		// if there was no intersection...
		if( !intersection )
			// create edge from the 2 vertices and add it to the local edges list
			edges.push_back( createEdge( h->v1, h->v2 ) );
	}


	// add the original vertices to the vertex list
	vertices.push_back( v[0] );
	vertices.push_back( v[1] );
	vertices.push_back( v[2] );

	// add the edges from the original triangle which circumfere the whole triangulation area
	edges.push_back( e[0] );
	edges.push_back( e[1] );
	edges.push_back( e[2] );


	// create triangles from our edges
	createTrianglesFromEdges( edges, normal );


	// done
}


//
// this method takes a list of already created edges and creates
// triangles for filling triangulated areas
//
//
void MeshEx::createTrianglesFromEdges( std::vector<Edge *> &edges, const math::Vec3f &normal )
{
	// find and create triangles from generated edges --------------------------------------------
	std::vector<Vertex *>    vertices; // vertices which are involved
	std::vector<Triangle *> triangles; // cached list of created triangles, for avoiding duplicates

	// collect outgoint edges for each vertes
	std::map<Vertex *, std::vector<Edge *> > vertexEdges;
	for( std::vector<Edge *>::iterator eit = edges.begin(); eit != edges.end(); ++eit)
	{
		Edge *edg = (*eit);

		vertexEdges[edg->v1].push_back( edg );
		vertexEdges[edg->v2].push_back( edg );

		vertices.push_back( edg->v1 );
		vertices.push_back( edg->v2 );
	}

	// remove duplicates in vertex list
	std::sort( vertices.begin(), vertices.end() );
	vertices.erase( std::unique( vertices.begin(), vertices.end() ), vertices.end()  );

	// now for each vertex
	for( std::vector<Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
	{
		// we look for 3 vertices which may build up a triangle
		Vertex *v1 = *it;
		Vertex *v2 = 0; // 2nd and 3rd vertices arent known yet
		Vertex *v3 = 0;

		// the 3 edges which build the triangle and contain v1, v2, v3 - arent known yet
		Edge  *e1 = 0;
		Edge  *e2 = 0;
		Edge  *e3 = 0;

		// the normalized direction vectors of the 3 edges
		math::Vec3f vec_edge1;
		math::Vec3f vec_edge2;
		math::Vec3f vec_edge3;

		// the first edge will together with the normal give the direction of the right vector
		// which helps building the triangle in the right vertex order (cw or ccw) - same goes for e2 in respect to e3
		math::Vec3f e1_right;
		math::Vec3f e2_right;

		// the side indicators are the edge directions projected onto the right vectors of the
		// previous edge - their sign tells whether the edge loop goes into the right direction
		// the crossproduct between the triangle normal and the edge candidates will indicate the direction
		// in which edges are selected to build the triangle - this is to make sure that the order of vertices
		// will not build triangles whos normal will point in the opposite direction as the normal of the
		// original triangle
		float e2_side;
		float e3_side;

		std::vector<Edge *> &v1_edgeList = vertexEdges[v1];

		// for each edge going out from the current vertex
		for( std::vector<Edge *>::iterator eit = v1_edgeList.begin(); eit != v1_edgeList.end(); ++eit )
		{
			e1 = (*eit); // our edge for e1

			// take the candidate for v2 to be the vertex on the other side of the edgecandidate e1
			v2 = e1->getOtherVertex(v1);

			vec_edge1 = math::normalize( v2->position - v1->position );

			// get the local right vector
			e1_right = math::crossProduct( vec_edge1, normal );

			std::vector<Edge *> &v2_edgeList = vertexEdges[v2];

			float minAngle = -1.0f;
			Edge *minEdge = 0;

			// find the outgoing edge with the minimal angle to the current edge of v1
			for( std::vector<Edge *>::iterator eit2 = v2_edgeList.begin(); eit2 != v2_edgeList.end(); ++eit2 )
			{
				Edge *ec2 = *eit2; // candidate for e2

				// skip if e1 and ec2 are the same (we dont want to test the same edge against itsself)
				if( e1 == ec2 )
					continue;

				math::Vec3f vec_edge2c = math::normalize( ec2->getOtherVertex(v2)->position - v2->position );

				// project the vector of edge 2 onto the vector which defines the direction in which
				// any next edge of the triangle may point
				e2_side = math::dotProduct( vec_edge2c, e1_right );


				// compute the enclosing angle between the edgecandidates e1 and e2
				float angle = math::dotProduct( vec_edge1, -vec_edge2c );

				// we are looking for the edge with the closest angle
				// by restricting e2_direction to be positive, we only look at edges which would build
				// the triangle in the right order
				if( (e2_side<0.0f) && ((angle > minAngle) || !minEdge) )
				{
					minEdge = ec2;
					minAngle = angle;
					vec_edge2 = vec_edge2c;
				}
			}

			// if we have not found a minimal edge (which only may be the case when there is no
			// edge pointing into the direction the triangle would have to be build)
			if( !minEdge )
				// we can skip looking for this edge
				continue;

			e2 = minEdge;

			// now we found the third vertex for our triangle candidate we look for the edge with
			// the closest angle to the edge coming from v2 and also restrict the edge to point
			// into the right direction. In addition of course the third edge must have contain
			// v1 and v2 so that the triangle is closed.
			v3 = e2->getOtherVertex( v2 );

			std::vector<Edge *> &v3_edgeList = vertexEdges[v3];

			// get the local right vector of edge 2
			e2_right = math::crossProduct( vec_edge2, normal );

			minAngle = -1.0f;
			minEdge = 0;

			// find the outgoing edge with the minimal angle to the current edge of v2
			for( std::vector<Edge *>::iterator eit3 = v3_edgeList.begin(); eit3 != v3_edgeList.end(); ++eit3 )
			{
				Edge *ec3 = *eit3;  // candidate for edge 3

				// skip if e2 and ec3 are the same (we dont want to test the same edge against itsself)
				// or if the edge does not contain v1 since then we woulnt get a triangle
				if( (e2 == ec3) || !ec3->contains(v1) )
					continue;

				vec_edge3 = math::normalize( v1->position - v3->position );

				// project the vector of edge 3 onto the vector which defines the direction in which
				// any next edge of the triangle may point
				e3_side = math::dotProduct( vec_edge3, e2_right );

				float angle = math::dotProduct( vec_edge2, -vec_edge3 );

				// by restricting e2_direction to be positive, we only look at edges which would build
				// the triangle in the right order
				if( (e3_side<0.0f) && ((angle > minAngle) || !minEdge) )
				{
					minEdge = ec3;
					minAngle = angle;
				}
			}

			// if we have not found a minimal edge (which only may be the case when there is no
			// edge pointing into the direction the triangle would have to be build)
			// we also go back to the next iteration when e3 does not contain v1 which means
			// that no path of 3 edges could be found to go from v1 to v2 to v3 back to v1
			if( !minEdge || !minEdge->contains(v1) )
				// we can skip looking for this edge
				continue;

			e3 = minEdge;

			// skip, if the triangle already exists
			bool alreadyExists = false;
			for( size_t t=0; t<triangles.size(); ++t )
				if( triangles[t]->contains(v1) && triangles[t]->contains(v2) && triangles[t]->contains(v3)  )
					alreadyExists = true;

			if( alreadyExists )
				continue;

			// now we have found our 3 edges and we can be sure by checking for
			// minimal edge-angles, that there are no edges/vertices within the area enclosed by
			// v1, v2 and v3
			// by checking for the direction of the edges we also made sure that the triangle
			// will have the right vertex ordering
			triangles.push_back( createTriangle( v1, v2, v3) );
			triangles.back()->normal = normal;
		}
	}
}

//
// MeshEx::Triangle -------------------------------------------------------------------------
//

//
// constructor
//
MeshEx::Triangle::Triangle()
{
	tag = false;
	tag2 = false;
}

//
// transforms the given coordinate in local 2d space into the global 3d space in which the triangle/element lies
//
//
math::Vec3f MeshEx::Triangle::convertFromLocalToGlobalSpace( const math::Vec2f &localSpace, bool isVector )
{
	// if the result has to be returned relative to the element origin...
	if( isVector )
		// ...then we dont add the center
		return localSpace.x*u + localSpace.y*v_v;
	else
		// the result has to be in absolute coordinates
		return v0->position + localSpace.x*u + localSpace.y*v_v;
}

//
// inverse operation of the above
//
// The argument isVector indicates whether the globalSpace argument is given relative to the element-center (isVector==true)
// or is given in absolut coordinates
//
math::Vec2f MeshEx::Triangle::convertFromGlobalToLocalSpace( const math::Vec3f &globalSpace, bool isVector )
{
	// substract the center
	math::Vec3f temp = globalSpace;

	// if the given position is given in absolute world coordinates
	if( !isVector )
		temp -= v0->position;

	// now project the vector onto u and v
	return math::Vec2f( math::dotProduct( temp, u ), math::dotProduct( temp, v_v ) );
}

//
// (re) computes normal
//
void MeshEx::Triangle::computeNormal()
{
	normal = math::normalize( math::crossProduct( v1->position - v0->position, v2->position - v0->position ) );
}


//
// computes normal, area, etc...
//
void MeshEx::Triangle::computeProperties()
{
	// compute element center
	for( size_t i=0; i<3; ++i )
		center += v[i]->position;
	center /= 3.0f;


	math::Vec3f u_unormalized = v[1]->position - v[0]->position;
	math::Vec3f         n3vec = v[2]->position - v[0]->position;

	// compute element normal
	u = math::normalize( u_unormalized );
	normal = math::normalize( math::crossProduct( u, n3vec ) );

	v_v = math::normalize( math::crossProduct( normal, u ) );

	// compute beta (barycentric matrix of the local coordinate frame)
	// the base vectors of beta are the positions of the m_vertices within the local coordinate frame

	// the first column is the position of the first node within the local coordinate frame which is the orign -> zero
	beta._11 = 0.0f;
	beta._12 = 0.0f;
	beta._13 = 1.0f;

	// the second column is the position of the second node within the local coordinate frame which is the
	// u-axis times the length of the original unormalized u vector
	beta._21 = u_unormalized.getLength();
	beta._22 = 0.0f;
	beta._23 = 1.0f;

	// the third column is the position of the third node within the local coordinate frame which is the
	// vector (e->n[2]->position - e->n[0]->position) projected onto u and v
	beta._31 = math::dotProduct( n3vec, u );
	beta._32 = math::dotProduct( n3vec, v_v );
	beta._33 = 1.0f;

	// invert beta
	// we dont get away with simply transposing it since betas basis is not orthonormal
	beta.invert();


	// compute area of the element
	float la = (v1->position - v0->position).getLength(); // compute lengths of the triangle sides
	float lb = (v2->position - v1->position).getLength();
	float lc = (v2->position - v0->position).getLength();
	float s = 0.5f*( la+lb+lc ); // compute the semiperimeter
	area = sqrt( s*(s-la)*(s-lb)*(s-lc) ); // compute the area	
}

//
// returns the third node which neither equals n1 nor n2
//
MeshEx::Vertex *MeshEx::Triangle::getOtherVertex( Vertex *vertex1, Vertex *vertex2 )
{
	for( size_t i=0; i<3; ++i )
		if( (v[i] != vertex1) && (v[i] != vertex2) )
			return v[i];

	printf( "error in MeshEx::Triangle::getOtherVertex() : could not find third node from given n1 and n2\n" );
	return 0;
}

//
// returns the third node which neither equals edge->v1 nor edge->v2
//
MeshEx::Vertex *MeshEx::Triangle::getOtherVertex( Edge *edge )
{
	for( size_t i=0; i<3; ++i )
		if( (v[i] != edge->v1) && (v[i] != edge->v2) )
			return v[i];

	printf( "error in MeshEx::Triangle::getOtherVertex() : could not find third node from given v1 and v2\n" );
	return 0;
}

//
// returns true if one of the element m_vertices equals the given node n
//
bool MeshEx::Triangle::contains( Vertex *vertex ) const
{
	if( (vertex == v0)||(vertex == v1)||(vertex == v2) )
		return true;
	return false;
}

//
// returns the local index of the node within the element
//
size_t MeshEx::Triangle::getVertexsLocalIndex( Vertex *vertex )
{
	for( size_t i=0; i<3; ++i )
		if( v[i] == vertex )
			return i;
	// node does not belong to the given element
	printf( "error in getVertexsLocalIndex() : node does not belong to element\n" );
	return 3;
}

//
// returns the edge which is shared with the given triangle or 0 if they dont share an edge
//
MeshEx::Edge *MeshEx::Triangle::getSharedEdge( Triangle *t )
{
	for( size_t i=0; i<3; ++i )
		if( e[i]->contains(t)  )
			return e[i];
	return 0;
}

//
// returns the edge of the triangle which contains the 2 given m_vertices
//
MeshEx::Edge *MeshEx::Triangle::getEdge( Vertex *vertex1, Vertex *vertex2 )
{
	for( int i=0; i<3; ++i )
		if( e[i]->contains( vertex1, vertex2 ) )
			return e[i];
	return false;
}

//
// returns the one edge of the triangle which does not contain the vertex v
//
MeshEx::Edge *MeshEx::Triangle::getOtherEdge( MeshEx::Vertex *v )
{
	// first we will make sure that the triangle actually contains v
	if( !contains(v) )
		// lets crash...
		return 0;

	for( size_t i=0; i<3; ++i )
		// only one edge within the triangle wont contain the vertex v
		if( !e[i]->contains(v) )
			return e[i];

	return 0;
}


//
// MeshEx::Vertex -------------------------------------------------------------------------
//

//
// Constructor
//
MeshEx::Vertex::Vertex()
{
}

//
// returns true if no element references this node
//
bool MeshEx::Vertex::isDesolate()
{
	return triangleRing.empty();
}

//
// adds the given element to the elementring list
//
void MeshEx::Vertex::registerTriangle( Triangle *t )
{
	triangleRing.push_back( t );
}

//
// remvoes the given element to the elementring list
//
void MeshEx::Vertex::unRegisterTriangle( Triangle *t )
{
	for( std::vector<Triangle*>::iterator eit=triangleRing.begin(); eit != triangleRing.end(); ++eit )
		if( *eit == t )
		{
			triangleRing.erase( eit );
			return;
		}


	printf( "Error Triangle was not member of the element ring of given node.\n" );
}

//
// Returns 1 if the node lies on the left side of the given plane (nodeplanedistance < 0.0f)
// or 2 otherwise.
// Returns 0 if point lies on the plane.
//
int MeshEx::Vertex::getPlaneSide( const math::Vec3f &normal, const float &distance )
{
	 float planeDistance = math::distancePointPlane( position, normal, distance );

	 // if the planeDistance is 0.0f then the point lies directly on the plane and can not be
	 // used to decide on which side the element lies
	 if( planeDistance == 0.0f )
		 return 0;

	 // signed plane distance
	 if( planeDistance < 0.0f )
		 // point lies on the left side
		 return 1;
	 else
		 // point lies on the right side
		 return 2;
}

//
// computes normal from the surrounding tris
// NOTE: works only if the surrounding tris have a valid normal
//
void MeshEx::Vertex::computeNormal()
{
	normal = math::Vec3f( 0.0f, 0.0f, 0.0f );
	for( std::vector<Triangle *>::iterator it=triangleRing.begin(); it != triangleRing.end(); ++it )
		normal += (*it)->normal;
	normal.normalize();
}


//
// MeshEx::Edge -------------------------------------------------------------------------
//

//
// constructor
//
MeshEx::Edge::Edge( Vertex *_v1, Vertex *_v2 ):left(0), right(0), v1(_v1), v2(_v2)
{
	tag = false;
};

//
// returns true if no m_triangles references this edge
//
bool MeshEx::Edge::isDesolate()
{
	if( left || right )
		return false;
	return true;
}

//
// convinience function to get the other node of the two m_vertices which are referenced by the edge
//
MeshEx::Vertex *MeshEx::Edge::getOtherVertex( Vertex *v )
{
	if( v == v1 )
		return v2;
	else
	if( v == v2 )
		return v1;
	else
		printf( "error in MeshEx::Edge::getOtherVertex() : given node n does not belong to the edge\n" );
	return 0;
}

//
// returns the opposite triangle
//
MeshEx::Triangle *MeshEx::Edge::getOtherTriangle( Triangle *t )
{
	if( left == t )
		return right;
	else
	if( right == t )
		return left;
	else
		return 0;
}

//
// returns true if the edge contains the given node
//
bool MeshEx::Edge::contains( Vertex *v )
{
	if( (v1 == v)||(v2 == v) )
		return true;
	return false;
}

//
// returns true if the edge contains the given m_vertices
//
bool MeshEx::Edge::contains( Vertex *vertex1, Vertex *vertex2 )
{
	if( ((vertex1 == v1)&&(vertex2 == v2))||((vertex2 == v1)&&(vertex1 == v2)) )
		return true;
	else
		return false;
}

//
// returns true if the edge contains the given m_triangles
//
bool MeshEx::Edge::contains( Triangle *t )
{
	if( (left==t)||(right==t) )
		return true;
	return false;
}

//
// returns true if the edge has only one element reference instead of two
//
bool MeshEx::Edge::isBoundaryEdge()
{
	// if both element references are valid...
	if( left && right )
		// then this is not a boundary edge
		return false;
	// boundary or desolate edge
	return true;
}

//
// convinience function to compute the intersection of the edge with a given plane
//
bool MeshEx::Edge::intersectsPlane( const math::Vec3f &normal, const float &distance, math::Vec3f &intersectionPoint )
{
	return math::intersectionRayPlane( math::Ray3d( v1->position, v2->position ), normal, distance, intersectionPoint );
}

//
// convinience function to compute the intersection of the edge with a given line
//
bool MeshEx::Edge::intersectsLine( const math::Vec3f &p1, const math::Vec3f &p2, math::Vec3f &intersectionPoint  )
{
	return math::intersectionRayRay( math::Ray3d( v1->position, v2->position ), math::Ray3d( p1, p2 ), intersectionPoint );	
}


//
// convinience function to compute the intersection
// of the edge with a given line without giving back the intersection point
//
bool MeshEx::Edge::intersectsLine( const math::Vec3f &p1, const math::Vec3f &p2 )
{
	math::Vec3f intersection;
	return math::intersectionRayRay( math::Ray3d( v1->position, v2->position ), math::Ray3d( p1, p2 ), intersection );	
}

//
// assigns the given element to the left or right wing of the edge (dependand on which wing is 0)
//
void MeshEx::Edge::registerTriangle( Triangle *t )
{
	if( !left )
	{
		left = t;
	}else
	if( !right )
	{
		right = t;

		if( left == right )
			printf( "error - edge must not reference the same triangle twice\n" );
	}else
	{
		printf( "error - edge must not belong to more than 2 m_triangles\n" );
		tag = true;
		left->tag = true;

		left->tag2 = true;
	}
}

//
// sets either the left or right reference to zero when it references the given element
//
void MeshEx::Edge::unregisterTriangle( Triangle *t )
{
	if( left == t )
	{
		left = 0;
	}else
	if( right == t )
	{
		right = 0;
	}else
	{
		printf( "error - edge did not reference the given element\n" );
	}
}