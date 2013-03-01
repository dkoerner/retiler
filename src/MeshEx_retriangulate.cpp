/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#include "MeshEx.h"
#include <list>

// CGAL -------------------------
typedef CGAL::Triangulation_vertex_base_2<K>                                                               Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K, CGAL::Triangulation_face_base_with_info_2<int, K> > Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                                       TDS;
typedef CGAL::Exact_predicates_tag                                                                       Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>                                Triangulation;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Face_iterator                 Face_iterator;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Face_handle                     Face_handle;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Vertex_handle                 Vertex_handle;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>::Point                                 Point;
typedef CGAL::Polygon_2<K, std::vector<Point> >                                                     Polygon_2;

math::Vec2f _get2d( math::Vec3f e1,  math::Vec3f e2, math::Vec3f point )
{
	return math::Vec2f( math::dotProduct( e1, point ), math::dotProduct( e2, point ) );
}

CGAL::Point_3<K> Point_3( const math::Vec3f &vec )
{
	return CGAL::Point_3<K>( vec.x, vec.y, vec.z );
}


///
/// \brief used as helper to connect to vertices together
///
struct DistanceHelper
{
	DistanceHelper( MeshEx::Vertex *_v1, MeshEx::Vertex *_v2 ) : v1(_v1), v2(_v2)
	{
		CGAL::Segment_3<K> seg( Point_3(v1->position), Point_3(v2->position) );
		sqDistance = seg.squared_length();
	}
	float   sqDistance;
	MeshEx::Vertex *v1;
	MeshEx::Vertex *v2;

	bool operator<( const DistanceHelper &rhs )
	{
		return sqDistance < rhs.sqDistance;
	}
};

///
/// \brief represents a simple polygon and falls back to CGAL for certain computations
///
struct Polygon
{
	void push_back( MeshEx::Vertex *v, math::Vec2f projectedPosition )
	{
		// insert point to cgal polygon and map vertex to index of the cgal-poly point
		indicees[ v ] = vertexPoints.size();
		vertices.push_back( v );
		vertexPoints.push_back( Point( projectedPosition.x, projectedPosition.y ) );
		//poly.push_back( Point( projectedPosition.x, projectedPosition.y ) );
	}

	void push_back( MeshEx::Vertex *v, Point projectedPosition )
	{
		// insert point to cgal polygon and map vertex to index of the cgal-poly point
		indicees[ v ] = vertexPoints.size();
		vertices.push_back( v );
		vertexPoints.push_back( projectedPosition );
		//poly.push_back( projectedPosition );
	}

	///
	/// \brief returns the winding of the polygon
	///
	CGAL::Orientation orientation()
	{
		return CGAL::orientation_2( vertexPoints.begin(), vertexPoints.end() );
	}

	bool contains( MeshEx::Vertex *v )
	{
		return( std::find( vertices.begin(), vertices.end(), v ) != vertices.end() );
	}

	bool isTriangle()
	{
		return (vertices.size() == 3);
	}

	MeshEx::Vertex *getNext( MeshEx::Vertex *current )
	{
		std::vector<MeshEx::Vertex *>::iterator it = std::find( vertices.begin(), vertices.end(), current );

		if( it == vertices.end() )
			return 0;

		if( ++it == vertices.end() )
			return vertices[0];
		else
			return *it;

		return 0;
	}

	///
	/// This method will setup 2 polyongs which will form the left and right polygon when
	/// it would be divided by the 2 given vertices v1 and v2.
	/// The polygon must contain v1 and v2 in order to work properly
	///
	void split( Polygon *&leftPoly, Polygon *&rightPoly, MeshEx::Vertex *v1, MeshEx::Vertex *v2 )
	{
		// start from v1 and go to v2
		leftPoly = new Polygon();

		MeshEx::Vertex *current = v1;
		do
		{
			leftPoly->push_back( current, vertexPoints[indicees[current]] );
			current = getNext(current);
		}while( current != v2 );
		leftPoly->push_back( current, vertexPoints[indicees[current]] );


		// setup right polygon
		// start from v2 and go to v1
		rightPoly = new Polygon();

		current = v2;
		do
		{
			rightPoly->push_back( current, vertexPoints[indicees[current]] );
			current = getNext(current);
		}while( current != v1 );
		rightPoly->push_back( current, vertexPoints[indicees[current]] );
	}

	///
	/// changes the orientation of the polygon
	///
	void flip()
	{
		std::reverse( vertices.begin(), vertices.end() );
		std::reverse( vertexPoints.begin(), vertexPoints.end() );

		int count = 0;
		for( std::vector<MeshEx::Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it )
			indicees[*it] = count++;
	}
/*
	void split( Polygon &leftPoly, Polygon &rightPoly, MeshEx::Vertex *v1, MeshEx::Vertex *v2 )
	{
		// start from v1 and go to v2

		MeshEx::Vertex *current = v1;
		do
		{
			leftPoly.push_back( current, vertexPoints[indicees[current]] );
			current = getNext(current);
		}while( current != v2 );
		leftPoly.push_back( current, vertexPoints[indicees[current]] );


		// setup right polygon
		// start from v2 and go to v1
		current = v2;
		do
		{
			rightPoly.push_back( current, vertexPoints[indicees[current]] );
			current = getNext(current);
		}while( current != v1 );
		rightPoly.push_back( current, vertexPoints[indicees[current]] );
	}
*/
	std::vector<MeshEx::Vertex *>                    vertices; // vertices in order of appearence
	std::vector<Point>                           vertexPoints;
	std::map<MeshEx::Vertex *, int>                  indicees; // associates each vertex in the polyon with its index into the vertexPoints list
};

///
/// \brief helper structure to hold certain information about the incoming and outgoing edges on a vertex
///
struct EdgeInfoHelper
{
	EdgeInfoHelper() : e1(0), e2(0), next(0), prev(0), vert(0)
	{
	}


	void swap()
	{
		std::swap(e1,e2);
		std::swap(next,prev);
	}

	void registerEdge( MeshEx::Edge *e )
	{
		if( !e1 )
		{
			e1 = e;
		}else
		if( !e2 )
		{
			e2 = e;
		}else
			printf( "error : there must not be more than 2 edges for a vertex on the vertex ring\n" );
	}

	void setVertex( MeshEx::Vertex *v )
	{
		vert = v;
		next = e2->getOtherVertex( v );
		prev = e1->getOtherVertex( v );
	}

	MeshEx::Edge          *e1; // incomming edge
	MeshEx::Edge          *e2; // outgoing edge
	MeshEx::Vertex      *next; // next vertex in the vertex ring
	MeshEx::Vertex      *prev; // previous vertex in the vertex ring
	MeshEx::Vertex      *vert; // the vertex for which this helper works
	math::Vec2f     projected; // 2d coordinates within the plane onto which the vertex has been projected
};

//
// removes the given vertex and retriangulates the hole which would be made
//
// The method returns true if the vertex has been removed and false if the vertex has been retained.
//
bool MeshEx::removeVertexAndReTriangulateNeighbourhood( Vertex *v )
{
	std::vector<MeshEx::Edge *>                                          boundaryEdges; // boundary edges which form the polygon which has to be triangulated
	std::vector<MeshEx::Edge *>                                          criticalEdges; // edges which are not only connected to other vertices of the vertexring through boundary edges
	std::map<MeshEx::Vertex *, EdgeInfoHelper>                            boundaryRing; // maps incomming and outgoing edge to each vertex of the boundary vertex-ring



	// gather and prepare information ------------------------------------------------
	for( std::vector<Triangle *>::iterator it = v->triangleRing.begin(); it != v->triangleRing.end(); ++it )
	{
		Triangle *t = (*it);

		Edge *boundaryEdge = t->getOtherEdge( v );

		boundaryRing[boundaryEdge->v1].registerEdge( boundaryEdge );
		boundaryRing[boundaryEdge->v2].registerEdge( boundaryEdge );

		boundaryEdges.push_back( boundaryEdge );
	}


	// align the edges so that for each vertex e1 is the incomming and e2 is the outgoing edge
	Vertex *first = boundaryRing.begin()->first;
	Vertex *current = first;
	Vertex *next = 0;
	Vertex *prev = 0;
	
	do
	{
		// we have to be sure that each boundaryRing Vertex has an incomming and outgoing edge
		if( !boundaryRing[current].e1 || !boundaryRing[current].e2 )
		{
			printf( "error : edge ring around vertex not closed boundary vertex has less than 2 edges (polygon hole?)\n" );
			//++g_iterations;
			return false;
		}

		boundaryRing[current].setVertex( current );

		next = boundaryRing[current].next;

		if( boundaryRing[next].e1 != boundaryRing[current].e2 )
			boundaryRing[next].swap();

		current = next;
	}while( current != first );

	// we have to collect all edges going out from the vertex-ring vertices which are connected to
	// other vertices of the vertex ring - these will later be used for consistency check #2 (see chapter 4.4 of the paper)
	// for each vertex of the vertex Ring V
	for( std::map<Vertex *, EdgeInfoHelper>::iterator it = boundaryRing.begin(); it != boundaryRing.end(); ++it )
	{
		Vertex *rv       = it->first; // current vertex of the vertex ring

		// for each triangle referencing rv...
		for( std::vector<Triangle *>::iterator tit = rv->triangleRing.begin(); tit != rv->triangleRing.end(); ++tit )
		{
			Triangle *rv_tri = *tit;

			// ... which doesnt belong to the triangleRing of v
			if( std::find( v->triangleRing.begin(), v->triangleRing.end(), rv_tri ) == v->triangleRing.end() )
			{
				// collect the edges containing rv...
				for( size_t i=0; i<3; ++i )
					if( rv_tri->e[i]->contains(rv) )
						// ...and dont belong to the edgeRing list
						if( std::find( boundaryEdges.begin(), boundaryEdges.end(), rv_tri->e[i] ) == boundaryEdges.end() )
							// store the edge if the node on the other side references another vertex of the vertex-ring
							if( boundaryRing.find( rv_tri->e[i]->getOtherVertex(rv) ) != boundaryRing.end() )
								criticalEdges.push_back( rv_tri->e[i] );

			}
		}
	}

	// remove duplicate entries
	std::sort( criticalEdges.begin(), criticalEdges.end() );
	criticalEdges.erase( std::unique(criticalEdges.begin(), criticalEdges.end()), criticalEdges.end() );

	// Now we will project the neighbourhood of v onto a plane so that we can employ the greedy
	// triangulation. For this we will have to find a projection of the neighbourhood of v so that
	// no edges intersect on that plane. If they would, then the greedy triangulation would introduce
	// folds which means that the topology would be destroyed.
	// Looking for a working projection may lead to a number of projection-trials. For the first try
	// we will take the plane to be the plane defined by the normal of v and (its distance to the origin).
	// This will do the job most of the time.
	math::Vec3f                           normal; // direction of the projection plane
	math::Vec3f                            base1; // base vector which builds the 2d-coordinate system
	math::Vec3f                            base2; // base vector which builds the 2d-coordinate system
	float                               distance; // distance of the plane to the origin
	CGAL::Orientation boundaryPolygonOrientation; // orientation (clockwise/counterclockwise of the boundary polyon)

	distance = -math::dotProduct( v->position, normal );

	size_t   trialCount = 13; // only one trial for now
	size_t currentTrial = 0;
	bool    success = true;
	do
	{
		// we assume that we will be successfull
		success = true;

		switch(currentTrial)
		{
		case 0:
			// first trial: we take the plane defined by the normal of v
			normal = v->normal;
			break;
		// for all other trials we will try one of 12 different directions
		case  1:normal = math::normalize( math::Vec3f( .8507f, .4472f, .2764f ));break;
		case  2:normal = math::normalize( math::Vec3f( -.8507f, .4472f, .2764f ));break;
		case  3:normal = math::normalize( math::Vec3f( .8507f, -.4472f, -.2764f ));break;
		case  4:normal = math::normalize( math::Vec3f( -.8507f, -.4472f, -.2764f ));break;
		case  5:normal = math::normalize( math::Vec3f( .5257f, -.4472f, .7236f ));break;
		case  6:normal = math::normalize( math::Vec3f( .5257f, .4472f, -.7236f ));break;
		case  7:normal = math::normalize( math::Vec3f( -.5257f, -.4472f, .7236f ));break;
		case  8:normal = math::normalize( math::Vec3f( -.5257f, .4472f, -.7236f ));break;
		case  9:normal = math::normalize( math::Vec3f( .0f, .4472f, .8944f ));break;
		case 10:normal = math::normalize( math::Vec3f( .0f, -1.0f, .0f ));break;
		case 11:normal = math::normalize( math::Vec3f( .0f, 1.0f, .0f ));break;
		case 12:normal = math::normalize( math::Vec3f( .0f, -.4472f, -.8944f ));break;
		};

		// compute the basis of the 2d-coordinate system of the plane defined by current normal
		base1 = math::normalize( math::projectPointOnPlane( normal, distance, boundaryRing.begin()->first->position ) - v->position );
		base2 = math::normalize( math::crossProduct( normal, base1 ) );


		// project neighbours into the given plane
		for( std::map<Vertex *, EdgeInfoHelper>::iterator it = boundaryRing.begin(); it != boundaryRing.end(); ++it )
			it->second.projected = _get2d( base1, base2, math::projectPointOnPlane( normal, distance, it->first->position ) );

		// topologic constistency check #1 (see chapter 4.4 of the paper)
		// now do the consistency check: test if all edges dont intersect and dont lie over each other
		// if any edge between the projected vertices intersect -> success = false
		// test each projected edge of the edge ring against each other edge
		for( std::vector<Edge *>::iterator it1 = boundaryEdges.begin(); it1 != boundaryEdges.end() - 1; ++it1 )
		{
			for( std::vector<Edge *>::iterator it2 = it1+1; it2 != boundaryEdges.end(); ++it2 )
			{
				Edge *e1 = *it1;
				Edge *e2 = *it2;
				math::Vec3f intersectionPoint;

				// skip if the 2 lines share a common vertex
				if( e2->contains(e1->v1) || e2->contains(e1->v2) )
					continue;

				math::Vec2f e1_v1_projected = boundaryRing[e1->v1].projected;
				math::Vec2f e1_v2_projected = boundaryRing[e1->v2].projected;
				math::Vec2f e2_v1_projected = boundaryRing[e2->v1].projected;
				math::Vec2f e2_v2_projected = boundaryRing[e2->v2].projected;

				CGAL::Segment_2<K> line_1( Point( e1_v1_projected.x, e1_v1_projected.y ), Point(e1_v2_projected.x, e1_v2_projected.y) );
				CGAL::Segment_2<K> line_2( Point( e2_v1_projected.x, e2_v1_projected.y ), Point(e2_v2_projected.x, e2_v2_projected.y) );

				if( CGAL::do_intersect( line_1, line_2 ) )
				{
					success = false;
					break;
				}
			}

			if( !success )
				break;
		}


		if( success )
		{
			//
			// now we do the topology consistency check #2 (see chapter 4.4 of the paper)
			// this is done by checking all critical edges whether they cross the interior of the
			// polygon boundary which is formed by the projected vertex ring vertices
			// to do this we create the 2 polygongs which can be made by using the edge as a divider
			// the line crosses the interior if both polgons have a counterclockwise orientation
			//
			for( std::vector<Edge *>::iterator it = criticalEdges.begin(); it != criticalEdges.end(); ++it )
			{
				Edge *criticalEdge = *it;


				// setup left polygon
				// start from v1 and go to v2
				Polygon leftPoly;

				prev = 0;
				current = criticalEdge->v1;
				do
				{
					leftPoly.push_back( current, boundaryRing[current].projected );
					prev = current;
					current = boundaryRing[current].next;

				}while( current != criticalEdge->v2 );
				leftPoly.push_back( current, boundaryRing[current].projected );


				// setup right polygon
				// start from v2 and go to v1
				Polygon rightPoly;

				prev = 0;
				current = criticalEdge->v2;
				do
				{
					rightPoly.push_back( current, boundaryRing[current].projected );
					prev = current;
					current = boundaryRing[current].next;

				}while( current != criticalEdge->v1 );
				rightPoly.push_back( current, boundaryRing[current].projected );


				// if orientations of both polygons are the same then the critical edge crosses the bounding polygon interior
				if( leftPoly.orientation() == rightPoly.orientation() )
				{
					printf( "info : unable to triangulate since critical edge(es) would cross the polygon interior - vertex is not removed\n" );

					success = false;
					break;
				}
			}
		}
	}while( (++currentTrial < trialCount)&&(!success) );

	
	if( !success )
	{
		// we dont remove the vertex since we didnt managed to find a planar
		// projection of the neighbourhood which would not destroy the topology
		printf( "info : unable to find planar projection for vertex remove and retriangulation - vertex is not removed\n" );
		//++g_iterations;
		return false;
	}

	boundaryPolygonOrientation = CGAL::orientation( Point_3( v->position ), Point_3( v->position + base1 ), Point_3( v->position + base2 ), Point_3( v->position + normal ) );



	// execute the result ------------------------------------------------------------

	// remove vertex v and the triangle ring around v
	removeVertex( v );


	// compute squared distances
	
	// compute the squared distances of all point pairs
	std::vector<DistanceHelper> sqDistances;  // sorted distances of each pair of points
	for( std::map<Vertex *, EdgeInfoHelper>::iterator it1 = boundaryRing.begin(); it1 != boundaryRing.end(); ++it1 )
		for( std::map<Vertex *, EdgeInfoHelper>::iterator it2 = it1; it2 != boundaryRing.end(); ++it2 )
		{
			Vertex *v1 = it1->first;
			Vertex *v2 = it2->first;

			if( v1 == v2 )
				continue;

			// we dont want to add pairs, which would build a boundary edge, since those
			// already exist -> simply check for next or prev vertex in boundaryRing
			if( (boundaryRing[v1].next == v2)||(boundaryRing[v1].prev == v2) )
				continue;


			sqDistances.push_back( DistanceHelper( it1->first, it2->first ) );
		}

	// sort
	std::sort( sqDistances.begin(), sqDistances.end() );


	//
	std::list<Polygon *>                                                polygons; // polygons which have to be triangulated
	std::vector<Edge *>        edges( boundaryEdges.begin(), boundaryEdges.end() ); // created edges

	// create polygon which is made by the boundary edges and put it on polygon list
	Polygon *boundaryPolygon = new Polygon();


	first = boundaryRing.begin()->first;
	current = first;
	prev = 0;

	do
	{
		boundaryPolygon->push_back( current, boundaryRing[current].projected );
		prev = current;
		current = boundaryRing[current].next;
	}while( current != first );

	if( boundaryPolygon->orientation() != boundaryPolygonOrientation )
		boundaryPolygon->flip();


	polygons.push_back( boundaryPolygon );

	// for each entry in the sqDistance list
	//     find the polygon on the polygon list which contains h->v1 and h->v2
	//     check if edge is outside the polygon by checking the orientations of the left and right sides (build left and right polygons)
	//     check if edge intersects with any other existing edge which dont have a point in h->v1 or h->v2
	//     if checks pass:
	//         create edge -> assign it to the left and right polygons | add it to the list of edges
	//         split the polygon by removing the polygon from polygon list and putting the left and right polygons on the list if it is not a triangle
	//     
	// if the any polygon remains on the polygon list, throw an error that the triangulation has created a hole in the mesh
	//std::vector<Edge *> edges( edgeRing.begin(), edgeRing.end() );

	// when the mother polygon is already a triangle, then we have sqDistances.size() == 0
	// and we have to handle the polygon
	if( boundaryPolygon->isTriangle() )
	{
		createTriangle( boundaryPolygon->vertices[0], boundaryPolygon->vertices[1], boundaryPolygon->vertices[2], edges );
		polygons.clear();
		delete boundaryPolygon;

		//++g_iterations;
		return true;

	}


	for( std::vector<DistanceHelper>::iterator it=sqDistances.begin(); it != sqDistances.end(); ++it )
	{
		DistanceHelper *h = &(*it);

		Polygon *poly = 0; // polygon which contains both vertices of h

		std::list<Polygon *>::iterator pit;

		// find the polygon which contains both vertices of h
		for( pit = polygons.begin(); pit != polygons.end(); ++pit )
		{
			Polygon *p = *pit;
			if( p->contains( h->v1 ) && p->contains( h->v2 ) )
			{
				poly = p;
				break;
			}
		}

		// if no polygon could be found which contains both vertices of h, then
		// we can skip this one
		if( !poly )
			continue;

		// test for intersection with all existing edges
		bool intersection = false;
		for( std::vector<Edge *>::iterator eit = edges.begin(); eit != edges.end(); ++eit)
		{
			Edge *e = (*eit);

			math::Vec3f intersectionPoint;

			// ---- test for intersection of the edges projected into 2d plane using cgal ----
			// skip if the 2 lines share a common vertex
			if( e->contains(h->v1) || e->contains(h->v2) )
				continue;

			CGAL::Segment_2<K> line_1( Point( boundaryRing[e->v1].projected.x, boundaryRing[e->v1].projected.y ), Point(boundaryRing[e->v2].projected.x, boundaryRing[e->v2].projected.y) );
			CGAL::Segment_2<K> line_2( Point( boundaryRing[h->v1].projected.x, boundaryRing[h->v1].projected.y ), Point(boundaryRing[h->v2].projected.x, boundaryRing[h->v2].projected.y) );

			if( CGAL::do_intersect( line_1, line_2 ) )
			{
				intersection = true;
				break;
			}
		}

		// if there was an intersection, then we can skip this one
		// since a non-compatible edge would be introduced
		if( intersection )
			continue;

		// no intersection, get the 2 polygons which would be created from dividing the mother polyon along the edge
		Polygon *left, *right;

		left = right = 0;

		poly->split( left, right, h->v1, h->v2 );

		// if the 2 polygons have the same orientation, then we can be sure that
		// the edge goes through the interior of the polygon

		if( left->orientation() == right->orientation() )
		{

			// remove the current polygon from the list of polygons
			polygons.erase( pit );
			delete poly;

			// create the edge h->v1 -> h->v2
			Edge *edge = createEdge( h->v1, h->v2 );

			edges.push_back( edge );

			// and add the left and righ polygon to the list if they are not triangles
			if( left->isTriangle() )
			{
				// add triangle to the mesh
				createTriangle( left->vertices[0], left->vertices[1], left->vertices[2], edges );

				// remove polygon helper structure
				delete left;
			}else
				polygons.push_back( left );

			if( right->isTriangle() )
			{
				// add triangle to the mesh
				createTriangle( right->vertices[0], right->vertices[1], right->vertices[2], edges );

				// remove polygon helper structure
				delete right;
			}else
				polygons.push_back( right );
		}else
		{
			delete left;
			delete right;
		}
	}

	if( !polygons.empty() )
		printf( "error : hole created during retriangulation\n" );

	for( std::list<Polygon *>::iterator pit = polygons.begin(); pit != polygons.end(); ++pit )
		delete *pit;
	polygons.clear();


	//++g_iterations;
	return true;
}

void findOutsideSegment( Triangulation &t, Face_handle fh, int currentSegment, int commingFromIndex, int &outsideSegment )
{

	// if the face is an infinite face
	// then we know the outside segment which is the
	// current segment
	if( t.is_infinite(fh) )
	{
		if( (outsideSegment != -1)&&(outsideSegment != currentSegment) )
			printf( "error : different outsideSegments detected during triangulation (edgeloop not closed?)\n" );
		outsideSegment = currentSegment;
		return;
	}

	// if there is already a segment identifier, then we know that
	// the face has already been visited
	if( fh->info() != -1 )
		return;

	fh->info() = currentSegment;

	for( int i=0; i<3; ++i )
	{
		// get edge associated with the index i
		std::pair<Face_handle, int> edge = std::make_pair( fh, i );

		if( i == commingFromIndex )
			continue;
		

		// if the edge is a constrained edge, then we know this is the border
		if(t.is_constrained(edge))
		{
			//...and we have to pass a madified currentSegment value
			findOutsideSegment( t, fh->neighbor(i), (currentSegment+1)%2, t.mirror_index( fh, i), outsideSegment );
		}else
			//...else we recurse and leave the currentSegment value untouched
			findOutsideSegment( t, fh->neighbor(i), currentSegment, t.mirror_index( fh, i), outsideSegment );
	}
}

//
// Retriangulates a hole within the mesh. The hole is specified through an edgeloop(closed sequence of edges).
// In addition an optional number of points can be specified which will be included in the triangulation.
//
void MeshEx::retriangulateHole( std::vector<MeshEx::Edge *> &boundaryEdges, std::map<MeshEx::Vertex *, math::Vec2f> &boundaryVertexProjections, std::vector<std::pair<math::Vec3f, math::Vec2f> > &interiorPoints )
{
	std::map<Vertex_handle, MeshEx::Vertex*>                                 vertexMap; // used to map cgal vertex_handles to vertices
	std::vector<MeshEx::Edge *>                                                  edges; // this vector will hold all edges which were involved (for faster edge search)


	// algorithm:
	// - prepare data
	//		- find all boundary vertices
	// - prepare CGAL constrained triangulation
	//		- insert boundary vertices into triangulation and build mapping from Triangulation vertices to MeshEx::Vertices
	//		- use the boundary edges as constrained edges
	//		- insert points into triangulation from interiorPoints and build mapping from Triangulation vertices to MeshEx::Vertices
	// - extract triangulation results
	//		- ?


	// prepare algorithm ----------------------------------------------------------

	/*
	// obsolete since we get the boundary vertices with the boundaryVertexProjections
	// find boundary vertices
	for( std::vector<MeshEx::Edge *>::iterator it = boundaryEdges.begin(); it != boundaryEdges.end(); ++it )
	{
		MeshEx::Edge *e = *it;
		boundaryVertices.push_back( e->v1 );
		boundaryVertices.push_back( e->v2 );
	}
	// remove duplicate entries
	std::sort( boundaryVertices.begin(), boundaryVertices.end() );
	boundaryVertices.erase( std::unique( boundaryVertices.begin(), boundaryVertices.end() ), boundaryVertices.end() );
	*/
	
	// algorithm ------------------------------------------------------------------
	Triangulation t;

	// constrain triangulation with the boundary edges

	// iterate over all boundary vertices
	for( std::map<MeshEx::Vertex *, math::Vec2f>::iterator it = boundaryVertexProjections.begin(); it != boundaryVertexProjections.end(); ++it )
	{
		MeshEx::Vertex *v = it->first;

		// add boundary vertex to triangulation
		Vertex_handle vh = t.insert( Point( it->second.x, it->second.y ) );

		// we dont need to create the vertex
		vertexMap[vh] = v;
	}


	// iterate over all boundary edges
	for( std::vector<MeshEx::Edge *>::iterator it = boundaryEdges.begin(); it != boundaryEdges.end(); ++it )
	{
		MeshEx::Edge *e = *it;

		Vertex_handle v1, v2;


		bool v1_found  = false;
		bool v2_found  = false;
		// find vertex_handles for the given edge vertices
		for( std::map<Vertex_handle, MeshEx::Vertex*>::iterator vmit = vertexMap.begin(); vmit != vertexMap.end(); ++vmit )
		{
			if( e->v1 == vmit->second )
			{
				v1 = vmit->first;
				v1_found = true;
			}
			if( e->v2 == vmit->second )
			{
				v2 = vmit->first;
				v2_found = true;
			}
		}

		// add constrainedge to the triangulation
		t.insert_constraint( v1, v2 );

		// add edge to the list of created/existing edges
		edges.push_back( e );
	}

	// add additional and optional interior points
	for( std::vector<std::pair<math::Vec3f, math::Vec2f> >::iterator it = interiorPoints.begin(); it != interiorPoints.end(); ++it )
	{
		// update triangulation
		Vertex_handle v = t.insert( Point( it->second.x, it->second.y ) );

		// insertion may return a vertex which already exists (when the position is the same)
		if( vertexMap.find( v ) == vertexMap.end() )
			// create according MeshEx::Vertex and keep mapping to the CGAL vertices
			vertexMap[v] = createVertex( it->first );
	}


	// extract results and create triangles ----------------------------------------

	// now we have the triangulation of the convex hull of the whole problem, now we
	// have to find the faces which are inside the polygon - we mark each face with a
	// segment(inside or outside) property and by finding a face which is adjacent to
	// a infinite face, we find the segment which is outside


	// we employ some floodfilling scheme
	for( Triangulation::Finite_faces_iterator it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it )
		// reset info to -1
		it->info() = -1;

	int outsideSegment = -1;
	findOutsideSegment( t, t.finite_faces_begin(), 0, -1, outsideSegment );

	if( outsideSegment == -1 )
		printf( "error : outsideSegment not found during triangulation\n" );

	for( Triangulation::Finite_faces_iterator it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it )
		if( (it->info() == -1) && (!t.is_infinite(it)) )
			printf( "triangle not touched!\n" );

	// iterate over all faces of the triangulation and create edges/triangles
	for( Triangulation::Finite_faces_iterator it = t.finite_faces_begin(); it != t.finite_faces_end(); ++it )
	{
		Face_handle fh = it;

		// we are only interested in interior triangles
		if( fh->info() == outsideSegment )
			continue;

		MeshEx::Vertex *v0, *v1, *v2;

		v0 = vertexMap[ fh->vertex(0) ];
		v1 = vertexMap[ fh->vertex(1) ];
		v2 = vertexMap[ fh->vertex(2) ];

		MeshEx::Edge *e0, *e1, *e2;

		e0 = e1 = e2 = 0;

		// look for the edges in the edge vector
		for( std::vector<MeshEx::Edge*>::iterator eit = edges.begin(); eit != edges.end(); ++eit )
		{
			MeshEx::Edge *e = *eit;

			if( e->contains(v0) && e->contains(v1) )
				e0 = e;
			else
			if( e->contains(v1) && e->contains(v2) )
				e1 = e;
			else
			if( e->contains(v2) && e->contains(v0) )
				e2 = e;
		}

		// create the edges which could not be found
		if( !e0 )
		{
			e0 = createEdge( v0, v1 );
			edges.push_back(e0);
		}
		if( !e1 )
		{
			e1 = createEdge( v1, v2 );
			edges.push_back(e1);
		}
		if( !e2 )
		{
			e2 = createEdge( v2, v0 );
			edges.push_back(e2);
		}

		// create triangle
		MeshEx::Triangle *tri = createTriangle( v0, v1, v2, e0, e1, e2 );
	}
}

//
//
//
void MeshEx::detectAndFillHoles()
{
	int count = 0;
	bool done;
	do
	{
		done = true;

		// detect holes
		for( std::vector<Edge *>::iterator it = m_edges.begin(); it != m_edges.end(); ++it )
		{
			Edge *e = *it;

			if( e->isBoundaryEdge() )
			{
				Edge *e_start = e;

				std::vector<Edge *> boundaryEdges;


				int count = 0;
				do
				{
					boundaryEdges.push_back( e );

					Vertex *v2 = e->v2;
					e = 0;

					// get next boundary edge
					// find all edges going out from v2
					for( std::vector<Triangle *>::iterator trit = v2->triangleRing.begin(); trit != v2->triangleRing.end(); trit++ )
					{
						Triangle *t = *trit;

						for( size_t i = 0; i<3; ++i )
							if( (t->e[i] != boundaryEdges.back())&&(t->e[i]->contains(v2))&&(t->e[i]->isBoundaryEdge()) )
							{
								e = t->e[i];
								break;
							}

						if( e )
							break;
					}
				}while( e && (++count < 10) && (e != e_start) );

				for( std::vector<Edge *>::iterator eit = boundaryEdges.begin(); eit != boundaryEdges.end(); ++eit )
				{
					Edge *e = *eit;
					//vectors.push_back( std::make_pair( e->v1->position, e->v2->position ) );
				}

				// project boundary vertices on a 2d plane
				std::map<MeshEx::Vertex *, math::Vec2f> boundaryVertexProjections;

				for( std::vector<Edge *>::iterator eit = boundaryEdges.begin(); eit != boundaryEdges.end(); ++eit )
				{
					Edge *e = *eit;

					boundaryVertexProjections[e->v1] = math::Vec2f();
					boundaryVertexProjections[e->v2] = math::Vec2f();
				}

				math::Vec3f normal = boundaryVertexProjections.begin()->first->normal;
				float distance = -math::dotProduct( boundaryVertexProjections.begin()->first->position, normal );

				// compute the basis of the 2d-coordinate system of the plane defined by current normal
				math::Vec3f base1 = math::normalize( math::projectPointOnPlane( normal, distance, boundaryEdges[0]->v2->position - boundaryEdges[0]->v1->position ) );
				math::Vec3f base2 = math::normalize( math::crossProduct( normal, base1 ) );

				for( std::map<MeshEx::Vertex *, math::Vec2f>::iterator vit = boundaryVertexProjections.begin(); vit != boundaryVertexProjections.end(); ++vit )
				{
					// project neighbours into the given plane
					vit->second = _get2d( base1, base2, math::projectPointOnPlane( normal, distance, vit->first->position ) );
				}

				printf( "filling detected hole\n" );
				retriangulateHole( boundaryEdges, boundaryVertexProjections );
				done = false;
				break;
			}
		}
	}while( (++count < 10)&&(!done) );
}