// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
#pragma once
#include <map>
#include <vector>


#include "math/Math.h"

// CGAL -------------------------
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Filtered_kernel< CGAL::Simple_cartesian<double> > K;

typedef CGAL::Segment_3<K>                              Segment_3;
typedef CGAL::Triangle_3<K>                            Triangle_3;
CGAL::Point_3<K>                Point_3( const math::Vec3f &vec );

#include "Mesh.h"

using namespace dk;


///
/// \brief winged edge data structure
///
class MeshEx
{
public:
	struct Triangle;

	///
	/// \brief correspondeces to the vertices of the input mesh
	///
	struct Vertex
	{
		Vertex();                                       // constructor
		void                                               registerTriangle( Triangle *t );  // adds the given element to the elementring list
		void                                             unRegisterTriangle( Triangle *t );  // remvoes the given element to the elementring list
		bool                                                                  isDesolate();  // returns true if no element references this node
		int               getPlaneSide( const math::Vec3f &normal, const float &distance );  // returns 1 if the node lies on the left side of the given plane (nodeplanedistance < 0.0f) or 2 otherwise - returns 0 if point lies on the plane
		void                                                               computeNormal();  // computes normal from the surrounding tris - works only if the surrounding tris have a valid normal

		math::Vec3f                                                               position;
		math::Vec3f                                                                 normal;
		math::Color                                                                  color;

		std::vector<Triangle *>                                               triangleRing;  // references to all elements which index into this node
	};

	///
	/// \brief We use the winged edge datastructure to ease mesh edit operations
	///
	struct Edge
	{
		Edge( Vertex *_v1, Vertex *_v2 );                                                                            // constructor

		void                                                                       registerTriangle( Triangle *t );  // assigns the given element to the left or right wing of the edge (dependand on which wing is 0)
		void                                                                     unregisterTriangle( Triangle *t );  // sets either the left or right reference to zero when it references the given element
		bool   intersectsPlane( const math::Vec3f &normal, const float &distance, math::Vec3f &intersectionPoint );  // convinience function to compute the intersection of the edge with a given plane
		bool        intersectsLine( const math::Vec3f &p1, const math::Vec3f &p2, math::Vec3f &intersectionPoint );  // convinience function to compute the intersection of the edge with a given line
		bool                                        intersectsLine( const math::Vec3f &p1, const math::Vec3f &p2 );  // convinience function to compute the intersection of the edge with a given line without giving back the intersection point
		Vertex                                                                        *getOtherVertex( Vertex *v );  // convinience function to get the other node of the two nodes which are referenced by the edge
		Triangle                                                                  *getOtherTriangle( Triangle *t );  // returns the opposite triangle
		bool                                                                                          isDesolate();  // returns true if no elements references this edge
		bool                                                                                 contains( Vertex *v );  // returns true if the edge contains the given vertices
		bool                                                          contains( Vertex *vertex1, Vertex *vertex2 );  // returns true if the edge contains the given vertices
		bool                                                                               contains( Triangle *t );  // returns true if the edge contains the given triangles
		bool                                                                                      isBoundaryEdge();  // returns true if the edge has only one element reference instead of two

		Vertex                                                                                                 *v1;  // Edge is made up of 2 nodes
		Vertex                                                                                                 *v2;
		Triangle                                                                                             *left;  // left neighbour - if it is 0, the edge is a borderedge
		Triangle                                                                                            *right;  // right neighbour - if it is 0, the edge is a borderedge

		bool tag;
	};

	///
	/// \brief correspondences to the triangles of the input mesh
	///
	struct Triangle
	{
		union
		{
			struct
			{
				Vertex *v0, *v1, *v2;
			};
			Vertex *v[3];
		};
		union
		{
			struct
			{
				Edge *e1, *e2, *e3;
			};
			Edge *e[3];
		};

		Triangle();                                                                                          // constructor
		math::Vec3f  convertFromLocalToGlobalSpace( const math::Vec2f &localSpace, bool isVector = false ); // transforms the given coordinate in local 2d space into the global 3d space in which the triangle/element lies
		math::Vec2f convertFromGlobalToLocalSpace( const math::Vec3f &globalSpace, bool isVector = false ); // inverse operation of the above
		void                                                                           computeProperties(); // computes normal, area, etc...
		Vertex                                         *getOtherVertex( Vertex *vertex1, Vertex *vertex2 ); // returns the third node which neither equals n1 nor n2
		Vertex                                                               *getOtherVertex( Edge *edge ); // returns the third node which neither equals n1 nor n2
		bool                                                              contains( Vertex *vertex ) const; // returns true if one of the triangles vertices equals the given vertex
		size_t                                                      getVertexsLocalIndex( Vertex *vertex ); // returns the local index of the node within the element
		Edge                                                                 *getSharedEdge( Triangle *t ); // returns the edge which is shared with the given triangle or 0 if they dont share an edge
		Edge                                                  *getEdge( Vertex *vertex1, Vertex *vertex2 ); // returns the edge of the triangle which contains the 2 given vertices
		Edge                                                                    *getOtherEdge( Vertex *v ); // returns the one edge of the triangle which does not contain the vertex v

		void                                                                               computeNormal(); // (re) computes normal

		math::Vec3f                                                                                 center; // barycentric center position of the triangle
		math::Vec3f                                                                                 normal;
		math::Vec3f                                                                                      u; // these 2 vectors (u and v) define the local 2dimensional coordinate
		math::Vec3f                                                                                    v_v; // frame of the element (although they are 3d vectors)

		math::Matrix33f                                                                               beta; // barycentric coordinate matrix -> base vectors are the position of the 3 nodes of the element in the local 2d cooridnate frame
																											// so beta transforms local coordinates into barycentric coordinates

		float                                                                                         area; // the area of the element

		bool tag;
		bool tag2;
	};


	MeshEx( Mesh *mesh );                                                                             ///< constructor which takes a mesh as input
	~MeshEx();                                                                                        ///< destructor
	Mesh                                                                                  *getMesh(); ///< creates a mesh which represents the MeshEx


	Vertex                                              *createVertex( const math::Vec3f &position ); ///< creates a Vertex and adds it to the node list with the given world position
	void                                                              removeVertex( Vertex *vertex ); ///< removes given node from the node list
	Triangle     *createTriangle( const size_t &index0, const size_t &index1, const size_t &index2 ); ///< creates a Triangle and adds it to the element list with the given node indices
	Triangle                                   *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2 ); ///< creates a Triangle and adds it to the triangle list from given vertices
	Triangle     *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2, Edge *e0, Edge *e1, Edge *e2 ); ///< creates a Triangle and adds it to the triangle list from given edges
	Triangle       *createTriangle( Vertex *v0, Vertex *v1, Vertex *v2, std::vector<Edge *> &edges ); ///< creates a Triangle and adds it to the triangle list -> edges are looked up within the given edges (for speedup)
	Edge                                                       *createEdge( Vertex *v1, Vertex *v2 ); ///< creates a new edge from 2 given m_vertices
	Edge                                                         *findEdge( Vertex *n1, Vertex *n2 ); ///< looks for an edge between the given nodes - returns 0 if it could not be found
	Vertex *                               splitEdge( Edge *edge, const math::Vec3f &splitPosition ); ///< splits the given edge into 2 and returns the splitnode which was created to seperate the edge
	void                                                                    removeEdge( Edge *edge ); ///< removes edge from the list of edges
	void    removeTriangle( Triangle *element, bool removeVertices = true, bool removeEdges = true ); ///< removes element from the list of elements
	void     replaceVertexWithinTriangle( Triangle *e, Vertex *originalVertex, Vertex *replacement ); ///< replaces all occurances of the original node within element with the replacementnode and changes the edges accordingly
	void                  reTriangulate( Triangle *t, const std::vector<math::Vec3f> &pointsWithin ); ///< retriangulates the given triangle taking the given points into acount (which have to be coplanar with the triangle)
	bool                                      removeVertexAndReTriangulateNeighbourhood( Vertex *v ); ///< removes the given vertex and retriangulate the hole which would be made
	void           createTrianglesFromEdges( std::vector<Edge *> &edges, const math::Vec3f &normal ); ///< this method takes a list of already created edges and creates triangles for filling triangulated areas

	bool findClosestIntersection( const math::Vec3f &position, const math::Vec3f &p1, const math::Vec3f &p2, math::Vec3f &intersection, math::Vec3f &normal, Triangle *&t );

	void retriangulateHole( std::vector<MeshEx::Edge *> &boundaryEdges, std::map<MeshEx::Vertex *, math::Vec2f> &boundaryVertexProjections, std::vector<std::pair<math::Vec3f, math::Vec2f> > &interiorPoints = std::vector<std::pair<math::Vec3f, math::Vec2f> >() );
	void                                                                        detectAndFillHoles();

	void                                                                                     clear();


	std::vector<Vertex*>                                                                  m_vertices;
	std::vector<Triangle*>                                                               m_triangles;
	std::vector<Edge*>                                                                       m_edges;

	bool     stop;
	math::Vec3f temp;
};



	///
	/// \brief holds a connection between 2 vertices (used by retiling algo)
	///
	struct helper
	{
		helper( MeshEx::Vertex *_v1, MeshEx::Vertex *_v2 ) : v1(_v1), v2(_v2)
		{
			sqDistance = (v2->position - v1->position).getSquaredLength();
		}
		float   sqDistance;
		MeshEx::Vertex *v1;
		MeshEx::Vertex *v2;

		bool operator<( const helper &rhs )
		{
			return sqDistance < rhs.sqDistance;
		}
	};

	///
	/// \brief This structure is another helper which is used to store the incoming and outgoing
	/// edge of an vertex of the vertex ring.
	///
	struct VertexRingEdges
	{
		VertexRingEdges() : m_e1(0), m_e2(0), triangleCount(0)
		{
		}



		/// \brief add edge
		/// Every ring-vertex has an incomming and outgoing edge. The VertexRingEdges structure is only
		/// interested in the vectors going out into the direction of the edges...
		///
		/// Registering the vectors isnt enough - to work properly, the normal of the plane has to be
		/// passed using setNormal
		///
		void registerEdgeVector( MeshEx::Edge *e, const math::Vec3f &vec )
		{
			if( !m_e1 )
			{
				m_vec1 = vec;
				m_e1 = e;
			}else
			if( !m_e2 )
			{
				m_vec2 = vec;
				m_e2 = e;
			}else
				printf( "error : There must NOT be more then 2 edges from the edgering\n" );
		}

		///
		/// Computes the clockwise angle between the given vector and the vector going out of m_v
		/// into the direction of m_e1
		///
		float computeClockwiseAngle( const math::Vec3f vec )
		{
			math::Vec3f cp = math::crossProduct( m_vec1, vec );
			float dp = math::dotProduct( m_vec1, vec );

			// if the crossproduct is zero, then m_vec1 and vec are colinear
			// which means the angle between those 2 vectors is either 0.0f or 180.0f
			// we can use the projection of the crossproduct onto the normal to decide
			// since the crossproduct is zero
			// but we know the 2 vecs are colinear, so we take the dotproduct between those 2
			// vectors to see whether they point in the same direction or not
			/*
			if( cp.getSquaredLength() < 0.000001f )
				if( dp < 0.0f )
					return MATH_PIf;
				else
					return 0.0f;
			*/

			// clamp the dotproduct to [-1,1] -> acos will be not happy otherwise
			if( dp < -1.0f )
				dp = -1.0f;
			if( dp > 1.0f )
				dp = 1.0f;

			float angle = acosf( dp );


			// if the crossproduct of the 2 vectors points in the same direction as the cross product
			if( math::dotProduct( m_normal, cp ) > 0.0f )
				// then the clockwise angle between these 2 vectors is just the dotproduct
				return angle;
			else
			{
				// if the crossproduct points in the opposite direction, then the angle
				// is 360° - angle
				//return math::degToRad( 360.0f ) - angle;
				//float t = math::degToRad( 360.0f ) - angle;
				float t = 6.283185307179586f - angle;
				//if( t >= 6.283185307179586f )
				//	t -= 6.283185307179586f;
				return t;
			}

			// we shouldnt come here
			return angle;
		}

		void setNormal( const math::Vec3f &normal )
		{
			m_normal = normal;

			// compute clockwise angle between the 2 edges
			m_cwAngle = computeClockwiseAngle( m_vec2 );
		}

		float getAngle()
		{
			return m_cwAngle;
		}

		void flipEdges()
		{
			std::swap( m_e1, m_e2 );
			std::swap( m_vec1, m_vec2 );
			// recompute clockwise angle between the 2 edge-vectors
			m_cwAngle = computeClockwiseAngle( m_vec2 );
		}

		math::Vec3f             m_normal;
		MeshEx::Edge               *m_e1;
		MeshEx::Edge               *m_e2;
		math::Vec3f               m_vec1;
		math::Vec3f               m_vec2;

		float                  m_cwAngle;

		size_t             triangleCount; // temp for debugging
	};
	
	
	
	