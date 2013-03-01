// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------



----------------------------------------------------------------------*/
#include "MeshEx.h"
#include <list>

//#include <CGAL/Polygon_2.h>
//#include <CGAL/Triangulation_2.h>
//#include <CGAL/Triangulation_face_base_with_info_2.h>



/*
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
typedef CGAL::Polygon_2<K, std::vector<Point> >                                                     Polygon_2;
*/

/*
math::Vec2f _get2d( math::Vec3f e1,  math::Vec3f e2, math::Vec3f point )
{
	return math::Vec2f( math::dotProduct( e1, point ), math::dotProduct( e2, point ) );
}
*/









//
//
//
bool MeshEx::findClosestIntersection( const math::Vec3f &position, const math::Vec3f &p1, const math::Vec3f &p2, math::Vec3f &intersection, math::Vec3f &normal, Triangle *&t )
{
	int count = 0;

	struct IntersectionHelper
	{
		IntersectionHelper( math::Vec3f intersection, Triangle *triangle, float sqDistance ) : m_intersection(intersection), m_triangle(triangle), m_sqDistance(sqDistance)
		{
		}

		// used for sorting
		bool operator<( const IntersectionHelper &rhs )
		{
			return m_sqDistance < rhs.m_sqDistance;
		}

		math::Vec3f m_intersection; // point of intersection
		Triangle       *m_triangle; // triangle in which the intersection point lies
		float         m_sqDistance; // distance to the near point
	};

	std::vector<IntersectionHelper> intersections;

	// test every triangle (TODO: do some spatial subdivision)
	for( std::vector<Triangle *>::iterator it = m_triangles.begin(); it != m_triangles.end(); ++it )
	{
		Triangle *tri = *it;

		Triangle_3 cgalTri( Point_3(tri->v0->position), Point_3(tri->v1->position), Point_3(tri->v2->position) );
		Segment_3 seg( Point_3(p1), Point_3(p2) );

		if( CGAL::do_intersect( cgalTri, seg ) )
		{
			CGAL::Point_3<K> tempIntersection;
			CGAL::Object o = CGAL::intersection(CGAL::Plane_3<K>( Point_3(tri->v0->position), Point_3(tri->v1->position), Point_3(tri->v2->position) ), seg);
			if( CGAL::assign(tempIntersection, o) )
			{
				float sd = CGAL::squared_distance( Point_3(position), tempIntersection );
				// intersection results in a point
				intersections.push_back( IntersectionHelper( math::Vec3f(tempIntersection.x(), tempIntersection.y(), tempIntersection.z()), tri, sd ) );
			}else
			{
				// segment!
				printf( "error : intersection during vertex movement is a line\n" );
			}
			//printf( "intersection\n" );
		}
	}



	// if we have found any intersection
	if( intersections.size() )
	{
		std::sort( intersections.begin(), intersections.end() );

		// find the closest one
		intersection = intersections[0].m_intersection;
		t = intersections[0].m_triangle;
		normal = t->normal;
		// and indicate success
		return true;
	}






	return false;
}