/*---------------------------------------------------------------------

ObjIO offers import and export support for obj files.

Obj files are very simple ascii based mesh description files which store
vertices, normals and texturecoordinates of n-polygonal meshes.
The format was originally developed by AliasWavefront a long time ago.
Despite being a little out of age it is very well known and supported
throught the industry.

----------------------------------------------------------------------*/
#include "ObjIO.h"

namespace dk
{
	namespace io
	{
		//
		// reads and creates a Mesh from the given *.obj file. returns 0 if it fails
		//
		Mesh *importFromObjFile( std::string fileName )
		{
			// temporary vertex and index buffer
			std::vector<math::Vec3f> vertices;
			std::vector<int>         indicees;


			// try to open a text file stream
			std::ifstream file;
			file.open( fileName.c_str(), std::ios::in );

			// file open?
			if( !file )
				// error file could not be opened
				return 0;

			// now read the file
			while( file )
			{
				std::string line;
				std::getline( file, line, '\n' );

				// seperate commandlet
				std::string command = util::getWord( line, " " );

				if( command == "#" )
				{
					// comment
				}
				if( command == "mtllib" )
				{
					// material library location
				}
				if( command == "v" )
				{
					// vertex

					// read in a vertex which has format %f %f %f
					std::vector<std::string> components;
					util::splitString( line, components, " " );

					if( components.size() == 3 )
					{
						math::Vec3f position;
						std::istringstream( components[0] ) >> position.x;
						std::istringstream( components[1] ) >> position.y;
						std::istringstream( components[2] ) >> position.z;
						vertices.push_back( position );
					}
				}

				if( command == "f" )
				{
					// face

					// read in a face which has format %f %f %f
					std::vector<std::string> verts;
					util::splitString( line, verts, " " );

					// handle triangles
					if( verts.size() == 3 )
					{
						int v[3]; // vertex indicees
						int n[3]; // normal indicees
						int t[3]; // texture indicees

						// for each vertex
						for( unsigned int i=0; i<verts.size(); ++i )
						{
							// split the string again, since each component(vertexindex, normalindex, textureindex) is
							// seperated with a '/'
							std::vector<std::string> components;
							util::splitString( verts[i], components, "/" );

							// first component is the vertex index
							std::istringstream( components[0] ) >> v[i];

							// second component is the texture index

							// third component is the normal index

						}

						// we have to substract one since the obj file has one as first index not zero!
						indicees.push_back( v[0] - 1 );
						indicees.push_back( v[1] - 1 );
						indicees.push_back( v[2] - 1 );
					}else
					// handle quads
					if( verts.size() == 4 )
					{
						int v[4]; // vertex indicees
						int n[4]; // normal indicees
						int t[4]; // texture indicees

						// for each vertex
						for( unsigned int i=0; i<verts.size(); ++i )
						{
							// split the string again, since each component(vertexindex, normalindex, textureindex) is
							// seperated with a '/'
							std::vector<std::string> components;
							util::splitString( verts[i], components, "/" );

							// first component is the vertex index
							std::istringstream( components[0] ) >> v[i];

							// second component is the texture index

							// third component is the normal index

						}

						// we have to substract one since the obj file has one as first index not zero!
						// triangulate the quad on the fly
						indicees.push_back( v[0] - 1 );
						indicees.push_back( v[1] - 1 );
						indicees.push_back( v[2] - 1 );
						indicees.push_back( v[0] - 1 );
						indicees.push_back( v[2] - 1 );
						indicees.push_back( v[3] - 1 );
					}

				}

				//printf("%s\n", command.c_str() );
			}


			file.close();

			return new Mesh( vertices, indicees );
		}

		//
		// creates and writes a Obj file from the given mesh
		//
		void exportToObjFile( Mesh *mesh, std::string fileName )
		{
			mesh->computeVertexIndicees();

			// create temporary index buffer
			std::vector<int> indexBuffer;
			indexBuffer.resize( mesh->triangles.size() * 3 );

			int currentTriangleIndex = 0;
			for( std::vector<Mesh::Triangle *>::iterator it=mesh->triangles.begin(); it != mesh->triangles.end(); ++it, ++currentTriangleIndex )
			{
				indexBuffer[currentTriangleIndex*3 + 0] = (*it)->v0->index;
				indexBuffer[currentTriangleIndex*3 + 1] = (*it)->v1->index;
				indexBuffer[currentTriangleIndex*3 + 2] = (*it)->v2->index;
			}



			// create the file
			std::ofstream file;
			file.open( fileName.c_str(), std::ios::out | std::ios::trunc );

			// create the timestamp
			char time_buffer[40];
			tm _tm;
			size_t len;
			time_t now;
			
			now = time ( NULL );
			localtime_s( &_tm, &now );

			len = strftime ( time_buffer, 40, "%d %B %Y %I:%M:%S %p", &_tm );

			// info header
			file << "# ObjIO.cpp (exportToObjFile()) - Obj Export\n";
			file << "# File created: "<< time_buffer << "\n";
			file << "# Coordinates exported from a left-handed coordinate system.\n";
			file << "\n";

			// vertices
			file << "#begin " << mesh->vertices.size() << " vertices\n";
			for( std::vector<Mesh::Vertex *>::iterator it=mesh->vertices.begin(); it != mesh->vertices.end(); ++it )
				file << "v " << (*it)->position.x << " " << (*it)->position.y << " " << (*it)->position.z << "\n";
			file << "#end " << mesh->vertices.size() << " vertices\n";

			file << "\n";
		/*
			// normals
			file << "#begin " << mesh->normalBuffer.size() << " normals\n";
			for( std::vector<math::Vec3f>::iterator it=mesh->normalBuffer.begin(); it != mesh->normalBuffer.end(); ++it )
				file << "n " << (*it).x << " " << (*it).y << " " << (*it).z << "\n";
			file << "#end " << mesh->normalBuffer.size() << " normals\n";

			file << "\n";
		*/
			// faces
			file << "#begin " << mesh->triangles.size() << " faces\n";
			for( unsigned int i=0; i<mesh->triangles.size(); ++i )
			{
				int i1,i2,i3;

				// be aware that the indexing starts from 1 instead of 0 !
				i1 = indexBuffer[i*3+0] + 1;
				i2 = indexBuffer[i*3+1] + 1;
				i3 = indexBuffer[i*3+2] + 1;

				file << "f ";
				// first
				//file << i1 << "/" << "" << "/" << i1 << " ";
				file << i1 << "/" << "" << "/" << "" << " ";
				// second
				//file << i2 << "/" << "" << "/" << i2 << " ";
				file << i2 << "/" << "" << "/" << "" << " ";
				// third
				//file << i3 << "/" << "" << "/" << i3 << "\n";
				file << i3 << "/" << "" << "/" << "" << "\n";
			}
			file << "#end " << mesh->triangles.size() << " faces\n";

			file << "\n";

			// finish by closing the file
			file.close();
		}
	} // namespace io
} // namespace dk