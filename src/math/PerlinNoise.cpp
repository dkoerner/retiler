//
//
// PerlinNoise generator (1d,2d,3d,4d)
//
//
#include "PerlinNoise.h"
#include <stdlib.h>
#include <math.h>

namespace math
{
	
	unsigned char PerlinNoise::g_permutationTable[256] =
	{
	225, 155, 210, 108, 175, 199, 221, 144, 203, 116, 70, 213, 69, 158, 33, 252,
	5, 82, 173, 133, 222, 139, 174, 27, 9, 71, 90, 246, 75, 130, 91, 191,
	169, 138, 2, 151, 194, 235, 81, 7, 25, 113, 228, 159, 205, 253, 134, 142,
	248, 65, 224, 217, 22, 121, 229, 63, 89, 103, 96, 104, 156, 17, 201, 129,
	36, 8, 165, 110, 237, 117, 231, 56, 132, 211, 152, 20, 181, 111, 239, 218,
	170, 163, 51, 172, 157, 47, 80, 212, 176, 250, 87, 49, 99, 242, 136, 189,
	162, 115, 44, 43, 124, 94, 150, 16, 141, 247, 32, 10, 198, 223, 255, 72,
	53, 131, 84, 57, 220, 197, 58, 50, 208, 11, 241, 28, 3, 192, 62, 202,
	18, 215, 153, 24, 76, 41, 15, 179, 39, 46, 55, 6, 128, 167, 23, 188,
	106, 34, 187, 140, 164, 73, 112, 182, 244, 195, 227, 13, 35, 77, 196, 185,
	26, 200, 226, 119, 31, 123, 168, 125, 249, 68, 183, 230, 177, 135, 160, 180,
	12, 1, 243, 148, 102, 166, 38, 238, 251, 37, 240, 126, 64, 74, 161, 40,
	184, 149, 171, 178, 101, 66, 29, 59, 146, 61, 254, 107, 42, 86, 154, 4,
	236, 232, 120, 21, 233, 209, 45, 98, 193, 114, 78, 19, 206, 14, 118, 127,
	48, 79, 147, 85, 30, 207, 219, 54, 88, 234, 190, 122, 95, 67, 143, 109,
	137, 214, 145, 93, 92, 100, 245, 0, 216, 186, 60, 83, 105, 97, 204, 52
	};

	float PerlinNoise::g_gradientTable[768];

	bool      PerlinNoise::g_tablesInitiated = false;

	//
	// constructor - initializes perlin noise parameters here
	//
	PerlinNoise::PerlinNoise()
	{
		m_amplitude = 0.2f;
        m_amplitudeRatio = .65f;
        m_frequency      = .05f;
		m_frequencyRatio = 2.0f;
        m_octaves = 8;
		m_inflection = false;

        if( !g_tablesInitiated )
			initGradientTable(357345);
	}

	//
	// 2D perlin noise function
	//
	float PerlinNoise::perlinNoise_2D( float s, float t )
	{
		int i1;

		float amplitudeFactor = m_amplitude;
		float frequencyFactor = m_frequency;

		float result = 0.0f;

		if( m_inflection )
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				result += fabs( interpolatedGradientNoise( s * frequencyFactor, t * frequencyFactor ) * amplitudeFactor );

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}

			return result * 2.0f - 1.0f;
		}else
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				result += interpolatedGradientNoise( s * frequencyFactor, t * frequencyFactor ) * amplitudeFactor;

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}

			return result;
		}

	}
	//
	// 3D perlin noise function
	//
	float PerlinNoise::perlinNoise_3D( float u, float v, float w  )
	{
		int i1;

		float amplitudeFactor = m_amplitude;
		float frequencyFactor = m_frequency;

		float result = 0.0f;



		if( m_inflection )
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				//result += interpolatedGradientNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor ) * amplitudeFactor;
				result += fabs(interpolatedGradientNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor ) * amplitudeFactor);

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}

			return result*2.0f - 1.0f;
		}else
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				//result += interpolatedNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor ) * amplitudeFactor;
				result += interpolatedGradientNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor ) * amplitudeFactor;

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}

			return result;
		}
	}
	//
	// 4D perlin noise function
	//
	float PerlinNoise::perlinNoise_4D( float u, float v, float w, float x  )
	{
		int i1;

		float amplitudeFactor = m_amplitude;
		float frequencyFactor = m_frequency;

		float result = 0.0f;

		if( m_inflection )
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				//result += interpolatedNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor, x * frequencyFactor ) * amplitudeFactor;
				result += fabs( interpolatedGradientNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor, x * frequencyFactor ) * amplitudeFactor );

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}
			return result*2.0f - 1.0f;
		}else
		{
			// perlin noise
			for( i1=0; i1<m_octaves;i1++ )
			{
				//result += interpolatedNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor, x * frequencyFactor ) * amplitudeFactor;
				result += interpolatedGradientNoise( u * frequencyFactor, v * frequencyFactor, w * frequencyFactor, x * frequencyFactor ) * amplitudeFactor;

				amplitudeFactor *= m_amplitudeRatio;
				frequencyFactor *= m_frequencyRatio;
			}
			return result;
		}


	}


	float PerlinNoise::getAmplitude( void )
	{
		return m_amplitude;
	}

	void PerlinNoise::setAmplitude( float amplitude )
	{
		m_amplitude = amplitude;
	}

	float PerlinNoise::getAmplitudeRatio( void )
	{
		return m_amplitudeRatio;
	}

	void PerlinNoise::setAmplitudeRatio( float amplitudeRatio )
	{
		m_amplitudeRatio = amplitudeRatio;
	}

	float PerlinNoise::getFrequency( void )
	{
		return m_frequency;
	}

	void PerlinNoise::setFrequency( float frequency )
	{
		m_frequency = frequency;
	}

	float PerlinNoise::getFrequencyRatio( void )
	{
		return m_frequencyRatio;
	}

	void PerlinNoise::setFrequencyRatio( float frequencyRatio )
	{
		m_frequencyRatio = frequencyRatio;
	}

	int PerlinNoise::getDepth( void )
	{
		return m_octaves;
	}
	void PerlinNoise::setDepth( int depth )
	{
		m_octaves = depth;
	}


	bool PerlinNoise::getInflection( void )
	{
		return m_inflection;
	}
	void PerlinNoise::setInflection( bool inflection )
	{
		m_inflection = inflection;
	}


	void PerlinNoise::initGradientTable( int seed )
	{
		g_tablesInitiated = true;

		float *f = g_gradientTable;

		float z,r,theta;
		const float contrast = 1.5f;

		srand( seed );

		for( unsigned int i=0; i<256; ++i )
		{
			z = (float(rand())/float(RAND_MAX)) * 2.0f - 1.0f;
			r = sqrt( 1.0f - z * z );
			theta = float(6.283185307179586476925286766559 * (float(rand()) / float(RAND_MAX)));
			*f++ = r * cos(theta) * contrast;
            *f++ = r * sin(theta) * contrast;
            *f++ = z * contrast;
		}
	}

	//
	// this is the noise interpolation function for 2dimensional perlin noise
	//
    float PerlinNoise::interpolatedNoise( float s, float t )
	{
		int x = (int)floor( s );
		float fracX = s - x;
		int y = (int)floor( t );
		float fracY = t - y;

		float v1,v2,v3,v4;

		// get the noise values at floor(s) and floor(s) + 1 and for t respectively
        v1 = noise(x, y);
		v2 = noise(x+1, y);
		v3 = noise(x,y+1);
		v4 = noise(x+1,y+1);

		// now perform interpolation between the for noise values dependent on the fractional parts of s and t
		// 2D -> bilinear

		// interpoliere die erste zeile (y)
		float f1 = interpolateEaseCurve( v1, v2, fracX );
		// interpoliere die zweite zeile (y)
		float f2 = interpolateEaseCurve( v3, v4, fracX );

		// interpoliere zwischen den beiden interpolationswerten beider x-zeilen in y
		return interpolateEaseCurve( f1, f2, fracY );
	}
	//
	// this is the noise interpolation function for 3dimensional perlin noise
	//
    float PerlinNoise::interpolatedNoise( float u, float v, float w )
	{
		int x = (int)floor( u );
		float fracX = u - x;
		int y = (int)floor( v );
		float fracY = v - y;
		int z = (int)floor( w );
		float fracZ = w - z;

		float v1,v2,v3,v4,v5,v6,v7,v8;


		// get the noise values at floor(u) and floor(u) + 1 and for t respectively
        v1 = noise( x,y,z );
		v2 = noise( (x+1),y,z );
		v3 = noise( x , (y+1),z );
		v4 = noise( (x+1),(y+1),z );
        v5 = noise( x,y,(z+1) );
		v6 = noise( (x+1),y,(z+1) );
		v7 = noise( x , (y+1),(z+1) );
		v8 = noise( (x+1),(y+1),(z+1) );


		// now perform interpolation between the for noise values dependent on the fractional parts of u,v and w
		// 3d->trilinear

		// interpoliere die erste zeile (y) bei z
		float f1 = interpolateEaseCurve( v1, v2, fracX );
		// interpoliere die zweite zeile (y+1) bei z
		float f2 = interpolateEaseCurve( v3, v4, fracX );

		// interpoliere die erste zeile (y) bei z+1
		float f3 = interpolateEaseCurve( v5, v6, fracX );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f4 = interpolateEaseCurve( v7, v8, fracX );

		// interpoliere die vorderen eckpunte bei z in y
		float f5 = interpolateEaseCurve( f1, f2, fracY );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f6 = interpolateEaseCurve( f3, f4, fracY );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		return interpolateEaseCurve( f5, f6, fracZ );
	}

	float PerlinNoise::interpolatedGradientNoise( float u, float v )
	{
		int x = (int)floor( u );
		float fracX0 = u - x;
		float fracX1 = fracX0 - 1;
		int y = (int)floor( v );
		float fracY0 = v - y;
		float fracY1 = fracY0 - 1;

		float v1,v2,v3,v4;


		// get the noise values at floor(u) and floor(u) + 1 and for t respectively
        v1 = gradientNoise(     x,   y, fracX0, fracY0  );
		v2 = gradientNoise(   x+1,   y, fracX1, fracY0 );
		v3 = gradientNoise(     x, y+1, fracX0, fracY1 );
		v4 = gradientNoise(   x+1, y+1, fracX1, fracY1 );

		// now perform interpolation between the for noise values dependent on the fractional parts of u,v and w
		// 3d->trilinear

		// interpoliere die erste zeile (y) bei z
		float f1 = interpolateEaseCurve( v1, v2, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z
		float f2 = interpolateEaseCurve( v3, v4, fracX0 );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		return interpolateEaseCurve( f1, f2, fracY0 );
	}
	//
	// this is the noise interpolation function for 3dimensional perlin gradient noise
	//
    float PerlinNoise::interpolatedGradientNoise( float u, float v, float w )
	{
		int x = (int)floor( u );
		float fracX0 = u - x;
		float fracX1 = fracX0 - 1;
		int y = (int)floor( v );
		float fracY0 = v - y;
		float fracY1 = fracY0 - 1;
		int z = (int)floor( w );
		float fracZ0 = w - z;
		float fracZ1 = fracZ0 - 1;

		float v1,v2,v3,v4,v5,v6,v7,v8;


		// get the noise values at floor(u) and floor(u) + 1 and for t respectively
        v1 = gradientNoise(     x,   y,   z, fracX0, fracY0, fracZ0  );
		v2 = gradientNoise(   x+1,   y,   z, fracX1, fracY0, fracZ0 );
		v3 = gradientNoise(     x, y+1,   z, fracX0, fracY1, fracZ0 );
		v4 = gradientNoise(   x+1, y+1,   z, fracX1, fracY1, fracZ0 );
        v5 = gradientNoise(     x,   y, z+1, fracX0, fracY0, fracZ1 );
		v6 = gradientNoise(   x+1,   y, z+1, fracX1, fracY0, fracZ1 );
		v7 = gradientNoise(     x, y+1, z+1, fracX0, fracY1, fracZ1 );
		v8 = gradientNoise(   x+1, y+1, z+1, fracX1, fracY1, fracZ1 );

		// now perform interpolation between the for noise values dependent on the fractional parts of u,v and w
		// 3d->trilinear

		// interpoliere die erste zeile (y) bei z
		float f1 = interpolateEaseCurve( v1, v2, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z
		float f2 = interpolateEaseCurve( v3, v4, fracX0 );

		// interpoliere die erste zeile (y) bei z+1
		float f3 = interpolateEaseCurve( v5, v6, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f4 = interpolateEaseCurve( v7, v8, fracX0 );

		// interpoliere die vorderen eckpunte bei z in y
		float f5 = interpolateEaseCurve( f1, f2, fracY0 );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f6 = interpolateEaseCurve( f3, f4, fracY0 );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		return interpolateEaseCurve( f5, f6, fracZ0 );
	}

	//
	// this is the noise interpolation function for 4dimensional perlin noise
	//
    float PerlinNoise::interpolatedGradientNoise( float _u, float _v, float _w, float _s )
	{
		int x = (int)floor( _u );
		float fracX0 = _u - x;
		float fracX1 = fracX0 - 1;
		int y = (int)floor( _v );
		float fracY0 = _v - y;
		float fracY1 = fracY0 - 1;
		int z = (int)floor( _w );
		float fracZ0 = _w - z;
		float fracZ1 = fracZ0 - 1;
		int w = (int)floor( _s );
		float fracW0 = _s - w;
		float fracW1 = fracW0 - 1;

		float v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16;


		// get the noise values at floor(u) and floor(u) + 1 and for t respectively
        v1  = gradientNoise(     x,   y,   z,   w, fracX0, fracY0, fracZ0, fracW0 );
		v2  = gradientNoise(   x+1,   y,   z,   w, fracX1, fracY0, fracZ0, fracW0 );
		v3  = gradientNoise(     x, y+1,   z,   w, fracX0, fracY1, fracZ0, fracW0 );
		v4  = gradientNoise(   x+1, y+1,   z,   w, fracX1, fracY1, fracZ0, fracW0 );
        v5  = gradientNoise(     x,   y, z+1,   w, fracX0, fracY0, fracZ1, fracW0 );
		v6  = gradientNoise(   x+1,   y, z+1,   w, fracX1, fracY0, fracZ1, fracW0 );
		v7  = gradientNoise(     x, y+1, z+1,   w, fracX0, fracY1, fracZ1, fracW0 );
		v8  = gradientNoise(   x+1, y+1, z+1,   w, fracX1, fracY1, fracZ1, fracW0 );
        v9  = gradientNoise(     x,   y,   z, w+1, fracX0, fracY0, fracZ0, fracW1 );
		v10 = gradientNoise(   x+1,   y,   z, w+1, fracX1, fracY0, fracZ0, fracW1 );
		v11 = gradientNoise(     x, y+1,   z, w+1, fracX0, fracY1, fracZ0, fracW1 );
		v12 = gradientNoise(   x+1, y+1,   z, w+1, fracX1, fracY1, fracZ0, fracW1 );
        v13 = gradientNoise(     x,   y, z+1, w+1, fracX0, fracY0, fracZ1, fracW1 );
		v14 = gradientNoise(   x+1,   y, z+1, w+1, fracX1, fracY0, fracZ1, fracW1 );
		v15 = gradientNoise(     x, y+1, z+1, w+1, fracX0, fracY1, fracZ1, fracW1 );
		v16 = gradientNoise(   x+1, y+1, z+1, w+1, fracX1, fracY1, fracZ1, fracW1 );


		// now perform interpolation between the for noise values dependent on the fractional parts of u,v and w
		// 3d->trilinear

		// interpoliere die erste zeile (y) bei z
		float f1  = interpolateEaseCurve( v1, v2, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z
		float f2  = interpolateEaseCurve( v3, v4, fracX0 );

		// interpoliere die erste zeile (y) bei z+1
		float f3  = interpolateEaseCurve( v5, v6, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f4  = interpolateEaseCurve( v7, v8, fracX0 );

		// interpoliere die vorderen eckpunte bei z in y
		float f5  = interpolateEaseCurve( f1, f2, fracY0 );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f6  = interpolateEaseCurve( f3, f4, fracY0 );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		float f7  = interpolateEaseCurve( f5, f6, fracZ0 );

		// interpoliere die erste zeile (y) bei z bei w
		float f8  = interpolateEaseCurve( v9, v10, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z
		float f9  = interpolateEaseCurve( v11, v12, fracX0 );

		// interpoliere die erste zeile (y) bei z+1
		float f10 = interpolateEaseCurve( v13, v14, fracX0 );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f11 = interpolateEaseCurve( v15, v16, fracX0 );

		// interpoliere die vorderen eckpunte bei z in y
		float f12 = interpolateEaseCurve( f8, f9, fracY0 );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f13 = interpolateEaseCurve( f10, f11, fracY0 );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		float f14 = interpolateEaseCurve( f12, f13, fracZ0 );

		return interpolateEaseCurve( f7, f14, fracW0 );
	}

	//
	// this is the noise interpolation function for 4dimensional perlin noise
	//
    float PerlinNoise::interpolatedNoise( float _u, float _v, float _w, float _s )
	{
		int x = (int)floor( _u );
		float fracX = _u - x;
		int y = (int)floor( _v );
		float fracY = _v - y;
		int z = (int)floor( _w );
		float fracZ = _w - z;
		int w = (int)floor( _s );
		float fracW = _s - w;

		float v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16;


		// get the noise values at floor(u) and floor(u) + 1 and for t respectively
        v1  = noise(     x,   y,   z,   w );
		v2  = noise(   x+1,   y,   z,   w );
		v3  = noise(     x, y+1,   z,   w );
		v4  = noise(   x+1, y+1,   z,   w );
        v5  = noise(     x,   y, z+1,   w );
		v6  = noise(   x+1,   y, z+1,   w );
		v7  = noise(     x, y+1, z+1,   w );
		v8  = noise(   x+1, y+1, z+1,   w );
        v9  = noise(     x,   y,   z, w+1 );
		v10 = noise(   x+1,   y,   z, w+1 );
		v11 = noise(     x, y+1,   z, w+1 );
		v12 = noise(   x+1, y+1,   z, w+1 );
        v13 = noise(     x,   y, z+1, w+1 );
		v14 = noise(   x+1,   y, z+1, w+1 );
		v15 = noise(     x, y+1, z+1, w+1 );
		v16 = noise(   x+1, y+1, z+1, w+1 );


		// now perform interpolation between the for noise values dependent on the fractional parts of u,v and w
		// 3d->trilinear

		// interpoliere die erste zeile (y) bei z
		float f1  = interpolateEaseCurve( v1, v2, fracX );
		// interpoliere die zweite zeile (y+1) bei z
		float f2  = interpolateEaseCurve( v3, v4, fracX );

		// interpoliere die erste zeile (y) bei z+1
		float f3  = interpolateEaseCurve( v5, v6, fracX );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f4  = interpolateEaseCurve( v7, v8, fracX );

		// interpoliere die vorderen eckpunte bei z in y
		float f5  = interpolateEaseCurve( f1, f2, fracY );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f6  = interpolateEaseCurve( f3, f4, fracY );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		float f7  = interpolateEaseCurve( f5, f6, fracZ );

		// interpoliere die erste zeile (y) bei z bei w
		float f8  = interpolateEaseCurve( v9, v10, fracX );
		// interpoliere die zweite zeile (y+1) bei z
		float f9  = interpolateEaseCurve( v11, v12, fracX );

		// interpoliere die erste zeile (y) bei z+1
		float f10 = interpolateEaseCurve( v13, v14, fracX );
		// interpoliere die zweite zeile (y+1) bei z+1
		float f11 = interpolateEaseCurve( v15, v16, fracX );

		// interpoliere die vorderen eckpunte bei z in y
		float f12 = interpolateEaseCurve( f8, f9, fracY );
		// interpoliere die hinteren eckpunte bei z+1 in y
		float f13 = interpolateEaseCurve( f10, f11, fracY );

		// interpoliere zwischen den beiden interpolationswerten beider y-zeilen in z
		float f14 = interpolateEaseCurve( f12, f13, fracZ );

		return interpolateEaseCurve( f7, f14, fracW );
	}
	//
	// this is a noise function which will return a pseudorandom number
	// based on the value x
	//
	float PerlinNoise::noise( int x )
    {
		int n = (x<<13) ^ x;
        float res = ( 1.0f - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff ) / 1073741824.0f);
        return res;
    }

	float PerlinNoise::noise( int x, int y )
    {
        int n = x + y * 57;
        n = (n<<13) ^ n;
        return ( 1.0f - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff ) / 1073741824.0f);
    }

	float PerlinNoise::noise( int x, int y, int z )
    {
        //int n = x + y * 57 + z*123;
		int i = x + y*57;
		int m = y + z*57;

		unsigned int n = i + m*57;
		n = (n<<13) ^ n;
        return ( 1.0f - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff ) / 1073741824.0f);
    }
	float PerlinNoise::noise( int x, int y, int z, int w )
    {
        int n = x + y * 57 + z*123 + w*389;
        n = (n<<13) ^ n;
        return ( 1.0f - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff ) / 1073741824.0f);
    }

	float PerlinNoise::gradientNoise(int x, int y, float fx, float fy )
    {
		float* g = g_gradientTable + 3 * g_permutationTable[ (x + g_permutationTable[ y & 255 ]) & 255 ];
        return g[0] * fx + g[1] * fy;
    }

	float PerlinNoise::gradientNoise(int x, int y, int z, float fx, float fy, float fz)
    {
		float* g = g_gradientTable + 3 * g_permutationTable[ (x + g_permutationTable[ (y + g_permutationTable[ z & 255 ]) & 255 ]) & 255 ];
        return g[0] * fx + g[1] * fy + g[2] * fz;
    }

	float PerlinNoise::gradientNoise(int x, int y, int z, int w, float fx, float fy, float fz, float fw )
    {
		float* g = g_gradientTable + 3 * g_permutationTable[ (x + g_permutationTable[ (y + g_permutationTable[ (z + g_permutationTable[w & 255]) & 255 ]) & 255 ]) & 255 ];
        return g[0] * fx + g[1] * fy + g[2] * fz + g[0]*fw;
    }

	//
	// simple linear interpolation between two values dependent on a interpolation value
	//
	float PerlinNoise::interpolateLinear( float v1, float v2, float t )
	{
		return (float)(v1 * ( 1.0 - t ) + v2 * t);
	}

	//
	// will interpolate between two values in a more harmonic and smooth manner
	//
    float PerlinNoise::interpolateEaseCurve( float v1, float v2, float t )
    {
        float fac1 = 3*pow(1-t, 2) - 2*pow(1-t,3);
        float fac2 = 3*pow(t, 2) - 2*pow(t, 3);
        return v1*fac1 + v2*fac2; //add the weighted factors
    }
}





















