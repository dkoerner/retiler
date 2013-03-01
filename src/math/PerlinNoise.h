// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
//
//
// PerlinNoise generator (1d,2d,3d,4d)
//
//
#pragma once


namespace math
{
	//
	// this class is a simple utilityclass implementing perlin noise
	//
	class PerlinNoise
	{
	public:

		PerlinNoise();     // constructor - initializes perlin noise parameters here

		float                    perlinNoise_2D( float s, float t ); // 2D perlin noise function
		float          perlinNoise_3D( float u, float v, float w  ); // 3D perlin noise function
		float perlinNoise_4D( float u, float v, float w, float x  ); // 4D perlin noise function

		float                                        getAmplitude( void );
		void                              setAmplitude( float amplitude );
		float                                   getAmplitudeRatio( void );
		void                    setAmplitudeRatio( float amplitudeRatio );
		float                                        getFrequency( void );
		void                              setFrequency( float frequency );
		float                                   getFrequencyRatio( void );
		void                    setFrequencyRatio( float frequencyRatio );
		int                                              getDepth( void );
		void                                        setDepth( int depth );
		bool                                        getInflection( void );
		void                             setInflection( bool inflection );

	private:
		// perlin noise parameters
		float                                                 m_amplitude;
		float                                            m_amplitudeRatio; // will control the influence of the higherfrequencies
		float                                                 m_frequency;
		float                                            m_frequencyRatio;

		int                                                     m_octaves; // number of passes - the more passes, the more higher frequencies
		bool                                                 m_inflection;

		static unsigned char                    g_permutationTable[ 256 ];
		static float                               g_gradientTable[ 768 ];
		static bool                                     g_tablesInitiated;

		static void                         initGradientTable( int seed );

        float                       interpolatedNoise( float s, float t ); // this is the noise interpolation function for 2dimensional perlin noise
        float              interpolatedNoise( float u, float v, float w ); // this is the noise interpolation function for 3dimensional perlin noise
        float interpolatedNoise( float _u, float _v, float _w, float _s ); // this is the noise interpolation function for 4dimensional perlin noise

		float               interpolatedGradientNoise( float u, float v );
		float      interpolatedGradientNoise( float u, float v, float w );
		float               interpolatedGradientNoise( float _u, float _v,
			                                         float _w, float _s );

		float                                              noise( int x ); // this is a noise function which will return a pseudorandom number based on the value x
		float                                       noise( int x, int y );
		float                                noise( int x, int y, int z );
		float                         noise( int x, int y, int z, int w );

		float            gradientNoise(int x, int y, float fx, float fy );
		float                           gradientNoise(int x, int y, int z,
			                                float fx, float fy, float fz);
		float                    gradientNoise(int x, int y, int z, int w,
			                     float fx, float fy, float fz, float fw );

		float            interpolateLinear( float v1, float v2, float t ); // simple linear interpolation between two values dependent on a interpolation value
        float         interpolateEaseCurve( float v1, float v2, float t ); // will interpolate between two values in a more harmonic and smooth manner
	};
}
