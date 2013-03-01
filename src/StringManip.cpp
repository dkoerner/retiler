// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This header offers some very usefull functions which you can use
for manipulating and working with std::strings.

----------------------------------------------------------------------*/
#include "StringManip.h"

namespace dk
{
	namespace util
	{
		inline std::string trim_right (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			std::string::size_type i (d.find_last_not_of (t));
			if (i == std::string::npos)
				return "";
			else
				return d.erase (d.find_last_not_of (t) + 1) ; 
		}  // end of trim_right

		inline std::string trim_left (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			return d.erase (0, s.find_first_not_of (t)) ; 
		}  // end of trim_left

		inline std::string trim (const std::string & s, const std::string & t)
		{ 
			std::string d (s); 
			return trim_left (trim_right (d, t), t) ; 
		}  // end of trim


		// split a line into the first word, and rest-of-the-line
		std::string getWord (std::string & s, const std::string delim,const bool trim_spaces)
		{
			// find delimiter  
			std::string::size_type i (s.find (delim));

			// split into before and after delimiter
			std::string w (s.substr (0, i));

			// if no delimiter, remainder is empty
			if (i == std::string::npos)
				s.erase ();
			else
				// erase up to the delimiter
				s.erase (0, i + delim.size ());

			// trim spaces if required
			if (trim_spaces)
			{
				w = trim (w);
				s = trim (s);
			}

			// return first word in line
			return w;
		} // end of getWord	

		// To be symmetric, we assume an empty string (after trimming spaces)
		// will give an empty vector.
		// However, a non-empty string (with no delimiter) will give one item
		// After that, you get an item per delimiter, plus 1.
		// eg.  ""      => empty
		//      "a"     => 1 item
		//      "a,b"   => 2 items
		//      "a,b,"  => 3 items (last one empty)
		void splitString( const std::string s, std::vector<std::string> & v, const std::string delim, const bool trim_spaces)
		{
			// start with initial string, trimmed of leading/trailing spaces if required
			std::string s1 (trim_spaces ? trim (s) : s);

			v.clear (); // ensure vector empty
			
			// no string? no elements
			if (s1.empty ())
				return;

			// add to vector while we have a delimiter
			while (!s1.empty () && s1.find (delim) != std::string::npos)
				v.push_back (getWord (s1, delim, trim_spaces));

			// add final element
			v.push_back (s1);
		} // end of splitString 







		// ------------------- FileName / Path Manipulation routines ------------------



		//
		// Returns the given filename with another extension.
		// if the extension parameter is an empty string, the extension is removed
		//
		std::string setExtension( const std::string &fileName, const std::string &extension )
		{
			std::string ext(extension);

			// is the extension not empty?
			if( extension.size() > 0 )
				// does the extension not has a leading period?
				if( extension[0] != '.' )
					// add a period
					ext = '.' + extension;

			//search rightmost period in file name
			std::string::size_type idx = fileName.rfind( '.' );
			if( idx == std::string::npos )
			{
				// filename does not contain any period
				return fileName + ext;
			}else
			{
				// is the rightmost period the first character within the filename?
				if( idx == 0 )
					// then the filename only consists of an extension
					return ext;

				return fileName.substr( 0, idx ) + ext;
			}
		}

		//
		// Returns the extension (including period) of the given path
		//
		std::string getExtension( const std::string &fileName )
		{
			//search rightmost period in file name
			std::string::size_type idx = fileName.rfind( '.' );
			if( idx == std::string::npos )
			{
				// filename does not contain any period -> so no extension
				return "";
			}else
			{
				return fileName.substr( idx, std::string::npos );
			}
		}


		// --------------------------- CHARACTER SET CONVERSIONS ----------------------------
		std::wstring toWString(const std::string& s)
		{
			std::wstring temp(s.length(),L' ');
			std::copy(s.begin(), s.end(), temp.begin());
			return temp;
		}

		std::string toString(const std::wstring& s)
		{
			std::string temp(s.length(), ' ');
			std::copy(s.begin(), s.end(), temp.begin());
			return temp;
		}
	} // namespace util
} // namespace dk