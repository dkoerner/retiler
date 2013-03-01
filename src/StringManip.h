// Copyright (c) David Koerner - https://github.com/dkoerner/retiler - see README.md for details
/*---------------------------------------------------------------------

This header offers some very usefull functions which you can use
for manipulating and working with std::strings.

----------------------------------------------------------------------*/
#pragma once
#include <string>
#include <vector>

namespace dk
{
	namespace util
	{
		inline std::string trim_right (const std::string & s, const std::string & t = " \t\r\n");

		inline std::string trim_left (const std::string & s, const std::string & t = " \t\r\n");

		inline std::string trim (const std::string & s, const std::string & t = " \t\r\n");


		/// split a line into the first word, and rest-of-the-line
		std::string getWord (std::string & s, const std::string delim = " ",const bool trim_spaces = true);

		/// split string through some tokens
		void splitString( const std::string s, std::vector<std::string> & v, const std::string delim = " ", const bool trim_spaces = true);


		// ------------------- FileName / Path Manipulation routines ------------------

		std::string setExtension( const std::string &fileName, const std::string &extension ); ///< Returns the given filename with another extension.
		std::string                               getExtension( const std::string &fileName ); ///< Returns the extension (including period) of the given path


		// -------------------- Character set conversions -----------------------------
		std::wstring                                          toWString(const std::string& s); ///< converts the given multibyte character string to unicode
		std::string                                           toString(const std::wstring& s); ///< converts the given widecharacter/unicode character string to multibyte version

	} // namespace util
} // namespace dk










