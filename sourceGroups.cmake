
#
# assigns sourcegroup from path
# the rootdir parameter specifies the part of the path which has to be striped of
# this function accepts an optional argument which allows to preappend some sourceGroup path
#
#
MACRO (sourceGroups sources rootdir)

	FOREACH (src ${${sources}})

		#delete source directory
		string(REGEX REPLACE ${rootdir} "" last_dir ${src})
		#delete last slash and filename
		string(REGEX REPLACE "[\\\\/][^\\\\/]*$" "" last_dir ${last_dir})
		#delete first slash
		string(REGEX REPLACE "^[\\\\/]" "" last_dir ${last_dir})
		#replace forward with backslash
		string(REGEX REPLACE "/" "\\\\" last_dir ${last_dir})
		#preappend optional argument - will be empty by default
		set( last_dir "${ARGN}\\${last_dir}" )
		# assign sourcegroup
		source_group(${last_dir} FILES ${src})
	ENDFOREACH (src)
ENDMACRO (sourceGroups)