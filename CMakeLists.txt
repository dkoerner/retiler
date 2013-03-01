project("retiler")
cmake_minimum_required(VERSION 2.8)

set( APP_NAME retiler )
MESSAGE( STATUS "APP_NAME  :         " ${APP_NAME} )


message( ERROR "need to add cgal in cmake script"


# global options ================================================================



# fill in path variables =========================================================
MESSAGE( STATUS "SOURCE_DIR   :         " ${CMAKE_CURRENT_SOURCE_DIR} )

# include some useful cmake scripts from base =========================================================
include( ${CMAKE_CURRENT_SOURCE_DIR}/sourceGroups.cmake )

# setup different things =========================================================
add_definitions(-DSRC_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# on windows
if( CMAKE_SYSTEM_NAME MATCHES Windows )
	add_definitions(-DNOMINMAX) # windows defines min and max in winbase.h which clashes with stl
	add_definitions(-DPA_USE_DS) # portaudio to use directx sound on windows

	# static runtime (http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F)
	foreach(flag_var
			CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
			CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
	   if(${flag_var} MATCHES "/MD")
		  string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
	   endif(${flag_var} MATCHES "/MD")
	endforeach(flag_var)
endif( CMAKE_SYSTEM_NAME MATCHES Windows )

# on linux
if( CMAKE_SYSTEM_NAME MATCHES Linux )
endif( CMAKE_SYSTEM_NAME MATCHES Linux )



# build app application from sources =========================================================

# specify app sources
file(GLOB_RECURSE app_sources src/*.c src/*.cpp src/*.h src/*.hpp src/*.glsl src/*.ui)
sourceGroups( app_sources ${CMAKE_CURRENT_SOURCE_DIR} )



# Create an executable
add_executable(${APP_NAME} ${app_sources} )

# set include directories =========================================================

#todo: cgal

# set link dependecies =========================================================

#todo: cgal

# on linux
if( CMAKE_SYSTEM_NAME MATCHES Linux )
	# X11
	find_package(X11)

	if(NOT X11_FOUND)
		message(FATAL_ERROR "Failed to find X11 which is required to build under linux")
	endif(NOT X11_FOUND)

	target_link_libraries(${APP_NAME} ${X11_LIBRARIES})

endif( CMAKE_SYSTEM_NAME MATCHES Linux )

