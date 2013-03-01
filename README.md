retiler
=======

![retiling in action](https://raw.github.com/dkoerner/retiler/master/retiling.jpg)

This is an implementation of Greg Turks paper "Re-Tiling Polygonal Surfaces". The paper and a more detailed explanation of what it is about can be found found on http://www.cc.gatech.edu/~turk/retile/retile.html.

###status
I wrote this code as part of my thesis some time ago and thought it might be useful for others. I havent tested it in its current state and it won't compile straight away. However it has been succesfully used on complex meshes which came out of a fluid sim - [see here](https://www.inf.tu-dresden.de/content/institutes/smt/cg/results/majorthesis/dkoerner/dkoerner.en.html).

###requirements
Retiler uses [CGAL](http://www.cgal.org) for computationally stable intersection and orientation tests. Only the basic "2D and 3D Linear Geometry Kernel" package of cgal is used. This particular package is licensed under LGPL. 

There is a cmake script which currently doesnt add the cgal dependency. Once this works CMake can be used to build the project.

###license
Copyright (c) David Koerner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

###notes

Directory Structure:

    src/        - source code
      main.cpp  - Shows a simple example for how to use the different classes to retile a mesh given as obj file on disk.
      Mesh.h    - very simple Mesh container for vertices and triangles
      ObjIO.h   - allows import and export of Mesh instances to disk. The obj file format is used.
      MeshEx.h  - takes a simple mesh and turns it into more complex winged edge datastructure. This is used by retiler.
      Retiler.h - Here the algorithm is actually implemented. See main.cpp how it is used. 
      math/     - simple math framework for vectors etc.
    
