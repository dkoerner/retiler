retiler
=======

![retiling in action](https://raw.github.com/dkoerner/retiler/master/retiling.jpg)

This is an implementation of Greg Turks paper "Re-Tiling Polygonal Surfaces". The paper and a more detailed explanation of what it is about can be found found on http://www.cc.gatech.edu/~turk/retile/retile.html.

status:
-------
I wrote this code as part of my thesis some time ago and thought it might be useful for others. I havent tested it in its current state and it won't compile straight away. However it has been succesfully used on complex meshes which came out of a fluid sim  [link](https://www.inf.tu-dresden.de/content/institutes/smt/cg/results/majorthesis/dkoerner/dkoerner.en.html).

requirements:
-------------
	CGAL (http://www.cgal.org) - retiler code uses only the basic "2D and 3D Linear Geometry Kernel" package of cgal for computationally stable intersection and orientation tests. The cgal package is licensed under LGPL. 
	
license:
--------
Copyright (c) David Koerner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

