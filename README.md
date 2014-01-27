GALA2D
======

Adaptive gradient-augmented level set method for a two-dimensional advection equation

1) Install

 git clone https://github.com/dkolom/GALA2D.git
 make

The code has only been tested with gcc.

2) Execute

 ./gala2d

All parameters are hardcoded. 
Most of them are in main.cpp.
The initial condition and the velocity field are defined in Fields.cpp.
Parameter eta is set in Solver.cpp (normally you should not change it).

3) Postprocess

Matlab scripts in post_scripts.tar.gz can be used to visualize the solution and the mesh. 

4) Read the theoretical manual

 acroread adaptive.pdf
