/*************************************************************************************************************************
 *
 *
 *                      Gradient-Augmented Level set Adaptive method for 2D advection equation
 *
 *                      Reference [1]: Dmitry Kolomenskiy, Jean-Christophe Nave and Kai Schneider
 *                      "Adaptive gradient-augmented levet set method with multiresolution error estimation"
 *
 *                      Notations are different in the code and in [1]; see comments
 *
 *
 *************************************************************************************************************************/

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include "QuadNode.h"
#include "QuadTree.h"
#include "Solver.h"
#include "Point.h"
using namespace std;

int main(){

    // This is a sample hi-level call sequence for a typical numerical simulation.
    // In this example, the swirl test described in section 3.2 in [1] is computed.
    
    // Create a new solver
    Solver* asolver = new Solver();

    // Output file for number of nodes vs time 
    ofstream fnumnodes;
    fnumnodes.open ("numnodes.dat", ios::out);

    // Define parameters of the swirl test
    // Algorithm in Section 2.5 in [1], steps 1a-1b
    char lmin = 1;            // m in [1]
    char lmax = 11;           // M in [1]
    double tolerance = 5e-3;  // epsilon in [1]
    double tmax = 10.0;       // T in [1]
    double size = 1.0;        // L in [1]
    // Set parameters
    asolver->setParams(lmin,lmax,size,tmax,tolerance);

    // Initialize tree data structure
    // Step 1c
    // Here we assume l_init = m = 1
    asolver->startupTree(0,0,0);
    asolver->getTree()->createListLeaves();

    // Initial multiresolution analysis that generates the initial mesh
    for (int imra=lmin; imra<lmax; imra++) {
       // List of grid points, step 1d
       asolver->createMesh();
       // Evaluate the initial condition, step 1e
       asolver->initPointData();
       // Remesh, step 1f
       asolver->remesh();
    }

    // Output for figures 3.4(a) and 3.5(a)
    // Save the initial mesh in a file
    asolver->saveMeshCoords();
    // Save the initial condition in a file
    //asolver->saveCart();  // uncomment (may be slow)
    system("mv meshfile.dat meshfile_t0.dat");
    system("mv pointfile.dat pointfile_t0.dat");
    //system("mv cartfile.dat cartfile_t0.dat");  // uncomment (may be slow)

    // Begin time stepping, step 2 in [1], Section 2.5
    cout << "begin iterations" << endl;
    // Initialize iteration counter
    int it = 0;
    // Iterate until final time is reached
    while (asolver->m_time < tmax) {

       // Create list of grid points, 
       // step 2a
       asolver->createMesh();

       // Find minimum, maximum and median levels in the mesh, for diagnostics
       unsigned char lleafmin,lleafmax,lleafmed;
       asolver->getTree()->levelStats(asolver->getTree()->getRoot(),lleafmin,lleafmax);
       lleafmed = (lleafmin+lleafmax)/2;

       // Find the number of the tree structure nodes. It gives the number of grid points shown in figure 3.7(a) in [1]
       unsigned int counter = 0;
       asolver->getTree()->countNodes(asolver->getTree()->getRoot(),counter);

       // Steps 2b-2i of the algorithm
       asolver->timeStep();

       // Remeshing, step 2j
       asolver->remesh();

       // Print diagnostics
       cout << "end it=" << ++it << " time=" << asolver->m_time << " lleafmin=" << int(lleafmin) 
            << " lleafmax=" << int(lleafmax) << " lleafmed=" << int(lleafmed) << endl;

       // Output at t=5, figures 3.4(b) and 3.5(b)
       if ( abs(asolver->m_time-5.0)<1e-12 ) {
          // Save mesh
          asolver->saveMeshCoords();
          // Save solution
          //asolver->saveCart(); // uncomment (may be slow)
          system("mv meshfile.dat meshfile_t5.dat");
          system("mv pointfile.dat pointfile_t5.dat");
          //system("mv cartfile.dat cartfile_t5.dat");  // uncomment (may be slow)
       }

       // Output for figure 3.7
       fnumnodes << it << " " << asolver->m_time << " " << counter << endl;
    }

    // Output at t=10, figures 3.4(c) and 3.5(c)
    // Save mesh
    asolver->saveMeshCoords();
    // Save solution
    //asolver->saveCart();  // uncomment (may be slow)
    system("mv meshfile.dat meshfile_t10.dat");
    system("mv pointfile.dat pointfile_t10.dat");
    //system("mv cartfile.dat cartfile_t10.dat");  // uncomment (may be slow)

    // Compute the error in L^\infty norm
    double linferr;
    linferr = asolver->linfError();
    cout << "Linfty error = " << linferr << endl;

    // Delete the solver
    delete asolver;

    // Close output file
    fnumnodes.close();

    // End computation
    return 0;
}
