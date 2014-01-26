/*************************************************************************************************************************
 *
 * Class Solver
 *
 * Functions for initial condition and velocity 
 *
 *************************************************************************************************************************/

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <ctime>
#include "QuadTree.h"
#include "QuadNode.h"
#include "Point.h"
#include "Solver.h" 
using namespace std;

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Initial condition
// For given x and y (correspond to x_1 and x_2 in [1])
// computes initial condition phi (corresponds to u in [1]),
// its first derivatives and its cross-derivative
void Solver::inco(double x, double y, double& phi, double& phix, double& phiy, double& phixy)
{

   // Periodized Gaussian, eq (3.2) in [1]
   // Initialize to zero
   double phi0 = 0.0;
   double phi0x = 0.0;
   double phi0y = 0.0;
   double phi0xy = 0.0; 
   // Sum up Gaussians
   for (int ix=-30; ix<=30; ix++) {
      for (int iy=-30; iy<=30; iy++) {
         double xi = x - double(ix);
         double yi = y - double(iy);
         phi0 = phi0 + exp( -10.0*((xi-0.5)*(xi-0.5)+(yi-0.75)*(yi-0.75)) );
         phi0x = phi0x - 20.0*(xi-0.5) * exp( -10.0*((xi-0.5)*(xi-0.5)+(yi-0.75)*(yi-0.75)) );
         phi0y = phi0y - 20.0*(yi-0.75) * exp( -10.0*((xi-0.5)*(xi-0.5)+(yi-0.75)*(yi-0.75)) );
         phi0xy = phi0xy + 400.0*(yi-0.75)*(xi-0.5) * exp( -10.0*((xi-0.5)*(xi-0.5)+(yi-0.75)*(yi-0.75)) );
   } }
   // Return values
   phi = phi0;
   phix = phi0x;
   phiy = phi0y;
   phixy = phi0xy;

}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Advection velocity
// For given x, y and time t
// compute the velocity components vx and vy (correspond to a_1 and a_2 in [1])
void Solver::velocity(double x, double y, double t, double& vx, double& vy)
{
   // Swirl test, eq (3.3) in [1]
   double T = 10.0;
   double xi = x;
   double yi = y;
   // TODO: M_PI may need to be redefined for double precision
   vx = cos(M_PI*t/T) * sin(M_PI*xi) * sin(M_PI*xi) * sin(2*M_PI*yi);
   vy = - cos(M_PI*t/T) * sin(2*M_PI*xi) * sin(M_PI*yi) * sin(M_PI*yi);
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

