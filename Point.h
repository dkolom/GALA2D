/*************************************************************************************************************************
 *
 * Class Point
 *
 * Data stored at a grid point
 *
 *************************************************************************************************************************/

#include <stdlib.h>
const double doubleinit = strtod("NAN",0);
//const double doubleinit = 0.0;

class Point
{
public:
   // Create point with user-set function values and derivatives
   // The point coordinates are defined as integer numbers px and py
   // that are indices of this point on the finest-scale grid (note that the scale cannot be larger than lmax)
   // ph is the finest-grid step size
   Point(unsigned int px, unsigned int py, double ph, double phi, double phix, double phiy, double phixy, Point* next)
   :  m_px(px),
      m_py(py),
      m_phi(phi),
      m_phix(phix),
      m_phiy(phiy),
      m_phixy(phixy),
      m_phi_upd(doubleinit),
      m_phix_upd(doubleinit),
      m_phiy_upd(doubleinit),
      m_phixy_upd(doubleinit),
      m_next(next)
   {
      // Compute physical coordinates using the finest-grid indices 
      m_x = px*ph;
      m_y = py*ph;
   }

   // Create point and set function values and derivatives to default values
   Point(unsigned int px, unsigned int py, double ph, Point* next)
   :  m_px(px),
      m_py(py),
      m_phi(doubleinit),
      m_phix(doubleinit),
      m_phiy(doubleinit),
      m_phixy(doubleinit),
      m_phi_upd(doubleinit),
      m_phix_upd(doubleinit),
      m_phiy_upd(doubleinit),
      m_phixy_upd(doubleinit),
      m_next(next)
   {
      m_x = px*ph;
      m_y = py*ph;
   }

   unsigned int m_px;   // Integer index in direction x (tree search algorithm operates integer values)
   unsigned int m_py;   // Integer index in direction y
   double m_x;          // x coordinate (corresponds to x_1 in [1])
   double m_y;          // y coordinate (corresponds to x_2 in [1])
   double m_phi;        // Function value (corresponds to u in [1])
   double m_phix;       // Derivative in x (corresponds to du/dx_1 in [1])
   double m_phiy;       // Derivative in y (corresponds to du/dx_2 in [1])
   double m_phixy;      // Cross-derivative
   double m_phi_upd;    // Function values 
   double m_phix_upd;   //   and derivatives
   double m_phiy_upd;   //   at the next time step,
   double m_phixy_upd;  //   as the time-stepping scheme requires that two layers in time are stored
   Point* m_next;       // Pointer to the next point in the mesh
};

