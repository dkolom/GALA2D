/*************************************************************************************************************************
 *
 * Class Solver
 *
 * Solves the advection equation on an adaptive grid. 
 *
 *************************************************************************************************************************/

class QuadTree;
class Point;

class Solver
{
private: 

public:
   unsigned char m_lmin;    // Minimum level (m in [1])
   unsigned char m_lmax;    // Maximum level (M in [1])
   unsigned int  m_jmax;    // Maximum number of grid points in one direction
   unsigned int  m_it;      // Iteration number
   double        m_size;    // Domain size (L in [1])
   double        m_tmax;    // Maximum time (T in [1])
   double        m_tol;     // Multiresolution thresholding tolerance (epsilon in [1])
   double        m_time;    // Current time (t in [1])
   double        m_dt;      // Global time step
   double        m_epsgrad; // Step for finite differences (eta in [1])
   QuadTree*     m_tree;    // Pointer to tree

   Solver();
   ~Solver();

   // Set solution parameters
   void setParams(unsigned char lmin, unsigned char lmax, double size, double tmax, double tol);
   // Initial condition for the scalar field evaluated at x,y
   void inco(double x, double y, double& phi, double& phix, double& phiy, double& phixy);
   // Velocity at point x,y and time t
   void velocity(double x, double y, double t, double& vx, double& vy);
   // Generate initial tree
   void startupTree(unsigned int njx, unsigned int njy, unsigned char nl);
   // Create mesh points
   void createMesh();
   // Advance solution in time
   void timeStep(); 
   // Time integration scheme
   void rk4Advect(double x0, double y0, double& xN, double& yN, double& Ex, double& Ey);
   // Remesh
   void remesh();
   // Evaluate initial condition at mesh points
   void initPointData();
   // Save mesh into file
   void saveMeshCoords();
   // Save solution into file
   void saveCart();
   // Returns pointer to tree
   QuadTree* getTree() {return m_tree;}
   // Delete solver
   void cleanup(); 
   // Hermite interpolation of function values, see Appendix in [1]
   void hermite2d(double phi[2][2], double dpdx[2][2], double dpdy[2][2], double dpdxy[2][2], double x0, double y0, double h, double xi, double yi, double& phi_interp);
   // Hermite interpolation of function values and derivatives at mid-points, required for multiresolution. See Appendix in [1]
   void hermite2d_mid(double phi[2][2], double dpdx[2][2], double dpdy[2][2], double dpdxy[2][2], double phi_int[3], double phix_int[3], double phiy_int[3], double phixy_int[3]);
   // Compute initial time step size
   void initTimeStep();
   // Compute L infinity error with respect to the initial condition
   double linfError();

   // Cubic Basis Functions
   inline long double f(long double x) {return (2.0*x*x*x-3.0*x*x+1.0);}
   inline long double g(long double x) {return (x*x*x-2.0*x*x+x);}
   inline long double df(long double x) {return (6.0*x*x-6.0*x);}
   inline long double dg(long double x) {return (3.0*x*x-4.0*x+1.0);}
   inline long double ddf(long double x) {return (12.0*x-6.0);}
   inline long double ddg(long double x) {return (6.0*x-4.0);}

   // Square
   inline long double sqr(long double a) {return a*a;}

};

