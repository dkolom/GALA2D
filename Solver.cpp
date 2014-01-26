/*************************************************************************************************************************
 *
 * Class Solver
 *
 * Solves the advection equation on an adaptive grid. 
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

Solver::Solver()
{
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

Solver::~Solver()
{
   cleanup();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Set model parameters
void Solver::setParams(unsigned char lmin, unsigned char lmax, double size, double tmax, double tol)
{
   m_lmin = lmin;
   m_lmax = lmax;
   m_size = size;
   m_tmax = tmax;
   m_tol = tol;
   m_time = 0;
   m_it = 0;
   m_epsgrad = 1e-7; 
   m_jmax = std::pow(2,int(lmax));

   m_tree = new QuadTree(this);
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Generate initial tree
void Solver::startupTree(unsigned int njx, unsigned int njy, unsigned char nl)
{
   if (nl < m_lmin) {
      // Divide and execute MRA for all sub-cells
      startupTree(2*njx,2*njy,nl+1);
      startupTree(2*njx+1,2*njy,nl+1);
      startupTree(2*njx,2*njy+1,nl+1);
      startupTree(2*njx+1,2*njy+1,nl+1);
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void Solver::createMesh()
{
   // Link root point
   m_tree->m_headpoint = m_tree->m_rootpoint;
   m_tree->m_headpoint->m_next = NULL;
   // Link all cell points
   m_tree->linkPoints(m_tree->getRoot());
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void Solver::timeStep()
{
   unsigned int clocx, clocy;
   unsigned int px, py;
   unsigned int dtflag;
   double time_error_norm;
   double x, y;
   double x_w, y_s, x_e, y_n;
   double x_foot, y_foot;
   double x_foot_sw, y_foot_sw, x_foot_se, y_foot_se, x_foot_nw, y_foot_nw, x_foot_ne, y_foot_ne;
   double Ex_foot, Ey_foot;
   double Ex_foot_sw, Ey_foot_sw, Ex_foot_se, Ey_foot_se, Ex_foot_nw, Ey_foot_nw, Ex_foot_ne, Ey_foot_ne;
   double x_cell, y_cell, h;
   double phi[2][2], dpdx[2][2], dpdy[2][2], dpdxy[2][2];
   double phi_sw, phi_se, phi_nw, phi_ne; 
   double dtmax, dtnew;
   QuadNode* cnode = NULL;
   Point* cpoint = NULL;

   // Time step size control parameters
   static double csmax = 2.0;
   static double csmin = 0.5;
   static double cs = 0.75;       
   dtmax = m_size;

   // Initialize time step size if this is zeroth iteration
   if (m_it == 0) initTimeStep();

   // Iterate until time error is smaller than tolerance
   dtflag = 1;
   while (dtflag) {

     // Adjust time step to stop exactly at t=T (see algorithm in Section 2.5 in [1], step 2b)
     if (m_time + m_dt > m_tmax) m_dt = m_tmax - m_time;

     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     // COMMENT THIS LINE IF NO OUTPUT IS REQUIRED AT t=5
     // TODO: Define a parameter that controls snapshot output frequency
     // Adjust time step for output at t=5
     if ((m_time < 0.5*m_tmax) && (m_time + m_dt > 0.5*m_tmax)) m_dt = 0.5*m_tmax - m_time;
     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     // Initialize time error estimate
     time_error_norm = 0.0;

     // Set first point
     cpoint = m_tree->m_headpoint;
     // Loop over all grid points (algorithm in Section 2.5, [1], step 2c)
     while (cpoint != NULL) {
        // Get point coordinates
        x = cpoint->m_x;
        y = cpoint->m_y;
        // Finite difference stencil
        x_w = x - m_epsgrad;
        x_e = x + m_epsgrad;
        y_s = y - m_epsgrad;
        y_n = y + m_epsgrad;
        // Advect backwards in time with given number or local time steps
        rk4Advect(x_w,y_s,x_foot_sw,y_foot_sw,Ex_foot_sw,Ey_foot_sw);
        rk4Advect(x_e,y_s,x_foot_se,y_foot_se,Ex_foot_se,Ey_foot_se); 
        rk4Advect(x_w,y_n,x_foot_nw,y_foot_nw,Ex_foot_nw,Ey_foot_nw);
        rk4Advect(x_e,y_n,x_foot_ne,y_foot_ne,Ex_foot_ne,Ey_foot_ne);
        // Centre point
        x_foot = 0.25 * (x_foot_sw+x_foot_se+x_foot_nw+x_foot_ne);
        y_foot = 0.25 * (y_foot_sw+y_foot_se+y_foot_nw+y_foot_ne);
        // Error estimate for centre point
        Ex_foot = 0.25 * (Ex_foot_sw+Ex_foot_se+Ex_foot_nw+Ex_foot_ne);
        Ey_foot = 0.25 * (Ey_foot_sw+Ey_foot_se+Ey_foot_nw+Ey_foot_ne);      
        // Ensure periodicity
        if (x_foot<0) {
           x_foot += m_size;
           x_foot_sw += m_size;
           x_foot_se += m_size;
           x_foot_nw += m_size;
           x_foot_ne += m_size;
        }
        if (x_foot>m_size) {
           x_foot -= m_size;
           x_foot_sw -= m_size;
           x_foot_se -= m_size;
           x_foot_nw -= m_size;
           x_foot_ne -= m_size;
        }
        if (y_foot<0) {
           y_foot += m_size;
           y_foot_sw += m_size;
           y_foot_se += m_size;
           y_foot_nw += m_size;
           y_foot_ne += m_size;
        }
        if (y_foot>m_size) {
           y_foot -= m_size;
           y_foot_sw -= m_size;
           y_foot_se -= m_size;
           y_foot_nw -= m_size;
           y_foot_ne -= m_size;
        }
        // Find cell that contains the foot point
        px = (unsigned int)(x_foot/m_size*m_jmax);
        py = (unsigned int)(y_foot/m_size*m_jmax);
        cnode = m_tree->findLeafXY_int(m_tree->getRoot(), px, py, 2*m_jmax);
        // Interpolation cell size and corner coordinates
        h = m_size/std::pow(2,int(cnode->m_l));
        x_cell = double(cnode->m_jx)*h;
        y_cell = double(cnode->m_jy)*h;
        // Values at cell corners
        for (clocx=0; clocx<2; clocx++)
        for (clocy=0; clocy<2; clocy++) {
          phi[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phi;
          dpdx[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phix;
          dpdy[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phiy;
          dpdxy[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phixy;
        }
        // Interpolate
        hermite2d(phi,dpdx,dpdy,dpdxy,x_cell,y_cell,h,x_foot_sw,y_foot_sw,phi_sw);
        hermite2d(phi,dpdx,dpdy,dpdxy,x_cell,y_cell,h,x_foot_se,y_foot_se,phi_se);
        hermite2d(phi,dpdx,dpdy,dpdxy,x_cell,y_cell,h,x_foot_nw,y_foot_nw,phi_nw);
        hermite2d(phi,dpdx,dpdy,dpdxy,x_cell,y_cell,h,x_foot_ne,y_foot_ne,phi_ne);
        // Update values
        cpoint->m_phi_upd   = ( phi_sw+phi_se+phi_nw+phi_ne)/(4.0);
        cpoint->m_phix_upd  = (-phi_sw+phi_se-phi_nw+phi_ne)/(4.0*m_epsgrad);
        cpoint->m_phiy_upd  = (-phi_sw-phi_se+phi_nw+phi_ne)/(4.0*m_epsgrad);
        cpoint->m_phixy_upd = ( phi_sw-phi_se-phi_nw+phi_ne)/(4.0*m_epsgrad*m_epsgrad);
        // Estimate time error
        time_error_norm = max( time_error_norm, abs( cpoint->m_phix_upd*Ex_foot + cpoint->m_phiy_upd*Ey_foot + cpoint->m_phixy_upd*Ex_foot*Ey_foot ) );
        // Next point
        cpoint = cpoint->m_next;
     }
     // New time step size
     dtnew = min(dtmax, m_dt * max(csmin,min(csmax,cs*std::pow(m_tol/time_error_norm,0.25))));
     // If error is large, adjust the time step size and iterate again
     if (time_error_norm <= m_tol) dtflag = 0; else m_dt = dtnew;
   }

   // Update time variable, next time step size and iteration number
   m_time += m_dt;
   m_dt = dtnew;
   m_it += 1;
   // Update points (Section 2.5 in [1], step 2i)
   // Set first point
   cpoint = m_tree->m_headpoint;
   // Loop over all grid points
   while (cpoint != NULL) {
      // Copy
      cpoint->m_phi = cpoint->m_phi_upd;
      cpoint->m_phix = cpoint->m_phix_upd;
      cpoint->m_phiy = cpoint->m_phiy_upd;
      cpoint->m_phixy = cpoint->m_phixy_upd;
      // Next point
      cpoint = cpoint->m_next;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void Solver::remesh()
{
   QuadNode* cnode = NULL;
   Point* cpoint = NULL;
   unsigned int j;
   unsigned char hold, ilev;
   int jj; 
   unsigned int clocx, clocy;
   double phi[2][2], phix[2][2], phiy[2][2], phixy[2][2];
   double phiM[3], phixM[3], phiyM[3], phixyM[3];
   double phiM_int[3], phixM_int[3], phiyM_int[3], phixyM_int[3];
   double r_phiM[3], r_phixM[3], r_phiyM[3], r_phixyM[3];
   double nloc, h, h2, dnorm, threshold;

   // The remeshing procedure is also explained in Section 2.5 in [1].

   // Analyze two levels to allow coarsening
   for (ilev=0; ilev<2; ilev++) {

      // Create new list of leaves
      m_tree->createListLeaves();

      // For all leaves, compute residuals and do thresholding
      for (j=0; j<m_tree->m_nleaf; j++) {
         // Current node
         cnode = m_tree->m_leaf[j];
         // Cell size
         nloc = exp(double(cnode->m_l)*log(2.0));
         h = m_size / nloc;
         // Cell size one level above (used for scaling)
         h2 = 0.5 * h;
         // Evaluate at finer scale
         // Four corner cell points
         // Scaled derivatives are calculated, as required for the stability of MRA
         for (clocx=0; clocx<2; clocx++)
         for (clocy=0; clocy<2; clocy++) {
            cpoint = cnode->getPoint(clocx,clocy);
            phi[clocx][clocy] = cpoint->m_phi;
            phix[clocx][clocy] = h2 * cpoint->m_phix;
            phiy[clocx][clocy] = h2 * cpoint->m_phiy;
            phixy[clocx][clocy] = h2*h2 * cpoint->m_phixy;
         }
         // Three points for horizontal, vertical and diagonal details, resp.
         // 0:M0=SE=1 1:0M=NW=2 2:MM=NE=3
         for (jj=0; jj<3; jj++) {
            cpoint = cnode->m_cpoint[jj];
            phiM[jj] = cpoint->m_phi;
            phixM[jj] = h2 * cpoint->m_phix;
            phiyM[jj] = h2 * cpoint->m_phiy;
            phixyM[jj] = h2*h2 * cpoint->m_phixy;
         }
         // Interpolate
         hermite2d_mid(phi,phix,phiy,phixy,phiM_int,phixM_int,phiyM_int,phixyM_int);
         // Compute scaled residuals
         for (jj=0; jj<3; jj++) {
            r_phiM[jj] = phiM_int[jj] - phiM[jj];
            r_phixM[jj] = phixM_int[jj] - phixM[jj];
            r_phiyM[jj] = phiyM_int[jj] - phiyM[jj];
            r_phixyM[jj] = phixyM_int[jj] - phixyM[jj];
         }
         // Compute max
         dnorm = 0.0;
         for (jj=0; jj<3; jj++) {
            if (abs(r_phiM[jj]) > dnorm) dnorm = abs(r_phiM[jj]);
            if (abs(r_phixM[jj]) > dnorm) dnorm = abs(r_phixM[jj]);
            if (abs(r_phiyM[jj]) > dnorm) dnorm = abs(r_phiyM[jj]);
            if (abs(r_phixyM[jj]) > dnorm) dnorm = abs(r_phixyM[jj]);
         }
         // Threshold
         threshold = m_tol;  // scale-independent thresholding
         // Thresholding: mark nodes for removal
         if ( ( (cnode->m_l > m_lmin) && (dnorm < threshold) ) || (cnode->m_l == m_lmax-1) ) {
           cnode->m_hold = 0;
         }
      }

      // Can only delete all 4 child nodes. Check if they are all marked. If not, unmark.
      for (j=0; j<m_tree->m_nleaf; j++) {
         // Current node
         cnode = m_tree->m_leaf[j];
         // Never mark the root node
         if (cnode->m_l == 0) {
            cnode->m_hold = 1;
         // Otherwise this node may be marked or unmarked
         } else {
            // 'hold' must be equal to zero to allow coarsening
            hold = 0;
            for (jj=0; jj<4; jj++) {
               hold += cnode->m_parent->m_child[jj]->m_hold;
            }
            // Unmark this leaf if cannot coarsen it
            if (hold > 0) cnode->m_hold = 1;
         }
      }

      // Make the tree graded and update it
      m_tree->balance();
      m_tree->updateTree(m_tree->getRoot());
      m_tree->holdTree(m_tree->getRoot());  // TODO: May be unnecessary !!!
   }

   // List of leaves
   m_tree->createListLeaves();

   // Add one layer
   m_tree->addLeaves();

   // List of leaves
   m_tree->createListLeaves();

   // Assign values to new leaves using interpolation
   // Go from root to top, level by level
   // Here, ilev is a level in the tree
   for (ilev=m_lmin; ilev<=m_lmax; ilev++) {
      for (j=0; j<m_tree->m_nleaf; j++) {
         // Current node
         cnode = m_tree->m_leaf[j];
         if (cnode->m_l == ilev) {
            // Cell size
            nloc = exp(double(cnode->m_l)*log(2.0));
            h = m_size / nloc;
            // Cell size one level above (used for scaling)
            h2 = 0.5 * h;
            // Evaluate at finer scale
            // Four corner cell points
            for (clocx=0; clocx<2; clocx++)
            for (clocy=0; clocy<2; clocy++) {
               cpoint = cnode->getPoint(clocx,clocy);
               phi[clocx][clocy] = cpoint->m_phi;
               phix[clocx][clocy] = h2 * cpoint->m_phix;
               phiy[clocx][clocy] = h2 * cpoint->m_phiy;
               phixy[clocx][clocy] = h2*h2 * cpoint->m_phixy;
            }
            // Interpolate
            hermite2d_mid(phi,phix,phiy,phixy,phiM_int,phixM_int,phiyM_int,phixyM_int);
            // Three points for horizontal, vertical and diagonal details, resp.
            // 0:M0=SE=1 1:0M=NW=2 2:MM=NE=3
            for (jj=0; jj<3; jj++) {
               cpoint = cnode->m_cpoint[jj];
               cpoint->m_phi = phiM_int[jj];
               cpoint->m_phix = phixM_int[jj] /h2;
               cpoint->m_phiy = phiyM_int[jj] /h2;
               cpoint->m_phixy = phixyM_int[jj] /h2/h2;
            }
         }
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Evaluate initial condition at mesh points
void Solver::initPointData()
{
   Point* cpoint = m_tree->m_headpoint;
   while (cpoint!=NULL) {
      inco(cpoint->m_x,cpoint->m_y,cpoint->m_phi,cpoint->m_phix,cpoint->m_phiy,cpoint->m_phixy);
      cpoint = cpoint->m_next;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Interplate solution on a uniform grid and save to file
void Solver::saveCart()
{
   // Grid size
   int nx = m_jmax;
   int ny = m_jmax;

   unsigned int clocx, clocy;
   double x, y;
   double x_cell, y_cell;
   double h;
   double phi[2][2], dpdx[2][2], dpdy[2][2], dpdxy[2][2];

   double phi_interp;
   QuadNode* cnode = NULL;

   // Open file
   ofstream cartfile;
   cartfile.open ("cartfile.dat", ios::out | ios::trunc);
   cartfile.precision(18);

   // Probe at Cartesian grid nodes
   for ( int iy=0; iy<ny; iy++ ) {
      y = m_size * double(iy) / double(ny);
      for ( int ix=0; ix<nx; ix++ ) {
         x = m_size * double(ix) / double(nx);
         // Find cell
         unsigned int px = (unsigned int)(x/m_size*m_jmax);
         unsigned int py = (unsigned int)(y/m_size*m_jmax);
         cnode = m_tree->findLeafXY_int(m_tree->getRoot(), px, py, 2*m_jmax);
         // Interpolate - get updated phi and derivatives
         x_cell = cnode->getPoint(0,0)->m_x;
         y_cell = cnode->getPoint(0,0)->m_y;
         h = cnode->getPoint(1,1)->m_x - x_cell;
         // Values at cell corners
         for (clocx=0; clocx<2; clocx++)
         for (clocy=0; clocy<2; clocy++) {
            phi[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phi;
            dpdx[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phix;
            dpdy[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phiy;
            dpdxy[clocx][clocy] = cnode->getPoint(clocx,clocy)->m_phixy;
         }
         // Interpolate
         hermite2d(phi,dpdx,dpdy,dpdxy,x_cell,y_cell,h,x,y,phi_interp);
         // Save to file
         cartfile << phi_interp << " ";
      }
      cartfile << "\n";
   }
   
   // Close file
   cartfile.close();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Save mesh to file
void Solver::saveMeshCoords()
{
   // Save cells
   m_tree->createListLeaves();
   ofstream meshfile;
   meshfile.open ("meshfile.dat", ios::out | ios::trunc);
   meshfile.precision(18);

   for (unsigned int jj=0; jj<m_tree->m_nleaf; jj++) {
      double h = m_size / exp(double(m_tree->m_leaf[jj]->m_l)*log(2.0));
      double x0 = double(m_tree->m_leaf[jj]->m_jx) * h;
      double x1 = x0 + h;
      double y0 = double(m_tree->m_leaf[jj]->m_jy) * h;
      double y1 = y0 + h;
      meshfile << x0 << " " << y0 << " " << x1 << " " << y1 << "\n";
   }

   meshfile.close();

   // Save points
   createMesh();
   ofstream pointfile;
   pointfile.open ("pointfile.dat", ios::out | ios::trunc);
   pointfile.precision(18);

   Point* cpoint = m_tree->m_headpoint;
   while (cpoint!=NULL) {
      pointfile << cpoint->m_x << " " << cpoint->m_y << "\n";
      cpoint = cpoint->m_next;
   }

   pointfile.close();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Clean exit
void Solver::cleanup()
{
   delete m_tree;
   m_tree = NULL;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Hermite interpolation at any point in a cell
// Only function value is interpolated
// See Appendix in [1]
void Solver::hermite2d(double phi[2][2], double dpdx[2][2], double dpdy[2][2], double dpdxy[2][2], double x0, double y0, double h, double xi, double yi, double& phi_interp)
{
    int i00,j00,i10,j10,i01,j01,i11,j11;
    double x,y,ap00,ap10,ap01,ap11,ax00,ax10,ax01,ax11,ay00,ay10,ay01,ay11,axy00,axy10,axy01,axy11;
    double dpdxy00,dpdxy10,dpdxy01,dpdxy11;

    // Input arrays are always 2x2
    i00=0; i10=1;i01=0;i11=1;
    j00=0; j10=0;j01=1;j11=1;

    dpdxy11= dpdxy[i11][j11];
    dpdxy00= dpdxy[i00][j00];
    dpdxy01= dpdxy[i01][j01];
    dpdxy10= dpdxy[i10][j10];

    x=(xi-x0)/h;
    y=(yi-y0)/h;

    // Weights
    ap00 =f(x)*f(y);ap10 = f(1.0-x)*f(y);ap01 = f(x)*f(1.0-y);ap11 = f(1.0-x)*f(1.0-y);
    ax00 =g(x)*f(y);ax10 =-g(1.0-x)*f(y);ax01 = g(x)*f(1.0-y);ax11 =-g(1.0-x)*f(1.0-y);
    ay00 =f(x)*g(y);ay10 = f(1.0-x)*g(y);ay01 =-f(x)*g(1.0-y);ay11 =-f(1.0-x)*g(1.0-y);
    axy00=g(x)*g(y);axy10=-g(1.0-x)*g(y);axy01=-g(x)*g(1.0-y);axy11= g(1.0-x)*g(1.0-y);

    // Function value
    phi_interp = (phi[i00][j00]*ap00+phi[i10][j10]*ap10+phi[i01][j01]*ap01+phi[i11][j11]*ap11)
      +h *(dpdx[i00][j00]*ax00+dpdx[i10][j10]*ax10+dpdx[i01][j01]*ax01+dpdx[i11][j11]*ax11)
      +h *(dpdy[i00][j00]*ay00+dpdy[i10][j10]*ay10+dpdy[i01][j01]*ay01+dpdy[i11][j11]*ay11)
      +h*h*(dpdxy00*axy00+dpdxy10*axy10+dpdxy01*axy01+dpdxy11*axy11);
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Hermite interpolation at midpoints using SCALED derivatives 
// If h is cell size, then  dpdx, dpdy, phix_int, phiy_int are the derivatives times h/2,
// dpdxy and phixy_int are the derivatives times (h/2)^2
// See Appendix in [1]
void Solver::hermite2d_mid(double phi[2][2], double dpdx[2][2], double dpdy[2][2], double dpdxy[2][2], double phi_int[3], double phix_int[3], double phiy_int[3], double phixy_int[3])
{
   // Phi
   phi_int[0] = 0.5*(phi[0][0]+phi[1][0]) + 0.25*(dpdx[0][0]-dpdx[1][0]);
   phi_int[1] = 0.5*(phi[0][0]+phi[0][1]) + 0.25*(dpdy[0][0]-dpdy[0][1]);
   phi_int[2] = 0.25*(phi[0][0]+phi[1][0]+phi[0][1]+phi[1][1]) + 0.125*((dpdx[0][0]-dpdx[1][0]+dpdx[0][1]-dpdx[1][1]) + (dpdy[0][0]+dpdy[1][0]-dpdy[0][1]-dpdy[1][1])) + 0.0625*(dpdxy[0][0]-dpdxy[1][0]-dpdxy[0][1]+dpdxy[1][1]);

   // Phix
   phix_int[0] = 0.75*(phi[1][0]-phi[0][0]) - 0.25*(dpdx[1][0]+dpdx[0][0]);
   phix_int[1] = 0.5*(dpdx[0][0]+dpdx[0][1]) + 0.25*(dpdxy[0][0]-dpdxy[0][1]);
   phix_int[2] = 0.375*(-phi[0][0]+phi[1][0]-phi[0][1]+phi[1][1]) - 0.125*(dpdx[0][0]+dpdx[1][0]+dpdx[0][1]+dpdx[1][1]) - 0.1875*(dpdy[0][0]-dpdy[1][0]-dpdy[0][1]+dpdy[1][1]) - 0.0625*(dpdxy[0][0]+dpdxy[1][0]-dpdxy[0][1]-dpdxy[1][1]);

   // Phiy
   phiy_int[0] = 0.5*(dpdy[0][0]+dpdy[1][0]) + 0.25*(dpdxy[0][0]-dpdxy[1][0]);
   phiy_int[1] = 0.75*(phi[0][1]-phi[0][0]) - 0.25*(dpdy[0][0]+dpdy[0][1]);
   phiy_int[2] = 0.375*(-phi[0][0]+phi[0][1]-phi[1][0]+phi[1][1]) - 0.125*(dpdy[0][0]+dpdy[0][1]+dpdy[1][0]+dpdy[1][1]) - 0.1875*(dpdx[0][0]-dpdx[0][1]-dpdx[1][0]+dpdx[1][1]) - 0.0625*(dpdxy[0][0]+dpdxy[0][1]-dpdxy[1][0]-dpdxy[1][1]);

   // Phixx and Phiyy are not required for error estimate

   //Phixy
   phixy_int[0] = 0.75*(dpdy[1][0]-dpdy[0][0]) - 0.25*(dpdxy[0][0]+dpdxy[1][0]);
   phixy_int[1] = 0.75*(dpdx[0][1]-dpdx[0][0]) - 0.25*(dpdxy[0][0]+dpdxy[0][1]);
   phixy_int[2] = 0.5625*(phi[0][0]-phi[1][0]-phi[0][1]+phi[1][1]) + 0.1875*( (dpdx[0][0]+dpdx[1][0]-dpdx[0][1]-dpdx[1][1]) + (dpdy[0][0]-dpdy[1][0]+dpdy[0][1]-dpdy[1][1]) ) + 0.0625*(dpdxy[0][0]+dpdxy[1][0]+dpdxy[0][1]+dpdxy[1][1]);

}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Embedded RK4(3) for advection
void Solver::rk4Advect(double x0, double y0, double& xN, double& yN, double& Ex, double& Ey)
{
   double vx,vy;
   double x1,y1,x2,y2,x3,y3;
   double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y;
   double endtime;
   double lambda;

   // Dormand-Prince parameter
   lambda = 0.1;

   // Final time
   endtime = m_time + m_dt;

   // RK4 scheme coefficients
   velocity(x0,y0,endtime,vx,vy);
   k1x = - m_dt * vx;
   k1y = - m_dt * vy;
   x1 = x0 + 0.5*k1x;
   y1 = y0 + 0.5*k1y;      

   velocity(x1,y1,endtime-0.5*m_dt,vx,vy);
   k2x = - m_dt * vx;
   k2y = - m_dt * vy;        
   x2 = x0 + 0.5*k2x;
   y2 = y0 + 0.5*k2y;      

   velocity(x2,y2,endtime-0.5*m_dt,vx,vy);
   k3x = - m_dt * vx;
   k3y = - m_dt * vy;
   x3 = x0 + k3x;
   y3 = y0 + k3y;       

   velocity(x3,y3,endtime-m_dt,vx,vy);
   k4x = - m_dt * vx;
   k4y = - m_dt * vy;

   // New point coordinates
   xN = x0 + (1.0/3.0)*( 0.5*k1x + k2x + k3x + 0.5*k4x );
   yN = y0 + (1.0/3.0)*( 0.5*k1y + k2y + k3y + 0.5*k4y );

   velocity(xN,yN,endtime-m_dt,vx,vy);
   k5x = - m_dt * vx;
   k5y = - m_dt * vy;
    
   // Error estimate
   Ex = lambda * (k4x - k5x);
   Ey = lambda * (k4y - k5y);
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Compute L infinity error with respect to the initial condition
// This is useful if the exact solution is periodic in time
double Solver::linfError()
{
   Point* cpoint;
   double errnorm, locerr;
   double inco_phi, inco_phix, inco_phiy, inco_phixy;

   // Make list of grid points
   createMesh();

   // Compute L infinity error
   errnorm = 0.0;
   cpoint = m_tree->m_headpoint;
   while (cpoint!=NULL) {
      // Evaluate initial condition at this point
      inco(cpoint->m_x,cpoint->m_y,inco_phi,inco_phix,inco_phiy,inco_phixy);
      // Error at this point
      locerr = abs(cpoint->m_phi - inco_phi);
      // Maximum error
      if (locerr > errnorm) errnorm = locerr;
      // Next point
      cpoint = cpoint->m_next;
   }

   return errnorm;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Compute initial time step size. 
void Solver::initTimeStep()
{
   unsigned int px,py;
   double x,y,vx,vy,h,dtopt,dtmin;
   QuadNode* cnode;
   Point* cpoint;

   // Initialize dt min
   dtmin = 1.0e6;
   // Set first point
   cpoint = m_tree->m_headpoint;
   // Loop over all grid points
   while (cpoint != NULL) {
      // Get point coordinates
      x = cpoint->m_x;
      y = cpoint->m_y;
      // Velocity
      velocity(x,y,m_time,vx,vy);
      // Find one of the nearest leaves to estimate cell size in the heighbourhood of this point
      px = (unsigned int)(x/m_size*m_jmax);
      py = (unsigned int)(y/m_size*m_jmax);
      cnode = m_tree->findLeafXY_int(m_tree->getRoot(), px, py, 2*m_jmax);
      // Nearest cell size 
      h = m_size/std::pow(2,int(cnode->m_l));
      // Locally optimal time step (cfl=1)
      dtopt = 2.0*h/sqrt((vx*vx+vy*vy)+1);
      // Minimum and maximum
      if (dtopt<dtmin) dtmin = dtopt;
      // Next point
      cpoint = cpoint->m_next;
   }
   // Set imitial time step size to minimum
   m_dt = dtmin;
}


