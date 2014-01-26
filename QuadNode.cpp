/*************************************************************************************************************************
 *
 * Class QuadNode
 *
 * Node of a tree
 *
 * Adapted from QuadtreeSimulator v1.21 by Frank Dockhorn
 *
 *************************************************************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>
#include <ctime>
#include "QuadNode.h" 
#include "QuadTree.h"
#include "Solver.h"
#include "Point.h"
using namespace std;

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode::QuadNode(QuadTree* ntree, QuadNode* nparent, unsigned int njx, unsigned int njy, unsigned char nl)
{
   m_tree = ntree;
   m_parent = nparent;
   m_child = NULL; 
   m_jx = njx;
   m_jy = njy;
   m_l = nl;

   // Initialize grid points
   assert(m_tree->m_solver->m_lmax-nl-1 >= 0);
   // Finest-level grid step size
   double ph = m_tree->m_solver->m_size / double(m_tree->m_solver->m_jmax);
   // Integer coordinates of the three points that correspond to the given cell (njx,njy,nl)
   unsigned int powl = std::pow(2,int(m_tree->m_solver->m_lmax-nl-1));
   unsigned int pxM0 = powl * (2*njx + 1);
   unsigned int pyM0 = powl * (2*njy);
   unsigned int px0M = powl * (2*njx);
   unsigned int py0M = powl * (2*njy + 1);
   unsigned int pxMM = powl * (2*njx + 1);
   unsigned int pyMM = powl * (2*njy + 1);
   // Create points
   m_cpoint[0] = new Point(pxM0,pyM0,ph,NULL);
   m_cpoint[1] = new Point(px0M,py0M,ph,NULL);
   m_cpoint[2] = new Point(pxMM,pyMM,ph,NULL);
   // Do not remove this point during update
   m_hold = 1;

//   printNode();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void QuadNode::allocateChildren()
{
   m_child = new QuadNode*[4];
   for (int jj=0; jj<4; jj++) {
      m_child[jj] = NULL;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void QuadNode::printNode()
{
   // Prind the position and the level
   cout << m_jx << " " << m_jy << " " << int(m_l) << endl;
}

/*-------------------------------------------------------------------------*\ 
\*-------------------------------------------------------------------------*/

// Destructor recursively removes all branches
QuadNode::~QuadNode()
{
   int i;

   // Recursively remove children
   for (i = 0; i < 4; i++) {
      if (m_child != NULL) {
         delete m_child[i];
         m_child[i] = NULL;
      } 
   }
   delete[] m_child;
   m_child = NULL;

   // Delete cell points
   for (i = 0; i < 3; i++) {
      delete m_cpoint[i];
      m_cpoint[i] = NULL;
   }

//   this = NULL;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadNode::findNode(unsigned int njx, unsigned int njy, unsigned char nl)
{
   // Initial guess
   QuadNode* fnode = this;

   // Check if the current node is sought
   if ( (m_jx!=njx) || (m_jy!=njy) || (m_l!=nl) ) 
   {
      // Search the branch
      location cloc;
      unsigned int cjx;
      unsigned int cjy;
      unsigned char cl;

      // Search begins from 'this' node 
      unsigned int njmax = (unsigned int)(floor(pow(double(2),double(nl-m_l))+0.5));
      for (cl = m_l+1; cl<=nl; cl++)
      {
         njmax /= 2;
         cjx = njx / njmax;
         cjy = njy / njmax;

         // Check jx and jy bounds
         unsigned int cjmax = (unsigned int)(floor(pow(double(2),double(cl))+0.5));
         assert ( (cjx<cjmax) && (cjy<cjmax) ); 

         // Find child position
         if (cjx % 2 == 0) {
            if (cjy % 2 == 0) {
               cloc = QuadNode::SW;
            } else {
               cloc = QuadNode::NW;
            }
         } else {
            if (cjy % 2 == 0) {
               cloc = QuadNode::SE;
            } else {
               cloc = QuadNode::NE;
            }
         }

         // Go to one step up the tree, if possible. 
         // Otherwise return NULL
         if (fnode->m_child != NULL) {
            fnode = fnode->m_child[cloc];
         } else {
            return NULL;
         }

         // Stop if cannot find
         assert(fnode != NULL);
      }
   }

   // Node found
   return fnode;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Returns a pointer to the point
// that is situated in a corner of this node,
// and its location is given by clocx, clocy
// clocx=0, clocy=0 is lower-left corner
// clocx=0, clocy=1 is upper-left corner etc.
Point* QuadNode::getPoint(unsigned int clocx, unsigned int clocy)
{
   unsigned int px, py;
   unsigned int powl = std::pow(2,int(m_tree->m_solver->m_lmax-m_l));
   // Point coordinates (integer number relative to the finest-level grid)
   px = ((m_jx+clocx) * powl) % m_tree->m_solver->m_jmax;
   py = ((m_jy+clocy) * powl) % m_tree->m_solver->m_jmax;
   // Find point if it's not the root point
//   if ( (px > 0) || (py > 0) ) return m_tree->getRoot()->findPoint(px,py,2*m_tree->m_solver->m_jmax);
   if ( (px > 0) || (py > 0) ) return m_tree->getRoot()->findPoint(px,py,m_tree->m_solver->m_jmax);
   // If px == 0 and py == 0, return root point
   return m_tree->m_rootpoint;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadNode::northNeighbor()
{
   QuadNode* cnode;
   // Check the parent cell's children
   if (m_parent == NULL) return this; // 'this' if periodic, 'NULL' if not periodic
   if (this == m_parent->m_child[QuadNode::SW]) return m_parent->m_child[QuadNode::NW];
   if (this == m_parent->m_child[QuadNode::SE]) return m_parent->m_child[QuadNode::NE];

   // Recursive search
   cnode = m_parent->northNeighbor();
   // Note that cnode cannot be NULL
   if (cnode->isLeaf()) {
      return cnode;
   } else {
      if (this == m_parent->m_child[QuadNode::NW]) {
         return cnode->m_child[QuadNode::SW];
      } else {
         return cnode->m_child[QuadNode::SE];
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadNode::southNeighbor()
{
   QuadNode* cnode;
   // Check the parent cell's children
   if (m_parent == NULL) return this;
   if (this == m_parent->m_child[QuadNode::NW]) return m_parent->m_child[QuadNode::SW];
   if (this == m_parent->m_child[QuadNode::NE]) return m_parent->m_child[QuadNode::SE];

   // Recursive search
   cnode = m_parent->southNeighbor();
   if (cnode->isLeaf()) {
      return cnode;
   } else {
      if (this == m_parent->m_child[QuadNode::SW]) {
         return cnode->m_child[QuadNode::NW];
      } else {
         return cnode->m_child[QuadNode::NE];
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadNode::westNeighbor()
{
   QuadNode* cnode;
   // Check the parent cell's children
   if (m_parent == NULL) return this;
   if (this == m_parent->m_child[QuadNode::SE]) return m_parent->m_child[QuadNode::SW];
   if (this == m_parent->m_child[QuadNode::NE]) return m_parent->m_child[QuadNode::NW];

   // Recursive search
   cnode = m_parent->westNeighbor();
   if (cnode->isLeaf()) {
      return cnode;
   } else {
      if (this == m_parent->m_child[QuadNode::SW]) {
         return cnode->m_child[QuadNode::SE];
      } else {
         return cnode->m_child[QuadNode::NE];
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadNode::eastNeighbor()
{
   QuadNode* cnode;
   // Check the parent cell's children
   if (m_parent == NULL) return this;
   if (this == m_parent->m_child[QuadNode::SW]) return m_parent->m_child[QuadNode::SE];
   if (this == m_parent->m_child[QuadNode::NW]) return m_parent->m_child[QuadNode::NE];

   // Recursive search
   cnode = m_parent->eastNeighbor();
   // Note that cnode cannot be NULL
   if (cnode->isLeaf()) {
      return cnode;
   } else {
      if (this == m_parent->m_child[QuadNode::SE]) {
         return cnode->m_child[QuadNode::SW];
      } else {
         return cnode->m_child[QuadNode::NW];
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

bool QuadNode::needToBeSplit()
{
   QuadNode* cnode;
   // Check north
   cnode = northNeighbor();
   if (!cnode->isLeaf())
      if ( (!cnode->m_child[SW]->isLeaf()) || (!cnode->m_child[SE]->isLeaf()) ) return true;
   // Check west
   cnode = westNeighbor();
   if (!cnode->isLeaf())
      if ( (!cnode->m_child[SE]->isLeaf()) || (!cnode->m_child[NE]->isLeaf()) ) return true;
   // Check south
   cnode = southNeighbor();
   if (!cnode->isLeaf())
      if ( (!cnode->m_child[NW]->isLeaf()) || (!cnode->m_child[NE]->isLeaf()) ) return true;
   // Check east
   cnode = eastNeighbor();
   if (!cnode->isLeaf())
      if ( (!cnode->m_child[NW]->isLeaf()) || (!cnode->m_child[SW]->isLeaf()) ) return true;
   return false;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Grid point search. Returns pointer to a node and indices plocx, ploy.
// Point 0: plocx = 1, plocy = 0; Point 1: plocx = 0, plocy = 1; Point 2: plocx = 1, plocy = 1.
Point* QuadNode::findPoint(unsigned int px, unsigned int py, unsigned int pow2)
{
   // The grid step that corresponds to the current cell is determined by its level,
   // and pow2 is the corresponding power of 2. However, in a recursive search algorithm ,
   // one should be very careful where to divide it by 2. Here this is done before computing the position indices
   assert(pow2 >= 2);
   unsigned int pow = pow2 / 2;
   unsigned int jpos = 2*( py >= m_jy * pow2 + pow ) + ( px >= m_jx * pow2 + pow );
   // Check if this cells contains the point
   for (unsigned int jj=0; jj<3; ++jj) if ( (m_cpoint[jj]->m_px == px) && (m_cpoint[jj]->m_py == py) ) return m_cpoint[jj];
   // Recursively search for cells that might contain that point
   assert(m_child != NULL);
   assert(m_child[jpos] != NULL);
   return m_child[jpos]->findPoint(px,py,pow);
}

