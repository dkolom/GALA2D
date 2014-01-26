/*************************************************************************************************************************
 *
 * Class QuadTree
 *
 * Data and functions that operate the tree (not just a node)
 *
 * Adapted from QuadtreeSimulator v1.21 by Frank Dockhorn
 *
 *************************************************************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <climits>
#include <forward_list>
#include "Point.h"
#include "QuadNode.h"
#include "QuadTree.h"
#include "Solver.h" 
using namespace std;



/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadTree::QuadTree(Solver* csolver)
{
   // Set pointer to the solver
   m_solver = csolver;
   // Create tree
   createTree();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadTree::~QuadTree()
{
   // Remove the tree
   cleanup();
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void QuadTree::cleanup()
{
   delete m_root;
   m_root = NULL;

   if (m_leaf != NULL) {
      delete[] m_leaf;
      m_leaf = NULL;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

void QuadTree::createTree()
{
   // Initialize root point
   double ph = m_solver->m_size / double(m_solver->m_jmax);
   unsigned int px = 0;
   unsigned int py = 0;
   m_rootpoint = new Point(px,py,ph,NULL);

   // Create root cell
   m_root = new QuadNode(this,NULL,0,0,0);
   m_nleaf = 0;
   m_leaf = NULL;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// the four sub regions of the quad tree
//+--------------+
//|  nw  |  ne   |
//c------c-------|
//|  sw  |  se   |
//+------c-------+
// c - m_cpoint[0:2]   in QuadNode
//
// Add a new node
// nl - level of the new node
// njx - horizontal location at level nl
// njy - vertical location at level nl
// If node already exists - nothing is changed
// Intermediate nodes are added, if necessary, to make a connected tree
void QuadTree::addBranch(unsigned int njx, unsigned int njy, unsigned char nl)
{
   QuadNode* cnode = m_root;
   QuadNode::location cloc;
   unsigned int cjx;
   unsigned int cjy;
   unsigned char cl = 0; 
   unsigned int njmax = (unsigned int)(floor(pow(double(2),double(nl))+0.5));

   while (cl<nl)
   {
      njmax /= 2;
      cjx = njx / njmax;
      cjy = njy / njmax;
      cl++;

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

      // If no children, create
      // Allocate pointers to children
      if (cnode->m_child == NULL) {
         cnode->allocateChildren();
      }
      // If child doesn't exist, create
      if (cnode->m_child[cloc] == NULL) {
         cnode->m_child[cloc] = new QuadNode(this,cnode,cjx,cjy,cl);
      }

      // Set current node to the new position
      cnode = cnode->m_child[cloc];

   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Add one layer of leaves
void QuadTree::addLeaves()
{
   // Only works if list of leaves exists
   assert(m_leaf!=NULL);

   // Loop for all leaves
   for (unsigned int jj=0; jj<m_nleaf; jj++) {
      // Ensure thin is indeed a leaf
      assert(m_leaf[jj]->m_child==NULL);
      // Allocate pointers to children
      m_leaf[jj]->allocateChildren();
      // Loop for all children
      for (int clocx=0; clocx<2; clocx++) 
      for (int clocy=0; clocy<2; clocy++) {
         unsigned char cl = m_leaf[jj]->m_l+1;
         unsigned int cjx = m_leaf[jj]->m_jx*2+clocx;
         unsigned int cjy = m_leaf[jj]->m_jy*2+clocy;
         m_leaf[jj]->m_child[clocx+2*clocy] = new QuadNode(this,m_leaf[jj],cjx,cjy,cl);
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// List of leaves
// ONLY INCLUDES LEAVES THAT ARE QuadNode::m_hold > 0
void QuadTree::listLeaves(QuadNode* nnode, unsigned int& nnleaf)
{
   // Assume this node is a leaf
   int thisIsLeaf = 1;

   // List all 4 children
   // OLD VERSION:  if (nnode->m_child != NULL) {
   if (!nnode->isLeaf()) {
      for (int jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         // This node is not a leaf
         thisIsLeaf = 0;
         // Recursively call all children etc.
         listLeaves(nnode->m_child[jj],nnleaf);
      }
   }

   // Put this node on the list if it is a leaf
   if (thisIsLeaf) {
      m_leaf[nnleaf] = nnode;
      nnleaf++;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Count leaves
// ONLY INCLUDES LEAVES THAT ARE QuadNode::m_hold > 0
void QuadTree::countLeaves(QuadNode* nnode, unsigned int& nnleaf)
{
   // Assume this node is a leaf
   int thisIsLeaf = 1;

   // List all 4 children
   // OLD VERSION:  if (nnode->m_child != NULL) {
   if (!nnode->isLeaf()) {
      for (int jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         // This node is not a leaf
         thisIsLeaf = 0;
         // Recursively call all children etc.
         countLeaves(nnode->m_child[jj],nnleaf);
      }
   }

   // Count this node if it is a leaf
   if (thisIsLeaf) nnleaf++;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Count and list leaves 
// ONLY INCLUDES LEAVES THAT ARE QuadNode::m_hold > 0
void QuadTree::createListLeaves()
{
   // Count leaves
   unsigned int nleaf = 0;
   countLeaves(m_root,nleaf);
   m_nleaf = nleaf;

   // Deallocate old list
   if (m_leaf != NULL) {
      delete[] m_leaf;
      m_leaf = NULL;
   }
   // Allocate list 
   m_leaf = new QuadNode*[m_nleaf];

   // Fill list
   nleaf = 0;
   listLeaves(m_root,nleaf);
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

QuadNode* QuadTree::findLeafXY_int(QuadNode* nnode, unsigned int px, unsigned int py, unsigned int pow2)
{
   QuadNode* snode = NULL;
   pow2 = pow2/2;
   int jpos = 2*( 2 * py >= (2 * nnode->m_jy + 1) * pow2 ) + ( 2 * px >= (2 * nnode->m_jx + 1) * pow2 );

   // This function is only called after a security zone has been added. It finds cells that were leaves before adding the security zone.
//   if (nnode->m_child[jpos]->m_child != NULL) {
   if (nnode->m_child != NULL) {
//   if (!nnode->isLeaf()) {
      snode = findLeafXY_int(nnode->m_child[jpos], px, py, pow2);
   } 
   else {
      snode = nnode;
   }

   return snode;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Find node with given coordinates in a tree
QuadNode* QuadTree::searchTree(unsigned int njx, unsigned int njy, unsigned char nl)
{
   QuadNode* fnode = this->getRoot()->findNode(njx,njy,nl);
   return fnode;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/

// Remove branch beginning from node nnode
void QuadTree::removeChildren(QuadNode* nnode)
{
   if (nnode->m_child != NULL) {
      for (int j = 0; j < 4; j++) {
         if (nnode->m_child[j] != NULL) removeChildren(nnode->m_child[j]);

         delete nnode->m_child[j];
         nnode->m_child[j] = NULL;
      }
      delete[] nnode->m_child;
      nnode->m_child = NULL;
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Print branch coordinates
void QuadTree::printBranch(QuadNode* nnode)
{
   // Print current node information
   nnode->printNode();

   // List all 4 children
   if (nnode->m_child != NULL) {
      for (int jj = 0; jj < 4; jj++) {
         if (nnode->m_child[jj] != NULL) {
            // Recursively print all children etc.
            printBranch(nnode->m_child[jj]);
         }
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Link points to a list
void QuadTree::linkPoints(QuadNode* nnode)
{
   // Link three points of this cell
   nnode->m_cpoint[0]->m_next = m_headpoint;
   nnode->m_cpoint[1]->m_next = nnode->m_cpoint[0];
   nnode->m_cpoint[2]->m_next = nnode->m_cpoint[1];
   m_headpoint = nnode->m_cpoint[2];

   // List all 4 children
   if (nnode->m_child != NULL) {
      for (int jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         // Recursively link all children's points
         linkPoints(nnode->m_child[jj]);
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Remove all cells marked for removal
void QuadTree::updateTree(QuadNode* nnode)
{
   int jj;

   // This node cannot be removed,
   // because its children may be removed
   nnode->m_hold = 1;

   // List all 4 children
   if (nnode->m_child != NULL) {
      unsigned char hold = 0;
      for (jj = 0; jj < 4; jj++) {
         if (nnode->m_child[jj] != NULL) if (nnode->m_child[jj]->m_hold == 1) hold = 1;
      }

      // Remove child nodes or consider them for coarsening at finer level
      if (hold == 0) {
         removeChildren(nnode); 
      } else {
         for (jj = 0; jj < 4; jj++) {
            if (nnode->m_child[jj] != NULL) {
               nnode->m_child[jj]->m_hold = 1;
               // Search tree to remove marked nodes
               updateTree(nnode->m_child[jj]);
            }
         }
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Balance to make a graded tree
void QuadTree::balance()
{
   unsigned int jj;
   QuadNode* cnode;
   QuadNode* snode;

   // Local list of pointers to nodes
   forward_list<QuadNode*> listleaf;
   // Update list of leaves (!!! MAY BE UNNECESSARY !!!)
   createListLeaves();
   // Copy leave pointers to the local list
   for (jj=0; jj<m_nleaf; jj++) listleaf.push_front(m_leaf[jj]);

   // Iterate over the list and subdivide leaves if required
   while (!listleaf.empty()) {
      // Pop the first element
      cnode = listleaf.front();
      listleaf.pop_front();
      // Subdivide if required
      if ( (cnode->isLeaf()) && (cnode->needToBeSplit()) ) {
         // Unmark children
         assert(cnode->m_child != NULL);
         for (jj=0; jj<4; jj++) cnode->m_child[jj]->m_hold = 1;
         // Insert the four new cells into the list
         for (jj=0; jj<4; jj++) listleaf.push_front(cnode->m_child[jj]);
         // Check if the four neighbors of 'cnode' need to be subdivided; if so, insert it/them into the list 'listleaf'
         // North neighbor:
         snode = cnode->northNeighbor();
         if (snode->needToBeSplit()) listleaf.push_front(snode);
         // West neighbor:
         snode = cnode->westNeighbor();
         if (snode->needToBeSplit()) listleaf.push_front(snode);
         // South neighbor:
         snode = cnode->southNeighbor();
         if (snode->needToBeSplit()) listleaf.push_front(snode);
         // East neighbor:
         snode = cnode->eastNeighbor();
         if (snode->needToBeSplit()) listleaf.push_front(snode);
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Unset removal flag
void QuadTree::holdTree(QuadNode* nnode)
{
   int jj;

   // This node cannot be removed,
   // because its children may be removed
   nnode->m_hold = 1;

   // List all 4 children
   if (nnode->m_child != NULL) {
      for (jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         holdTree(nnode->m_child[jj]);
      }
   }
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Count nodes in a tree
void QuadTree::countNodes(QuadNode* nnode, unsigned int& counter)
{
   // List all 4 children
   if (nnode->m_child != NULL) {
      for (int jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         // Recursively call all children etc.
         countNodes(nnode->m_child[jj],counter);
      }
   }

   // Count this node 
   counter++;
}

/*-------------------------------------------------------------------------*\
\*-------------------------------------------------------------------------*/
// Minimum and maximum level of leaves in the tree
void QuadTree::levelStats(QuadNode* nnode, unsigned char& lleafmin, unsigned char& lleafmax)
{
   // Initialize
   if (nnode == m_root) {
      lleafmax = 0;
      lleafmin = m_solver->m_lmax+1;
   }

   // List all 4 children
   if (nnode->m_child != NULL) {
      // Recursively search if this is not a leaf
      for (int jj = 0; jj < 4; jj++) {
         assert(nnode->m_child[jj] != NULL);
         // Recursively call all children etc.
         levelStats(nnode->m_child[jj],lleafmin,lleafmax);
      }
   } else {
     // Check if max and min change
     if (nnode->m_l < lleafmin) lleafmin = nnode->m_l;
     if (nnode->m_l > lleafmax) lleafmax = nnode->m_l;
   }
}



