/*************************************************************************************************************************
 *
 * Class QuadTree
 *
 * Data and functions that operate the tree (not just a node)
 *
 * Adapted from QuadtreeSimulator v1.21 by Frank Dockhorn
 *
 *************************************************************************************************************************/

class QuadNode;
class Point;
class Solver;

class QuadTree
{
private:   
   QuadNode*     m_root;      // Pointer to the root	

public:
   Solver*       m_solver;    // Pointer to the solver that uses this tree
   QuadNode**    m_leaf;      // Pointers to leaves
   unsigned int  m_nleaf;     // Number of leaves
   Point*        m_rootpoint; // Lower-left corner point that does not belong to any cell
   Point*        m_headpoint; // Pointer to list of grid points

   QuadTree(Solver* csolver);
   ~QuadTree();

   // Create new tree
   void      createTree(); 
   // Delete tree
   void      cleanup();
   // Add a new node at a given position and add all intermediary nodes to make a connected tree
   void      addBranch(unsigned int njx, unsigned int njy, unsigned char nl);
   // Remove all child branches
   void      removeChildren(QuadNode* nnode);
   // Print out all nodes of a branch
   void      printBranch(QuadNode* nnode);
   // Link all points in a list
   void      linkPoints(QuadNode* nnode);
   // Remove all cells marked for removal
   void      updateTree(QuadNode* nnode);
   // Generate a list of leaves
   void      listLeaves(QuadNode* nnode, unsigned int& nnleaf); 
   // Add one layer of leaves
   void      addLeaves(); 
   // Count leaves
   void      countLeaves(QuadNode* nnode, unsigned int& nnleaf);
   // Count and list leaves
   void      createListLeaves();
   // Modify removal flags 'QuadNode::m_hold' to ensure that the tree is graded
   void      balance();
   // Set all removal flags 'QuadNode::m_hold' to 1 
   void      holdTree(QuadNode* nnode);
   // Count nodes in the tree
   void      countNodes(QuadNode* nnode, unsigned int& counter); 
   // Minimum and maximum leaf level
   void      levelStats(QuadNode* nnode, unsigned char& lleafmin, unsigned char& lleafmax);
   // Returns pointer to root
   QuadNode* getRoot() {return m_root;}
   // Search for a node with given coordinates
   QuadNode* searchTree(unsigned int njx, unsigned int njy, unsigned char nl);
   // Find leaf containing point px,py
   QuadNode* findLeafXY_int(QuadNode* nnode, unsigned int px, unsigned int py, unsigned int pow2);
   // Find grid point px,py
   QuadNode* findPoint_CB(QuadNode* nnode, unsigned int px, unsigned int py, unsigned int pow2, unsigned int& plocx, unsigned int& plocy);
};

