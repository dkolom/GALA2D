/*************************************************************************************************************************
 *
 * Class QuadNode
 *
 * Node of a tree
 *
 * Adapted from QuadtreeSimulator v1.21 by Frank Dockhorn
 *
 *************************************************************************************************************************/

class Point;
class QuadTree;

class QuadNode
{
public:
   enum location // type used to enumerate 4 children of a node
   {
      SW = 0,
      SE = 1,
      NW = 2,
      NE = 3
   };

   QuadTree*      m_tree;      // Hosting tree 
   QuadNode*      m_parent;    // Parent
   QuadNode**     m_child;     // Children
   Point*         m_cpoint[3]; // Midpoints. 0:M0(SE) 1:0M(NW) 2:MM(NE). See figure 2.3 in [1]
   unsigned int   m_jx;        // Horizontal coordinate at level m_l
   unsigned int   m_jy;        // Vertical coordinate at level m_l
   unsigned char  m_l;         // Level
   unsigned char  m_hold;      // Flag that indicates if to remove this node at remeshing. 0: remove; 1: do not remove

   QuadNode(QuadTree* ntree, QuadNode* nparent, unsigned int njx, unsigned int njy, unsigned char nl);       
   ~QuadNode();      
   // Find pointer to the node at a given location
   QuadNode*      findNode(unsigned int njx, unsigned int njy, unsigned char nl); 
   // Get a corner point of this node (note that it belongs another node). clocx and cloc y take values 0 or 1
   Point*         getPoint(unsigned int clocx, unsigned int clocy); 
   // Find grid point px,py
   Point*         findPoint(unsigned int px, unsigned int py, unsigned int pow2);
   // Find north neighbour
   QuadNode*      northNeighbor();
   // Find south neighbour
   QuadNode*      southNeighbor();
   // Find west neighbour
   QuadNode*      westNeighbor();
   // Find east neighbour
   QuadNode*      eastNeighbor();
   // Check if the node needs to be split for balancing
   bool           needToBeSplit();
   // Print out node information
   void           printNode();
   // Allocate memory for pointers to children
   void           allocateChildren();
   // Check if this is a leaf
   inline bool    isLeaf() {
       if (m_child == NULL) return true; 
       for (unsigned int jj = 0; jj<4; jj++) 
       if (m_child[jj]->m_hold == 0) return true; 
       return false;
   }
};

