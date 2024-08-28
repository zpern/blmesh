#ifndef __tri_mesh_h__
#define __tri_mesh_h__
#include <vector>
#include <cstdio>
#include <cstring>
#include <set>
namespace TriMesh  {

enum {MAX_NODE_EDGE_SIZE = 128};
enum {MAX_EDGE_FACE_SIZE = 128};



typedef struct BoundFace {
	int form[3];	/* 3个形成节点 */
//	int edge[3];	/* 3条边 */
	int neig[3];	/* 该面的3个邻居 */
} BoundFace;
	
typedef struct BoundEdge {
	int form[2];	/* 起点和终点 */
	int face;		/* 包含该边的面片 */
	std::set<int> faces;
 	int hashNxt;	/* 用于边哈希表，指向下一条哈希边 */
} BoundEdge;

class StaticMesh {
public:
	StaticMesh();
	~StaticMesh();

	int initialize(
		int nbp,			/* number of boundary points */
		double *lcoord,     /* coord. of boundary points , size = 3 * nbp */
		int nbf,		    /* number of boundary surface */
		int *ibnd			/* topu. of boundary entities, size = 3 * nbf (i1, i2, i3) */);

	bool isManifold();       

	int edgeColoring(int faceColor[], int *numOfColors);
	
	int natureOrient(int faceSign[]);	

	int faceNodes(int faceIdx, int faceNodes[]);
	int nodeCords(int nodeIdx, double nodeCords[]);
	
	int nodeSize() {return m_vecBoundNodes.size()/3 - 1;}
	int edgeSize() {return m_vecBoundEdges.size() - 1;}
	int faceSize() {return m_vecBoundFaces.size() - 1;}
protected:
	/* 边界点周围边编号 */
	int nodeEdges(int nodeIdx, int nodeEdges[], int *nodeEdgeSize);
	int edgeFaces(int edgeIdx, int edgeFaces[], int *edgeFaceSize);
protected:
	std::vector<double>    m_vecBoundNodes; /* 基于1 */
	std::vector<BoundEdge> m_vecBoundEdges; /* 基于1 */
	std::vector<BoundFace> m_vecBoundFaces; /* 基于1 */
	std::vector<int>	   m_vecRefEdgeHash;	/* 包含节点的第一条边, 基于1 */
};
};

#endif
