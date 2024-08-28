#include "trimesh.h"
#include <spdlog/spdlog.h> 
 #include <assert.h>
#include <stack>
#include<stdexcept>
#include <iostream>
#include <stdexcept>
using namespace std;
namespace TriMesh
{
StaticMesh::StaticMesh() {
}

StaticMesh::~StaticMesh() {
	m_vecBoundNodes.clear();
	m_vecBoundEdges.clear();
	m_vecBoundFaces.clear();
	m_vecRefEdgeHash.clear();
}

int StaticMesh::initialize( int nbp, double *lcoord, int nbf, int *ibnd) {
	int faceIdx, edgeHdIdx, edgeAdIdx, edNdIdx1, edNdIdx2, minedNdIdx, maxedNdIdx;
	int faceIt, faceHd, faceTl;
	int nodeEdges[MAX_NODE_EDGE_SIZE];
	int nodeEdgeSize;
	BoundEdge edgeAd, zeroEdge;
	BoundFace zeroFace;
	int i, j, k;
	static int fe[3][2] = {{1, 2}, {2, 0}, {0, 1}};
	
	/* 初始化面片的节点坐标 */
	m_vecBoundNodes.resize(3*nbp + 3);
	for (i = 0; i < 3*nbp; i++) {
		m_vecBoundNodes[i+3] = lcoord[i];
	}

	/* 初始化面片的形成节点 */
	memset(&zeroFace, 0, sizeof(BoundFace));
	m_vecBoundFaces.resize(nbf + 1, zeroFace);
 	for (i = 0; i < nbf; i++) {
		m_vecBoundFaces[i+1].form[0] = ibnd[3*i + 0] + 1;
		m_vecBoundFaces[i+1].form[1] = ibnd[3*i + 1] + 1;
		m_vecBoundFaces[i+1].form[2] = ibnd[3*i + 2] + 1;
		m_vecBoundFaces[i].neig[0] = -1;
		m_vecBoundFaces[i].neig[1] = -1;
		m_vecBoundFaces[i].neig[2] = -1;
 	}
 	
	/* 构建面片的拓扑关系 */
	m_vecRefEdgeHash.resize(nbp + 1, 0);
	zeroEdge.face = 0;
	zeroEdge.hashNxt = 0;
	m_vecBoundEdges.push_back(zeroEdge);
	for (faceIdx = 1; faceIdx < m_vecBoundFaces.size(); faceIdx++)
		for (i = 0; i < 3; i++) {
			edNdIdx1 = m_vecBoundFaces[faceIdx].form[fe[i][0]];
			edNdIdx2 = m_vecBoundFaces[faceIdx].form[fe[i][1]];
			
			minedNdIdx = edNdIdx1 > edNdIdx2 ? edNdIdx2 : edNdIdx1;
			maxedNdIdx = edNdIdx1 < edNdIdx2 ? edNdIdx2 : edNdIdx1;
			edNdIdx1 = minedNdIdx;
			edNdIdx2 = maxedNdIdx;
			
			/* 更新节点edNdIdx1对应的边链表 */
			this->nodeEdges(edNdIdx1, nodeEdges, &nodeEdgeSize);
    		for (j = 0; j < nodeEdgeSize; j++)
				if (m_vecBoundEdges[nodeEdges[j]].form[0] == edNdIdx2 ||
					m_vecBoundEdges[nodeEdges[j]].form[1] == edNdIdx2) {
					/* 找到这条边 */
					break;
				}
			if (j >= nodeEdgeSize) {
				/* 没有找到这条边，创建这条边 */
				edgeAd.form[0] = edNdIdx1;
				edgeAd.form[1] = edNdIdx2;
				edgeAd.face = faceIdx;
				edgeAdIdx = m_vecBoundEdges.size();
				
				/* 更新edNdIdx1对应的边链表 */
				edgeHdIdx = m_vecRefEdgeHash[edNdIdx1];
				edgeAd.hashNxt = edgeHdIdx > 0 ? edgeHdIdx : 0;
				m_vecRefEdgeHash[edNdIdx1] = edgeAdIdx;

				m_vecBoundEdges.push_back(edgeAd);
				m_vecBoundEdges.rbegin()->faces.insert(faceIdx);
				m_vecBoundFaces[faceIdx].neig[i] = faceIdx;
				//m_vecBoundFaces[faceIdx].edge[i] = edgeAdIdx;
			}
			else {
				/* 将faceIdx增加到包含边的面片链表中 */
				faceIt = faceHd = m_vecBoundEdges[nodeEdges[j]].face;
				faceTl = -1;
				if (faceHd > 0) {
					do { 
						faceTl = faceIt;
						for (k = 0; k < 3; k++)
							if (m_vecBoundFaces[faceIt].form[k] != edNdIdx1 &&
								m_vecBoundFaces[faceIt].form[k] != edNdIdx2) {
								faceIt = m_vecBoundFaces[faceIt].neig[k];
								break;
							}
						assert(k < 3);
					} while (faceIt != faceHd);
				}
				assert(faceTl > 0);
				m_vecBoundEdges[nodeEdges[j]].faces.insert(faceIdx);
				m_vecBoundFaces[faceTl].neig[k] = faceIdx;
				m_vecBoundFaces[faceIdx].neig[i] = faceHd;

// 				for (k = 0; k < 3; k++)
// 					if (m_vecBoundFaces[faceIt].neig[k] == faceHd) {
// 						
// 						break;
// 					}
			}
		}
	
	for (auto i : m_vecBoundEdges) {
		if (i.faces.size() > 2) {
			for (auto j : i.faces) {
				for (auto k : i.faces) {
					for (int l = 0; l < 3; l++) {
						if (m_vecBoundFaces[j].neig[l] == k) {
							m_vecBoundFaces[j].neig[l] = -1;
						}
					}
				}
			}
		}
	}
	return 0; /* S_OK */
}

bool StaticMesh::isManifold() {
	int edgeIdx, edgeFaces[MAX_EDGE_FACE_SIZE], edgeFaceSize = 0;

	for (edgeIdx = 1; edgeIdx < m_vecBoundEdges.size(); edgeIdx++) {
		this->edgeFaces(edgeIdx, edgeFaces, &edgeFaceSize);
 		if (edgeFaceSize != 2) {
			printf("edgeface size=%d\n",edgeFaceSize);
			printf("Nonmanifold edge (%d-->%d): No. of faces around it is %d.",
				m_vecBoundEdges[edgeIdx].form[0] - 1, 
				m_vecBoundEdges[edgeIdx].form[1] - 1,
				edgeFaceSize);
			cout << endl;
			for (int j = 0; j < edgeFaceSize; j++) {
				cout << edgeFaces[j] << " ";
			}
			cout<<endl;
			throw std::runtime_error("Nonmanifold edge.");
			return false;
		}
	}

	return true;
}

/* 依照边相邻关系对面片进行染色 */
int StaticMesh::edgeColoring(int faceColor[], int *numOfColors) {
	int faceIdx = 0, color = 0, k, neigIdx = 0, seedIdx;
	std::stack<int> stackSrch;
	
	memset(faceColor, 0, sizeof(int)*faceSize());
	for (faceIdx = 1; faceIdx < m_vecBoundFaces.size(); faceIdx++) {
		if (faceColor[faceIdx - 1] > 0)
			continue;
		color++;
		stackSrch.push(faceIdx);
		faceColor[faceIdx - 1] = color;
		while (!stackSrch.empty()) {
			seedIdx = stackSrch.top();
			stackSrch.pop();
			for (k = 0; k < 3; k++) {
				neigIdx = m_vecBoundFaces[seedIdx].neig[k];
				if (neigIdx > 0 && faceColor[neigIdx - 1] <= 0) {
					stackSrch.push(neigIdx);
					faceColor[neigIdx - 1] = color;
				}
			}
		}
	}

	*numOfColors = color;

	return 0;
}

/* 使得边在共享其的两个面中的绕向相反，仅针对二边流形 
 * faceSign[i]: 1, 正向；-1：反向；0：不确定 */
/*yhf尝试将其改成适应非二边流型的形式，其实对于总体流通的形状，只需绕过非二边流型即可 2021/7/10*/
int StaticMesh::natureOrient(int faceSign[]) {
	int faceIdx = 0, neigIdx = 0, seedIdx, sign;
	int n1, n2, k, m;
	std::stack<int> stackSrch;
	
	memset(faceSign, 0, sizeof(int)*faceSize());
	//zlj 2011.12.10
// 	faceIdx = 1;
// 	stackSrch.push(faceIdx);
//     do 
//     {
// 	  seedIdx = stackSrch.top();
// 	  stackSrch.pop();
// 	  faceSign[seedIdx - 1] = 1;
// 	  for (k = 0; k < 3; k++)
// 	  {
// 		  n1 = m_vecBoundFaces[seedIdx].form[(k+1) % 3];
// 		  n2 = m_vecBoundFaces[seedIdx].form[(k+2) % 3];
// 		  neigIdx =m_vecBoundFaces[seedIdx].neig[k];
// 		  stackSrch.push(neigIdx);
// 		  if (neigIdx > 0 && faceSign[neigIdx - 1] != 0)
// 		  {
// 			  stackSrch.pop();
// 			  for (m = 0; m < 3; m++)
// 			  {
// 				  if (m_vecBoundFaces[neigIdx].form[m] != n1 &&
// 					  m_vecBoundFaces[neigIdx].form[m] != n2)
// 				  {
// 					  	faceSign[seedIdx - 1] = (m_vecBoundFaces[neigIdx].form[(m+1)%3] == n2) ? faceSign[neigIdx - 1] : - faceSign[neigIdx - 1];
// 						break;
// 
// 				  }
// 				  assert(m < 3);
// 			  }
// 		  }
// 	  }
// 	 
//     } while (!stackSrch.empty());
	for (faceIdx = 1; faceIdx < m_vecBoundFaces.size(); faceIdx++) {
		if (faceSign[faceIdx - 1] != 0)
			continue;
		stackSrch.push(faceIdx);
		faceSign[faceIdx - 1] = 1;
		while (!stackSrch.empty()) {
			seedIdx = stackSrch.top();
			stackSrch.pop();
			sign = faceSign[seedIdx - 1];
			assert(sign != 0);
			for (k = 0; k < 3; k++) {
				n1 = m_vecBoundFaces[seedIdx].form[(k+1)%3];
				n2 = m_vecBoundFaces[seedIdx].form[(k+2)%3];
				neigIdx = m_vecBoundFaces[seedIdx].neig[k];
				if (neigIdx > 0 && faceSign[neigIdx - 1] == 0) {
					stackSrch.push(neigIdx);
					for (m = 0; m < 3; m++)
						if (m_vecBoundFaces[neigIdx].form[m] != n1 &&
							m_vecBoundFaces[neigIdx].form[m] != n2)
						{
							faceSign[neigIdx - 1] = (m_vecBoundFaces[neigIdx].form[(m+1)%3] == n2) ? sign : -sign;
							break;
						}
					assert(m < 3);
				}
			}
		}
	}
	
	return 0;
}

/* 面片的形成节点编号 */
int StaticMesh::faceNodes(int faceIdx, int faceNodes[]) {
	faceNodes[0] = m_vecBoundFaces[faceIdx + 1].form[0] - 1;
	faceNodes[1] = m_vecBoundFaces[faceIdx + 1].form[1] - 1;
	faceNodes[2] = m_vecBoundFaces[faceIdx + 1].form[2] - 1;
	return 0; /* S_OK */
}

/* 节点坐标 */
int StaticMesh::nodeCords(int nodeIdx, double nodeCords[]) {
	nodeCords[0] = m_vecBoundNodes[3*(nodeIdx + 1) + 0];
	nodeCords[1] = m_vecBoundNodes[3*(nodeIdx + 1) + 1];
	nodeCords[2] = m_vecBoundNodes[3*(nodeIdx + 1) + 2];
	return 0; /* S_OK */
}

/* 边界点周围边编号 */
int StaticMesh::nodeEdges(int nodeIdx, int nodeEdges[], int *nodeEdgeSize) {
	int edgeHd, edgeIt;
	int i = 0;
	edgeIt = edgeHd = m_vecRefEdgeHash[nodeIdx];
	if (edgeIt > 0) {
		do {
			nodeEdges[i++] = edgeIt;
			edgeIt = m_vecBoundEdges[edgeIt].hashNxt;
		} while (edgeIt > 0);
	}
	*nodeEdgeSize = i;
	return 0; /* S_OK */
}

int StaticMesh::edgeFaces(int edgeIdx, int edgeFaces[], int *edgeFaceSize) {
	/* 边界边周围面编号 */
	int faceHd, faceIt, n1, n2;
	int i = 0, k;
	faceIt = faceHd = m_vecBoundEdges[edgeIdx].face;
	if (faceHd > 0) {
		n1 = m_vecBoundEdges[edgeIdx].form[0];
		n2 = m_vecBoundEdges[edgeIdx].form[1];
		do {
			edgeFaces[i++] = faceIt;
			for (k = 0; k < 3; k++)
				if (m_vecBoundFaces[faceIt].form[k] != n1 &&
					m_vecBoundFaces[faceIt].form[k] != n2) {
					faceIt = m_vecBoundFaces[faceIt].neig[k];
					break;
				}
				assert(k < 3);
		} while (faceIt != faceHd);
	}
	*edgeFaceSize = i;
	
	return 0; /* S_OK */
}

}; /* namespace */
