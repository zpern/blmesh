#include <assert.h>
#include <spdlog/spdlog.h> 
 #include "stacksrch.h"
#include "tetramesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
TetraMesh::TetraMesh() {
	m_nMaxNodeIdx = -1;
}

TetraMesh::~TetraMesh() {
}

int TetraMesh::initialize(REAL nodes[], int ndSize, int clNodes[], int clSize, bool autoNeig) {
	int i, j;

	m_nMaxNodeIdx = -1e10;
	for (i = 0; i < ndSize; i++) 
		for (j = 0; j < 3; j++)
			m_vecNodeCords.push_back(nodes[3*i + j]);
	for (i = 0; i < clSize; i++) 
		for (j = 0; j < 4; j++) {
			m_vecCellNodes.push_back(clNodes[4*i + j]);
			if (clNodes[4*i + j] > m_nMaxNodeIdx)
				m_nMaxNodeIdx = clNodes[4*i + j];
		}
	
	initialNodeCellHash(m_nMaxNodeIdx + 1);

	if (autoNeig) {
		setupCellNeig(m_nMaxNodeIdx + 1, clSize);
	}

	return 1;
}

int TetraMesh::initialize(REAL nodes[], int ndSize, int clNodes[], int clNeigs[], int clSize) {
	int i, j;
	
	m_nMaxNodeIdx = -1e10;
	for (i = 0; i < ndSize; i++) 
		for (j = 0; j < 3; j++)
			m_vecNodeCords.push_back(nodes[3*i + j]);	
	for (i = 0; i < clSize; i++) 
		for (j = 0; j < 4; j++) {
			m_vecCellNodes.push_back(clNodes[4*i + j]);
			m_vecCellNeigs.push_back(clNeigs[4*i + j]);
			if (clNodes[4*i + j] > m_nMaxNodeIdx)
				m_nMaxNodeIdx = clNodes[4*i + j];
		}
	
	initialNodeCellHash(m_nMaxNodeIdx + 1);
	return 1;
}

int TetraMesh::nodeSize() {
	return m_vecNodeCords.size() / 3;
}

int TetraMesh::nodeCoords(int iNode, REAL coords[]) {
	coords[0] = m_vecNodeCords[3*iNode    ];
	coords[1] = m_vecNodeCords[3*iNode + 1];
	coords[2] = m_vecNodeCords[3*iNode + 2];

	return 0;
}

// int TetraMesh::addNode(double coords[]) {
// 
// }

int TetraMesh::cellSize() {
	return m_vecCellNodes.size() / 4;
}

int TetraMesh::cellNodes(int iCell, int clNodes[]) {
	clNodes[0] = m_vecCellNodes[4*iCell    ];
	clNodes[1] = m_vecCellNodes[4*iCell + 1];
	clNodes[2] = m_vecCellNodes[4*iCell + 2];
	clNodes[3] = m_vecCellNodes[4*iCell + 3];

	return 0;

}

int TetraMesh::cellNeigs(int iCell, int clNeigs[]) {
	clNeigs[0] = m_vecCellNeigs[4*iCell    ];
	clNeigs[1] = m_vecCellNeigs[4*iCell + 1];
	clNeigs[2] = m_vecCellNeigs[4*iCell + 2];
	clNeigs[3] = m_vecCellNeigs[4*iCell + 3];
	
	return 0;

}

int TetraMesh::cellNeig(int iCell, int iCode) {
	return m_vecCellNeigs[4*iCell + iCode];
}

int TetraMesh::nodeCells(int node, int cells[], int *size) {
	int iCell, iNeig;
	int i = 0, m;
	INTEGER cellHd;
	StackSrch stackSrch; 
	int *visited = nullptr;
	int clSize = 0;

	clSize = cellSize();
	visited = new int[clSize];
	memset(visited, 0, sizeof(int)*clSize);

	cellHd = m_vecNodClHash[node];
	assert(cellHd >= 0);
	cells[i++] = cellHd;
	stackSrch.push(cellHd);
	visited[cellHd] = 1;
	while (!stackSrch.empty()) {
		iCell = stackSrch.top();
		stackSrch.pop();
		
		for (m = 0; m <= 3; m++) {
			iNeig = m_vecCellNeigs[4*iCell + m];
			if (iNeig >= 0 && visited[iNeig] == 0 && isCellNode(iNeig, node)) {
				if (i >= MAX_NODE_CELL_SIZE)
				{
					printf("Increase MAX_NODE_CELL_SIZE(%d), please.\n", MAX_NODE_CELL_SIZE);
					exit(1);
				}
				cells[i++] = iNeig;
				stackSrch.push(iNeig);
			}
			if (iNeig >= 0)
				visited[iNeig] = 1;
		}
	}
	
	*size = i;

	if (visited)
		delete[] visited;
	
	return 0; /* S_OK */
}

int TetraMesh::edgeCells(int node1, int node2, int cells[], int *size) {
	int nodeCls[MAX_NODE_CELL_SIZE], nodeClSize = 0;
	int i, j;

	nodeCells(node1, nodeCls, &nodeClSize);
	for (i = 0, j = 0; i < nodeClSize; i++) {
		if (isCellNode(nodeCls[i], node2))
			cells[j++] = nodeCls[i];
	}	
	*size = j;

	return 0; /* S_OK */
}

int TetraMesh::faceCells(int node1, int node2, int node3, int *lftCell, int *rgtCell){
	int code1, code2, code3;
	int nodeCls[MAX_NODE_CELL_SIZE], faceCls[MAX_NODE_CELL_SIZE], nodeClSize = 0;
	int i, j, m;

	nodeCells(node1, nodeCls, &nodeClSize);
	for (i = 0, j = 0; i < nodeClSize; i++) {
		if (isCellNode(nodeCls[i], node2) && isCellNode(nodeCls[i], node3))
			faceCls[j++] = nodeCls[i];
	}	
	
	assert(j <= 2);
	
	*lftCell = -1;
	*rgtCell = -1;
	for (i = 0; i < j; i++) {
		for (m = 0; m < 4; m++) {
			if (m_vecCellNodes[4 * faceCls[i] + m] == node1)
				code1 = m;
			else if (m_vecCellNodes[4 * faceCls[i] + m] == node2)
				code2 = m;
			else if (m_vecCellNodes[4 * faceCls[i] + m] == node3)
				code3 = m;
		}
		if (isLeftCell(code1, code2, code3))
			*lftCell = faceCls[i];
		else
			*rgtCell = faceCls[i];
	}
	
	return 0; /* S_OK */
}

int TetraMesh::initialNodeCellHash(int ndSize) {
	int i;
	m_vecNodClHash.resize(ndSize);
	for (i = 0; i < m_vecNodClHash.size(); i++)
		m_vecNodClHash[i] = -1;
	for (i = 0; i < m_vecCellNodes.size(); i++)
		m_vecNodClHash[m_vecCellNodes[i]] = i/4;

	return 0;
}

int TetraMesh::setupCellNeig(int ndSize, int clSize) {
	std::vector<int> vecRefIntFHash;
	std::vector<InterFace> vecInterFaces;
	int cellIdx, faceAdIdx, lftCell, rgtCell;
	int facNdIdx1, facNdIdx2, facNdIdx3, minFacNdIdx;
	int ndIFaces[MAX_NODE_INTERIOR_FACE_SIZE];
	int ndIFaceSize = 0;
	int i, j, k, nCommon;
	InterFace faceAd;
	const int cf[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
	
	vecRefIntFHash.resize(ndSize);
	m_vecCellNeigs.resize(4*clSize);
	for (i = 0; i < ndSize; i++)
		vecRefIntFHash[i] = -1;
	for (i = 0; i < 4*clSize; i++)
		m_vecCellNeigs[i] = -1;

	for (cellIdx = 0; cellIdx < m_vecCellNodes.size()/4; cellIdx++)
		for (i = 0; i < 4; i++) {
			facNdIdx1 = m_vecCellNodes[4*cellIdx + cf[i][0]];
			facNdIdx2 = m_vecCellNodes[4*cellIdx + cf[i][1]];
			facNdIdx3 = m_vecCellNodes[4*cellIdx + cf[i][2]];
			minFacNdIdx = facNdIdx1 < facNdIdx2 ? facNdIdx1 : facNdIdx2;
			minFacNdIdx = minFacNdIdx < facNdIdx3 ? minFacNdIdx : facNdIdx3;
			nodeInterFace(minFacNdIdx, ndIFaces, &ndIFaceSize, vecInterFaces, vecRefIntFHash);
			
			nCommon = 0;
			for (j = 0;	j < ndIFaceSize; j++) { 
				nCommon = 0;
				for (k = 0; k < 3; k++)
					if (vecInterFaces[ndIFaces[j]].form[k] == facNdIdx1 ||
						vecInterFaces[ndIFaces[j]].form[k] == facNdIdx2 ||
						vecInterFaces[ndIFaces[j]].form[k] == facNdIdx3)
						nCommon++;
				if (nCommon >= 3)
					break;
			}
			
			if (nCommon < 3) {
				/* 没有找到 */
				faceAd.form[0] = facNdIdx1;
				faceAd.form[1] = facNdIdx2;
				faceAd.form[2] = facNdIdx3;
				faceAd.lftCell = cellIdx;
				faceAd.rgtCell = -1;
				
				faceAdIdx = vecInterFaces.size();
				faceAd.hashNxt = vecRefIntFHash[minFacNdIdx];
				vecRefIntFHash[minFacNdIdx] = faceAdIdx;
				vecInterFaces.push_back(faceAd);
			}
			else {
				assert(vecInterFaces[ndIFaces[j]].lftCell >= 0 &&
					vecInterFaces[ndIFaces[j]].rgtCell < 0);
				vecInterFaces[ndIFaces[j]].rgtCell = cellIdx;
				
				lftCell = vecInterFaces[ndIFaces[j]].lftCell;
				rgtCell = vecInterFaces[ndIFaces[j]].rgtCell;		
				for (k = 0; k < 4; k++)
					if (m_vecCellNodes[4*lftCell + k] != facNdIdx1 &&
						m_vecCellNodes[4*lftCell + k] != facNdIdx2 &&
						m_vecCellNodes[4*lftCell + k] != facNdIdx3)
						break;
				assert(k < 4);
				m_vecCellNeigs[4*rgtCell + i] = lftCell;
				m_vecCellNeigs[4*lftCell + k] = rgtCell;
			}
		}
		
		
	return 0; /* S_OK */
}

int TetraMesh::nodeInterFace(int minFacNdIdx, int ndIFaces[], int *ndIFaceSize, 
						 std::vector<InterFace> &vecInterFaces, std::vector<int>& vecRefIntFHash) {
	int faceIt;
	int i;

	i = 0;
	faceIt = vecRefIntFHash[minFacNdIdx];
	while (faceIt >= 0) {
		if (i >= MAX_NODE_INTERIOR_FACE_SIZE)
		{
			spdlog::info("Increase MAX_NODE_INTERIOR_FACE_SIZE ({}).\n", MAX_NODE_INTERIOR_FACE_SIZE);
			exit(1);
		}
		ndIFaces[i++] = faceIt;
		faceIt = vecInterFaces[faceIt].hashNxt;
	}
	*ndIFaceSize = i;
	
	return 0; /* S_OK */
}

int TetraMesh::isCellNode(int iCell, int iNode) {
	return m_vecCellNodes[4*iCell    ] == iNode ||
		   m_vecCellNodes[4*iCell + 1] == iNode ||
		   m_vecCellNodes[4*iCell + 2] == iNode ||
		   m_vecCellNodes[4*iCell + 3] == iNode;
}

bool TetraMesh::isLeftCell(int code1, int code2, int code3) {
	/* 0 left; 1 right; 2 incorrect */
	const char leftCellFlag[4][4][4] = {
		{{2, 2, 2, 2}, {2, 2, 0, 1}, {2, 1, 2, 0}, {2, 0, 1, 2}},
		{{2, 2, 1, 0}, {2, 2, 2, 2}, {0, 2, 2, 1}, {1, 2, 0, 2}},
		{{2, 0, 2, 1}, {1, 2, 2, 0}, {2, 2, 2, 2}, {0, 1, 2, 2}},
		{{2, 1, 0, 2}, {0, 2, 1, 2}, {1, 0, 2, 2}, {2, 2, 2, 2}}
	};
	assert(leftCellFlag[code1][code2][code3] < 2);

	return leftCellFlag[code1][code2][code3] == 0;
}