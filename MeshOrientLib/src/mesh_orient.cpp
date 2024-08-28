#include "mesh_orient.h"
#include "trimesh.h"
#include "vector.h"
#include <spdlog/spdlog.h> 
 #include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <map>
namespace MeshOrient {
int removeBox(
	int nbp,			/* number of boundary points */
	double *lcoord,     /* coord. of boundary points , size = 3 * nbp */
	int &nbf,		    /* number of boundary surface */
	int *ibnd          /* topu. of boundary entities, size = 3 * nbf (i1, i2, i3) */) 
{
	TriMesh::StaticMesh mesh;
	int *faceColor = NULL, *faceSign = NULL, numOfColors = 0;
	double shellVolume[MAX_SHELL_SIZE], maxVolume, nodeCords[3][3];
	int shellSign[MAX_SHELL_SIZE], maxShell, shellIdx;
	int faceIdx, faceNodes[3];
	int k, t;
	int nRet = 0;
	BLVector nodeVect[3];

	assert(nbf > 0);
	mesh.initialize(nbp, lcoord, nbf, ibnd);
	if (!mesh.isManifold()) {
		nRet = 2;
	}

	if (!(faceColor = new int[nbf])) {
		throw(std::string("not enough memory"));
	}
	mesh.edgeColoring(faceColor, &numOfColors);

	faceSign = new int[nbf];
	mesh.natureOrient(faceSign);

	memset(shellVolume, 0.0, sizeof(double)*MAX_SHELL_SIZE);
	for (faceIdx = 0; faceIdx < mesh.faceSize(); faceIdx++) {
		shellIdx = faceColor[faceIdx];
		mesh.faceNodes(faceIdx, faceNodes);
		if (faceSign[faceIdx] == 1) {
			t = faceNodes[1];
			faceNodes[1] = faceNodes[2];
			faceNodes[2] = t;
		}
		for (k = 0; k < 3; k++) {
			mesh.nodeCords(faceNodes[k], nodeCords[k]);
			nodeVect[k].x = nodeCords[k][0];
			nodeVect[k].y = nodeCords[k][1];
			nodeVect[k].z = nodeCords[k][2];
		}
		shellVolume[shellIdx] += (nodeVect[1] ^ nodeVect[0])*nodeVect[2];
	}


	maxShell = 0;
	maxVolume = 0.0;
	for (shellIdx = 1; shellIdx <= numOfColors; shellIdx++) {
		if (fabs(shellVolume[shellIdx]) >= maxVolume) {
			maxVolume = fabs(shellVolume[shellIdx]);
			maxShell = shellIdx;
		}
	}
	for (shellIdx = 1; shellIdx <= numOfColors; shellIdx++) {
		if (shellVolume[shellIdx] == 0.0)
		{
			shellSign[shellIdx] = 0;
			nRet = 1; /* S_FALSE */
		}
		else
			shellSign[shellIdx] =
			(((shellIdx != maxShell) && shellVolume[shellIdx] < 0.0) ||
			((shellIdx == maxShell) && shellVolume[shellIdx] > 0.0)) ? 1 : -1;
	}

	int count = 0;
	for (int i = 0; i < nbf; i++) {
		if (faceColor[i] != maxShell) {
			ibnd[3 * count] = ibnd[3 * i];
			ibnd[3 * count+1] = ibnd[3 * i+1];
			ibnd[3 * count+2] = ibnd[3 * i+2];
			count++;
		}
	}
	nbf = count;

	return 0;
	
	
}
//#pragma optimize( "", off )
int orientTriangularSurface(
	int nbp,			/* number of boundary points */
	double *lcoord,     /* coord. of boundary points , size = 3 * nbp */
	int nbf,		    /* number of boundary surface */
	int *ibnd,          /* topu. of boundary entities, size = 3 * nbf (i1, i2, i3) */
	int *sign,			/* 1: positive; -1: negative; 0: undetermined */
	double *volume,		/* volume */
	std::vector<double>& holes
) 
{
	TriMesh::StaticMesh mesh;
	int *faceColor = NULL, *faceSign = NULL, numOfColors = 0;
	std::map<int,double> shellVolume;
	double maxVolume, nodeCords[3][3];
	std::map<int,int> shellSign;
	int maxShell, shellIdx;
	int faceIdx, faceNodes[3];
	int k, t;
	int nRet = 0;
	BLVector nodeVect[3];


	assert(nbf > 0);
	mesh.initialize(nbp, lcoord, nbf, ibnd);
	//try { mesh.isManifold(); }
	//catch (std::exception e) {
	//	spdlog::info("no manifold mesh");
	//}


	faceColor = new int[nbf];
	mesh.edgeColoring(faceColor, &numOfColors);

	
	faceSign = new int[nbf];

	mesh.natureOrient(faceSign);
	//There we use averagy coordinate to define the hole,not robust
	holes.resize(numOfColors*3);
	std::vector<int> count(numOfColors);
	//memset(shellVolume, 0.0, sizeof(double)*MAX_SHELL_SIZE);
	for (faceIdx = 0; faceIdx < mesh.faceSize(); faceIdx++) {
		shellIdx = faceColor[faceIdx];
		mesh.faceNodes(faceIdx, faceNodes);
		if (faceSign[faceIdx] == 1) {
			t = faceNodes[1];
			faceNodes[1] = faceNodes[2];
			faceNodes[2] = t;
		}
		for (k = 0; k < 3; k++) {
			mesh.nodeCords(faceNodes[k], nodeCords[k]);
			nodeVect[k].x = nodeCords[k][0];
			nodeVect[k].y = nodeCords[k][1];
			nodeVect[k].z = nodeCords[k][2];
		}
		holes[shellIdx * 3 - 3 + 0] += nodeVect[1].x;
		holes[shellIdx * 3 - 3 + 1] += nodeVect[1].y;
		holes[shellIdx * 3 - 3 + 2] += nodeVect[1].z;
		count[shellIdx - 1]++;
		//spdlog::info("shellIdx=" << shellIdx);
		shellVolume[shellIdx] += (nodeVect[1]^nodeVect[0])*nodeVect[2];
	}
	
	for (int i = 0; i < holes.size()/3; i++) {
		holes[i * 3+ 0] /= count[i];
		holes[i * 3  + 1] /= count[i];
		holes[i * 3  + 2] /= count[i];
	}
	maxShell = 0;
	maxVolume = 0.0;
	for (shellIdx = 1; shellIdx <= numOfColors; shellIdx++) {
		if (fabs(shellVolume[shellIdx]) >= maxVolume) {
			maxVolume = fabs(shellVolume[shellIdx]);
			maxShell = shellIdx;
		}
	}
	for (shellIdx = 1; shellIdx <= numOfColors; shellIdx++) {
		if (shellVolume[shellIdx] == 0.0)
		{
			shellSign[shellIdx] = 0;
			nRet = 1; /* S_FALSE */
		}
		else 
			shellSign[shellIdx] = 
				(((shellIdx != maxShell) && shellVolume[shellIdx] < 0.0) ||
				((shellIdx == maxShell) && shellVolume[shellIdx] > 0.0)) ? 1 : -1; 	
	}


	for (faceIdx = 0; faceIdx < nbf; faceIdx++){
		shellIdx = faceColor[faceIdx];
		sign[faceIdx] = shellSign[shellIdx] * faceSign[faceIdx];
	}
	holes.erase(holes.begin() + 3 * (maxShell - 1) + 2);
	holes.erase(holes.begin() + 3 * (maxShell - 1) + 1);
	holes.erase(holes.begin() + 3 * (maxShell - 1) + 0);

	if (volume) {
		*volume = fabs(maxVolume);
		for (shellIdx = 1; shellIdx <= numOfColors; shellIdx++) {
			if (shellIdx != maxShell) {
				*volume -= fabs(shellVolume[shellIdx]);
			}
		}
		*volume /= 6.0;
	}

	for (faceIdx = 0; faceIdx < nbf; faceIdx++){
		if(sign[faceIdx] == 1)
		{
			int temp = ibnd[faceIdx*3];
			ibnd[faceIdx*3] = ibnd[faceIdx*3+1];
			ibnd[faceIdx*3+1] = temp;
		}
	}
	
END:
	if (faceColor) delete[] faceColor;
	if (faceSign) delete[] faceSign;

	return nRet;
}
//#pragma optimize( "", on)
};