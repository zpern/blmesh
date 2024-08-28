#ifndef __iso3d_multidomain_h__
#define __iso3d_multidomain_h__

//add by zyj
int splitMultiDomain(
	int bndPtNum, double bndPts[],
	int bndFctNum, int bndFcts[],
	int *domainNum,
	int *domainBndPtNum[], double* domainBndPts[],			//DIM of *domainBndPts = sum(*domainBndPtNum[i]) * 3
	int *domainBndFctNum[], int* domainBndFcts[],			//DIM of *domainBndFcts = sum(*domainBndFctNum[i]) * 3
	int *ltowBndPtMap[]									//DIM of *ltowBndPtMap = domainNum,  DIM of *ltowBndPtMap[domainIndex] = domainBndPtNum[domainIndex]
);

#endif
