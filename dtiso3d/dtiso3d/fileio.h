#ifndef _FILEIO_H_
#define _FILEIO_H_

#include <fstream>
#include <string>
#include "iso3d_utility.h"
#include "solvezone.h"
#include "tools/EdgeGrapher.h"

inline int readPLS(const char *fname, 
			int *bndPtNum,
			double **bndPts, 
			int *bndFctNum,
			int **bndFcts)
{
	FILE* fp = NULL;
	int iTok = 0, iPatch;
	int i, j;
	char buf[512];
	int nEdges;
	//int nConnect;
	MYPOINT minW, maxW;
	INTEGER iMemory, nAllocNodes;
	int nNodes, nFaces;
	int errCode = 0;

	if (!fname || !bndPtNum || !bndFctNum || !bndPts || !bndFcts)
	{
		errCode = 2;
		goto FAIL;
	}
	*bndPtNum = *bndFctNum = 0;
	*bndPts = NULL;
	*bndFcts = NULL;

	fp = fopen(fname, "r");
	if (!fp)
	{
		printf("cannot read FILE %s\n", fname);
		errCode = 2;
		goto FAIL;
	}

	fgets(buf, 512, fp);
	sscanf(buf, "%d%d", &nFaces, &nNodes);

	printf("SURFACE NODES: %d\tSURFACE FACETS: %d\n", nNodes, nFaces);
	if (nFaces <= 0 || nNodes <= 0)
	{
		printf("Invalid PLS file %s\n", fname);
		errCode = 2;
		goto FAIL;
	}
	
	*bndPts = (double*) malloc(nNodes * sizeof(double)*3);
	*bndFcts = (int*) malloc(nFaces * sizeof(int)*3);

	if (NULL == *bndPts || NULL == *bndFcts)
	{
		printf("Not enough memory.\n", fname);
		errCode = -1;
		goto FAIL;
	}

	for (i = 0, j = 0; i < nNodes; i++, j++)
	{
		fscanf(fp, "%d%lf%lf%lf", &iTok, &(*bndPts)[3*j], &(*bndPts)[3*j+1], &(*bndPts)[3*j+2]);
	} 

	fgets(buf, 512, fp);

	for (i = 0, j = 0; i < nFaces; i++, j++)
	{
		fgets(buf, 512, fp);
#if 1
		sscanf(buf, "%d%d%d%d%d", &iTok, &(*bndFcts)[3*j], &(*bndFcts)[3*j+1], &(*bndFcts)[3*j+2], &iPatch);
			
		if ((*bndFcts)[3*j] > nNodes ||
			(*bndFcts)[3*j + 1] > nNodes ||
			(*bndFcts)[3*j + 2] > nNodes)
		{
			printf("Invalid file format, facet %d with too big node index\n", i + 1);
			errCode = -1;
			goto FAIL;
		}
		(*bndFcts)[3*j]--;
		(*bndFcts)[3*j + 1]--;
		(*bndFcts)[3*j + 2]--;
#endif
	}

	*bndPtNum = nNodes;
	*bndFctNum = nFaces;
	goto SUCC;

FAIL:
	if (*bndPts)
	{
		free(*bndPts);
		*bndPts = NULL;
	}
	if (*bndFcts)
	{
		free(*bndFcts);
		*bndFcts = NULL;
	}
	*bndPtNum = *bndFctNum = 0;
SUCC:
	
	if (fp)
		fclose(fp);

	return errCode;
}

inline int writePL3(const char *fname,
			 int mshPtNum, double *mshPts, 
			 int volElmNum, int *volElms)
{
	int elemSize = volElmNum, nodeSize = mshPtNum;
	FILE *fp = NULL;
	INTEGER i, j;
	MYPOINT pnt;

	if (NULL == (fp = fopen(fname, "w")))
	{
		printf("cannot write FILE %s\n", fname);
		return 2;
	}

	printf("\nWrite resulting file %s, please wait...\n", fname);

	fprintf(fp, "%d %d %d\n", elemSize, nodeSize, 0);

	for (i = 0, j = 0; i < nodeSize; i++, j++)
	{
		fprintf(fp, "%d %24.20g %24.20g %24.20g\n", i+1, 
			mshPts[3*j], mshPts[3*j+1], mshPts[3*j+2]);
	}

	for (i = 0, j = 0; i < elemSize; i++, j++)
	{
		fprintf(fp, "%d %d %d %d %d\n", i+1, 
			volElms[4*j] + 1, volElms[4*j + 1] + 1,
			volElms[4*j + 2] + 1, volElms[4*j + 3] + 1);
	}

	if (fp)
		fclose(fp);

	return 0;
}

inline int readVtk(const char *fname,
	int *bndPtNum,
	double **bndPts,
	int *bndFctNum,
	int **bndFcts,
	ZoneAttachInfo* attachInfo = nullptr) {

	using namespace std;
	std::ifstream fin(fname);

	if (!fin) return -1;

	string tempstr;

	//skip head
	for (int i = 0; i < 4; i++) {
		getline(fin, tempstr);
	}
	
	fin >> tempstr;
	if (tempstr != "POINTS") {
		return -2;
	}
	
	//get number of points
	if (!(fin >> *bndPtNum))
		return -2;
	*bndPts = new double[3 * (*bndPtNum)];
	fin >> tempstr;

	for (int i = 0; i < *bndPtNum; i++) {
		if (!(fin >> (*bndPts)[i * 3] >> (*bndPts)[i * 3 + 1] >> (*bndPts)[i * 3 + 2]))
			return -2;
	}

	fin >> tempstr;
	if (tempstr != "CELLS") {
		return -3;
	}
	//get number of faces
	fin >> *bndFctNum;
	*bndFcts = new int[3 * (*bndFctNum)];
	fin >> tempstr;

	for (int i = 0; i < *bndFctNum; i++){
		int nindices;
		fin >> nindices;
		if (nindices != 3) {
			return -4;
		}
		fin >> (*bndFcts)[i * 3] >> (*bndFcts)[i * 3 + 1] >> (*bndFcts)[i * 3 + 2];
		(*bndFcts)[i * 3];
		(*bndFcts)[i * 3 + 1];
		(*bndFcts)[i * 3 + 2];
	}

	fin >> tempstr;
	if (tempstr != "CELL_TYPES") {
		return -5;
	}
	int ntypes;
	fin >> ntypes;
	if (ntypes != *bndFctNum) {
		return -6;
	}

	for (int i = 0; i < *bndFctNum; i++) {
		int type_;
		fin >> type_;
		if (type_ != 5) {
			return -7;
		}
	}

	if (attachInfo) {
		fin >> tempstr;
		if (tempstr != "CELL_DATA" || !fin) {
			std::cout << "No attached information in VTK file!" << '\n';
			attachInfo->valid = false;
			return 0;
		}
		fin >> tempstr;
		while (true) {
			fin >> tempstr >> tempstr;
			if (!fin) break;
			if (tempstr == "tri_elm_to_org_face_tags") {
				std::string s1, s2, s3;
				fin >> s1 >> s2 >> s3;
				if (s1 != "int" || s2 != "LOOKUP_TABLE" || s3 != "default") return -8;
				int faceCount = 0;
				attachInfo->tri2FaceOffsets.push_back(0);
				for (int i = 0; i < *bndFctNum; i++) {
					int nface;
					int domain;
					fin >> domain;
					attachInfo->tri2Domain.push_back(domain);
					fin >> nface;
					for (int j = 0; j < nface; j++) {
						int iface;
						fin >> iface;
						attachInfo->tri2Faces.push_back(iface);
					}
					faceCount += nface;
					attachInfo->tri2FaceOffsets.push_back(faceCount);
				}
			}
			else if (tempstr == "body_to_face_tags") {
				int nblock;
				fin >> nblock;
				attachInfo->block2FaceOffsets.push_back(0);
				int faceCount = 0;
				for (int i = 0; i < nblock; i++) {
					std::string name, s_int;
					fin >> name;
					int fix1, nface;
					fin >> fix1 >> nface;
					fin >> s_int;
					if (s_int != "int" || fix1 != 1) {
						return -9;
					}
					for (int j = 0; j < nface; j++) {
						int iface;
						fin >> iface;
						attachInfo->block2Faces.push_back(iface);
					}
					faceCount += nface;
					attachInfo->block2FaceOffsets.push_back(faceCount);
				}
				attachInfo->blockSize.clear();
				attachInfo->blockSize.resize(nblock, -1);
			}
			else if (tempstr == "pec_body_tags") {
				int fix1_1,fix1_2;
				std::string s1, s2;
				int npec;
				fin >> fix1_1 >> s1 >> fix1_2 >> npec >> s2;
				if (fix1_1 != 1 || fix1_2 != 1 || s1 != "pec_bodies" || s2 != "int")
					return -10;
				for (int i = 0; i < npec; i++) {
					int ipec;
					fin >> ipec;
					attachInfo->pecBlocks.push_back(ipec);
				}
			}
			else if (tempstr == "body_sizes" || tempstr == "body_size") {
				int fix1, fix2;
				std::string s1, s2;
				int npair;
				fin >> fix1 >> s1 >> fix2 >> npair >> s2;
				if (fix1 != 1 || (s1 != "body_sizes" && s1 != "body_size") || fix2 != 2 || s2 != "double")
					return -11;
				if (npair != 0 && attachInfo->blockSize.empty()) return -12;
				for (int i = 0;  i < npair; i++) {
					int iblock;
					double size;
					fin >> iblock >> size;
					attachInfo->blockSize[iblock] = size;
				}
			}
			else {
				std::cout << "Unknown words " << tempstr << " in VTK file!" << '\n';
				return -14;
			}
		}
		attachInfo->valid = true;
	}
	return 0;
}

inline int writeVtk(const char *fname,
	int mshPtNum, double *mshPts,
	int volElmNum, int *volElms,
	int surfTriNum = 0, int * triElems = nullptr, bool tagTri = false, 
	int *triToFace = nullptr, int *triFaceNums = nullptr, int* triToDomain = nullptr,
	int *elemToSubdomain = nullptr, int subdomainNum=0, int *subdomainToBlock=nullptr, int *subdomainBlockNums=nullptr)
{
	ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0" << '\n';
	fout << fname << '\n';
	fout << "ASCII" << '\n';
	fout << "DATASET UNSTRUCTURED_GRID" << '\n';
	fout << "POINTS " << mshPtNum << " double" << '\n';

	fout.precision(16);

	for (int i = 0; i < mshPtNum; i++) {
		fout << mshPts[3 * i] << '\t' << mshPts[3 * i + 1] << '\t' << mshPts[3 * i + 2] << '\n';
	}

	fout << "CELLS " << volElmNum + surfTriNum << ' ' << volElmNum * 5 + surfTriNum * 4 << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 3 << '\t' << triElems[3 * i] << '\t' << triElems[3 * i + 1] << '\t' << triElems[3 * i + 2];
		if (tagTri) fout << '\t' << i;
		fout << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 4 << '\t' << volElms[4 * i] << '\t' << volElms[4 * i + 1] << '\t' << volElms[4 * i + 2] << '\t' << volElms[4 * i + 3];
		fout << '\n';
	}

	fout << "CELL_TYPES " << volElmNum + surfTriNum << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 5 << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 10 << '\n';
	}

	fout << "CELL_DATA " << volElmNum + surfTriNum << '\n';
	fout << "LOOKUP_TABLE default" << '\n';

	if (triToFace || elemToSubdomain)
	{
		fout << "SCALARS tri_elm_to_org_face_tags int" << '\n';
		fout << "LOOKUP_TABLE default" << '\n';
	}
	if (triToFace) {
		int offset = 0;
		for (int i = 0; i < surfTriNum; i++) {
			fout << triToDomain[i] << '\t' << triFaceNums[i];
			for (int j = 0; j < triFaceNums[i]; j++) {
				fout << '\t' << triToFace[offset + j];
			}
			fout << '\n';
			offset += triFaceNums[i];
		}
	}
	if (elemToSubdomain) {
		for (int i = 0; i < volElmNum; i++) {
			fout << elemToSubdomain[i] << '\n';
		}
	}

	if (subdomainNum) {
		fout << "FIELD subdomain_to_org_body_tags " << subdomainNum << '\n';
		int offset = 0;
		for (int i = 0; i < subdomainNum; i++) {
			fout << "subdomain_body" << i << " 1 " << subdomainBlockNums[i] << " int" << '\n';
			for (int j = 0; j < subdomainBlockNums[i]; j++) {
				if (j != 0)
					fout << " ";
				fout << subdomainToBlock[offset + j];
			}
			fout << '\n';
			offset += subdomainBlockNums[i];
		}
	}
	return 0;
}

inline int writeSurfVtk(const char *fname,
	int mshPtNum, double *mshPts,
	int surfElmNum, int *surfElms) {
	ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0" << '\n';
	fout << fname << '\n';
	fout << "ASCII" << '\n';
	fout << "DATASET UNSTRUCTURED_GRID" << '\n';
	fout << "POINTS " << mshPtNum << " double" << '\n';

	fout.precision(16);

	for (int i = 0; i < mshPtNum; i++) {
		fout << mshPts[3 * i] << '\t' << mshPts[3 * i + 1] << '\t' << mshPts[3 * i + 2] << '\n';
	}

	fout << "CELLS " << surfElmNum << ' ' << surfElmNum * 4 << '\n';

	for (int i = 0; i < surfElmNum; i++) {
		fout << 3 << '\t' << surfElms[3 * i] << '\t' << surfElms[3 * i + 1] << '\t' << surfElms[3 * i + 2];
		fout << '\n';
	}

	fout << "CELL_TYPES " << surfElmNum << '\n';

	for (int i = 0; i < surfElmNum; i++) {
		fout << 5 << '\n';
	}

	fout << "CELL_DATA " << surfElmNum << '\n';
	return 0;
}

//额外输出中点信息
inline int writeVtk_CenterPoint(const char *fname,
	int mshPtNum, double *mshPts,
	int volElmNum, int *volElms, 
	int surfTriNum = 0, int * triElems = nullptr, bool tagTri = false, int *triToFace = nullptr, int *triFaceNums = nullptr, int* triToDomain = nullptr,
	int *elemToSubdomain = nullptr, int subdomainNum = 0, int *subdomainToBlock = nullptr, int *subdomainBlockNums = nullptr) {

	ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0" << '\n';
	fout << fname << '\n';
	fout << "ASCII" << '\n';
	fout << "DATASET UNSTRUCTURED_GRID" << '\n';

	// 计算高阶中点信息
	std::vector<std::array<double, 3>> cps;
	const int iedge[6][2] = {
		{0,1},
		{1,2},
		{2,0},
		{0,3},
		{1,3},
		{2,3}};
	EdgeGrapher<int> eg;
	for (int i = 0; i < volElmNum; i++) {
		for (int e = 0; e < 6; e++) {
			int n1 = volElms[4 * i + iedge[e][0]];
			int n2 = volElms[4 * i + iedge[e][1]];
			if (eg.hasEdge(n1, n2)) continue;
			std::array<double, 3> cp;
			for (int d = 0; d < 3; d++) {
				cp[d] = (mshPts[3 * n1 + d] + mshPts[3 * n2 + d]) / 2;
			}
			eg.addEdgeOverride(n1, n2, cps.size());
			cps.push_back(cp);
		}
	}

	fout << "POINTS " << mshPtNum + cps.size() << " double" << '\n';

	fout.precision(16);

	for (int i = 0; i < mshPtNum; i++) {
		fout << mshPts[3 * i] << '\t' << mshPts[3 * i + 1] << '\t' << mshPts[3 * i + 2] << '\n';
	}

	for (const auto& cp : cps) {
		fout << cp[0] << '\t' << cp[1] << '\t' << cp[2] << '\n';
	}

	fout << "CELLS " << volElmNum + surfTriNum << ' ' << volElmNum * 11 + surfTriNum * 4 << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 3 << '\t' << triElems[3 * i] << '\t' << triElems[3 * i + 1] << '\t' << triElems[3 * i + 2];
		if (tagTri) fout << '\t' << i;
		fout << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 10 << '\t' << volElms[4 * i] << '\t' << volElms[4 * i + 1] << '\t' << volElms[4 * i + 2] << '\t' << volElms[4 * i + 3];
		for (int e = 0; e < 6; e++) {
			int n1 = volElms[4 * i + iedge[e][0]];
			int n2 = volElms[4 * i + iedge[e][1]];
			int icp;
			fout << '\t' << eg.getEdgeVal(n1, n2) + mshPtNum;
		}
		fout << '\n';
	}

	fout << "CELL_TYPES " << volElmNum + surfTriNum << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 5 << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 24 << '\n';
	}

	fout << "CELL_DATA " << volElmNum + surfTriNum << '\n';
	if (triToFace || elemToSubdomain)
	{
		fout << "SCALARS tri_elm_to_org_face_tags int" << '\n';
		fout << "LOOKUP_TABLE default" << '\n';
	}
	if(triToFace){
		int offset = 0;
		for (int i = 0; i < surfTriNum; i++) {
			fout << triToDomain[i] << '\t' << triFaceNums[i];
			for (int j = 0; j < triFaceNums[i]; j++) {
				fout << '\t' << triToFace[offset + j];
			}
			fout << '\n';
			offset += triFaceNums[i];
		}
	}
	if (elemToSubdomain) {
		for (int i = 0; i < volElmNum; i++) {
			fout << elemToSubdomain[i] << '\n';
		}
	}

	if (subdomainNum) {
		fout << "FIELD subdomain_to_org_body_tags " << subdomainNum << '\n';
		int offset = 0;
		for (int i = 0; i < subdomainNum; i++) {
			fout << "subdomain_body" << i << " 1 " << subdomainBlockNums[i] << " int" << '\n';
			for (int j = 0; j < subdomainBlockNums[i]; j++) {
				if (j != 0)
					fout << " ";
				fout << subdomainToBlock[offset + j];
			}
			fout << '\n';
			offset += subdomainBlockNums[i];
		}
	}
	return 0;
}
inline int writeVtk_CenterPoint_new(const char* fname,
	int mshPtNum, double* mshPts,
	int volElmNum, int* volElms,
	int surfTriNum = 0, int* triElems = nullptr, bool tagTri = false, int* triToFace = nullptr, int* triFaceNums = nullptr, int* triToDomain = nullptr,
	int* elemToSubdomain = nullptr, int subdomainNum = 0, int* subdomainToBlock = nullptr, int* subdomainBlockNums = nullptr) {

	ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0" << '\n';
	fout << fname << '\n';
	fout << "ASCII" << '\n';
	fout << "DATASET UNSTRUCTURED_GRID" << '\n';
	fout << "POINTS " << mshPtNum << " double" << '\n';

	fout.precision(16);

	for (int i = 0; i < mshPtNum; i++) {
		fout << mshPts[3 * i] << '\t' << mshPts[3 * i + 1] << '\t' << mshPts[3 * i + 2] << '\n';
	}
	fout << "CELLS " << volElmNum + surfTriNum << ' ' << volElmNum * 11 + surfTriNum * 4 << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 3 << '\t' << triElems[3 * i] << '\t' << triElems[3 * i + 1] << '\t' << triElems[3 * i + 2];
		if (tagTri) fout << '\t' << i;
		fout << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 10 << '\t' << volElms[10 * i] << '\t' << volElms[10 * i + 1] << '\t' << volElms[10 * i + 2] << '\t' << volElms[10 * i + 3]
			<< '\t' << volElms[10 * i + 4] << '\t' << volElms[10 * i + 5] << '\t' << volElms[10 * i + 6] << '\t' << volElms[10 * i + 7]
			<< '\t' << volElms[10 * i + 8] << '\t' << volElms[10 * i + 9] << '\n';
	}

	fout << "CELL_TYPES " << volElmNum + surfTriNum << '\n';

	for (int i = 0; i < surfTriNum; i++) {
		fout << 5 << '\n';
	}

	for (int i = 0; i < volElmNum; i++) {
		fout << 24 << '\n';
	}

	fout << "CELL_DATA " << volElmNum + surfTriNum << '\n';
	if (triToFace || elemToSubdomain)
	{
		fout << "SCALARS tri_elm_to_org_face_tags int" << '\n';
		fout << "LOOKUP_TABLE default" << '\n';
	}
	if (triToFace) {
		int offset = 0;
		for (int i = 0; i < surfTriNum; i++) {
			fout << triToDomain[i] << '\t' << triFaceNums[i];
			for (int j = 0; j < triFaceNums[i]; j++) {
				fout << '\t' << triToFace[offset + j];
			}
			fout << '\n';
			offset += triFaceNums[i];
		}
	}
	if (elemToSubdomain) {
		for (int i = 0; i < volElmNum; i++) {
			fout << elemToSubdomain[i] << '\n';
		}
	}

	if (subdomainNum) {
		fout << "FIELD subdomain_to_org_body_tags " << subdomainNum << '\n';
		int offset = 0;
		for (int i = 0; i < subdomainNum; i++) {
			fout << "subdomain_body" << i << " 1 " << subdomainBlockNums[i] << " int" << '\n';
			for (int j = 0; j < subdomainBlockNums[i]; j++) {
				if (j != 0)
					fout << " ";
				fout << subdomainToBlock[offset + j];
			}
			fout << '\n';
			offset += subdomainBlockNums[i];
		}
	}
	return 0;
}
#endif _FILEIO_H_