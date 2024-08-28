#ifndef _WRITE_BLMESH_H_
#define _WRITE_BLMESH_H_
#include <fstream>
#include <string>
/*
enum EntityTopology
{
	NODE = 5,
	LINE = 6,
	POLYGON = 7,
	TRIANGLE = 8,
	QUADRILATERAL = 9,
	POLYHEDRON = 10,
	TETRAHEDRON = 11,
	HEXAHEDRON = 12,
	PRISM = 13,
	PYRAMID = 14,
	SEPTAHEDRON = 15,
	MIXED = 16,
};
*/
static void writeBLVTK(
	std::string filename,
	double	   *ppdMNC,		/* 体网格节点坐标 */
	int         pnMN,			/* 体网格节点数目 */
	int         *ppnMEFm,		/* 体网格单元节点编号 */
	int         *ppnMETp,		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 */
	int         pnME			/* 体网格单元数目 */
) {
	int i, j, sidx;
	FILE *fout = nullptr;

	int	npt = pnMN;
	int nelm = pnME;


	fout = fopen(filename.c_str(), "w");


	//vtk
	fprintf(fout, "# vtk DataFile Version 2.0\n");
	fprintf(fout, "boundary layer mesh\n");
	fprintf(fout, "ASCII\n");
	fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fout, "POINTS %d double\n", npt);

	for (i = 0; i < npt; i++)
	{
		fprintf(fout, "%.15lf %.15lf %.15lf\n", ppdMNC[3 * i + 0], ppdMNC[3 * i + 1], ppdMNC[3 * i + 2]);
	}

	int nPrism = 0; int nTet = 0; int nPyramid = 0;
	for (i = 0; i < pnME; i++) {
		if (ppnMETp[i] == 13)
			nPrism++;
		if (ppnMETp[i] == 14)
			nPyramid++;
		if (ppnMETp[i] == 11)
			nTet++;
	}

	fprintf(fout, "CELLS %d %d\n", nelm, nPrism * 7 + nTet * 5 + nPyramid * 6);



	int count = 0;

	for (i = 0; i < pnME; i++)
	{
		if (ppnMETp[i] == 13) {
			fprintf(fout, "6 %d %d %d %d %d %d\n", ppnMEFm[count],
				ppnMEFm[count+1], ppnMEFm[count+2],
				ppnMEFm[count+3], ppnMEFm[count+4], ppnMEFm[count+5]);
			count += 6;
		}
		if (ppnMETp[i] == 14) {
			fprintf(fout, "5 %d %d %d %d %d\n", ppnMEFm[count],
				ppnMEFm[count + 1], ppnMEFm[count + 2],
				ppnMEFm[count + 3], ppnMEFm[count + 4]);
			count += 5;
		}
		if (ppnMETp[i] == 11) {
			fprintf(fout, "4 %d %d %d %d\n", ppnMEFm[count],
				ppnMEFm[count + 1], ppnMEFm[count + 2],
				ppnMEFm[count + 3]);
			count += 4;
		}
	}

	fprintf(fout, "CELL_TYPES %d\n", nelm);


	for (i = 0; i < pnME; i++)
	{
		if (ppnMETp[i] == 13) {
			fprintf(fout, "%d\n", 13);
		}
		if (ppnMETp[i] == 14) {
			fprintf(fout, "%d\n", 14);
		}
		if (ppnMETp[i] == 11) {
			fprintf(fout, "%d\n", 10);
		}
	}
	return ;
}
static void writeBC(std::string filename, int *boundary_face) {

}

static void writeBoundary(
	std::string filename,
	int         num_boundary_face, /*边界面网格数量 */
	int         *boundary_mesh, /*边界面网格，注意每四个为一组，而且注意如果为三角形，最后一项为-1，id从0开始 */
	int         *boundary_face /*边界面网格对应的面id，长度为num_boundary_face */
	
	) {

	int i, j, nelm, cnt;
	FILE *fout = nullptr;


	fout = fopen(filename.c_str(), "w");

	cnt = 0;
	//only output volume mesh
	int nTri = 0;
	int nQuad = 0;
	for (int i = 0; i < num_boundary_face; i++) {
		if (boundary_mesh[4 * i + 3] == -1)
			nTri++;
		else
			nQuad++;
	}

	nelm = nTri + nQuad;

	fprintf(fout, "%d %d %d\n", nelm, nTri, nQuad);

	for (i = 0; i < num_boundary_face; i++)
	{
		if (boundary_mesh[4 * i + 3] == -1)
			fprintf(fout, "%d 3 %d %d %d %d\n", cnt + 1, boundary_mesh[4 * i + 0] + 1,
				boundary_mesh[4 * i + 1] + 1, boundary_mesh[4 * i + 2] + 1, boundary_face[i]);

		else
		{
			fprintf(fout, "%d 4 %d  %d %d %d %d\n", cnt + 1, boundary_mesh[4 * i + 0] + 1,
				boundary_mesh[4 * i + 1] + 1, boundary_mesh[4 * i + 2] + 1, boundary_mesh[4 * i + 3] + 1, boundary_face[i]);

		}
	}

	fclose(fout);
	fout = nullptr;

	return ;
}
#endif //! _WRITE_BLMESH_H_
