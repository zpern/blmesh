#ifndef __mesh_orient_h__
#define __mesh_orient_h__
#include <stdio.h>

namespace MeshOrient {

	enum { MAX_SHELL_SIZE = 128 };
	/*删除包围盒部分的网格 yhf20191226*/
	int removeBox(int nbp, double * lcoord, int &nbf, int * ibnd);
	/* --------------------------------------------------------
	* 基于以下基本假设对三角形网格的绕向进行排序
	* 1. 只有一个外环，多个内环
	* 2. 联通度（Genus)为1
	* 3. 所有面片法向指向区域外部
	* ------------------------------------------------------*/



	int orientTriangularSurface(
		int nbp,				/* number of boundary points */
		double *lcoord,			/* coord. of boundary points , size = 3 * nbp */
		int nbf,				/* number of boundary surface */
		int *ibnd,				/* topu. of boundary entities, size = 3 * nbf (i1, i2, i3) */
		int *sign,				/* 1: positive; -1: negative; 0: undetermined */
		double *volume,	/* volume */
		std::vector<double>& holes
	);
};

#endif /* __mesh_orient_h__ */