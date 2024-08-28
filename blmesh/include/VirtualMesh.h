#ifndef _VM_H_
#define _VM_H_
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <vector>
#include <array>
#include <set>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using std::cout;


struct VM {
	VM() :ppdMNC(nullptr), pnMN(0), ppnMEFm(nullptr), ppnMETp(nullptr), pnME(0), ppdSNC0(nullptr), nSN0(0), pnSNO(0)
		, ppnSFTpO(nullptr), ppnSFFmO(nullptr) {}

	double	   *ppdMNC;		/* 体网格节点坐标 **/
	int        pnMN;			/* 体网格节点数目 **/
	int        *ppnMEFm;		/* 体网格单元节点编号 **/
	int        *ppnMETp;		/* 体网格单元类型。当前支持单元类型：三棱柱、金字塔、四面体 **/
	int        pnME;			/* 体网格单元数目 **/


	double		*ppdSNC0;			/* 顶面网格节点坐标，coord[3*i], coord[3*i+1], coord[3*i+2]为第i个节点x,y,z坐标 **/
	int			nSN0;		/* 顶面网格边界节点数目 **/
	int			pnSNO;		/* 顶面网格单元数目 **/
	int			*ppnSFTpO;		/* 顶面网格单元类型。当前仅支持三角形单元 **/
	int			*ppnSFFmO;		/* 顶面网格单元节点编号 **/
    double      *sizing; /* 顶面网格目标尺寸 */

	int      num_boundary_face;     /*边界面单元数目**/
	int *    boundary_mesh;         /*边界面网格，4个为一个单位，如果是三角形则第四个为-1**/
	int *    boundary_face;         /*边界面网格所在的面id**/
	
	int* l2g; /* map local point index into global **/

	bool checkOri();
	bool checkManifold();
	void convert2Manifold();
	void convert2NonManifold();
	void rmoveFreeNodeInTop();

	void saveVTK(std::string filename);

};
#endif
