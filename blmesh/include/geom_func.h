/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * �����������ɵļ��κ����� (�汾�ţ�0.1)
 * liwgeom Function Library for Mesh Generation (Version 0.1)
 *
 * �½��� �й� �㽭��ѧ�������ѧ�����о�����
 * ��Ȩ����	  2012��11��28��
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 28/11/2012
 * 
 * ��ϵ��ʽ
 *   �绰��+86-571-87951883
 *   ���棺+86-571-87953167
 *   ���䣺chenjj@zju.edu.cn
 * For further information, please conctact
 *  Tel: +86-571-87951883
 *  Fax: +86-571-87953167
 * Mail: chenjj@zju.edu.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/
#ifndef __geom_func_h__
#define __geom_func_h__

#include <stdio.h>
//#include<math.h>
//#include<stdlib.h>

namespace GEOM_FUNC
{
// extern double orient2d(double *pa, double *pb, double *pc);
// extern double orient3d(double *pa, double *pb, double *pc, double *pd);
// extern void exactinit();

int cmpTwoTriDir(double *,double *,double *,double *,double *,double *);
/* �����ཻ�ж� */
int lnFacInt(double ln1[], double ln2[], 
			 double fac1[], double fac2[], double fac3[],  
			 double pnt[], int *val,
			 int edgFac1 = -1, int edgFac2 = -1, int edgFac3 = -1, 
			 double *pt1 = NULL, double *pt2 = NULL, double *pt3 = NULL);


#ifdef __cplusplus
extern "C" {
#endif
/* function prototype */

/* extern the predicates */
extern double macheps;
double exactinit();
double fixedSplitPoint(double s1, double s2, double pnt1, double pnt2);
double orient2d(double *pa, double *pb, double *pc);
double incircle(double*, double*, double*, double*);
double orient3d(double *pa, double *pb, double *pc, double *pd);
double orient3dexact(double *pa, double *pb, double *pc, double *pd);
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);

void norm_3p(double *p1 , double *p2 ,double *p3, double *normal);
extern int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], 
			    double p2[3], double q2[3], double r2[3]);
extern int one_node_same_tri_tri_overlap_3d(double p1[3], double q1[3], double r1[3], 								 
			    double p2[3], double q2[3], double r2[3]);  //zhaodawei add 2010-07-16
int lin_tri_intersect3d_check(double linep[2][3], double facep[3][3]);
/* -----------------------------------------------------------------------------------
 * ����һ���ж�1���߶κ��������Ƿ��ཻ�Ĵ��룬�����ؽ���λ�� 
 * -----------------------------------------------------------------------------------*/
enum LIN_TRI_INT_TYPE 
{
	LTI_INTERSECT_NUL = 0,
	LTI_INTERSECT_NOD,
	LTI_INTERSECT_EDG,
	LTI_INTERSECT_FAC,
	LTI_INTERSECT_INS,			/* ֱ����ȫ��ƽ���� */
	LTI_INTERSECT_DEG_FACE		/* ������ʱ���õ�ö�ٱ����������ܻᱻ�滻��
								   ��ʾ�������ƬΪ�˻���Ƭ��Ŀǰ�ݲ��ж�һ��
								   �˻�����Ƭ�Ƿ�����ཻ��û���������*/
};

/* -----------------------------------------------------------------------------------
 * ��������������p1/p2/p3/p4���棬��p1/p2/p3/p4��ʹ��p1p2p3ͶӰ�����������ƽ��ͶӰ
 * ���p1p2p3���˻���Ƭ������0�����򣬷���1
 * -----------------------------------------------------------------------------------*/
int proj_four_coplanr_points(double p1[3], double p2[3], double p3[3], double p4[3],
		double proj1[2], double proj2[2], double proj3[2], double proj4[2]);

/* -----------------------------------------------------------------------------------
 * ����һ���ж�1���߶κ��������Ƿ��ཻ�Ĵ���(2D)�������ؽ���λ��
 * �ཻ���ͣ�
 * PNT �ཻ��1���㣨intCod���ص�ı��0~2��i������ĵ�i�����㣩��
 * EDG �ཻ��1���ߣ�intCod���رߵı��0~2��i����(i,(i+1)%3�γɵı�)
 * FAC �ཻ��1����
 * intPnt���ؽ����ֵ
 * ע�⣺�ڹ��������£�һ���߿��ܻ��һ����������߶��ཻ����ʱ��
 * ���ݵ��ø��㷨�ı߽�ָ�������������ֻ�������������Ǹ�����
 * ��������������ڱ�������ཻ������ʱ��������Ҫ���ݾ���������и���
 * -----------------------------------------------------------------------------------*/
int lin_tri_intersect2d(double linep[2][2], double facep[3][2], int *intTyp, int *intCod, double intPnt[2]);

/* -----------------------------------------------------------------------------------
 * ����һ���ж�1��λ�����ϵ�3ά�����Ĺ�ϵ
 * �ཻ���ͣ�
 * PNT �ཻ��1���㣨intCod���ص�ı��0~2��i������ĵ�i�����㣩��
 * EDG �ཻ��1���ߣ�intCod���رߵı��0~2��i����(i,(i+1)%3�γɵı�)
 * FAC �ཻ��1����
 * intPnt���ؽ����ֵ
 * -----------------------------------------------------------------------------------*/
extern int pnt_tri_intersect3d(double poinp[3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3]);

extern int lin_tri_intersect3d(double linep[2][3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3], bool bEpsilon = true);

extern int lin_tri_intersect3d_idx(int iline[2], int iface[3], double linep[2][3], double facep[3][3], int *intTyp, int *intCod, double intPnt[3], bool bEpsilon = true);

extern int isintersect_oneSharePoint(int iface[], int iline[], double facep[][3], double linep[][3], int *intTyp, int *intCod, double intPnt[3]);

/* ------------------------------------------------------------------------------------
 * �����ཻ�ж�
 * ����һ�����ɺ��������ǳ�����lin_tri_intersect3dȥʵ���������������
 * ---------------------------------------------------------------------------------- */
int lnFacInt2(double ln1[], double ln2[], 
			 double fac1[], double fac2[], double fac3[],  
			 double pnt[], int *val,
			 int edgFac1 = -1, int edgFac2 = -1, int edgFac3 = -1, 
			 double *pt1 = NULL, double *pt2 = NULL, double *pt3 = NULL, bool bEpsilon = true);

/* -----------------------------------------------------------------------------------
 * ����һ���ж�2��������Ƭ�Ƿ��ཻ�Ĵ��룬������small polyhedron reconnection�㷨�е���
 * -----------------------------------------------------------------------------------*/
extern int tri_tri_intersect3d(int facei1[3], int facei2[3], double facep1[3][3], double facep2[3][3]);
extern int tri_tri_intersect3d_fast(int facei1[3], int facei2[3], double *facep1[3], double *facep2[3]);

/* ����һ�������嵥Ԫ����״����ֵ */
extern double tetrahedron_gamma(double p1[3], double p2[3], double p3[3], double p4[3]);

#ifdef __cplusplus
}
#endif

}
#endif