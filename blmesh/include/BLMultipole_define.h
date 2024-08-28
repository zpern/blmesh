
#ifndef _BLMultipole_Define_
#define _BLMultipole_Define_
#include "algebra.h"

#define TOLERANCE     1/10000000000.0
#define PROXIMITY     1/1000000.0
#ifndef PI
	#define PI            3.14159265358979323846264338327950
#endif
#define ONE_PI        3.14159265358979323846264338327950
#define ONEby4PI      0.079577471546

typedef struct double3D
{
	double x, y, z;

	void init(double xx, double yy, double zz)
	{
		x = xx;
		y = yy;
		z = zz;
	}

	struct double3D operator=(struct double3D op)
	{
		x=op.x;	y=op.y;	z=op.z;
		return *this;
	}
	friend struct double3D operator+(struct double3D op1, struct double3D op2)
	{
		struct double3D temp;
		temp.x=op1.x+op2.x;  temp.y=op1.y+op2.y;  temp.z=op1.z+op2.z;
		return temp;
	}
	friend struct double3D operator-(struct double3D op1, struct double3D op2)
	{
		struct double3D temp;
		temp.x=op1.x-op2.x;  temp.y=op1.y-op2.y;  temp.z=op1.z-op2.z;
		return temp;
	}
	friend struct double3D operator*(struct double3D op1, double op2)
	{
		struct double3D temp;
		temp.x=op1.x*op2;  temp.y=op1.y*op2;  temp.z=op1.z*op2;
		return temp;
	}
	friend struct double3D operator/(struct double3D op1, double op2)
	{
		struct double3D temp;
		temp.x=op1.x/op2;  temp.y=op1.y/op2;  temp.z=op1.z/op2;
		return temp;
	}
	friend double operator *(struct double3D op1, struct double3D op2)
	{
		return op1.x*op2.x+op1.y*op2.y+op1.z*op2.z;
	}
} D3POINT, D3NORMAL, D3VECTOR;

typedef struct CellSpacePoint
{
	double s, t;//xi(sai), eta

	CellSpacePoint()
	{
		s=t=0.0;
	}
	CellSpacePoint(double s0, double t0)
	{
		s=s0; t=t0;
	}

	CellSpacePoint operator=(CellSpacePoint op)
	{
		s=op.s;	t=op.t;
		return *this;
	}
	friend struct CellSpacePoint operator+(struct CellSpacePoint op1, struct CellSpacePoint op2)
	{
		struct CellSpacePoint temp;
		temp.s=op1.s+op2.s;  temp.t=op1.t+op2.t;
		return temp;
	}
	friend struct CellSpacePoint operator-(struct CellSpacePoint op1, struct CellSpacePoint op2)
	{
		struct CellSpacePoint temp;
		temp.s=op1.s-op2.s;  temp.t=op1.t-op2.t;
		return temp;
	}
	friend struct CellSpacePoint operator*(struct CellSpacePoint op1, double op2)
	{
		struct CellSpacePoint temp;
		temp.s=op1.s*op2;  temp.t=op1.t*op2;
		return temp;
	}
	friend struct CellSpacePoint operator/(struct CellSpacePoint op1, double op2)
	{
		struct CellSpacePoint temp;
		temp.s=op1.s/op2;  temp.t=op1.t/op2;
		return temp;
	}
	friend double operator *(struct CellSpacePoint op1, struct CellSpacePoint op2)
	{
		return op1.s*op2.s+op1.t*op2.t;
	}
} CEPOINT;

#define TRIANGLE_GAUSS_NUM_4
typedef struct TrianglePatch
{
	double x1, y1, x2, y2, x3, y3;

	TrianglePatch()
	{
		x1=y1=x2=y2=x3=y3=0.0;
	}
	~TrianglePatch(){};

	void getElemtCoordAndJb_Full(CEPOINT *vPt, double *vJb)
	{
		///double S=(x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1;
		////(pt1.s-pt3.s)*(pt2.t-pt3.t)-(pt2.s-pt3.s)*(pt1.t-pt3.t);
		double S=(x1-x3)*(y2-y3)-(x2-x3)*(y1-y3);
		if(S < 0.0) S=0.0-S;
		S=S*0.5;

#ifdef TRIANGLE_GAUSS_NUM_7
		for(int i=0; i<7; i++)
		{
			vPt[i].s=x1*gsTrgl1_g7[i]+x2*gsTrgl2_g7[i]+x3*gsTrgl3_g7[i];
			vPt[i].t=y1*gsTrgl1_g7[i]+y2*gsTrgl2_g7[i]+y3*gsTrgl3_g7[i];
			vJb[i]=gsTrigl_w7[i]*S;
		}
#endif
#ifdef TRIANGLE_GAUSS_NUM_4
		for(int i=0; i<4; i++)
		{
			vPt[i].s=x1*gsTrgl1_g4[i]+x2*gsTrgl2_g4[i]+x3*gsTrgl3_g4[i];
			vPt[i].t=y1*gsTrgl1_g4[i]+y2*gsTrgl2_g4[i]+y3*gsTrgl3_g4[i];
			vJb[i]=gsTrigl_w4[i]*S;
		}
#endif
#ifdef TRIANGLE_GAUSS_NUM_3
		for(int i=0; i<3; i++)
		{
			vPt[i].s=x1*gsTrgl1_g3[i]+x2*gsTrgl2_g3[i]+x3*gsTrgl3_g3[i];
			vPt[i].t=y1*gsTrgl1_g3[i]+y2*gsTrgl2_g3[i]+y3*gsTrgl3_g3[i];
			vJb[i]=gsTrigl_w3[i]*S;
		}
#endif
#ifdef TRIANGLE_GAUSS_NUM_1
		vPt[0].s=x1/3.0+x2/3.0+x3/3.0;
		vPt[0].t=y1/3.0+y2/3.0+y3/3.0;
		vJb[0]=S;
#endif
	}
	///////////////////////////////////////////////////////////////////
	void getElemtCoordAndJb_Rank(CEPOINT *vPt, double *vJb)
	{
		double S=(x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1;
		if(S < 0.0) S=0.0-S;

#ifdef RANK_TRIANGLE_GAUSS_NUM_7
		for(int i=0; i<7; i++)
		{
			vPt[i].s=x1*gsTrgl1_g7[i]+x2*gsTrgl2_g7[i]+x3*gsTrgl3_g7[i];
			vPt[i].t=y1*gsTrgl1_g7[i]+y2*gsTrgl2_g7[i]+y3*gsTrgl3_g7[i];
			vJb[i]=gsTrigl_w7[i]*S;
		}
#endif
#ifdef RANK_TRIANGLE_GAUSS_NUM_4
		for(int i=0; i<4; i++)
		{
			vPt[i].s=x1*gsTrgl1_g4[i]+x2*gsTrgl2_g4[i]+x3*gsTrgl3_g4[i];
			vPt[i].t=y1*gsTrgl1_g4[i]+y2*gsTrgl2_g4[i]+y3*gsTrgl3_g4[i];
			vJb[i]=gsTrigl_w4[i]*S;
		}
#endif
#ifdef RANK_TRIANGLE_GAUSS_NUM_3
		for(int i=0; i<3; i++)
		{
			vPt[i].s=x1*gsTrgl1_g3[i]+x2*gsTrgl2_g3[i]+x3*gsTrgl3_g3[i];
			vPt[i].t=y1*gsTrgl1_g3[i]+y2*gsTrgl2_g3[i]+y3*gsTrgl3_g3[i];
			vJb[i]=gsTrigl_w3[i]*S;
		}
#endif
#ifdef RANK_TRIANGLE_GAUSS_NUM_1
		vPt[0].s=x1/3.0+x2/3.0+x3/3.0;
		vPt[0].t=y1/3.0+y2/3.0+y3/3.0;
		vJb[0]=S;
#endif
	}

}CTrglPatch;



///////////////////////////////////////////////////////////////
struct double3D CrossProduct(struct double3D op1, struct double3D op2);
void Normalize(struct double3D *op);
double Distance(struct double3D pt1, struct double3D pt2);
double Distance(struct double3D *pt1, struct double3D *pt2);
double DistSquare(struct double3D, struct double3D);
double DistSquare(struct double3D*, struct double3D*);
double ScalarProduct(struct double3D vt1, struct double3D vt2, struct double3D vt3);
struct double3D VectorProduct(struct double3D v1,struct double3D v2,struct double3D v3);
void getAnOrthogonalVector(struct double3D v0,struct double3D *orthogonal);
//through a point outside a plane
D3POINT getVecticalDotPlane(D3POINT out, D3POINT pt0, D3POINT pt1, D3POINT pt2);
//through a point outside a line
D3POINT getVecticalDotLine(D3POINT out, D3POINT pt0, D3POINT pt1);
double getAngleParameter(double theta1, double Dtheta, double theta);

//立方体区域结构
struct TREECUBE
{
	D3POINT center;
	double halfD;
};

//节点列表结构
struct NODELIST
{
	long nodeNo_;
	NODELIST* next;
};

//八叉树节点结构
struct CUBENODELIST;

struct FMMTREENODE
{
	FMMTREENODE* father;
	FMMTREENODE* child[8];

	TREECUBE cube;
	NODELIST* listNodes;
	long tlNode; //number of nodes included within this treenode
	int level;
	CUBENODELIST* listNeighbors;
	CUBENODELIST* listInteract;

	complex<double>*** Mmn;//远端系数
	complex<double>*** Lmn;//近端系数

	int tlLocalNode;//for leaf treenode only
	int* sequence;
	//	double** localU;
	//	double** localQ;
	double** M;
	//	complex*** INT_RmnU;
	//	complex*** INT_RmnQ;
};

//邻居树节点列表结构
struct CUBENODELIST
{
	FMMTREENODE* cubeNode;
	CUBENODELIST* next;
};
#endif
