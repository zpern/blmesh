#include <spdlog/spdlog.h> 
 #include "BLMultipole_define.h"

//inline 
struct double3D CrossProduct(struct double3D op1, struct double3D op2)
{
	struct double3D temp;
	temp.x=op1.y*op2.z-op2.y*op1.z;
	temp.y=op1.z*op2.x-op2.z*op1.x;
	temp.z=op1.x*op2.y-op2.x*op1.y;
	return temp;
}
//inline 
void Normalize(struct double3D *op)
{
	double L=sqrt((*op)*(*op));
	*op=*op/L;
}

//inline 
double Distance(struct double3D pt1, struct double3D pt2)
{
	return sqrt((pt1-pt2)*(pt1-pt2));
}
double Distance(struct double3D *pt1, struct double3D *pt2)
{
	return sqrt((*pt1-*pt2)*(*pt1-*pt2));
}

//inline 
double ScalarProduct(struct double3D vt1, struct double3D vt2, struct double3D vt3)
{
	return vt1.x*vt2.y*vt3.z+vt1.y*vt2.z*vt3.x+vt1.z*vt2.x*vt3.y
		-vt1.z*vt2.y*vt3.x+vt1.y*vt2.x*vt3.z+vt1.x*vt2.z*vt3.y;
}
//inline struct 
double3D VectorProduct(struct double3D v1,struct double3D v2,struct double3D v3)
{
	struct double3D temp;
	double a=v1*v3;
	double b=v2*v3;
	temp.x=a*v2.x-b*v1.x;	temp.y=a*v2.y-b*v1.y;	temp.y=a*v2.y-b*v1.y;
	return temp;
}
//inline 
double DistSquare(struct double3D pt1, struct double3D pt2)
{
	return (pt1-pt2)*(pt1-pt2);
}
double DistSquare(struct double3D *pt1, struct double3D *pt2)
{
	return (*pt1-*pt2)*(*pt1-*pt2);
}
//////////////////////////////////////////////////////////////////////////////
void getAnOrthogonalVector(struct double3D v0,struct double3D *orthogonal)
{
	D3VECTOR v1, v2;
	double L=sqrt(v0*v0);
	if(fabs(v0.z) > L*0.2)
	{
		v1.x=v0.x; v1.y=-v0.z; v1.z=v0.y;
	}
	else
	{
		v1.x=-v0.y; v1.y=v0.x; v1.z=v0.z;
	}
	v2=CrossProduct(v0, v1);
	//ASSERT(v2*v2 > L*1.0e-15);
	*orthogonal=CrossProduct(v2, v0);
	Normalize(orthogonal);
}
//through a point outside a plane
D3POINT getVecticalDotPlane(D3POINT out, D3POINT pt0, D3POINT pt1, D3POINT pt2)
{
	D3VECTOR v1, v2, vo;
	double t1, t2,j, a11, a12, a22, c1, c2;
	v1=pt1-pt0;
	v2=pt2-pt0;
	vo=out-pt0; 
	a11=v1*v1;
	a22=v2*v2;
	a12=v1*v2;
	c1=v1*vo;
	c2=v2*vo;
	j=a11*a22-a12*a12;
	t1=(c1*a22-c2*a12)/j;
	t2=(c2*a11-c1*a12)/j;
	return pt0+v1*t1+v2*t2;
}
//through a point outside a line
D3POINT getVecticalDotLine(D3POINT out, D3POINT pt0, D3POINT pt1)
{
	D3VECTOR v1, vo;

	v1=pt1-pt0; 
	vo=out-pt0; 
	double t=(v1*vo)/(v1*v1);
	return pt0+v1*t;
}

double getAngleParameter(double theta1, double Dtheta, double theta)
{
	double t;
	t=(theta-theta1)/Dtheta;
	if(t>=0 && t<=1.0) return t;
	t=(theta-theta1+2.0*ONE_PI)/Dtheta;
	if(t>=0 && t<=1.0) return t;
	t=(theta-theta1-2.0*ONE_PI)/Dtheta;
	if(t>=0 && t<=1.0) return t;
	return -1.0;
}
