#ifndef SCENDATA_H
#define SCENDATA_H
#include "dataio.h"
#include "kdtree.h"
class scenData {
public:
    scenData();
    kdTree *kdt;
    TNode* head;

    void rayTriIntersect(Ray r, Triangle t, bool& b, BLVector& interPos, bool &inBou, float &tt, float &u, float &v);//使用重心坐标法
    void rayTriFirstIntersect(Ray r,bool& b,BLVector& interPos,Triangle* &t,bool & inBou,float &u,float &v);//最耗时函数
    bool planeAABBIntersect(hyperPlanType tp,AABB b,BLVector min,BLVector max);
	void raySphereFirstIntersect(Ray y, bool &b, LightSou ls, int &t);
};
#endif // SCENDATA_H
