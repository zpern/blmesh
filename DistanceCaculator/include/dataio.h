#pragma once
#ifndef _DATAIO_H
#define _DATAIO_H
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <climits>
#include <cmath>
using namespace std;
#include "BLVector.h"
enum hyperPlanType{
    X=0,
    Y,
    Z
};
enum LSType{
    sphere=0,
    quad
};
static inline double max(double num1,double num2){
    return num1>num2?num1:num2;
}
static inline double min(double num1,double num2){
    return num1>num2?num2:num1;
}
map<int, double> GetDistanceRatio(map<int, pair<BLVector, BLVector> > &point, vector<vector<int> > &graph);
class LightSou
{
public:
	BLVector source;
    double size;
    LSType type;
	BLVector le;
    BLVector nor;//for quad only
};
struct AABB{//axis aligned bounding box
    BLVector max;
    BLVector min;
};
struct Triangle{
  BLVector pos[3];
  AABB box;
  unsigned int index;
};
struct Ray{
   BLVector startPos;//起始点
   BLVector oritation;//方向
   double Ix,Iy,Iz;
};

struct Camera{
    BLVector d1;
    BLVector d2;
    BLVector pos;
    BLVector head;
    BLVector min;//这两个变量定义了屏幕所在的位置
    BLVector max;//相当于定义了center（lookat）和fovy
};
bool rayAABBIntersect(Ray r,AABB b, double length);
#endif
