#ifndef KDTREE_H
#define KDTREE_H
#include "dataio.h"

struct TNode{
    TNode* left;
    TNode* right;
    double plane;
    int dim;
    BLVector min;
    BLVector max;
    std::vector<Triangle*> currentTriangle;
    bool isLeafNode(){
        return !currentTriangle.empty();
    }
	inline bool PointInTNode(BLVector p) {
		for (int k = 0; k < 3; k++) {
			if (p[k]<min[k] || p[k]>max[k])
				return false;
		}
		return true;
	}
    TNode():left(nullptr),right(nullptr){}
};
class kdTree
{
public:
    kdTree();
    static TNode *creatTree(BLVector min, BLVector max);
    bool buildTree(TNode* t,std::vector<Triangle*> triArr,int level);
    void print(TNode* t);
    std::vector<Triangle*> getTriangleMayInterWithRay(TNode* head,Ray r, double length);
private:
    inline bool isInLeft(Triangle *t, double plane, int d);
    inline bool isInRight(Triangle *t, double plane, int d);
};

#endif // KDTREE_H
