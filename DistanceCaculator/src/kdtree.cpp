#include <spdlog/spdlog.h> 
 #include "../include/kdtree.h"

kdTree::kdTree()
{

}

TNode *kdTree::creatTree(BLVector min,BLVector max)
{
    TNode* t=new TNode;
    t->min=min;
    t->max=max;
    return t;
}
inline bool kdTree::isInLeft(Triangle* t,double plane,int d){
    return !(t->box.min[d]>plane);
}
inline bool kdTree::isInRight(Triangle* t,double plane,int d){
    return !(t->box.max[d]<plane);
}
bool kdTree::buildTree(TNode* t,std::vector<Triangle*> triArr,int level)
{
	//cout << "level=" << level << "triSize" << triArr.size() << endl;;
    t->currentTriangle=triArr;
    if(triArr.size()<=10)
        return false;
    int c=0;//cost
    int cMin=INT_MAX;//minCost
    int dMin,iMin;
    std::vector<Triangle*> Tleft;
    std::vector<Triangle*> Tright;
	int dim = 10;
    for (int d=0;d<3;d++) {
      for (int i=0;i<dim-1;i++) {
		  int leftc=0;
		  int rightc=0;
          double m1=t->min[d];
          double m2=t->max[d];
          double m3=(i+1)*1.0/ dim *(m2-m1)+m1;
          for(int j=0;j<triArr.size();j++){
              if(isInLeft(triArr[j],m3,d)){
				  leftc++;
              }
              if(isInRight(triArr[j],m3,d)){
				  rightc++;
              }
          }
          if(abs(leftc - rightc)+ leftc + rightc<cMin){
              cMin= abs(leftc - rightc)+ leftc + rightc;
              dMin=d;
              iMin=i;
          }
      }
    }
    if(level>15||triArr.size()<200){
        return false;
    }
    Tleft.clear();
    Tright.clear();
    double m1=t->min[dMin];
    double m2=t->max[dMin];
    double m3=(iMin+1)*(1.0/ dim) *(m2-m1)+m1;
    for(int j=0;j<triArr.size();j++){
        if(isInLeft(triArr[j],m3,dMin)){
            Tleft.push_back(triArr[j]);
        }
        if(isInRight(triArr[j],m3,dMin)){
            Tright.push_back(triArr[j]);
        }
    }
    t->left=new TNode;
    t->left->min=t->min;
    t->left->max=t->max;
    t->left->max[dMin]=m3;
    t->right=new TNode;
    t->right->min=t->min;
    t->right->max=t->max;
    t->right->min[dMin]=m3;
    if(Tright.empty()||Tleft.empty()||!buildTree(t->left,Tleft,level+1)||!buildTree(t->right,Tright,level+1)){
        delete t->left;
        t->left=nullptr;
        delete t->right;
        t->right=nullptr;  
		//cout<< "level=" << level << "triSize" << triArr.size();
    }
    else{
        t->currentTriangle.clear();
    }
    return true;

}

void kdTree::print(TNode *t)
{
     
     if(t->left==nullptr||t->right==nullptr){

     }

    else{
        print(t->left);
        print(t->right);
    }
}

std::vector<Triangle*> kdTree::getTriangleMayInterWithRay(TNode* head,Ray r, double length)
{
    AABB box;
    box.max=head->max;
    box.min=head->min;
    if(rayAABBIntersect(r,box,length)){
        if(!head->currentTriangle.empty()){
            return head->currentTriangle;
        }
        else {
            std::vector<Triangle*> v1;
            if(head->left&&head->right){
               v1= getTriangleMayInterWithRay(head->left,r,length);


               std::vector<Triangle*> v2=getTriangleMayInterWithRay(head->right,r,length);
               for(int i=0;i<v2.size();i++){
                      v1.push_back(v2[i]);
               }
            }
           return v1;
        }

    }
    else {
        return std::vector<Triangle*>();
    }
}
