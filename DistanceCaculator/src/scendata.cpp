#include "scendata.h"



void  scenData::rayTriIntersect(Ray r,Triangle t,bool& b,BLVector& interPos,bool& inBou,float &tt,float &u,float &v){
//    if(!rayAABBIntersect(r,t.box)){
    //    b=false;
    //    return ;
    //}
    BLVector e1=t.pos[1]-t.pos[0];
    BLVector e2=t.pos[2]-t.pos[0];
    BLVector q=BLVector::crossProduct(r.oritation,e2);
    float a=BLVector::dotProduct(q,e1);
    if(abs(a)<float(1e-9)){
        b=false;
        return;
    }
    float f=float(1.0)/a;
    BLVector s=r.startPos-t.pos[0];
    u=f*BLVector::dotProduct(q,s);
    if(u<float(0.0)){
        b=false;
        return;
    }
    BLVector rp=BLVector::crossProduct(s,e1);
    v=f*BLVector::dotProduct(rp,r.oritation);
    if(v<float(0.0f)||u+v>float(1.0f)){
        b=false;
        return;
    }
    tt=f*BLVector::dotProduct(rp,e2);
    if(tt<0){
        b=false;
        return;
    }
    b=true;
    interPos=r.startPos+tt*r.oritation;
    inBou=false;
    /*for test*/
    if(u<1e-2||v<1e-2||u>1-1e-2||v>1-1e-2||u+v>1-1e-2)
        inBou=true;
    return;
}

scenData::scenData()
{
   kdt=new kdTree;
}


bool  scenData::planeAABBIntersect(hyperPlanType tp, AABB b, BLVector min, BLVector max)
{
    if(min[tp]>b.max[tp]||max[tp]<b.min[tp])
        return false;
    int id1,id2;
    bool cont1=false,cont2=false;
    if(tp==0)
        id1=1;id2=2;
    if(tp==1)
        id1=0;id2=2;
    if(tp==2)
        id1=0;id2=1;
    if(min[id1]>b.min[id1]&&min[id1]<b.max[id1]){
        cont1=true;
    }
    if(min[id2]>b.min[id2]&&min[id2]<b.max[id2]){
        cont2=true;
    }
    if(max[id1]>b.min[id1]&&max[id1]<b.max[id1]){
        cont1=true;
    }
    if(max[id2]>b.min[id2]&&max[id2]<b.max[id2]){
        cont2=true;
    }
    if(min[id1]<b.min[id1]&&max[id1]>b.max[id1]){
        cont1=true;
    }
    if(min[id2]<b.min[id2]&&max[id2]>b.max[id2]){
        cont2=true;
    }
    return cont1&&cont2;

}
void  scenData::raySphereFirstIntersect(Ray y,bool &b,LightSou ls,int &t){

}
void  scenData::rayTriFirstIntersect(Ray r, bool &b, BLVector &interPos, Triangle* &t, bool &inBou, float &u, float &v)
{
    bool success;
    bool inBouIn=false;
    b=false;
    float tt,um,vm;
    BLVector pos;
    Triangle* bestTri;
    BLVector bestPos;
    float min=float(1e10);
    std::vector<Triangle*> vt;
   // vt=kdt->getTriangleMayInterWithRay(head,r);
    std::map<int,int> iToi;
    //qDebug()<<vt.size();
    for(unsigned int i=0;i< vt.size();i++){
        if(iToi.find(vt[i]->index)!=iToi.end())
            continue;
        else {
            iToi[vt[i]->index]=1;
        }
        rayTriIntersect(r, *vt[i],success,pos,inBouIn,tt,um,vm);
        inBou=inBou||inBouIn;
        if(success){
            if((pos-r.startPos).length()<min){
                min=(pos-r.startPos).length();
                bestTri=vt[i];
                bestPos=pos;
                u=um;
                v=vm;
                b=true;
            }
        }
    }
    t=bestTri;
    interPos=bestPos;
    return;
}
