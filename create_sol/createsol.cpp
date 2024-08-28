#include "createsol.h"
#include <spdlog/spdlog.h> 
 #include "bkmeshimpl.h"
#include <Eigen/Eigen>
#include <exception>
#include <math.h>
using namespace std;
namespace CREATESOL
{
bkMesh* bkMesh::createFromBoundData(int NumPoints, double bndPts[], double bndPtSizes[], int NumTris, int bndFcts[], double beta)
{
    if(NumPoints<=0||NumTris<=0||beta<=0) return nullptr;
	//throw std::exception("123");
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd Psize;

    //reading points part
    V.resize(NumPoints,3);
    for(int i=0;i<NumPoints;i++)
    {
        V(i,0)=bndPts[3*i+0];
        V(i,1)=bndPts[3*i+1];
        V(i,2)=bndPts[3*i+2];
    }
    //end points

    //reading facets
    pair<string,bool>trianglePair;
    F.resize(NumTris,3);
    for(int i=0;i<NumTris;i++)
    {
        F(i,0)=bndFcts[3*i+0];
        F(i,1)=bndFcts[3*i+1];
        F(i,2)=bndFcts[3*i+2];
    }

    Psize.resize(NumPoints,1);
    for(int i=0;i<NumPoints;i++)
    {
        Psize(i,0)=bndPtSizes[i];
    }

    auto ret = new bkMesh;
    ret->impl = new bkMeshImpl(V,F,Psize,beta);
    return ret;
}

void bkMesh::destroy_Data()
{
	//GEO::terminate();

}

double bkMesh::getSize(double x, double y, double z)
{
    return impl->getSize(x,y,z);
}

double bkMesh::getSizeSpatial(double x, double y, double z)
{
    return impl->getSizeSpatial(x,y,z);
}

bkMesh::~bkMesh()
{
	if (impl) {
		delete impl;
		impl = nullptr;
	}
	destroy_Data();
}
}
