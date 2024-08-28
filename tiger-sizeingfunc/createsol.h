#ifndef CREATESOL_H
#define CREATESOL_H

#include <iostream>
#include <Eigen/Eigen>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include "vector_own.h"
#include <string>

using namespace  std;
//#include <QString>
namespace AABBSER
{
void readfile(std::string &filename,Eigen::MatrixXd &V,Eigen::MatrixXi&F,Eigen::VectorXd&Psize);
double getTriSquare(Vector &pa, Vector &pb, Vector &pc);
void writeSol(std::string filename, Eigen::MatrixXd&V, Eigen::VectorXd&size);
void writeTriVTK(std::string filename, Eigen::MatrixXd &VI, Eigen::MatrixXi &FI, Eigen::VectorXd &markI);
void readEpsFromVTK(string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int> > &mark, Eigen::VectorXd &eps);
void rm_duplicate_facets(Eigen::MatrixXi &F);
void writeLocalParameterForMMGS(string&filename,Eigen::MatrixXd &V, Eigen::MatrixXi &F,Eigen::VectorXd &edgeSize,Eigen::VectorXd &hausd);
//void read_mesh_version2(QString &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F,Eigen::MatrixXi &T);
//void read_mesh_version2(QString &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F,Eigen::MatrixXi &T);
class bkMesh
{
public:
    bkMesh(Eigen::MatrixXd &V,Eigen::MatrixXi&F,Eigen::VectorXd&Psize,double beta=1.2);
    ~bkMesh();
    double getSize(double x,double y,double z);
    double getSizeSpatial(double x,double y,double z);
	double getSpanSize(double x, double y, double z);
public:
    GEO::Mesh supportMesh;
    Eigen::VectorXd PSize;
    GEO::index_t nearest_hint;
    GEO::MeshFacetsAABB *geo_sf_tree;
    double grad;
    GEO::Attribute<double> pointSize;
};


}






#endif // CREATESOL_H
