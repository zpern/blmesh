#ifndef BKMESHIMPL_H
#define BKMESHIMPL_H
#include <unordered_map>
#include "vector_own.h"
#define SLIB
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <geogram/basic/geometry_nd.h>
#include <memory>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include "vector_own.h"
#include <iostream>
#include <string>
#ifndef SLIB
#include <QString>
#include <QFile>
#include <QTextStream>
#include <functional>
#endif
#include <Eigen/Dense>
#include <Eigen/Core>

namespace CREATESOL {
class bkMeshImpl
{
public:
    /**
     * @brief getSize 输入坐标，返回该点在边界上投影点的尺寸值
     */
    double getSize(double x,double y,double z);

    /**
     * @brief getSizeSpatial 查询尺寸值
     * @return 返回该点在空间中的的尺寸值，失败则返回负数
     */
    double getSizeSpatial(double x,double y,double z);

    bkMeshImpl(const Eigen::MatrixXd &V,const Eigen::MatrixXi&F,const Eigen::VectorXd&Psize,double beta=1.2);
private:

    GEO::Mesh supportMesh;
    Eigen::VectorXd PSize;
    GEO::index_t nearest_hint;
    std::shared_ptr<GEO::MeshFacetsAABB> geo_sf_tree;
    double grad;
    GEO::Attribute<double> pointSize;
    double querytimesum;
    int searchtimes;
};
}
#endif // BKMESHIMPL_H
