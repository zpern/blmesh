#include <spdlog/spdlog.h> 
 #include "bkmeshimpl.h"

using namespace std;

namespace CREATESOL {
bkMeshImpl::bkMeshImpl(const Eigen::MatrixXd &V,const Eigen::MatrixXi&F,const Eigen::VectorXd&Psize,double beta)
{
    cout<<"begining..."<<endl;
    PSize=Psize;
    nearest_hint=GEO::NO_FACET;
    grad=beta;

    GEO::initialize();

   // getchar();


    supportMesh.vertices.clear();
    supportMesh.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < V.rows(); i++) {
        GEO::vec3 &p = supportMesh.vertices.point(i);
        for (int j = 0; j < 3; j++)
            p[j] = V(i, j);
    }
  //  getchar();

    cout<<"vertex finished!!!..."<<endl;

    pointSize.create_vector_attribute(supportMesh.vertices.attributes(),"size",1);

    cout<<"binding finished!!!"<<endl;

    int attrsizze=pointSize.size();
   // pointSize.data();
    for(int i=0;i<Psize.rows();i++)
    {
        pointSize[i]=Psize[i];
    }

  //  getchar();
    supportMesh.facets.clear();
    supportMesh.facets.create_triangles((int) F.rows());
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            supportMesh.facets.set_vertex(i, j, F(i, j));
    }
  //  getchar();

 //   GEO::MeshFacetsAABB *geo_sf_treeMid=new GEO::MeshFacetsAABB(supportMesh,true);
    cout<<"facet fineished!!!"<<endl;

    geo_sf_tree=std::shared_ptr<GEO::MeshFacetsAABB>(new GEO::MeshFacetsAABB(supportMesh));
    cout<<"construct finished!!!"<<endl;

   // getchar();

   // double
    cout<<"size of points: (kB): "<<(double)sizeof(supportMesh.vertices)/1024<<endl;
    cout<<"size of facets: (kB): "<<(double)sizeof(supportMesh.facets)/1024<<endl;
    cout<<"size of backgroundmesh : (kB): "<<(double)sizeof(supportMesh)/1024<<endl;
    cout<<"size of AABBtree : (kB): "<<(double)sizeof(*geo_sf_tree)/1024<<endl;

    querytimesum=0;
    searchtimes=0;
}

double getTriSquare(Vector &pa, Vector &pb, Vector &pc)
{
    double volume=0.0;
    double ab=pa.getDistance(pb);
    double ac=pa.getDistance(pc);
    double bc=pb.getDistance(pc);

    double c=(ab+ac+bc)/2;

    volume=sqrt(c*(c-ab)*(c-ac)*(c-bc));

    if(std::isnan(volume))
    {
    //    cout<<"Nan"<<endl;
        volume=0;
      //  return volume;
    }

    //volume=sqrt(c*(c-ab)*(c-ac)*(c-bc));

    return volume;
}

double bkMeshImpl::getSizeSpatial(double x, double y, double z)
{
    searchtimes++;
    double result=0;

    GEO::vec3 nearest_point;
    double sq_dist = std::numeric_limits<double>::max();
    GEO::index_t prev_facet = GEO::NO_FACET;
   // nearest_hint=GEO::NO_FACET;
    int cnt = 0;

    GEO::vec3 current_point;
    current_point[0]=x;
    current_point[1]=y;
    current_point[2]=z;

    GEO::index_t facetID=-1;
    double dis=0;

  //  facetID=nearest_hint;

    facetID=nearest_hint;

    clock_t start = clock();//开始计时

    facetID=geo_sf_tree->nearest_facet(current_point,nearest_point,dis);

    double duration = (clock() - start)*1000.0/CLOCKS_PER_SEC;
    querytimesum+=duration;

   // geo_sf_tree->nearest_facet_with_hint(current_point,facetID,nearest_point,dis);

    nearest_hint=facetID;
    //geo_sf_tree.get_nearest_facet_hint();


    int id1=supportMesh.facets.vertex(facetID,0);
    int id2=supportMesh.facets.vertex(facetID,1);
    int id3=supportMesh.facets.vertex(facetID,2);

    GEO::vec3 point1=supportMesh.vertices.point(id1);
    GEO::vec3 point2=supportMesh.vertices.point(id2);
    GEO::vec3 point3=supportMesh.vertices.point(id3);
    Vector p1,p2,p3,p_current;
    p1.x=point1[0];
    p1.y=point1[1];
    p1.z=point1[2];
    p2.x=point2[0];
    p2.y=point2[1];
    p2.z=point2[2];
    p3.x=point3[0];
    p3.y=point3[1];
    p3.z=point3[2];
    p_current.x=nearest_point[0];
    p_current.y=nearest_point[1];
    p_current.z=nearest_point[2];

    double size1=pointSize[id1];
    double size2=pointSize[id2];
    double size3=pointSize[id3];

    double volume1=getTriSquare(p1,p2,p_current);
    double volume2=getTriSquare(p1,p3,p_current);
    double volume3=getTriSquare(p2,p3,p_current);

    double sum=volume1+volume2+volume3;

    if(sum==0)
    {
        result=(size1+size2+size3)/3;
    }
    else
    {
        result=size3*volume1/sum+size2*volume2/sum+size1*volume3/sum;
    }

    if(std::isnan(sum))
    {
        return -1;
    }
//    cout<<"squared dis:"<<dis<<endl;



    //calculate the size

    //geo_sf_tree.get_nearest_facet_hint(current_point,facetID,nearest_point,dis);

    dis=sqrt(dis);
  //  return dis;

    static int count=0;
    //spdlog::info(count<<": "<<dis<<" "<<grad<<" "<<result);

//	cout << "dis=" << dis << " orisize=" << result << endl;
    result=std::max(0.0,dis-1.0*result)*(grad-1)+result;

    return result;
}

double bkMeshImpl::getSize(double x, double y, double z)
{
    searchtimes++;
    double result=0;

    GEO::vec3 nearest_point;
    double sq_dist = std::numeric_limits<double>::max();
    GEO::index_t prev_facet = GEO::NO_FACET;
   // nearest_hint=GEO::NO_FACET;
    int cnt = 0;

    GEO::vec3 current_point;
    current_point[0]=x;
    current_point[1]=y;
    current_point[2]=z;

    GEO::index_t facetID=-1;
    double dis=0;

  //  facetID=nearest_hint;

    facetID=nearest_hint;

    clock_t start = clock();//开始计时
    facetID=geo_sf_tree->nearest_facet(current_point,nearest_point,dis);

    double duration = (clock() - start)*1000.0/CLOCKS_PER_SEC;
    querytimesum+=duration;
   // geo_sf_tree->nearest_facet_with_hint(current_point,facetID,nearest_point,dis);

    nearest_hint=facetID;
    //geo_sf_tree.get_nearest_facet_hint();


    int id1=supportMesh.facets.vertex(facetID,0);
    int id2=supportMesh.facets.vertex(facetID,1);
    int id3=supportMesh.facets.vertex(facetID,2);

    GEO::vec3 point1=supportMesh.vertices.point(id1);
    GEO::vec3 point2=supportMesh.vertices.point(id2);
    GEO::vec3 point3=supportMesh.vertices.point(id3);
    Vector p1,p2,p3,p_current;
    p1.x=point1[0];
    p1.y=point1[1];
    p1.z=point1[2];
    p2.x=point2[0];
    p2.y=point2[1];
    p2.z=point2[2];
    p3.x=point3[0];
    p3.y=point3[1];
    p3.z=point3[2];
    p_current.x=nearest_point[0];
    p_current.y=nearest_point[1];
    p_current.z=nearest_point[2];

    double size1=pointSize[id1];
    double size2=pointSize[id2];
    double size3=pointSize[id3];

    double volume1=getTriSquare(p1,p2,p_current);
    double volume2=getTriSquare(p1,p3,p_current);
    double volume3=getTriSquare(p2,p3,p_current);

    double sum=volume1+volume2+volume3;

    if(sum==0)
    {
        cout<<"Degenerated!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        result=size1+size2+size3;
        result=result/3;
        return result;
    }

    result=size3*volume1/sum+size2*volume2/sum+size1*volume3/sum;


    if(std::isnan(sum))
    {
        cout<<"Degenerated!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        exit(0);
    }
//    cout<<"squared dis:"<<dis<<endl;



    //calculate the size

    //geo_sf_tree.get_nearest_facet_hint(current_point,facetID,nearest_point,dis);

    return result;
}
}

