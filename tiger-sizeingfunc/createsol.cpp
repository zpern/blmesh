#include "createsol.h"
#include <spdlog/spdlog.h> 
 #include <geogram/mesh/mesh_AABB.h>
#include <unordered_map>
using namespace std;
namespace AABBSER
{
void writeTriVTK(std::string filename, Eigen::MatrixXd &VI, Eigen::MatrixXi &FI, Eigen::VectorXd &markI)
{
      int ngridpts=VI.rows();

       int index=-1;
       const char*file=filename.data();
       FILE *fout = fopen(file, "w");
       if (!fout)
       {
           spdlog::info("Cannot open file %s.", filename.c_str());
           return  ;
       }

       fprintf(fout,"# vtk DataFile Version 2.0\n");
       fprintf(fout,"Discrete surface mesh\n");
       fprintf(fout,"ASCII\n");
       fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
       fprintf(fout,"POINTS %d double\n", ngridpts);
       //output background mesh node //start from Zero
       for (int i = 0; i < ngridpts; i++)
       {
           int index = i + 1;
           double x = VI(i,0);
           double y = VI(i,1);
           double z = VI(i,2);
           fprintf(fout,"%.5f %.5f %.5f\n",x, y, z);
       }

       int NoVertex = 0;
   //    fprintf(fout,"VERTICES %d %d\n", NoVertex,NoVertex*2);

       int NoLines=0;
       int p1=-1;
       int p2=-1;

     //  fprintf(fout,"LINES %d %d\n", NoLines,NoLines*3);


       int ngridels=FI.rows();
       int p3;
       fprintf(fout,"CELLS %d %d\n",ngridels,ngridels*4);
       int marksize=0;
       for(int i=0;i<ngridels;i++)
       {
           p1=FI(i,0);
           p2=FI(i,1);
           p3=FI(i,2);

           fprintf(fout,"%d %d %d %d\n",3,p1,p2,p3);

       }
       fprintf(fout,"CELL_TYPES %d\n",ngridels);

       for(int i=0;i<ngridels;i++)
       {


           fprintf(fout,"5\n");

       }

       int NoCells=NoVertex+NoLines+ngridels;
       //cout<<NoCells<<" "<<NoLines<<" "<<ngridels<<endl;
       fprintf(fout,"POINT_DATA %d\n",ngridpts);
       fprintf(fout,"SCALARS fixed float\n");
       fprintf(fout,"LOOKUP_TABLE default\n");
       int counter=0;

       for(int i=0;i<ngridpts;i++)
       {

           double attribute=markI(i,0);
           fprintf(fout,"%lf ",attribute);

           fprintf(fout,"\n");


       }
       fclose(fout);
       return ;
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
/*
void read_mesh_version2(QString &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F,Eigen::MatrixXi &T)
{
    QFile file;

    file.setFileName(filename);
    QByteArray t;
    if(file.open(QIODevice::ReadOnly |QIODevice::Text))
    {

        while(!file.atEnd())
        {
            QByteArray line=file.readLine();

            QString str(line);
            str=str.trimmed();
           // qDebug()<<str;

         //   exit(0);
            QStringList charlist=str.split(" ");
            QString matchWord=charlist[0];
            if(matchWord=="Vertices")
            {
                QByteArray tempLine=file.readLine();
                QString strtemp(tempLine);
                strtemp=strtemp.trimmed();
                QStringList liststr=strtemp.split(" ");
                QString Vnum=liststr[0];

                int num=Vnum.toInt();
                V.resize(num,3);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    liststr=strtemp.split(" ");
                    V(i,0)=liststr[0].toDouble();
                    V(i,1)=liststr[1].toDouble();
                    V(i,2)=liststr[2].toDouble();
                }
            }
            if(matchWord=="Triangles")
            {
                QByteArray tempLine=file.readLine();
                QString strtemp(tempLine);
                strtemp=strtemp.trimmed();
                QStringList liststr=strtemp.split(" ");
                QString Fnum=liststr[0];

                int num=Fnum.toInt();
                F.resize(num,3);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    liststr=strtemp.split(" ");
                    F(i,0)=liststr[0].toInt()-1;
                    F(i,1)=liststr[1].toInt()-1;
                    F(i,2)=liststr[2].toInt()-1;
                }

            }
            if(matchWord=="Tetrahedra")
            {
                QByteArray tempLine=file.readLine();
                QString strtemp(tempLine);
                strtemp=strtemp.trimmed();
                QStringList liststr=strtemp.split(" ");
                QString Fnum=liststr[0];

                int num=Fnum.toInt();
                T.resize(num,4);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    liststr=strtemp.split(" ");
                    T(i,0)=liststr[0].toInt()-1;
                    T(i,1)=liststr[1].toInt()-1;
                    T(i,2)=liststr[2].toInt()-1;
                    T(i,3)=liststr[3].toInt()-1;
                }

            }
         //   auto iter=find(specificLine.begin(),specificLine.end(),matchWord)
        }
    }

    file.close();
    return;
}

*/
void writeSol(std::string filename, Eigen::MatrixXd &V, Eigen::VectorXd &size)
{
    const char*file=filename.data();
    FILE *fout = fopen(file, "w");
    if (!fout)
    {
        spdlog::info("Cannot open file %s.", filename.c_str());
        return  ;
    }

    fprintf(fout,"MeshVersionFormatted 2\n");
    fprintf(fout,"\n");
    fprintf(fout,"Dimension 3\n");
    fprintf(fout,"\n");
    fprintf(fout,"SolAtVertices\n");
    fprintf(fout,"%d\n",size.rows());
    fprintf(fout,"1 1\n");
    fprintf(fout,"\n");
    for(int i=0;i<size.rows();i++)
    {
        double temp=size(i,0);
        temp=1/(temp*temp);
        fprintf(fout,"%lf\n",size(i,0));
        //fprintf(fout,"%lf 0 0 %lf 0 %lf\n",temp,temp,temp);
    }

    fprintf(fout,"End\n");

    fclose(fout);

    //output background mesh node //start from Zero


    return ;
}
void readfile(std::string &filename,Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXd &Psize)
{
    int NumPoints;
    int NumVertices;
    int NumLines;
    int NumPolygons;
    double x,y,z;

    const char*file=filename.data();
    FILE *fin = fopen(file, "r");
    if (fin == NULL)
    {
        fprintf(stderr,"Error in reading vtk file:%s\n",filename.c_str());
        return;
    }
    char ignore[256];
    char ignore2[256];
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    //ignore the first four lines

    //reading points part
    fscanf(fin, "%*s%d%s",&NumPoints,ignore);
    V.resize(NumPoints,3);
    for(int i=0;i<NumPoints;i++)
    {
        fscanf(fin, "%lf%lf%lf", &x, &y,&z);
        V(i,0)=x;
        V(i,1)=y;
        V(i,2)=z;
    }
    //end points

    //reading  vertex
    int pointIndex;

    int count=0;
    int pointID1,pointID2;


    //reading facets
    count=0;
    int point[3];
    //unordered_map<string,bool>::iterator triangleIter;
    pair<string,bool>trianglePair;
    fscanf(fin, "%s %d %[^\n]", ignore,&NumPolygons,ignore2);
   // vtkinput.NumTraingles=NumPolygons;
   // vtkinput.triangles=new TriangleInput[NumPolygons];
    F.resize(NumPolygons,3);
//    mark.resize(NumPolygons);
    for(int i=0;i<NumPolygons;i++)
    {
        fscanf(fin, "%s %d %d %d", ignore,&point[0],&point[1],&point[2]);

        F(i,0)=point[0];
        F(i,1)=point[1];
        F(i,2)=point[2];
    }
    fscanf(fin, "\n%[^\n]", ignore);
    for(int i=0;i<NumPolygons;i++)
    {
         fscanf(fin, "\n%[^\n]", ignore);
    }
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    Psize.resize(NumPoints,1);
    double size=0;
    for(int i=0;i<NumPoints;i++)
    {
        fscanf(fin, "%lf",&size );

        Psize(i,0)=size;
    }

    fclose(fin);
    cout<<"Reading VTK file successfully!"<<endl;
    return ;

}
void writeLocalParameterForMMGS(string &filename,Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXd &edgeSize, Eigen::VectorXd &hausd)
{
    const char*file=filename.data();
    FILE *fout = fopen(file, "w");
    if (!fout)
    {
        spdlog::info("Cannot open file %s.", filename.c_str());
        return  ;
    }
    int numTri=F.rows();

    fprintf(fout,"parameters\n");
    fprintf(fout,"%d\n",numTri);
    for(int i=0;i<numTri;i++)
    {
        int pID[3];
        pID[0]=F(i,0);
        pID[1]=F(i,1);
        pID[2]=F(i,2);
        double minSize=-1;
        double maxSize=-1;
        double minHausd=-1;
        vector<double> psize;
        psize.resize(3);
        vector<double>allEps;
        allEps.resize(3);
        for(int j=0;j<3;j++)
        {
            psize[j]=edgeSize(pID[j],0);
            allEps[j]=hausd(pID[j],0);
        }
        sort(psize.begin(),psize.end());
        sort(allEps.begin(),allEps.end());
        minSize=psize[0];
        maxSize=psize[2];
        minHausd=allEps[0];
        fprintf(fout,"%d Triangles %lf %lf %lf\n",i+1,minSize,maxSize,minHausd);


    }
 /*   parameters
    3
    38 Triangles 1.8     2.2   0.01
    36 Triangles 0.098   0.12  0.1
    37 Triangles 4.8     5.2   1
   */

    fclose(fout);
    return;
}
void rm_duplicate_facets(Eigen::MatrixXi &F)
{
    int fsize=F.rows();
    set<set<int>>facets;
    for(int i=0;i<fsize;i++)
    {
        set<int>temp;
        temp.insert(F(i,0));
        temp.insert(F(i,1));
        temp.insert(F(i,2));
        facets.insert(temp);
    }
    assert(fsize>=facets.size());
    F.resize(facets.size(),3);
    int count=0;
    for(auto it=facets.begin();it!=facets.end();it++)
    {
        set<int> tri=*it;
        int subcount=0;
        for(auto k=tri.begin();k!=tri.end();k++)
        {
            F(count,subcount)=*k;
            subcount++;
        }
        count++;
    }

    return;
}
using namespace std;
void readEpsFromVTK(string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int> > &mark, Eigen::VectorXd  &eps)
{
    int NumPoints;
    int NumVertices;
    int NumLines;
    int NumPolygons;
    double x,y,z;

    const char*file=filename.data();
    FILE *fin = fopen(file, "r");
    if (fin == NULL)
    {
        //fprintf(stderr,"Error in reading vtk file:%s\n",filename);
        cout<<"Error in reading vtk file:"<<endl;
        return;
    }
    char ignore[256];
    char ignore2[256];
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);
    //ignore the first four lines

    //reading points part
    fscanf(fin, "%*s%d%s",&NumPoints,ignore);
    V.resize(NumPoints,3);
    for(int i=0;i<NumPoints;i++)
    {
        fscanf(fin, "%lf%lf%lf%", &x, &y,&z);
        V(i,0)=x;
        V(i,1)=y;
        V(i,2)=z;
    }
    cout<<"initial points: "<<NumPoints<<endl;
    //end points

    //reading  vertex
    int pointIndex;;
    fscanf(fin, "%s %d %s%[^\n]", ignore,&NumVertices,ignore2);
    int count=0;
    int pointID1,pointID2;
    fscanf(fin, "%s %d %s%[^\n]", ignore,&NumLines,ignore2);

    cout<<"initial line: "<<NumLines<<endl;

    //reading facets
    count=0;
    int point[3];
    unordered_map<string,bool>::iterator triangleIter;
    pair<string,bool>trianglePair;
    fscanf(fin, "%s %d %s%[^\n]", ignore,&NumPolygons,ignore2);
    cout<<"initial facets: "<<NumPolygons<<endl;
   // vtkinput.NumTraingles=NumPolygons;
   // vtkinput.triangles=new TriangleInput[NumPolygons];
    F.resize(NumPolygons,3);
    mark.resize(NumPolygons);
    for(int i=0;i<NumPolygons;i++)
    {
        fscanf(fin, "%s %d %d %d", ignore,&point[0],&point[1],&point[2]);
     //   sort(point, point + 3);
        //vtkinput.triangles[i].points[0];
       // vtkinput.triangles[i].points[0]=point[0];
       // vtkinput.triangles[i].points[1]=point[1];
       // vtkinput.triangles[i].points[2]=point[2];
        F(i,0)=point[0];
        F(i,1)=point[1];
        F(i,2)=point[2];
    }
    cout<<"initial facets: "<<NumPolygons<<endl;
    //end facets
    int epsSize=0;
    fscanf(fin, "%s %d", ignore,&epsSize);
    fscanf(fin, "\n%[^\n]", ignore);
    fscanf(fin, "\n%[^\n]", ignore);


    eps.resize(epsSize,1);
    for(int i=0;i<epsSize;i++)
    {
        double outsize=0;
        fscanf(fin,"%lf",&outsize);
        eps(i,0)=outsize;
    }


    //reading attributes
    int NumCells;
    int NumAttributes=0;
    int attributes=-1;
    fscanf(fin, "%s %d", ignore,&NumCells);
    fscanf(fin, "%s %*s %*s", ignore);
    fscanf(fin, "%s %d %d%*s", ignore,&NumAttributes,&NumCells);


    assert(NumCells==(NumVertices+NumLines+NumPolygons));

    for(int i=0;i<NumPolygons;i++)
    {
        for(int j=0;j<NumAttributes;j++)
        {

            fscanf(fin,"%d",&attributes);
            if(attributes!=-1)
            {
                mark[i].insert(attributes);
            }
        }
        // discretsolid->discreteFacets[i].NumGeometry=discretsolid->discreteFacets[i].geometryID.size();
    }

    //end attributes

    fclose(fin);
    cout<<"Reading VTK file successfully!"<<endl;
    return ;
}
bkMesh::~bkMesh(){
    delete geo_sf_tree;
}
bkMesh::bkMesh(Eigen::MatrixXd &V, Eigen::MatrixXi&F, Eigen::VectorXd&Psizem, double beta)
{
    cout<<"begining..."<<endl;
    PSize=Psizem;
    nearest_hint=GEO::NO_FACET;
    grad=beta;

    GEO::initialize();

    supportMesh.vertices.clear();
    supportMesh.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < V.rows(); i++) {
        GEO::vec3 &p = supportMesh.vertices.point(i);
        for (int j = 0; j < 3; j++)
            p[j] = V(i, j);
    }

    cout<<"vertex finished!!!..."<<endl;

    pointSize.create_vector_attribute(supportMesh.vertices.attributes(),"size",1);

    cout<<"binding finished!!!"<<endl;

    int attrsizze=pointSize.size();
   // pointSize.data();
    for(int i=0;i<Psizem.rows();i++)
    {
        pointSize[i]=Psizem[i];
    }

    supportMesh.facets.clear();
    supportMesh.facets.create_triangles((int) F.rows());
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            supportMesh.facets.set_vertex(i, j, F(i, j));
    }

 //   GEO::MeshFacetsAABB *geo_sf_treeMid=new GEO::MeshFacetsAABB(supportMesh,true);
    cout<<"facet fineished!!!"<<endl;

    geo_sf_tree=new GEO::MeshFacetsAABB(supportMesh);
    cout<<"construct finished!!!"<<endl;

}
double bkMesh::getSizeSpatial(double x, double y, double z)
{
    {
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

        facetID=nearest_hint;


        facetID=geo_sf_tree->nearest_facet(current_point,nearest_point,dis);
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



//        double size1=PSize(id1,0);
//        double size2=PSize(id2,0);
//        double size3=PSize(id3,0);
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
           // cout<<"Degenerated!!!!!!!!!!!!!!!!!!!!!!"<<endl;
           // exit(0);
        }
        else
        {

            result=size3*volume1/sum+size2*volume2/sum+size1*volume3/sum;
        }

        if(std::isnan(sum))
        {
            cout<<"Degenerated!!!!!!!!!!!!!!!!!!!!!!"<<endl;
            exit(0);
        }
    //    cout<<"squared dis:"<<dis<<endl;



        //calculate the size

        //geo_sf_tree.get_nearest_facet_hint(current_point,facetID,nearest_point,dis);

        dis=sqrt(dis);
      //  return dis;
        result=dis*(grad-1)+result;

        return result;
    }

}
double bkMesh::getSpanSize(double x, double y, double z)
{

	return 1.;
}
double bkMesh::getSize(double x, double y, double z)
{
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


    facetID=geo_sf_tree->nearest_facet(current_point,nearest_point,dis);
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



//    double size1=PSize(id1,0);
//    double size2=PSize(id2,0);
//    double size3=PSize(id3,0);
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
        exit(0);
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
