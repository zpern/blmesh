#include "createsol.h"
#include <spdlog/spdlog.h> 
 #include <unordered_map>
#include <algorithm>
#include <iostream>
#include <QFile>
#include <QTextStream>
#include <igl/readMESH.h>
#include <igl/readSTL.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
using namespace std;
namespace AABBSER
{
class ActionAdd
{
public:
    ActionAdd()
    {
        ids.clear();
    }
    bool operator ()(int a)
    {
        ids.push_back(a);
    }
public:
    vector<int>ids;
};
void boxAABBcompute(Eigen::MatrixXd &V, Eigen::MatrixXi &F,vector<int>&speFacets,double eps,vector<double>&vlocaleps);
void readTriFromVTK(string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int> > &mark, map<int, int> &face_body)
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

    return;
    //end facets

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

            if(attributes==120)
            {
                /*int tempIndex=F(i,1);
                F(i,1)=F(i,2);
                F(i,2)=tempIndex;*/
            }

        }
        // discretsolid->discreteFacets[i].NumGeometry=discretsolid->discreteFacets[i].geometryID.size();
    }

    fscanf(fin, "%s %*s %*s", ignore);
    fscanf(fin, "%s %d %d%*s", ignore,&NumAttributes,&NumCells);
    int bodynum=-1;
    vector<int>bodyvec;
    for(int i=0;i<NumPolygons;i++)
    {
        fscanf(fin,"%d",&bodynum);
        bodyvec.push_back(bodynum);
    }

    for(int i=0;i<mark.size();i++)
    {
        set<int> temp=mark[i];
        for(auto it=temp.begin();it!=temp.end();it++)
        {
            int faceID=*it;
            face_body[faceID]=bodyvec[i];
        }

    }

    //end attributes

    fclose(fin);
    cout<<"Reading VTK file successfully!"<<endl;
    return ;
}

void testBoxIntersection(Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int> > &mark, map<int, int> &face_body,vector<double>&vlocaleps)
{
    set<int>bodies;


    map<int,vector<int>>body_face;
    for(auto it=face_body.begin();it!=face_body.end();it++)
    {
        bodies.insert(it->second);
    }
    map<int,vector<int>> body_facets;
    for(auto it=bodies.begin();it!=bodies.end();it++)
    {
        vector<int>facets,face;
        body_facets[*it]=facets;
        body_face[*it]=face;
    }

    for(auto it=face_body.begin();it!=face_body.end();it++)
    {
        //bodies.insert(it->second);
        int faceID=it->first;
        int bodyID=it->second;
        body_face[bodyID].push_back(faceID);

    }


    for(int i=0;i<mark.size();i++)
    {
        set<int>temp=mark[i];
        for(auto it=temp.begin();it!=temp.end();it++)
        {
            int tempMark=*it;
            int bodyMark=face_body[tempMark];
            body_facets[bodyMark].push_back(i);
        }
    }
    map<double,int >bodyEPS;
    bodyEPS[0.02]=0;
    bodyEPS[0.01]=1;

    boxAABBcompute(V,F,body_facets[0],0.02,vlocaleps);


    return;
}
/*
int main()
{
    string filename="cc_intial_dis.vtk";

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_temp;
    Eigen::MatrixXi IV, _;

    Eigen::MatrixXi F_temp;
    vector<set<int>>mark;
    map<int,int>face_body;

    readTriFromVTK(filename,V_temp,F_temp,mark,face_body);


    igl::remove_duplicate_vertices(V_temp, F_temp, 1e-4, V, IV, _, F);

    vector<double>vlocaleps;

    testBoxIntersection(V,F,mark,face_body,vlocaleps);


    return 0;
}
*/
void readpointsID(string &fname, Eigen::VectorXi &points,Eigen::VectorXi&triIDs)
{
    QFile file;

    string mid=fname;
    QString filename=QString::fromStdString(mid);
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
            if(matchWord=="vtkOriginalPointIds")
            {
                QString Vnum=charlist[2];
                QByteArray tempLine;

                QString strtemp;
            //    strtemp=strtemp.trimmed();
           //     strtemp.replace(QRegExp("[\\s]+"), " ");
                QStringList liststr;


                int num=Vnum.toInt();
                points.resize(num);
                int cnt=0;
               // int currcnt=0;
                while(cnt != num)
                {
                   // file.re
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    strtemp.replace(QRegExp("[\\s]+"), " ");
                    liststr=strtemp.split(" ");
                    int nsize=liststr.size();
                 //   cnt+=nsize;
                    for(int i=0;i<liststr.size();i++)
                    {
                        points[cnt]=liststr[i].toInt();
                        cnt++;
                    }
                }
            }
            if(matchWord=="vtkOriginalCellIds")
            {
                QString Fnum=charlist[2];
                QByteArray tempLine;

                QString strtemp;
            //    strtemp=strtemp.trimmed();
           //     strtemp.replace(QRegExp("[\\s]+"), " ");
                QStringList liststr;


                int num=Fnum.toInt();
                triIDs.resize(num);
                int cnt=0;
               // int currcnt=0;
                while(cnt != num)
                {
                   // file.re
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    strtemp.replace(QRegExp("[\\s]+"), " ");
                    liststr=strtemp.split(" ");
                    int nsize=liststr.size();
                 //   cnt+=nsize;
                    for(int i=0;i<liststr.size();i++)
                    {
                        triIDs[cnt]=liststr[i].toInt();
                        cnt++;
                    }
                }
            }
        }
    }
    file.close();
    return;
}
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
                strtemp.replace(QRegExp("[\\s]+"), " ");
                QStringList liststr=strtemp.split(" ");
                QString Vnum=liststr[0];

                int num=Vnum.toInt();
                V.resize(num,3);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    strtemp.replace(QRegExp("[\\s]+"), " ");
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
                strtemp.replace(QRegExp("[\\s]+"), " ");
                QStringList liststr=strtemp.split(" ");
                QString Fnum=liststr[0];

                int num=Fnum.toInt();
                F.resize(num,3);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    strtemp.replace(QRegExp("[\\s]+"), " ");
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
                 strtemp.replace(QRegExp("[\\s]+"), " ");
                QStringList liststr=strtemp.split(" ");
                QString Fnum=liststr[0];

                int num=Fnum.toInt();
                T.resize(num,4);
                for(int i=0;i<num;i++)
                {
                    tempLine=file.readLine();
                    strtemp=tempLine;
                    strtemp=strtemp.trimmed();
                    strtemp.replace(QRegExp("[\\s]+"), " ");
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

void boxAABBcompute(Eigen::MatrixXd &V, Eigen::MatrixXi &F,vector<int>&speFacets,double eps,vector<double>&vlocaleps)
{
    GEO::Mesh geo_sf_mesh;

    Eigen::MatrixXd VI;
    Eigen::MatrixXi FI;
    Eigen::VectorXd size;
    set<int>vids;
    map<int,int>new_old_ids;
    map<int,int>old_new_ids;
    for(int i=0;i<speFacets.size();i++)
    {
        int tempID=speFacets[i];
        for(int j=0;j<3;j++)
        {
            int verticeID=F(tempID,j);
            vids.insert(verticeID);
        }
    }
    int count=0;
    cout<<"zhvliu"<<endl;
    cout<<vids.size()<<endl;
    for(auto it=vids.begin();it!=vids.end();it++)
    {
        int verID=*it;
        cout<<verID<<" ";
        new_old_ids[count]=verID;
        old_new_ids[verID]=count;
        count++;
    }
    //  exit(0);

    VI.resize(vids.size(),3);
    FI.resize(speFacets.size(),3);


    geo_sf_mesh.vertices.clear();
    geo_sf_mesh.vertices.create_vertices((int) vids.size());
    count=0;
    for (auto  it = vids.begin(); it!=vids.end(); it++)
    {

        GEO::vec3 &p = geo_sf_mesh.vertices.point(count);
        for (int j = 0; j < 3; j++)
        {
            VI(count,j)=V(new_old_ids[count], j);
            p[j] = V(new_old_ids[count], j);
        }
        count++;
    }

    geo_sf_mesh.facets.clear();
    geo_sf_mesh.facets.create_triangles((int) speFacets.size());
    for (int i = 0; i < speFacets.size(); i++)
    {

        int tempID=speFacets[i];
        for (int j = 0; j < 3; j++)
        {

            int oldID=F(tempID, j);
            int newID=old_new_ids[oldID];
            FI(i,j)=newID;
            geo_sf_mesh.facets.set_vertex(i, j,newID );
        }
    }


    size.resize(VI.rows(),1);
    writeTriVTK("testname.vtk",VI,FI,size);

   // exit(0);

    GEO::MeshFacetsAABB geo_sf_tree(geo_sf_mesh,false);

    vector<int>intersectTri;

    for(int i=0;i<F.rows();i++)
    {
        auto it= find(speFacets.begin(),speFacets.end(),i);
        if(it!=speFacets.end())
        {
            continue;
        }
        else
        {

            int pid[3];
            for(int j=0;j<3;j++)
            {
                pid[j]=F(i,j);
            }
            vector<double >xvec,yvec,zvec;
            for(int j=0;j<3;j++)
            {
                xvec.push_back(V(pid[j],0));
                yvec.push_back(V(pid[j],1));
                zvec.push_back(V(pid[j],2));
            }
            sort(xvec.begin(),xvec.end());
            sort(yvec.begin(),yvec.end());
            sort(zvec.begin(),zvec.end());


            //            double minx,miny,minz;
            //            minx=miny=minz=DBL_MAX;
            //            double maxx,maxy,maxz;
            //            maxx=maxy=maxz=DBL_MIN;


            GEO::Box triBox;
            triBox.xyz_max[0]=xvec[2]+eps;
            triBox.xyz_max[1]=yvec[2]+eps;
            triBox.xyz_max[2]=zvec[2]+eps;
            triBox.xyz_min[0]=xvec[0]-eps;
            triBox.xyz_min[1]=yvec[0]-eps;
            triBox.xyz_min[2]=zvec[0]-eps;
            ActionAdd action;
            geo_sf_tree.compute_bbox_facet_bbox_intersections(triBox,action);
            if(action.ids.size()!=0)
            {
                cout<<"intersection within epsilion!!!"<<endl;
                intersectTri.push_back(i);
            }

        }
    }
    Eigen::MatrixXi interF;
    interF.resize(intersectTri.size(),3);
    for(int i=0;i<intersectTri.size();i++)
    {
        interF(i,0)=F(intersectTri[i],0);
        interF(i,1)=F(intersectTri[i],1);
        interF(i,2)=F(intersectTri[i],2);
    }

    size.resize(V.rows(),1);
    writeTriVTK("intersection.vtk",V,interF,size);

}
}
