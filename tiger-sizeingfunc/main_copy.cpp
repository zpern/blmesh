//#include <QCoreApplication>
#include <spdlog/spdlog.h> 
 #include "createsol.h"
#include <igl/readMESH.h>
#include <igl/readSTL.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/writePLY.h>
#include <igl/writeSTL.h>
#include <QString>
#define PROJECTSIZE 1
//using namespace Eigen;
using namespace AABBSER;
#ifndef SLIB

#if PROJECTSIZE
int main(int argc, char *argv[])
{
    // QCoreApplication a(argc, argv);

    //return a.exec();
    char*temp=nullptr;
    char*sizefile=nullptr;
    char*epsfile=nullptr;
    // char*sizefile=argv[]

    if(argc==2+1)//means only size field
    {
        temp=argv[1];
        sizefile=argv[2];

        if(1)
        {
            std::string filename;
            filename=temp;
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::VectorXd Psize;
            readfile(filename,V,F,Psize);

            bkMesh meshtemp(V,F,Psize);

            double x,y,z;
            x=30.9621;
            y=32.1197;
            z=9.75;
            //30.9621
          //  32.1197
         //   9.75
            double size=meshtemp.getSize(x,y,z);
            cout<<size<<endl;
            return 0;
        }

    }
    if(argc==3+1)//means size and eps
    {
        temp=argv[1];
        sizefile=argv[3];
        epsfile=argv[2];

    }

    bool biflag=true;
    if(biflag==true)
    {

        std::string intialm=sizefile;
        Eigen::MatrixXd mV,mN,iV;
        Eigen::MatrixXi mT;
        Eigen::MatrixXi mF,iF;
        Eigen::VectorXi IV;
        Eigen::VectorXi _;
        //igl::readMESH(intialm,mV,mT,mF);
#if 1
        QString mesh2=QString::fromStdString(intialm);
        read_mesh_version2(mesh2,mV,mF,mT);
#endif
        //   igl::read_triangle_mesh(intialm,iV,iF);
        //igl::readMESH();

        // igl::remove_duplicate_vertices(iV, iF, 1e-6, mV, IV, _, mF);
        std::string filename;
        filename=temp;
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd Psize;
        readfile(filename,V,F,Psize);

        bkMesh meshtemp(V,F,Psize);

    //    getchar();

//        for(int i=0;i<1;i++)
//        {
//            double x,y,z;
//            x=0;
//            y=0;
//            z=0;
//            double sizetemp=meshtemp.getSize(x,y,z);
//            sizetemp=meshtemp.getSizeSpatial(x,y,z);
//            cout<<"zhvliu"<<endl;
//        }



        Eigen::VectorXd size;
        size.resize(mV.rows(),1);
        for(int i=0;i<mV.rows();i++)
        {
            double x,y,z;
            x=mV(i,0);
            y=mV(i,1);
            z=mV(i,2);
            size(i,0)=meshtemp.getSize(x,y,z);
        }
        cout<<"end size sol"<<endl;

        cout<<"search times : "<<meshtemp.searchtimes<<endl;
        cout<<"query time: "<<(double)(meshtemp.querytimesum)/1000<<endl;

        std::string sizename=intialm+"_size.vtk";
        writeTriVTK(sizename,mV,mF,size);
        std::string tempstr=intialm+"size.sol";
        writeSol(tempstr,mV,size);

        if(argc==3+1)
        {
            Eigen::MatrixXd Veps;
            Eigen::MatrixXi Feps;
            vector<set<int>>Fmark;
            Eigen::VectorXd Peps;
            string epsvtk=epsfile;
            readEpsFromVTK(epsvtk,Veps,Feps,Fmark,Peps);
            rm_duplicate_facets(Feps);
            bkMesh epsMesh(Veps,Feps,Peps);
            Eigen::VectorXd eps;
            eps.resize(mV.rows(),1);
            for(int i=0;i<mV.rows();i++)
            {
                double x,y,z;
                x=mV(i,0);
                y=mV(i,1);
                z=mV(i,2);
                eps(i,0)=epsMesh.getSize(x,y,z);
            }
            string epsname=intialm+"_eps.vtk";
            writeTriVTK(epsname,mV,mF,eps);
            string localParaname="DEFAULT.mmgs";
            writeLocalParameterForMMGS(localParaname,mV,mF,size,eps);
            cout<<"end eps sol"<<endl;
        }

        //  meshtemp.getSize(4300,0,215.98);


    }

    return 0;
}
#else
int main(int argc, char *argv[])
{
    bool transPLY=false;

    char*temp=nullptr;
    char*sizefile=nullptr;
    char*epsfile=nullptr;

    char *pointsFile=nullptr;
    if(transPLY==true)
    {
        sizefile=argv[1];
        std::string intialm=sizefile;
        Eigen::MatrixXd mV,mN,iV;
        Eigen::MatrixXi mT;
        Eigen::MatrixXi mF,iF;
        Eigen::VectorXi IV;
        Eigen::VectorXi _;
        //igl::readMESH(intialm,mV,mT,mF);
        QString mesh2=QString::fromStdString(intialm);
        read_mesh_version2(mesh2,mV,mF,mT);
        string plyname=intialm+".ply";
        igl::writePLY(plyname,mV,mF);
        return 0;
    }
    bool patchMesh=false;
    if(patchMesh==true)
    {
        sizefile=argv[1];
        std::string intialm=sizefile;
        Eigen::MatrixXd mV,mN,iV;
        Eigen::MatrixXi mT;
        Eigen::MatrixXi mF,iF,mE;
        Eigen::VectorXi IV,Fmark;
        Eigen::VectorXi _;
        //igl::readMESH(intialm,mV,mT,mF);
        QString mesh2=QString::fromStdString(intialm);
        read_mesh_version2(mesh2,mV,mF,mT);

        pointsFile=argv[2];
        string pfile=pointsFile;
        Eigen::VectorXi pointsIDs,triIDs;
        Fmark.setZero(mF.rows());
        readpointsID(pfile,pointsIDs,triIDs);
        for(int i=0;i<triIDs.rows();i++)
        {
            int fid=triIDs[i];
            Fmark[fid]=1;
        }

        string plyname=intialm+".ply";
        igl::writePLY(plyname,mV,mF);
        string medit=intialm.substr(0,intialm.length()-5);
        medit+="_patch.mesh";
        writeMeditMesh(medit,mV,mE,mF,Fmark,0);
   //     return 0;
        return 0;
    }

    if(argc==2+1+1)//means only size field
    {
        temp=argv[1];
        sizefile=argv[2];
        pointsFile=argv[3];
    }

    bool biflag=true;
    if(biflag==true)
    {

        std::string intialm=sizefile;
        Eigen::MatrixXd mV,mN,iV;
        Eigen::MatrixXi mT;
        Eigen::MatrixXi mF,iF;
        Eigen::VectorXi IV;
        Eigen::VectorXi _;
        //igl::readMESH(intialm,mV,mT,mF);
#if 1
        QString mesh2=QString::fromStdString(intialm);
        read_mesh_version2(mesh2,mV,mF,mT);
#endif

        std::string filename;
        filename=temp;
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd Psize;

        vector<set<int>>Fmark;
        map<int,int>face_body;



        readTriFromVTK(filename,V,F,Fmark,face_body);

     //   igl::writePLY("vtktest.ply",V,F);

        Psize.setZero(V.rows());
        bkMesh meshtemp(V,F,Psize);

        string pfile=pointsFile;
        Eigen::VectorXi pointsIDs,triIDs;
        readpointsID(pfile,pointsIDs,triIDs);

      //  bool projection=true;
        cout<<"projectSize: "<<pointsIDs.size()<<endl;
        for(int i=0;i<pointsIDs.rows();i++)
        {
            int pid=pointsIDs[i];
            GEO::vec3 current_point;
            current_point[0]=mV(pid,0);
            current_point[1]=mV(pid,1);
            current_point[2]=mV(pid,2);
            GEO::vec3 pro_point;
            double dis=0;
            meshtemp.geo_sf_tree->nearest_facet(current_point,pro_point,dis);
            mV(pid,0)=pro_point[0];
            mV(pid,1)=pro_point[1];
            mV(pid,2)=pro_point[2];
        }
        std::string sizename=intialm+"_pro.ply";
        igl::writePLY(sizename,mV,mF);
        sizename=intialm+"_pro.stl";
        igl::writeSTL(sizename,mV,mF);
    //    writeTriVTK(sizename,mV,mF,size);
   //     std::string tempstr=intialm+"size.sol";
    //    writeSol(tempstr,mV,size);

    }

    return 0;
}



#endif
#endif
