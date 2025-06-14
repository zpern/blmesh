#include "inputdata.h"
#include <spdlog/spdlog.h> 
 #include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <cstring>
bool operator<( const vertex & _v0,  const vertex & _v1 )
{
    const float eps_=1e-6;
    if (fabs(_v0.x - _v1.x) <= eps_)
    {
        if (fabs(_v0.y - _v1.y) <= eps_)
        {
            return (_v0.z < _v1.z - eps_);
        }
        else return (_v0.y < _v1.y - eps_);
    }
    else return (_v0.x < _v1.x - eps_);
}
InputData::InputData():source_t(SOUCETYPE::MESH_SF_SOUCE) {};
InputData::InputData(string g,string s,string t):gm3file(g),Tran_ConfigFile(t),stlplsfile(s), source_t(SOUCETYPE::MESH_SF_SOUCE)
{

}

double InputData::getSolidSouceSize(double p[3])
{
	Vector point(p);
	point = point - coord_zero;

	double length = 0.0;
	if (source_t == SOUCETYPE::CUBOID_SOLID_SOUCE) {
		double dot;
		for (int k = 0; k < 3; k++) {
			dot = abs(point * coord[k].normalized()) - coord[k].magnitude();
			if (dot > 0) {
				length += dot *dot;
			}
		}
		length = sqrt(length);
	}
	else if(source_t==SOUCETYPE::ELLIOPSOID_SOLID_SOUCE){
		double dot=0.;
		for (int k = 0; k < 3; k++) {
			double exp_length = point * coord[k]/ (coord[k].magnitude()*coord[k].magnitude());
			dot += exp_length* exp_length;
		}
		if (dot > 1) {
			dot = (sqrt(dot) - 1)/sqrt(dot);
			length = dot * point.magnitude();
		}

	}
	//cout << "length=" << length << endl;
	return basic_source+(volmGrowthRatio-1)*length;

}

//int cvNums; int cvPtNums[]; double cvPts[];
//int fcNums; int fcPtNums[]; double fcPts[];
//int lpCvNums[]; int lpCvs[];
bool InputData::readGm3() {
    int lpCvsIdx=0;
    char idumy[25];
    int nvertices_G,ncurves_G,nfaces_G,nloops_G,ndomains_G;
	//cout<<"gm3:" << gm3file << endl;
    FILE *fin = fopen(gm3file.c_str(), "rb");
    if (fin == NULL)
    {
        fprintf(stderr,"Error in reading gm3 file\n");
        return 0;
    }
    fscanf(fin, "\n%[^\n]", idumy);
    if (fscanf(fin, "%d %d %d\n",&nvertices_G, &ncurves_G, &nfaces_G) == EOF)
    {
        fprintf( stderr, "ERROR MESSAGE: ncurves and nfaces.\n");
        return 0;
    }
    cvNums=ncurves_G;

    //point label
    fscanf(fin, "\n%[^\n]", idumy);
    gm3V.resize(nvertices_G);
    for (int i = 0; i < nvertices_G; i++)
    {
        float x,y,z;
        int idx;
        fscanf(fin, "%d%e%e%e",&idx, &x, &y,&z);
        gm3V[i].x = x;
        gm3V[i].y = y;
        gm3V[i].z = z;
        gm3V[i].idx = idx;
    }//end for adding points
    //curve label
    fscanf(fin, "\n%[^\n]", idumy);
    gm3C.resize(ncurves_G);

    cvPtNums.resize(cvNums);
    for (int i = 0; i < ncurves_G; i++)
    {
        int idx,type,ncp;
        fscanf(fin,"%d %d %d",&idx, &type, &ncp);
        gm3C[i].idx=idx;
        gm3C[i].type=type;
        gm3C[i].ncp=ncp;

        cvPtNums[i]=ncp;
        for(int j=0;j<ncp;j++) {
            float x,y,z;int idx;
            fscanf(fin, "%e%e%e%d",&x, &y, &z,&idx);
            gm3C[i].v.push_back(idx);
        }
    }
    cvPtNumsAll=0;
    for(int i=0;i<cvNums;i++) {
        cvPtNumsAll+=cvPtNums[i];
    }
    cvPts.resize(cvPtNumsAll*3);
    int cvPtsIdx=0;
    for(int i=0;i<cvNums;i++) {
        for(int j=0;j<cvPtNums[i];j++) {
            cvPts[cvPtsIdx*3]=gm3V[gm3C[i].v[j]-1].x;
            cvPts[cvPtsIdx*3+1]=gm3V[gm3C[i].v[j]-1].y;
            cvPts[cvPtsIdx*3+2]=gm3V[gm3C[i].v[j]-1].z;
            cvPtsIdx++;
        }
    }
    // Support_Surfaces label
    fscanf(fin, "\n%[^\n]", idumy);
    gm3S.resize(nfaces_G);
    fcNums=nfaces_G;
    fcPtNums.resize(2*nfaces_G);
    for(int i=0;i<nfaces_G;i++) {
        int idx,type,nsp1,nsp2;
        fscanf(fin, "%d %d %d %d\n", &idx, &type, &nsp1, &nsp2);
        gm3S[i].idx=idx;
        gm3S[i].type=type;
        gm3S[i].nsp1=nsp1;
        gm3S[i].nsp2=nsp2;
        gm3S[i].v.resize(nsp2);

        fcPtNums[i*2]=nsp1;
        fcPtNums[i*2+1]=nsp2;

        for(int j=0;j<nsp2;j++) {
            gm3S[i].v[j].resize(nsp1);
            for(int k=0;k<nsp1;k++) {
                int idx;float x,y,z;
                fscanf(fin, "%e%e%e%d",&x, &y, &z,&idx);
                gm3S[i].v[j][k]=idx;
            }
        }
    }
    fcPtsNumsAll=0;
    for(int i=0;i<fcNums;i++) {
        fcPtsNumsAll+=fcPtNums[i*2]*fcPtNums[i*2+1];
    }
    ptNumsAll=fcPtsNumsAll+cvPtNumsAll;
    int fcPtsId_local=0;
    int lpCvsIdx_local=0;
    fcPts.resize(3*fcPtsNumsAll);
    for(int i=0;i<fcNums;i++) {
        for(int j=0;j<gm3S[i].nsp2;j++) {
            for(int k=0;k<gm3S[i].nsp1;k++) {
                fcPts[fcPtsId_local*3]=gm3V[gm3S[i].v[j][k]-1].x;
                fcPts[fcPtsId_local*3+1]=gm3V[gm3S[i].v[j][k]-1].y;
                fcPts[fcPtsId_local*3+2]=gm3V[gm3S[i].v[j][k]-1].z;
                fcPtsId_local++;
            }
        }
    }
    /* Boundary segments*/
    if (fscanf(fin, "\n%[^\n]",idumy) == EOF)
        goto close_file;
    int idumy1;
    if (fscanf(fin, "%d %d %d", &idumy1, &nloops_G,&ndomains_G) == EOF)
        goto close_file;
    //"Segments in curves" label
    if (fscanf(fin, "\n%[^\n]",idumy) == EOF)
        goto close_file;
    for (int j = 0; j < ncurves_G; j++)
    {
        if (fscanf(fin, "\n%[^\n]",idumy) == EOF)
            goto close_file;
    }
    //"Regions on curves" label
    if (fscanf(fin, "\n%[^\n]",&idumy) == EOF)
        goto close_file;
    gm3L.resize(nloops_G);
    lpCvNums.resize(fcNums);
    for(int i=0;i<nloops_G;i++) {
        int idx1,idx2,idumy3;
        if (fscanf(fin, "%d %d %d", &idx1, &idx2, &idumy3) == EOF)
            goto close_file;
        gm3L[i].idx=idx1;
        gm3L[i].sfidx=idx2;
        int nlc;
        if (fscanf(fin, "%d", &nlc) == EOF)
            goto close_file;
        gm3L[i].nlc=nlc;
        gm3L[i].lc.resize(nlc);
        lpCvNums[i]=nlc;
        for(int j=0;j<nlc;j++) {
            int tmp;
            fscanf(fin, "%d", &tmp);
            gm3L[i].lc[j]=tmp;
        }
    }
    lpCvNumsAll=0;
    for(int i=0;i<fcNums;i++) {
        lpCvNumsAll+=lpCvNums[i];
    }
    lpCvs.resize(lpCvNumsAll);
    for(int i=0;i<fcNums;i++) {
        for(int j=0;j<lpCvNums[i];j++) {
            lpCvs[lpCvsIdx_local++]=gm3L[i].lc[j];
        }
    }
close_file:
    fclose(fin);
    return true;
}
//int bndPtNum; double* bndPts;
//int bndFctNum; int* bndFcts;
bool InputData::readStlPls() {
    FILE *fin = NULL;
	//cout <<"stl:" << stlplsfile << endl;
    fin = fopen(stlplsfile.c_str(), "r");
    if (fin == NULL)
    {
        spdlog::info("Error: cann't open file %s\n", stlplsfile.c_str());
        return 0;
    }
    char strn[128];
    int ntri,npt;
    fscanf(fin, "%d %d%[^\n]", &ntri, &npt, &strn);
    stlplsE.resize(ntri);
    stlplsV.resize(npt);
    bndPtNum=npt;
    bndFctNum=ntri;
    bndPts.resize(npt*3);
    bndFcts.resize(ntri*3);
    fctFaces.resize(ntri);
    for(int i=0;i<npt;i++) {
        int idx;
        fscanf(fin, "%d %e %e %e",&idx, &stlplsV[i].x, &stlplsV[i].y, &stlplsV[i].z);
        bndPts[i*3]=stlplsV[i].x;
        bndPts[i*3+1]=stlplsV[i].y;
        bndPts[i*3+2]=stlplsV[i].z;
    }
    for(int i=0;i<ntri;i++) {
        int idx,pt0,pt1,pt2,pt3;
        fscanf(fin, "%d %d %d %d %d%[^\n]", &idx, &pt0, &pt1, &pt2, &pt3, &strn);
        stlplsE[i].idx=idx;
        stlplsE[i].a=pt0;
        stlplsE[i].b=pt1;
        stlplsE[i].c=pt2;
        stlplsE[i].sfidx=pt3;
        bndFcts[i*3]=pt0;
        bndFcts[i*3+1]=pt1;
        bndFcts[i*3+2]=pt2;
        fctFaces[i]=pt3;
    }
    fclose(fin);
    fin = NULL;
    return 0;
}
bool InputData::readTF() {
    tf.pntSetNums=spSettingNum[0];
    tf.cvSetNums=spSettingNum[1];
    tf.fcSetNums=spSettingNum[2];

	tf.GlSettings.hmax=spGlSettings[0];
	tf.GlSettings.hmin=spGlSettings[1];
	tf.GlSettings.beta=spGlSettings[2];
	tf.GlSettings.theta=spGlSettings[3];
	tf.GlSettings.nsegPt=spGlSettings[4];
	tf.GlSettings.nsegCv=spGlSettings[5];
	tf.GlSettings.nsegFc=spGlSettings[6];

    for(int i=0;i<tf.pntSetNums;i++) {
		TigerConfig::pointSettings tmp;
		tmp.id=spPtSettings[i*7+0];
		tmp.hmax=spPtSettings[i*7+1];
		tmp.hmin=spPtSettings[i*7+2];
		tmp.theta=spPtSettings[i*7+3];
		tmp.nsegPt=spPtSettings[i*7+4];
		tmp.nsegCv=spPtSettings[i*7+5];
		tmp.nsegFc=spPtSettings[i*7+6];
		tf.ptSettings.push_back(tmp);
        //tf.ptSettings.push_back({spPtSettings[0],spPtSettings[1],spPtSettings[2],spPtSettings[3],spPtSettings[4],spPtSettings[5],spPtSettings[6]});
    }
    for(int i=0;i<tf.cvSetNums;i++) {
		TigerConfig::curveSettings tmp={spCvSettings[i*6+0],spCvSettings[i*6+1],spCvSettings[i*6+2],spCvSettings[i*6+3],spCvSettings[i*6+4],spCvSettings[i*6+5]};
		tf.cvSettings.push_back(tmp);
		
		//tf.cvSettings.push_back({spCvSettings[0],spCvSettings[1],spCvSettings[2],spCvSettings[3],spCvSettings[4],spCvSettings[5]});
    }
    for(int i=0;i<tf.fcSetNums;i++) {
		TigerConfig::faceSettings tmp={spFcSettings[i*6+0],spFcSettings[i*6+1],spFcSettings[i*6+2],spFcSettings[i*6+3],spFcSettings[i*6+4],spFcSettings[i*6+5]};
		tf.fcSettings.push_back(tmp);
		//tf.fcSettings.push_back({spFcSettings[0],spFcSettings[1],spFcSettings[2],spFcSettings[3],spFcSettings[4],spFcSettings[5]});
    }
    tf.volumeGrowthRatio=volmGrowthRatio;
    return true;
}
bool InputData::readTranConfig() {
    char word[100];      /* variable name from config file */
    int value = 0;         /* variable value from config file */
    double fvalue = 0.0; /* variable for floats in config file */
    char str[100];
    char fstr[100];
    char temp[100];

    /* temp stuff for line reading */
    int nbytes = 1024;
    char line[1024];
    int numassigned;
	//cout << "config:" << Tran_ConfigFile
	//	<< endl;
    FILE *fp = fopen(Tran_ConfigFile.c_str(), "r");
    if (!fp)
    {
        spdlog::info("Config File does not exist.\n");
        exit(-1);
    }
	tf.pntSetNums=0;
	tf.cvSetNums=0;
	tf.fcSetNums=0;
    while (fgets(line, nbytes, fp) != NULL)
    {
        /* attempt to fetch a variable name and value from the config file */
        numassigned = sscanf(line, "%s %s", word, temp);
        fvalue = atof(fstr);
        /* check if this is a comment */
        if (word[0] == '#' || word[0] == '\n' || numassigned < 2) continue;

        if (!strcmp(word, "filename")) {
            sscanf(line, "%s %s", word, str);
            strcpy(cf.filenam, str);
            Gm3AndStlFileName=string(str);
        }
        if (!strcmp(word, "step")) {
            sscanf(line, "%s %d", word, &value);
            cf.step = value;
        }
        // face
        if (!strcmp(word, "face")) {
            int iface;
            sscanf(line, "%s %d %lf", word, &iface, &fvalue);
            cf.fcsize.insert(std::make_pair(iface-1, fvalue));
            tf.fcSetNums++;
            TigerConfig::faceSettings t={double(iface),fvalue,0,0,0,0};
            tf.fcSettings.push_back(t);
        }
        // point
        if (!strcmp(word, "point")) {
            int ipt;
            sscanf(line, "%s %d %lf", word, &ipt, &fvalue);
            cf.ptsize.insert(std::make_pair(ipt-1, fvalue));
            tf.pntSetNums++;
            TigerConfig::pointSettings t={double(ipt),fvalue,0,0,0,0,0};
            tf.ptSettings.push_back(t);
        }
        // curve
        if (!strcmp(word, "curve")) {
            int icv;
            sscanf(line, "%s %d %lf", word, &icv, &fvalue);
            cf.lnsize.insert(std::make_pair(icv-1, fvalue));
            tf.cvSetNums++;
            TigerConfig::curveSettings t={double(icv),fvalue,0,0,0,0};
            tf.cvSettings.push_back(t);
        }
        // global hmax
        if (!strcmp(word, "global_spacing")) {
            sscanf(line, "%s %lf", word, &fvalue);
            cf.global_space = fvalue;
            tf.GlSettings.hmax=fvalue;
        }
        // global hmin
        if (!strcmp(word, "min_spacing")) {
            sscanf(line, "%s %lf", word, &fvalue);
            cf.min_space = fvalue;
            tf.GlSettings.hmin=fvalue;
        }
        // global theta
        if (!strcmp(word, "curvature_angle")) {
            sscanf(line, "%s %lf", word, &fvalue);
            cf.curvature_angle = fvalue;
            tf.GlSettings.theta=fvalue;
        }
        // global beta
        if (!strcmp(word, "expand_ratio")) {
            sscanf(line, "%s %lf", word, &fvalue);
            cf.expand_ratio = fvalue;
            tf.GlSettings.beta=fvalue;
        }
        // global nsegPt,Cv,Fc ?
        if (!strcmp(word, "proximity_num")) {
            sscanf(line, "%s %d", word, &value);
            cf.proximity_num = value;
            tf.GlSettings.nsegCv=value;
            tf.GlSettings.nsegFc=value;
            tf.GlSettings.nsegPt=value;
        }
    }
    spSettingNum.push_back(tf.pntSetNums);
    spSettingNum.push_back(tf.cvSetNums);
    spSettingNum.push_back(tf.fcSetNums);
    spGlSettings.resize(7);
    spGlSettings[0]=tf.GlSettings.hmax;
    spGlSettings[1]=tf.GlSettings.hmin;
    spGlSettings[2]=tf.GlSettings.beta;
    spGlSettings[3]=tf.GlSettings.theta;
    spGlSettings[4]=tf.GlSettings.nsegPt;
    spGlSettings[5]=tf.GlSettings.nsegCv;
    spGlSettings[6]=tf.GlSettings.nsegFc;
    spPtSettings.resize(tf.pntSetNums*7);
    for(int i=0;i<tf.pntSetNums;i++) {
        spPtSettings[i*7+0]=tf.ptSettings[i].id;
        spPtSettings[i*7+1]=tf.ptSettings[i].hmax;
        for(int j=2;j<7;j++)spPtSettings[i*7+j]=0;
//        spPtSettings[i*7+2]=tf.ptSettings[i].hmin;
//        spPtSettings[i*7+3]=tf.ptSettings[i].theta;
//        spPtSettings[i*7+4]=tf.ptSettings[i].nsegPt;
//        spPtSettings[i*7+5]=tf.ptSettings[i].nsegCv;
//        spPtSettings[i*7+6]=tf.ptSettings[i].nsegFc;
    }
    spCvSettings.resize(tf.cvSetNums*6);
    for(int i=0;i<tf.cvSetNums;i++) {
        spCvSettings[i*6+0]=tf.cvSettings[i].id;
        spCvSettings[i*6+1]=tf.cvSettings[i].hmax;
        for(int j=2;j<6;j++)spCvSettings[i*6+j]=0;
    }
    spFcSettings.resize(tf.fcSetNums*6);
    for(int i=0;i<tf.fcSetNums;i++) {
        spFcSettings[i*6+0]=tf.fcSettings[i].id;
        spFcSettings[i*6+1]=tf.fcSettings[i].hmax;
        for(int j=2;j<6;j++)spFcSettings[i*6+j]=0;
    }
    volmGrowthRatio=0;
    fclose(fp);
    return 0;

}
bool InputData::writeGm3String(string str) {
	return true;
}
bool InputData::writeGm3(string fileName) {
    map<vertex,int> mpv;
    int cvPtsIdx=0;
    int pntIdx=0;
    TigerC.resize(cvNums);
    for(int i=0;i<cvNums;i++) {
        TigerC[i].idx=i;
        TigerC[i].ncp=cvPtNums[i];
        for(int j=0;j<cvPtNums[i];j++) {
            vertex tmp;
            tmp.x=cvPts[cvPtsIdx*3];
            tmp.y=cvPts[cvPtsIdx*3+1];
            tmp.z=cvPts[cvPtsIdx*3+2];
            cvPtsIdx++;
            if(mpv.find(tmp)==mpv.end()) {
                tmp.idx=pntIdx;
                TigerC[i].v.push_back(pntIdx);
                mpv[tmp]=pntIdx++;
                TigerV.push_back(tmp);

            }else {
                TigerC[i].v.push_back(mpv[tmp]);
            }
        }
    }
    TigerS.resize(fcNums);
    int fcPtsIdx=0;
    for(int i=0;i<fcNums;i++) {
        TigerS[i].nsp1=fcPtNums[i*2];
        TigerS[i].nsp2=fcPtNums[i*2+1];
        TigerS[i].v.resize(TigerS[i].nsp2);
        for(int j=0;j<fcPtNums[i*2+1];j++) {
            for(int k=0;k<fcPtNums[i*2];k++) {
                vertex tmp;
                tmp.x=fcPts[fcPtsIdx*3];
                tmp.y=fcPts[fcPtsIdx*3+1];
                tmp.z=fcPts[fcPtsIdx*3+2];
                fcPtsIdx++;
                tmp.idx=pntIdx;
                TigerS[i].v[j].push_back(pntIdx);
                pntIdx++;
                TigerV.push_back(tmp);
            }
        }
    }
    ofstream outfile;
    outfile.open(fileName);
    outfile <<" 1.- liwgeom Definition"<<endl;
    outfile <<"\t"<<TigerV.size()<<"\t"<<cvNums<<"\t"<<fcNums<<endl;
    outfile <<"  Points"<<endl;
    for(int i=0;i<TigerV.size();i++) {
        outfile <<scientific<<"\t"<<i+1<<"\t"<<TigerV[i].x<<"\t"<<TigerV[i].y<<"\t"<<TigerV[i].z<<"\n";
    }
    outfile <<"  Curves"<<endl;
    for(int i=0;i<cvNums;i++) {
        outfile <<"\t"<<i+1<<'\t'<<1<<endl
               <<"\t"<<cvPtNums[i]<<endl;
        for(int j=0;j<TigerC[i].ncp;j++) {
            outfile<<scientific<<"\t"<<TigerV[TigerC[i].v[j]].x<<"\t"<<TigerV[TigerC[i].v[j]].y<<"\t"<<TigerV[TigerC[i].v[j]].z<<"\t"<<TigerC[i].v[j]+1<<endl;
        }
    }
    outfile <<"  Support_Surfaces"<<endl;
    for(int i=0;i<fcNums;i++) {
        outfile <<"\t"<<i+1<<"\t"<<1<<endl
               <<"\t"<<TigerS[i].nsp1<<"\t"<<TigerS[i].nsp2<<endl;
        for(int j=0;j<TigerS[i].nsp2;j++) {
            for(int k=0;k<TigerS[i].nsp1;k++)
                outfile <<scientific<<"\t"<<TigerV[TigerS[i].v[j][k]].x<<"\t"<<TigerV[TigerS[i].v[j][k]].y<<"\t"<<TigerV[TigerS[i].v[j][k]].z<<"\t"<<TigerS[i].v[j][k]+1<<endl;
        }

    }
    outfile <<"  2.- Mesh Generation"<<endl;
    outfile <<"\t"<<cvNums<<"\t"<<fcNums<<"\t"<<0<<endl;
    outfile <<" Segments in curves"<<endl;
    for(int i=0;i<cvNums;i++) {
        outfile<< "\t"<<i+1<<"\t"<<i+1<<"\t"<<1<<endl;
    }
    outfile<<" Regions on curves"<<endl;
    int lpCvIdx=0;
    for(int i=0;i<fcNums;i++) {
        outfile<< "\t"<<i+1<<"\t"<<i+1<<"\t"<<1<<endl;
        outfile<<"\t"<<lpCvNums[i]<<endl;
        for(int j=0;j<lpCvNums[i];j++) {
            outfile <<"\t"<<lpCvs[lpCvIdx++];
        }
        outfile<<endl;
    }
    outfile.close();
    return true;
}

bool InputData::writeTigerGm3(string fileName) {
    ofstream outfile;
    outfile.open(fileName);
    outfile << "cvNums:"<<endl;
    outfile << cvNums<<endl;
    outfile << "cvPtNums:"<<endl;
    for(int i=0;i<cvNums;i++) {
        outfile <<cvPtNums[i]<<" ";
    }
    outfile <<endl;
    int cvPtsIdx=0;
    for(int i=0;i<cvNums;i++) {
        outfile <<"cv idx:"<<i<<" ";
        for(int j=0;j<cvPtNums[i];j++) {
            outfile <<scientific<<cvPts[cvPtsIdx*3]<<" "<<cvPts[cvPtsIdx*3+1]<<" "<<cvPts[cvPtsIdx*3+2]<<endl;
            cvPtsIdx++;
        }
    }
    outfile << "fcNums:"<<endl;
    outfile << fcNums<<endl;
    outfile << "fcPtNums:"<<endl;
    for(int i=0;i<fcNums;i++) {
        outfile <<fcPtNums[i*2]<<" "<<fcPtNums[i*2+1]<<" / ";
    }
    outfile << "fcPts:"<<endl;
    int fcPtsIdx=0;
    for(int i=0;i<fcNums;i++) {
        outfile <<"fc idx:"<<i<<endl;
        for(int j=0;j<fcPtNums[i*2]*fcPtNums[i*2+1];j++) {
            outfile <<fcPts[fcPtsIdx*3]<<" "<<fcPts[fcPtsIdx*3+1]<<" "<<fcPts[fcPtsIdx*3+2]<<" * ";
            fcPtsIdx++;
        }
    }
    outfile <<"lpCvNums:"<<endl;
    for(int i=0;i<fcNums;i++) {
        outfile<<lpCvNums[i]<<" ";
    }
    outfile <<"lpCvs:"<<endl;
    int lpIdx=0;
    for(int i=0;i<fcNums;i++) {
        outfile <<"lp idx:"<<i<<endl;
        for(int j=0;j<lpCvNums[i];j++) {
            outfile <<lpCvs[lpIdx]<<" ";
            lpIdx++;
        }
        outfile<<"\n";
    }
    outfile.close();
    return true;
}

bool InputData::writeTigerStlPls(string fileName) {
    ofstream outfile;
    outfile.open(fileName);
    outfile<<"bndPtNum:"<<bndPtNum<<endl;
    for(int i=0;i<bndPtNum;i++) {
        outfile<<"idx:"<<i+1<<"\t"<<bndPts[i*3]<<"\t"<<bndPts[i*3+1]<<"\t"<<bndPts[i*3+2]<<endl;
    }
    outfile <<"bndFctNum:"<<bndFctNum<<endl;
    for(int i=0;i<bndFctNum;i++) {
        outfile<<"idx:"<<i+1<<"\t"<<bndFcts[i*3]<<"\t"<<bndFcts[i*3+1]<<"\t"<<bndFcts[i*3+2]<<endl;
    }
    outfile.close();
    return true;
}
//HJD 缺少面索引
bool InputData::writeStlPls(string fileName) {
    ofstream outfile;
    outfile.open(fileName);
    outfile<<"\t"<<bndFctNum<<"\t"<<bndPtNum<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
    for(int i=0;i<bndPtNum;i++) {
        outfile <<"\t"<<i+1<<scientific<<"\t"<<bndPts[i*3]<<"\t"<<bndPts[i*3+1]<<"\t"<<bndPts[i*3+2]<<endl;
    }
    for(int i=0;i<bndFctNum;i++) {
        outfile <<"\t"<<i+1<<scientific<<"\t"<<bndFcts[i*3]<<"\t"<<bndFcts[i*3+1]<<"\t"<<bndFcts[i*3+2]<<"\t"<<fctFaces[i]<<endl;
    }
    outfile.close();
    return true;
}
bool InputData::writeTranConfigFromCF(string fileName) {
    ofstream outfile;
    outfile.open(fileName);
    outfile <<"#parameters for adasize"<<endl;
    outfile <<"filename "<<cf.filenam<<endl;
    if(cf.step!=-1)
        outfile <<"step "<<cf.step<<endl;
    for(auto& x:cf.fcsize) {
        outfile <<"face "<<x.first+1<<" "<<x.second<<endl;
    }
    for(auto& x:cf.ptsize) {
        outfile <<"point "<<x.first+1<<" "<<x.second<<endl;
    }
    for(auto& x:cf.lnsize) {
        outfile <<"curve "<<x.first+1<<" "<<x.second<<endl;
    }
    if(cf.global_space!=-1)
        outfile <<"global_spacing "<<cf.global_space<<endl;
    if(cf.min_space!=-1)
        outfile <<"min_spacing "<<cf.min_space<<endl;
    if(cf.curvature_angle!=-1)
        outfile <<"curvature_angle "<<cf.curvature_angle<<endl;
    if(cf.expand_ratio!=-1)
        outfile <<"expand_ratio "<<cf.expand_ratio<<endl;
    if(cf.proximity_num!=-1)
        outfile <<"proximity_num "<<cf.proximity_num<<endl;
    outfile.close();
    return true;
}
//int* spSettingNum;
//double* spGlSettings; double* spPtSettings;
//double* spCvSettings; double* spFcSettings;
//double volmGrowthRatio;
/**
 * @brief InputData::writeTranConfig
 * @param fileName
 * @return
 */
bool InputData::writeTranConfig(string fileName) {
//    char filenam[1024];		// the project name
//    int step;
//    // hmin
//    double min_space;
//    // hmax
//    double global_space;
//    // theta
//    double curvature_angle;
//    // beta
//    double expand_ratio;
//    // 间隔 nsegPt,nsegCv,nsegFc
//    int proximity_num;
//    // 点的hmax
//    std::map<int, double> ptsize;
//    // 点的hmax
//    std::map<int, double> lnsize;
//    // 点的hmax
//    std::map<int, double> fcsize;
    ConfigArgc tigerCF;
    strcpy(tigerCF.filenam,Gm3AndStlFileName.c_str());
    tigerCF.global_space=tf.GlSettings.hmax;
    tigerCF.min_space=tf.GlSettings.hmin;
    tigerCF.expand_ratio=tf.GlSettings.beta;
    // 2020/06/14 HJD 与comsol 对标测试

    tigerCF.curvature_angle=tf.GlSettings.theta;
    // 目前GlSettings.的三个seg为同样的值都为proximity_num;
    tigerCF.proximity_num=tf.GlSettings.nsegPt;
    for(auto& x:tf.ptSettings) {
        tigerCF.ptsize[x.id]=x.hmax;
    }
    for(auto& x:tf.cvSettings) {
        tigerCF.lnsize[x.id]=x.hmax;
    }
    for(auto& x:tf.fcSettings) {
        tigerCF.fcsize[x.id]=x.hmax;
    }
    ofstream outfile;
    outfile.open(fileName);
    outfile <<"#parameters for adasize"<<endl;
    outfile <<"filename "<<tigerCF.filenam<<endl;
    if(tigerCF.step!=-1)
        outfile <<"step "<<tigerCF.step<<endl;
    for(auto& x:tigerCF.fcsize) {
        outfile <<"face "<<x.first<<" "<<x.second<<endl;
    }
    for(auto& x:tigerCF.ptsize) {
        outfile <<"point "<<x.first<<" "<<x.second<<endl;
    }
    for(auto& x:tigerCF.lnsize) {
        outfile <<"curve "<<x.first<<" "<<x.second<<endl;
    }
    if(tigerCF.global_space!=-1)
        outfile <<"global_spacing "<<tigerCF.global_space<<endl;
    if(tigerCF.min_space!=-1)
        outfile <<"min_spacing "<<tigerCF.min_space<<endl;
    if(tigerCF.curvature_angle!=-1)
        outfile <<"curvature_angle "<<tigerCF.curvature_angle<<endl;
    if(tigerCF.expand_ratio!=-1)
        outfile <<"expand_ratio "<<tigerCF.expand_ratio<<endl;
    if(tigerCF.proximity_num!=-1)
        outfile <<"proximity_num "<<tigerCF.proximity_num<<endl;
//    outfile<<"method pardiso"<<endl;
//        outfile<<"method ma86"<<endl;
    outfile<<"step 500"<<endl;
    outfile.close();
    return true;
}
string InputData::getFileName(const string& s) {

	char sep = '/';

#ifdef _WIN32
	sep = '\\';
#endif

	size_t i = s.rfind(sep, s.length());
	if (i != string::npos) {
		return(s.substr(i + 1, s.length() - i));
	}

	return("");
}
string InputData::getPathName(const string& s) {

	char sep = '/';

#ifdef _WIN32
	sep = '\\';
#endif

	size_t i = s.rfind(sep, s.length());
	if (i != string::npos) {
		return(s.substr(0, i));
	}

	return("");
}
