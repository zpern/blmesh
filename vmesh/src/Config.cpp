#include <spdlog/spdlog.h> 
 #include "Config.h"

using namespace std;
using namespace tiger;

Config::Config(void)
{
    SetDefaultConfig();
}

Config::~Config(void)
{

}

int Config::SetDefaultConfig()
{
    /* basic members */
    sConfigName = "TRAN_CONFIG";  // config name
    iWrtLvl = 0; // 0: all menbers; 1: geometry members; 2: mesh members
    dElemNum = 500000;

    /* geometry members */
    sKwisotropic_stop="isotropic_stop";
    sCmtisotropic_stop="#isotropic_stop";
    isotropic_stop=0;


    sGeoFileName = "";	// the geometry file name
    sKwGeoFileName = "filename";
    sCmtGeoFileName= "#filename,as it is"; // the comment

    iGeoType = 1;
    sKwGeoType = "geometry_type";
    sCmtGeoType = "# liwgeom type. 1: STEP";

    dSTolerance = 1.0e-6;				//set sTolerance, default value 1.0e-8
    sKwSTolerance = "sTolerance";
    sCmtSTolerance = "#sTolerance: Default is 1.0e-6.";

    dGTolerance = 1.0e-6;				//set gTolerance, default value 1.0e-8
    sKwGTolerance = "gTolerance";
    sCmtGTolerance = "#gTolerance: Default is 1.0e-6.";

    bDelSmallEdges = false;			//set del_small_edges, default value false
    sKwDelSmallEdges = "del_small_edges";
    sCmtDelSmallEdges = "#";

    /* parameters for mesh size */
    dMinSize = 50;				//the minimize mesh size
    sKwMinSize = "min_spacing";
    sCmtMinSize = "#min_spacing";

    dGlblSize = 5;				//the global mesh size
    sKwGlblSize = "global_spacing";
    sCmtGlblSize = "#global_spacing";

    dExpandRatio = 1.2;
    sKwExpandRatio = "expand_ratio";
    sCmtExpandRatio = "#expand_ratio";

    mapPtSize.clear();
    sKwPtSize = "point";
    sCmtPtSize = "#point size";

    mapCvSize.clear();
    sKwCvSize = "curve";
    sCmtCvSize = "#curve size";

    mapFcSize.clear();
    sKwFcSize = "face";
    sCmtFcSize = "#face size";

    /* parameters for adaptive */
    iAdaptive = 1;
    sKwAdaptive = "adaptive";
    sCmtAdaptive = "#set adaptive";

    iProximityNum = 1;
    sKwProximityNum = "proximity_num";
    sCmtProximityNum = "#proximity_num";

    dCurvatureAngle = 10.0;
    sKwCurvatureAngle = "curvature_angle";
    sCmtCurvatureAngle = "#curvature_angle";

    mapFcBku.clear();
    sKwFcBku = "bkfaceu";
    sCmtFcBku = "#bkfaceu";

    mapFcBkv.clear();
    sKwFcBkv = "bkfacev";
    sCmtFcBkv = "#bkfacev";

    iuNum = 4;
    sKwuNum = "unum";
    sCmtuNum = "#unum";

    ivNum = 4;
    sKwvNum = "vnum";
    sCmtvNum = "#vnum";

    /* layer mesh */
    nLayerNum = 0;
    sKwLayerNum = "layer_num";
    sCmtLayerNum ="#layer_num";

    vecBoxFc.clear();
    sKwBoxFc = "box_face";
    sCmtBoxFc = "#box_face";

	vecSymmFc.clear();
	sKwSymmFc = "sym_face";
	sCmtBoxFc = "#sym_face";

    dStepLen = 0.01;
    sKwStepLen = "step_len";
    sCmtStepLen = "#initial step length";

    dStepLenRatio = 1.35;
    sKwStepLenRatio = "ratio";
    sCmtStepLenRatio = "#ratio";



    dStepLenInitRatio = 1.15;
    //sKwStepLenInitRatio;
    //sCmtStepLenInitRatio;

    nInitLayerNum = 5;
    //sKwInitLayerNum;
    //sCmtInitLayerNum;

    /* mesh improvement */
    nSmoothLoop = 10;
    sKwSmoothLoop = "smooth_attempt";
    sCmtSmoothLoop = "#smooth_attempt";

    return 0;
}

int Config::ReadConfigFile(const char* filename)
{
    char word[100];      /* variable name from config file */
    int value=0;         /* variable value from config file */
    double fvalue, fvalue2; /* variable for floats in config file */
    char str[100];
    char fstr[100];
    char temp[100];

    /* temp stuff for line reading */
    int nbytes = 1000;
    char line[1000]="";
    int numassigned;

    int ipoint, icurve, iface;

    FILE *fp = fopen(filename,"r");
    if(!fp)
    {
        spdlog::info("Config File does not exist.\n");
        //exit(ERRORCONFIG);
        exit(-1);
    }

    while (fgets(line, nbytes, fp) != NULL)
    {
        /* attempt to fetch a variable name and value from the config file */
        numassigned = sscanf(line, "%s %s", word, temp);
        fvalue = atof(fstr);
        /* check if this is a comment */
        if (word[0] == '#' || word[0] == '\n' || numassigned < 2) continue;

        /* geometry members */
        if (!strcmp(word, sKwGeoType.c_str())) {
            sscanf(line, "%s %d", word, &value);
            iGeoType = value;
        }
        else if (!strcmp(word, sKwGeoFileName.c_str())) {
            sscanf(line, "%s %s", word, str);
            sGeoFileName = string(str);
        }
        else if (!strcmp(word, sKwSTolerance.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dSTolerance = fvalue;
        }
        else if (!strcmp(word, sKwGTolerance.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dGTolerance = fvalue;
        }
        else if (!strcmp(word, sKwDelSmallEdges.c_str())) {
            sscanf(line, "%s %d", word, &value);
            bDelSmallEdges = bool(value);
        }
        /* parameters for mesh size */
        else if (!strcmp(word, sKwGlblSize.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dGlblSize = fvalue;
        }
        else if (!strcmp(word, sKwMinSize.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dMinSize = fvalue;
        }
        else if (!strcmp(word, sKwExpandRatio.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dExpandRatio = fvalue;
        }
        else if (!strcmp(word, sKwPtSize.c_str())) {
            sscanf(line, "%s %d %lf", word, &ipoint, &fvalue);
            mapPtSize.insert(std::make_pair(ipoint, fvalue));
        }
        else if (!strcmp(word, sKwCvSize.c_str())) {
            sscanf(line, "%s %d %lf", word, &icurve, &fvalue);
            mapCvSize.insert(std::make_pair(icurve, fvalue));
        }
        else if (!strcmp(word, sKwFcSize.c_str())) {
            sscanf(line, "%s %d %lf", word, &iface, &fvalue);
            mapFcSize.insert(std::make_pair(iface, fvalue));
        }
        /* parameters for adaptive */
        else if (!strcmp(word, sKwAdaptive.c_str())) {
            sscanf(line, "%s %d", word, &value);
            iAdaptive = value;
        }
        else if (!strcmp(word, sKwProximityNum.c_str())) {
            sscanf(line, "%s %d", word, &value);
            iProximityNum = value;
        }
        else if (!strcmp(word, sKwisotropic_stop.c_str())) {
            sscanf(line, "%s %d", word, &value);
            isotropic_stop = value;
        }


        else if (!strcmp(word, sKwCurvatureAngle.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dCurvatureAngle = fvalue;
        }
        else if (!strcmp(word, sKwuNum.c_str())) {
            sscanf(line, "%s %d", word, &value);
            iuNum = value;
        }
        else if (!strcmp(word, sKwvNum.c_str())) {
            sscanf(line, "%s %d", word, &value);
            ivNum = value;
        }
        else if (!strcmp(word, sKwFcBku.c_str())) {
            sscanf(line, "%s %d %d", word, &iface, &value);
            mapFcBku.insert(std::make_pair(iface-1, value));
        }
        else if (!strcmp(word, sKwFcBkv.c_str())) {
            sscanf(line, "%s %d %d", word, &iface, &value);
            mapFcBkv.insert(std::make_pair(iface-1, value));
        }
        else if (!strcmp(word, sKwLayerNum.c_str())) {
            sscanf(line, "%s %d", word, &value);
            nLayerNum = value;
        }
        else if (!strcmp(word, sKwBoxFc.c_str())) {
            sscanf(line, "%s %d", word, &value);
            vecBoxFc.push_back(value);
        }
        else if (!strcmp(word, sKwStepLen.c_str())) {
            sscanf(line, "%s %lf", word, &fvalue);
            dStepLen = fvalue;
        }
        else if (!strcmp(word, sKwStepLenRatio.c_str())) {
            sscanf(line, "%s %d %lf %lf", word, &value, &fvalue, &fvalue2);
            nInitLayerNum = value;
            dStepLenInitRatio = fvalue;
            dStepLenRatio = fvalue2;
        }
        else if (!strcmp(word, sKwSmoothLoop.c_str())) {
            sscanf(line, "%s %d", word, &value);
            nSmoothLoop = value;
        }

        //else if (!strcmp(word, "filefin")) {
        //	sscanf(line, "%s %s", word, str);
        //	strcpy(cf.filfin, str);
        //}
    }
    dGTolerance = dSTolerance;
    if (fp)
        fclose(fp);
    return 0;//INPUTCONFIG;
}

/* 根据iLv输出，可选取部分参数输出（未实现） */
int Config::WriteConfigFile(const int iLv,std::string filename)
{
    FILE *fp;

    if(filename.empty()){
        fp=fopen(sConfigName.c_str(),"w");
    }
    else {
        fp=fopen(filename.c_str(),"w");
    }
    if(!fp)
    {
        fprintf(stderr,"Can not write Config File.\n");
        exit(-1);
    }

    fprintf(fp, "#parameters for tran2fli adasize blmesh\n\n");

    WriteConfigGeomPart(fp);

    WriteConfigMeshSizePart(fp);

    WriteConfigAdapPart(fp);

    WriteConfigLayerMeshPart(fp);

    WriteConfigMeshOptPart(fp);

    if(fp)
        fclose(fp);

    return 0;
}

int Config::WriteConfigGeomPart(FILE *fp)
{
    if(!fp)
        return 1;

    fprintf(fp, "%s\n", sCmtGeoType.c_str());
    fprintf(fp, "%s %d\n\n", sKwGeoType.c_str(), iGeoType);

    fprintf(fp, "%s\n", sCmtGeoFileName.c_str());
    fprintf(fp, "%s %s\n\n", sKwGeoFileName.c_str(), sGeoFileName.c_str());

    fprintf(fp, "%s\n", sCmtSTolerance.c_str());
    fprintf(fp, "%s %e\n\n", sKwSTolerance.c_str(), dSTolerance);

    fprintf(fp, "\n");

    return 0;
}

int Config::WriteConfigMeshSizePart(FILE *fp)
{
    map<int,double>::iterator itMap;
    if(!fp)
        return 1;

    fprintf(fp, "%s\n", sCmtGlblSize.c_str());
    fprintf(fp, "%s %lf\n\n", sKwGlblSize.c_str(), dGlblSize);

    fprintf(fp, "%s\n", sCmtMinSize.c_str());
    fprintf(fp, "%s %lf\n\n", sKwMinSize.c_str(), dMinSize);

    fprintf(fp, "%s\n", sCmtExpandRatio.c_str());
    fprintf(fp, "%s %lf\n\n", sKwExpandRatio.c_str(), dExpandRatio);

    if (mapPtSize.size() > 0)
    {
        fprintf(fp, "%s\n", sCmtPtSize.c_str());
        for (itMap = mapPtSize.begin(); itMap != mapPtSize.end(); itMap++)
        {
            if((itMap->second-dGlblSize)>0)
            fprintf(fp, "%s %d %lf\n", sKwPtSize.c_str(), itMap->first, itMap->second);
        }
        fprintf(fp, "\n");
    }

    if (mapCvSize.size() > 0)
    {
        fprintf(fp, "%s\n", sCmtCvSize.c_str());
        for (itMap = mapCvSize.begin(); itMap != mapCvSize.end(); itMap++)
        {
            if((itMap->second-dGlblSize)>0)
            fprintf(fp, "%s %d %lf\n", sKwCvSize.c_str(), itMap->first, itMap->second);
        }
        fprintf(fp, "\n");
    }

    if (mapFcSize.size() > 0)
    {
        fprintf(fp, "%s\n", sCmtFcSize.c_str());
        for (itMap = mapFcSize.begin(); itMap != mapFcSize.end(); itMap++)
        {
            if((itMap->second-dGlblSize)>0)
            fprintf(fp, "%s %d %lf\n", sKwFcSize.c_str(), itMap->first, itMap->second);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");

    return 0;
}

int Config::WriteConfigAdapPart(FILE *fp)
{
    map<int,int>::iterator itMap;

    if(!fp)
        return 1;

    fprintf(fp, "%s\n", sCmtAdaptive.c_str());
    fprintf(fp, "%s %d\n\n", sKwAdaptive.c_str(), iAdaptive);

    fprintf(fp, "%s\n", sCmtProximityNum.c_str());
    fprintf(fp, "%s %d\n\n", sKwProximityNum.c_str(), iProximityNum);

    fprintf(fp, "%s\n", sCmtCurvatureAngle.c_str());
    fprintf(fp, "%s %lf\n\n", sKwCurvatureAngle.c_str(), dCurvatureAngle);

    if (iuNum != 4)
    {
        fprintf(fp, "%s\n", sCmtuNum.c_str());
        fprintf(fp, "%s %d\n\n", sKwuNum.c_str(), iuNum);
    }

    if (ivNum != 4)
    {
        fprintf(fp, "%s\n", sCmtvNum.c_str());
        fprintf(fp, "%s %d\n\n", sKwvNum.c_str(), ivNum);
    }

    if (mapFcBku.size() > 0)
    {
        fprintf(fp, "%s\n", sCmtFcBku.c_str());
        for (itMap = mapFcBku.begin(); itMap != mapFcBku.end(); itMap++)
        {
            fprintf(fp, "%s %d %d\n", sKwFcBku.c_str(), itMap->first+1, itMap->second);
        }
    }

    if (mapFcBkv.size() > 0)
    {
        fprintf(fp, "%s\n", sCmtFcBkv.c_str());
        for (itMap = mapFcBkv.begin(); itMap != mapFcBkv.end(); itMap++)
        {
            fprintf(fp, "%s %d %d\n", sKwFcBkv.c_str(), itMap->first+1, itMap->second);
        }
    }

    fprintf(fp, "\n");
    return 0;
}

int Config::WriteConfigLayerMeshPart(FILE *fp)
{
    if(!fp)
        return 1;


    fprintf(fp, "%s\n", sCmtLayerNum.c_str());
    fprintf(fp, "%s %d\n\n", sKwLayerNum.c_str(), nLayerNum);

    fprintf(fp, "%s\n", sCmtBoxFc.c_str());
    for (int i = 0; i < vecBoxFc.size(); i++)
    {
        fprintf(fp, "%s %d\n", sKwBoxFc.c_str(), vecBoxFc[i]);
    }
    fprintf(fp, "\n");

	fprintf(fp, "%s\n", sCmtSymmFc.c_str());
	for (int i = 0; i < vecSymmFc.size(); i++)
	{
		fprintf(fp, "%s %d\n", sKwSymmFc.c_str(), vecSymmFc[i]);
	}
	fprintf(fp, "\n");

    fprintf(fp, "%s\n", sCmtStepLen.c_str());
    fprintf(fp, "%s %lf\n\n", sKwStepLen.c_str(), dStepLen);

    fprintf(fp, "%s\n", sCmtisotropic_stop.c_str());
    fprintf(fp, "%s %d\n\n", sKwisotropic_stop.c_str(),isotropic_stop);


    fprintf(fp, "%s\n", sCmtStepLenRatio.c_str());
    fprintf(fp, "%s %d %lf %lf\n\n", sKwStepLenRatio.c_str(), nInitLayerNum, dStepLenInitRatio, dStepLenRatio);

    fprintf(fp, "\n");

    return 0;
}

int Config::WriteConfigMeshOptPart(FILE *fp)
{
    if(!fp)
        return 1;

    fprintf(fp, "%s\n", sCmtSmoothLoop.c_str());
    fprintf(fp, "%s %d\n", sKwSmoothLoop.c_str(), nSmoothLoop);

    fprintf(fp, "\n");

    return 0;
}


int Config::CopyConfig(const Config *cf)
{
    /* basic members */
    sConfigName = cf->sConfigName;
    iWrtLvl = cf->iWrtLvl;
    dElemNum = cf->dElemNum;

    /* geometry members */
    sGeoFileName = cf->sGeoFileName;
    sKwGeoFileName = cf->sKwGeoFileName;
    sCmtGeoFileName= cf->sCmtGeoFileName;

    iGeoType = cf->iGeoType;

    dSTolerance = cf->dSTolerance;

    dGTolerance = cf->dGTolerance;

    bDelSmallEdges = cf->bDelSmallEdges;

    /* parameters for mesh size */
    dMinSize = cf->dMinSize;				//the minimize mesh size

    dGlblSize = cf->dGlblSize;			//the global mesh size

    dExpandRatio = cf->dExpandRatio;

    mapPtSize = cf->mapPtSize;

    mapCvSize = cf->mapCvSize;

    mapFcSize = cf->mapFcSize;

    /* parameters for adaptive */
    iAdaptive = cf->iAdaptive;

    iProximityNum = cf->iProximityNum;

    dCurvatureAngle = cf->dCurvatureAngle;

    mapFcBku = cf->mapFcBku;

    mapFcBkv = cf->mapFcBkv;

    iuNum = cf->iuNum;

    ivNum = cf->ivNum;

    /* layer mesh */
    nLayerNum = cf->nLayerNum;

    vecBoxFc = cf->vecBoxFc;

    dStepLen = cf->dStepLen;

    dStepLenRatio = cf->dStepLenRatio;

    dStepLenInitRatio = cf->dStepLenInitRatio;
    nInitLayerNum = cf->nInitLayerNum;

    /* mesh improvement */
    nSmoothLoop = cf->nSmoothLoop;

    return 0;
}
