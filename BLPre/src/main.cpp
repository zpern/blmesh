
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <spdlog/spdlog.h> 
 #include "MeshInfo.h"
#include "../multiple_normal/include/MNormalMesh.h"
#include "blpre.h"
#include<iomanip>
#include <sstream>
#include <array>
using namespace std;
int *pntidx = NULL, *pntelm = NULL;

#define min(a,b)    (((a) < (b)) ? (a) : (b))
namespace PRE {
	void setpntelm(int npt, int nElm, int *pElm)
	{
		int i, j, start, *idx, *pelem;
		int npoints, nelm;

		nelm = nElm;
		npoints = npt;

		idx = new int[npoints + 1];
		pelem = new int[nelm * 2];
		pntidx = idx;
		pntelm = pelem;

		for (i = 0; i <= npoints; i++)
			idx[i] = 0;

		/*----------------------------------------------------------------------------
		| Count the number of elements for each point:
		----------------------------------------------------------------------------*/
		for (i = 0; i < nelm; i++)
		{
			for (j = 0; j < 2; j++)
				idx[pElm[3 * i + j]]++;
		}

		/*----------------------------------------------------------------------------
		| From the numbers of elements compute the startindex for each point.
		| The elements of point 'i' are then stored in pelem[idx[i]] to
		| pelem[idx[i+1]]-1.
		----------------------------------------------------------------------------*/
		start = 0;
		for (i = 0; i <= npoints; i++)
		{
			int count = idx[i];
			if (count == 0 && i < npoints)
			{
				spdlog::info("Node {} is not connected with any element!\n", i);
				return;
			}

			idx[i] = start;
			start += count;
		}

		/*----------------------------------------------------------------------------
		| Store the elements for each point. Thereby the 'idx' field is modified
		| and has to be restored afterwards.
		----------------------------------------------------------------------------*/
		for (i = 0; i < nelm; i++)
		{
			for (j = 0; j < 2; j++)
			{
				int pnt = pElm[3 * i + j];
				pelem[idx[pnt]] = i;
				idx[pnt]++;
			}
		}

		/*----------------------------------------------------------------------------
		| Restore 'idx' field:
		----------------------------------------------------------------------------*/
		start = 0;
		for (i = 0; i < npoints; i++)
		{
			int nstart = idx[i];
			idx[i] = start;
			start = nstart;
		}
	}

	void orientloop(int nPt, int nElm, double *pPt, int *pElm, int **elmN, int **ortEl)
	{
		int *ortElm = NULL, *pElmN = NULL, i, j, ptidxs, ptidxe;
		int beg, end, eidxp, eidx, k, fstidx, cntl = 0, ortidx = 0;
		bool *flag, newloop = false;

		int loopidx[100];
		double looparea[100];

		ortElm = new int[nElm];
		pElmN = new int[2 * nElm];
		flag = new bool[nElm];

		*ortEl = ortElm;
		*elmN = pElmN;

		for (i = 0; i < nElm; i++)
			flag[i] = false;

		for (i = 0; i < nElm; i++)
		{
			if (flag[i])
				continue;

			loopidx[cntl++] = ortidx;

			ortElm[ortidx] = i;
			fstidx = i;
			ptidxs = pElm[i * 3 + 1];
			eidxp = i;

			pElmN[ortidx * 2 + 0] = pElm[i * 3 + 0];
			pElmN[ortidx * 2 + 1] = pElm[i * 3 + 1];

			newloop = false;
			flag[i] = true;
			ortidx++;

			while (true)
			{
				beg = pntidx[ptidxs];
				end = pntidx[ptidxs + 1];

				eidxp = ortElm[ortidx - 1];

				for (j = beg; j < end; j++)
				{
					eidx = pntelm[j];

					if (eidx == eidxp)
						continue;
					else if (eidx == fstidx)
					{
						newloop = true;
						break;
					}

					ortElm[ortidx] = eidx;
					flag[eidx] = true;

					for (k = 0; k < 2; k++)
					{
						int idx = pElm[3 * eidx + k];
						if (idx != ptidxs)
						{
							ptidxe = idx;
							pElmN[2 * ortidx + 0] = ptidxs;
							pElmN[2 * ortidx + 1] = ptidxe;

							ptidxs = ptidxe;
							ortidx++;

							break;
						}
					}
				}

				if (newloop)
					break;

			}//while true

		}

		loopidx[cntl] = nElm;

		spdlog::info("loop number: {}\n", cntl);

		//orient each loop according to its area
		double max = 0.0;
		int maxi;
		for (i = 0; i < cntl; i++)
		{
			double t = 0.0, x1, y1, x2, y2;
			int idx1, idx2;
			beg = loopidx[i], end = loopidx[i + 1];
			for (j = beg; j < end; j++)
			{
				idx1 = pElmN[2 * j + 0];
				idx2 = pElmN[2 * j + 1];

				x1 = pPt[2 * idx1 + 0];
				y1 = pPt[2 * idx1 + 1];
				x2 = pPt[2 * idx2 + 0];
				y2 = pPt[2 * idx2 + 1];

				t += x1 * y2 - x2 * y1;
			}

			looparea[i] = t;
			spdlog::info("i{}: %lf\n", i, t);
			if (fabsf(t) > max)
			{
				max = fabsf(t);
				maxi = i;
			}
		}

		//orient
		for (i = 0; i < cntl; i++)
		{
			bool borient = false;
			if (i == maxi)
			{
				if (looparea[i] < 0.0)
				{
					borient = true;
				}
			}
			else
			{
				if (looparea[i] > 0.0)
				{
					borient = true;
				}
			}

			if (borient)
			{
				beg = loopidx[i], end = loopidx[i + 1];
				int tmp;
				for (j = beg; j < end; j++)
				{
					tmp = pElmN[2 * j + 0];
					pElmN[2 * j + 0] = pElmN[2 * j + 1];
					pElmN[2 * j + 1] = tmp;
				}
			}
		}
	}




#if 1
#define FILENAMESIZE 120
#define DEFAULTCFNAME "TRAN_CONFIG"

	char config[FILENAMESIZE];

#define ERRORCONFIG -1
#define DEFAULTCONFIG 1
#define INPUTCONFIG 0
#define NODEFAULTCONFIG -2

#define N0_LICENSE_FILE -1
#define INCORRECT_LICENSE_FILE 0
#define CORRECT_LICENSE_FILE 1
#define EXPIRED_LICENSE 2
#define INCORRECT_SYSTEM_TIME 3

	struct ConfigArgc {
		char filenam[1024];		// the project name
		int layer_num;
		double step_len;
		double ratio;
		std::vector<int> matchfc;
		std::vector<int> boxfc;
		std::vector<int> symfc;
		std::vector<int> adjacent;
		std::map<int, double> lnsize;
		std::map<int, double> fcsize;
	};

	void SetDefaultConfig(ConfigArgc &cf)
	{
		strcpy(cf.filenam, "automesh");
		cf.layer_num = 10;
		cf.step_len = 13;
		cf.ratio = 1.25;
	}

	int GetConfigName(int argc, char ** argv)
	{
		int i = argc;
		if (i > 2)
		{
			spdlog::info("Error : Please type in the config filename.\n");
			exit(ERRORCONFIG);
		}
		else if (i == 1)
		{
			//spdlog::info("Note : Use the default config.\n");
			strcpy(config, DEFAULTCFNAME);
			return DEFAULTCONFIG;
		}
		else
		{
			strcpy(config, argv[1]);
			return INPUTCONFIG;
		}
	}

	int parsecommand(int argc, char ** argv, struct ConfigArgc &cf)
	{
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

		GetConfigName(argc, argv);
		FILE *fp = fopen(config, "r");
		if (!fp)
		{
			spdlog::info("Config File does not exist.\n");
			throw(std::string("Config File does not exist."));
		}
		SetDefaultConfig(cf);

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
			}
			if (!strcmp(word, "layer_num")) {
				sscanf(line, "%s %d", word, &value);
				cf.layer_num = value;
			}
			if (!strcmp(word, "box_face")) {
				int iface;
				sscanf(line, "%s %d", word, &iface);
				cf.boxfc.push_back(iface);
			}
			if (!strcmp(word, "match_face")) {
				int iface;
				sscanf(line, "%s %d", word, &iface);
				cf.matchfc.push_back(iface);
			}
			if (!strcmp(word, "sym_face")) {
				int iface;
				sscanf(line, "%s %d", word, &iface);
				cf.symfc.push_back(iface);
			}
			if (!strcmp(word, "step_len")) {
				sscanf(line, "%s %lf", word, &fvalue);
				cf.step_len = fvalue;
			}
			if (!strcmp(word, "ratio")) {
				sscanf(line, "%s %lf", word, &fvalue);
				cf.ratio = fvalue;
			}

		}
		fclose(fp);
		return 0;
	}

	int printConfig(const struct ConfigArgc *cf)
	{
		char filenam[1024];		// the project name
		int layer_num;
		double step_len;
		double ratio;
		std::vector<int> boxfc;

		spdlog::info("------------- input parameters -------------");
		spdlog::info("file name:          {}", cf->filenam);
		spdlog::info("total layer number: {}", cf->layer_num);
		spdlog::info("initial step length:{}", cf->step_len);
		spdlog::info("growth rate:        {}", cf->ratio);
		spdlog::info("\n");

		return 0;
	}

	
#define _CONFIG_FILE
	double *pt = NULL;
	int *elm = NULL, npt, nelm, *idx, *pelem;

	std::tuple<std::string, double*, int*, int*, std::vector<double>> blpre(
		std::string& f,
		blpreConfig blcf,
	std::vector < std::array< double, 3> > points ,ControlVolume& cv1,ControlVolume& cv2) {
	vector<int> symm = blcf.symm;
    vector<int> wall = blcf.wall;
	vector<int> box = blcf.box;
	vector<int> match = blcf.match;
	vector<int> per_face;
	vector<int> adjacent_face = blcf.adjacent;
    bool fast_intersection = blcf.fast_intersection;
    bool exist_prism = blcf.preMultiple;
	for (auto i : blcf.per) {
		per_face.push_back(i);
		symm.push_back(i);
	}
    int faceCount = symm.size() + wall.size() + box.size() + match.size() + adjacent_face.size();
	int n = blcf.n;
	double len = blcf.len;
	double Ro = blcf.Ro;
	int argc = 1;
	bool use_multiple_normals = blcf.use_multiple_normals;
	char* argv[10];
	string mfile = "";
    if (use_multiple_normals) {
        // --- Step 1: split points by faceID ---
        std::vector<std::array<double, 3>> points_multiply, points_nonwall;
        std::string f_multiply, f_nonwall;

        splite_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall, wall);

        spdlog::info("finish splite.");

        // for fallback
        bool multiplySuccess = true;

        // --- Step 2: map point → length ---
        std::map<std::array<double, 3>, double> point_to_length;
        if (!blcf.length_vec.empty()) {
            for (int i = 0; i < points.size(); i++) {
                point_to_length[points[i]] = blcf.length_vec[i];
            }
        }

        // --- Step 3: prism pre-layer processing ---
        if (exist_prism) {
            ChamferBehavior behavior;
            MNormalMesh chamfer;

            chamfer.number_of_layer = 1;
            chamfer.step_of_length = blcf.multiple_steplength;
            chamfer.point_to_length = point_to_length;
            chamfer.fast_intersection = fast_intersection;

            chamfer.SetBehavior(behavior);
            chamfer.ReadPlsBuf(f_multiply, points_multiply);

            chamfer.CalculateMultiNormal();
            chamfer.BuildTopo(faceCount);

            spdlog::info("Handling output mesh!");

            chamfer.pre_WriteVol(cv2.v, cv2.f, cv2.lower_point_num, cv2.add_point_num);

            if (chamfer.multiplySuccess) {
                chamfer.WriteMesh(f_multiply, points_multiply, blcf.len);
                spdlog::info("PreJob Finished.");
            } else {
                multiplySuccess = chamfer.multiplySuccess;
            }
        }

        // --- Step 4: actual multilayer extrusion ---
        if (blcf.multiple_numlayer > 0 && multiplySuccess) {
            ChamferBehavior behavior;
            MNormalMesh chamfer;

            chamfer.number_of_layer = blcf.multiple_numlayer;
            chamfer.step_of_length = blcf.multiple_steplength;
            chamfer.point_to_length = point_to_length;
            chamfer.fast_intersection = fast_intersection;
            chamfer.exist_prism = exist_prism;

            chamfer.SetBehavior(behavior);
            chamfer.ReadPlsBuf(f_multiply, points_multiply);

            spdlog::info("Done!");

            chamfer.CalculateMultiNormal();
            chamfer.BuildTopo(faceCount);

            spdlog::info("Handling output mesh!");
            chamfer.WriteVol(cv1.v, cv1.f, cv1.lower_point_num, cv1.add_point_num);

            if (chamfer.multiplySuccess) {
                chamfer.WriteMesh(f_multiply, points_multiply, blcf.len);
            } 
        } else {
            spdlog::info("Skipping Chamfer processing because number_of_layer is 0");
        }

        // --- Step 5: merge back into original mesh ---
        combine_by_faceID(points, points_multiply, points_nonwall, f, f_multiply, f_nonwall);

        spdlog::info("Job Finished.");
    }


	if (mfile.size() == 0) {
		spdlog::info("use single-pass method");
		
		mfile = f;
	}
	spdlog::info("************************************************");
	stringstream file(mfile);
	stringstream fout;
	int tmp;
	char bdryfile[256], blmfile[256], outfile[256];
	double /**pt = NULL, */*bc = NULL, *symcoord = NULL;
	int /**elm = NULL, */npt1, npt2, nelm1, nelm2, nsym, *symaxis = NULL, *wbc = NULL;
	int i,j, npt, nelm, nbou, ibeg, iend, *bptratio = NULL;
	//FILE *fin1 = NULL, *fin2 = NULL, *fout = NULL, tmp;
	ConfigArgc cf;
#ifndef _CONFIG_FILE
	memset(bdryfile, 0, sizeof(bdryfile));
	strcpy(bdryfile, argv[1]);
	strcat(bdryfile, ".pls");

	memset(blmfile, 0, sizeof(blmfile));
	strcpy(blmfile, argv[1]);
	strcat(blmfile, "_box.fr2");

	memset(outfile, 0, sizeof(outfile));
	strcpy(outfile, argv[1]);
	strcat(outfile, "_o.pls");
#else
	//parsecommand(argc, argv, cf);
	cf.ratio = Ro;
	cf.layer_num = n;
	cf.step_len = len;
	cf.boxfc = box;
	cf.symfc = symm;
	cf.adjacent = adjacent_face;

	cf.matchfc = match;
	//cf.
	//ofstream fout1("1.txt");
	//fout1 << *f;
	strcpy(cf.filenam , "virtualmesh");
	//printConfig(&cf);

	memset(bdryfile, 0, sizeof(bdryfile));
	strcpy(bdryfile, cf.filenam);
	strcat(bdryfile, ".pls");

	memset(blmfile, 0, sizeof(blmfile));
	strcpy(blmfile, cf.filenam);
	strcat(blmfile, "_box.fr2");

	memset(outfile, 0, sizeof(outfile));
	strcpy(outfile, cf.filenam);
	strcat(outfile, "_o.pls");
#endif
	spdlog::info("begin to covert...");
//	ofstream	fout(string(outfile));
	string trash;
	string line;
	
	getline(file,trash);
	stringstream ss(trash);
	ss >> nelm1 >> npt1 >> nbou >> trash >> trash >> trash;

//	fscanf(fin1,"%d %d %d  0  0  0\n", &nelm1, &npt1, &nbou);

	npt = npt1;
	nelm = nelm1;
	
	pt = new double[3*npt];
	//elm = new int[3*nelm];
	elm = new int[4*nelm];
	bc = new double[2*nelm];
	wbc = new int[nelm];
	
	nsym = cf.symfc.size();
	if (nsym) {
		cout << "symm face";
		for (auto i : cf.symfc) {
			cout << i << " ";
		}
		cout << endl;
	}
	symcoord = new double[nsym];	//stores the coordinate of each symmetry plane
	symaxis = new int[nsym];	//stores the axis of each symmetry plane

	double x, y, z;
	for (i=0; i<npt1; i++)
	{
		if(points.empty())
			file >> tmp >> x >> y >> z;
		else {
			x = points[i][0];
			y = points[i][1];
			z = points[i][2];
		}
		//fscanf(fin1,"%d %lf %lf %lf\n", &tmp, &x, &y, &z);
		pt[i*3+0] = x;
		pt[i*3+1] = y;
		pt[i*3+2] = z;
	}

	std::set<int> nmpts;
	std::set<int>::iterator nmptit;

	std::multimap<int, int> mapFaceElms;

	int p1, p2, p3, fi, ci;
	for (i=0; i<nelm1; i++)
	{
		file >> tmp >> p1 >> p2 >> p3 >> fi;
	//	fscanf(fin1,"%d %d %d %d %d\n", &tmp, &p1, &p2, &p3, &fi);

		elm[4*i+0] = p1-1;
		elm[4*i+2] = p2-1;
		elm[4*i+1] = p3-1;
		elm[4*i+3] = fi;

		int nnmlt = 0;
#if 0
 		if(fi == 15 || fi == 18 || fi == 21 || fi == 33 || fi == 35 || 
			fi == 36 || fi == 40 || fi == 48)
 		{
 			nmpts.insert(p1);
 			nmpts.insert(p2);
 			nmpts.insert(p3);
 		}

		nnmlt = nmpts.size();
#endif
	
#ifndef _CONFIG_FILE
		//wbc[i] = 0(wall)\1(farfield)\2(symmetry)
		//if ((fi >= 16 && fi <= 21))
		//if ((fi >= 354 && fi <= 359))
		//if ((fi >= 2 && fi <= 6))
		if(fi != 5 /*&& fi != 6*/)
		{
			bc[2*i+0] = 1;
			bc[2*i+1] = 0.1;
			wbc[i] = 0;
		}
		else /*if(fi == 1)*/
		{
			bc[2*i+0] = 2;
			bc[2*i+1] = 0;
			wbc[i] = 2;
		}
		/*else
		{
			bc[2*i+0] = 1;
			bc[2*i+1] = 0.1;
			wbc[i] = 0;
		}*/
#else
		std::vector<int>::iterator iterBox, iterSym,iterMatch,iterPer,iterAdjacent;
		iterBox = find(cf.boxfc.begin(), cf.boxfc.end(), fi);
		iterSym = find(cf.symfc.begin(), cf.symfc.end(), fi);
		iterMatch= find(cf.matchfc.begin(), cf.matchfc.end(), fi);
		iterPer = find(per_face.begin(), per_face.end(), fi);
		iterAdjacent = find(cf.adjacent.begin(),cf.adjacent.end(),fi);
		if(iterBox!=cf.boxfc.end())
		{
			bc[2*i+0] = 1;
			bc[2*i+1] = 0;
			wbc[i] = 1;
		}
		else if(iterSym!=cf.symfc.end())
		{
			bc[2*i+0] = 2;
			bc[2*i+1] = 0;
			wbc[i] = 2;
			mapFaceElms.insert(std::make_pair(fi, i));
			if (iterPer != per_face.end()) {
				bc[2 * i + 0] = 5;
				bc[2 * i + 1] = 0;
				wbc[i] = 5;
			}
		}
		else if (iterMatch != cf.matchfc.end()) {
			bc[2 * i + 0] = 4;
			bc[2 * i + 1] = 0;
			wbc[i] = 4;
		}
		else if (iterAdjacent != cf.adjacent.end()) {
			bc[2 * i + 0] = 6;
			bc[2 * i + 1] = 0;
			wbc[i] = 6;
		}
		else
		{
			bc[2*i+0] = 1;
			bc[2*i+1] = 0.1;
			wbc[i] = 0;
		}
#endif
	}
	
	//获取边界点
	int nbdryPt, *bpt, nbdry, *bdry;
	MeshInfo mshInfo(MeshType::MESH_3D);
	mshInfo.Initialize(pt, npt, nelm, elm);
	npt = mshInfo.GetNumPoint();
	mshInfo.GetBdryPt(&nbdryPt, &bpt);

	//获取边界边
	mshInfo.GetBdryElm(&nbdry, &bdry);

	// determines the symmetry planes and coordinate of each symmetry plane
	int *nSymBdry = NULL, **pSymBdry = NULL, nTtlSymBdry = 0, *pSymFidx = NULL;
	nSymBdry = new int[nsym];
	pSymBdry = new int *[nsym];
	pSymFidx = new int[nsym];
	std::vector<int>::iterator itbeg, itend;
	itbeg = cf.symfc.begin();
	itend = cf.symfc.end();
	int fcnt = 0;
	while (itbeg != itend)
	{
		double xcoord=0, ycoord, zcoord, eps = 1.0e-5;
		double xttl = 0.0, yttl = 0.0, zttl = 0.0;
		double xavg, yavg, zavg;
		double xeps = 0, yeps = 0, zeps = 0;
		int fidx, icnt = 0;
		bool bx = true, by = true, bz = true;
		fidx = *itbeg;

		pSymFidx[fcnt] = fidx;
		std::multimap<int, int>::iterator mapbeg, mapend;
		mapbeg = mapFaceElms.lower_bound(fidx);
		mapend = mapFaceElms.upper_bound(fidx);
		while (mapbeg != mapend)
		{
			int eidx = mapbeg->second;
			int pidx = elm[4*eidx+0];

			xcoord = pt[pidx*3+0];
			ycoord = pt[pidx*3+1];
			zcoord = pt[pidx*3+2];

			xttl += pt[pidx*3+0];
			yttl += pt[pidx*3+1];
			zttl += pt[pidx*3+2];

			xavg = xttl/(icnt+1);
			yavg = yttl/(icnt+1);
			zavg = zttl/(icnt+1);
			if (icnt != 0)
			{
				if (fabs(xavg - xcoord) > xeps)
					xeps = fabs(xavg - xcoord);

				if (fabs(yavg - ycoord) > yeps)
					yeps = fabs(yavg - ycoord);

				if (fabs(zavg - zcoord) > zeps)
					zeps = fabs(zavg - zcoord);
			}
			
			++icnt;
			++mapbeg;
		}

		double mineps = min(xeps, min(yeps, zeps));

		if(fabs(xeps-mineps)<eps)
		{
			symaxis[fcnt] = 0;
			symcoord[fcnt] = xcoord;
		}
		else if(fabs(yeps-mineps)<eps)
		{
			symaxis[fcnt] = 1;
			symcoord[fcnt] = ycoord;
		}
		else if (fabs(zeps-mineps)<eps)
		{
			symaxis[fcnt] = 2;
			symcoord[fcnt] = zcoord;
		}
		else
		{
			symaxis[fcnt] = 0;
			symcoord[fcnt] = xcoord;
		}


		mshInfo.GetBdryElm(&nSymBdry[fcnt], &pSymBdry[fcnt], fidx);

		nTtlSymBdry += nSymBdry[fcnt];

		++fcnt;
		++itbeg;
	}

	int nboupt = 0;
#if 0
	//collect points on curves
	std::set<int> lpts;
	std::set<int>::iterator ptit;
	for (i=0; i<nbou; i++)
	{
		fscanf(fin1, "%d %d %d %d %d\n", &tmp, &p1, &p2, &ci, &tmp);

		if(ci == 60 || ci == 62 ||ci == 71)
		{
			lpts.insert(p1);
			lpts.insert(p2);
		}
	}
	nboupt = lpts.size();
#endif


#if 0
	ptit = lpts.begin();
	while (ptit != lpts.end())
	{
		int pidx = *ptit;
		ibeg = idx[pidx];
		iend = idx[pidx+1];
		for (i=ibeg; i<iend; i++)
		{
			bc[2*i+1] = 0.15;
		}
		++ptit;
	}
#endif

	
	std::vector<std::array<double, 3>> point_array;
	std::vector<double> sym_coord;
	fout << nelm << " " << npt << " "<<nsym << " " << nTtlSymBdry << " " << nbdryPt << endl;

	nmptit = nmpts.begin();
	i=0;
	while (nmptit != nmpts.end())
	{
		++i;
		int pidx = *nmptit;
		//fprintf(fout, "%d %d %d\n", i, pidx, 2);
		fout << i << " " << pidx << " " << 2 << endl;
		++nmptit;
	}
	
	for (i=0; i<nsym; i++)
	{
		//fprintf(fout, "%d %d %d %lf\n", nSymBdry[i], symaxis[i], pSymFidx[i], symcoord[i]);
		sym_coord.push_back(symcoord[i]);
		fout << nSymBdry[i]<<" "<< symaxis[i]<<" "<< pSymFidx[i] << endl;
		for (j=0; j<nSymBdry[i]; j++)
			//fprintf(fout, "%d %d %d\n", j+1, pSymBdry[i][j*2+0]+1, pSymBdry[i][j*2+1]+1);
			fout << j + 1 << " " << pSymBdry[i][j * 2 + 0] + 1<<" "<<pSymBdry[i][j * 2 + 1] + 1 << endl;
	}

	for (i=0; i<nbdryPt; i++)
	{
		//fprintf(fout, "%d ", bpt[i]+1);
		fout << bpt[i] + 1 << " ";
			if ((i + 1) % 10 == 0)
				//fprintf(fout, "\n");
				fout<<endl;
	}
	//fprintf(fout, "\n");
	fout << endl;

	//fclose(fout);
	//fout = NULL;

	spdlog::info("finished converting files!");

	delete []symcoord;
	delete []symaxis;
	if (bdry) {
		delete[]bdry;
		bdry = nullptr;
	}
	if (bpt) {
		delete[]bpt;
		bpt = nullptr;
	}

	if (pSymBdry) {
		for (int i = 0; i < nsym; i++)
			delete[] pSymBdry[i];
		delete pSymBdry;
		pSymBdry = nullptr;
	}
	if (bc) {
		delete[]bc;
		bc = nullptr;
	}

	if (nSymBdry) {
		delete[]nSymBdry;
		nSymBdry = nullptr;
	}
	if (pSymFidx) {
		delete[]pSymFidx;
		pSymFidx = nullptr;
	}
	
	return std::make_tuple(fout.str(),pt,elm,wbc,sym_coord);
}

void setpntelm()
{
	int i, j, start/*, *idx, *pelem*/;
	int npoints, nelm;
	
	nelm = nelm;
	npoints = npt;
	
	idx = new int[npoints + 1];
	pelem = new int[nelm*3];
	
	for(i = 0; i <= npoints; i++)
		idx[i] = 0;
	
	/*----------------------------------------------------------------------------
	| Count the number of elements for each point:
	----------------------------------------------------------------------------*/
	int *pelm = elm;
	for (i=0; i<nelm; i++)
	{
		for(j=0; j<3; j++)
			idx[pelm[3*i + j]]++;
	}
	
	/*----------------------------------------------------------------------------
	| From the numbers of elements compute the startindex for each point.
	| The elements of point 'i' are then stored in pelem[idx[i]] to
	| pelem[idx[i+1]]-1.
	----------------------------------------------------------------------------*/
	start = 0;
	for(i = 0; i <= npoints; i++)
	{
		int count = idx[i];
		if(count == 0 && i < npoints)
		{
			spdlog::info("Node {} is not connected with any element!\n", i);
			return ;
		}
		
		idx[i] = start;
		start += count;
	}
	
	/*----------------------------------------------------------------------------
	| Store the elements for each point. Thereby the 'idx' field is modified
	| and has to be restored afterwards.
	----------------------------------------------------------------------------*/
	for(i = 0; i < nelm; i++)
	{
		for(j = 0; j < 3; j++)
		{
			int pnt = pelm[3*i + j];
			pelem[idx[pnt]] = i;
			idx[pnt]++;
		}
	}
	
	/*----------------------------------------------------------------------------
	| Restore 'idx' field:
	----------------------------------------------------------------------------*/
	start = 0;
	for(i = 0; i < npoints; i++)
	{
		int nstart = idx[i];
		idx[i] = start;
		start = nstart;
	}
}

#endif

bool isnear(int pidx, int eidx, double *pt, int npt, int *elm)
{
	int pt0, pt1;
	double dis, base = 0.05;

	pt0 = elm[2*eidx+0];
	pt1 = elm[2*eidx+1];

	if(pidx < 0 || pidx > npt)
		return false;

	dis = sqrt((pt[pidx*2+0]-pt[pt0*2+0])*(pt[pidx*2+0]-pt[pt0*2+0]) + (pt[pidx*2+1]-pt[pt0*2+1])*(pt[pidx*2+1]-pt[pt0*2+1]));
	if(dis < base)
		return true;

	dis = sqrt((pt[pidx*2+0]-pt[pt1*2+0])*(pt[pidx*2+0]-pt[pt1*2+0]) + (pt[pidx*2+1]-pt[pt1*2+1])*(pt[pidx*2+1]-pt[pt1*2+1]));
	if(dis < base)
		return true;

	return false;
}




}
