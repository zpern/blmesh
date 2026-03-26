
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
#include "../MNormal/include/MNormalMesh.h"
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

    for (i = 0; i <= npoints; i++) {
        idx[i] = 0;
    }

    /*----------------------------------------------------------------------------
    | Count the number of elements for each point:
    ----------------------------------------------------------------------------*/
    for (i = 0; i < nelm; i++) {
        for (j = 0; j < 2; j++) {
            idx[pElm[3 * i + j]]++;
        }
    }

    /*----------------------------------------------------------------------------
    | From the numbers of elements compute the startindex for each point.
    | The elements of point 'i' are then stored in pelem[idx[i]] to
    | pelem[idx[i+1]]-1.
    ----------------------------------------------------------------------------*/
    start = 0;
    for (i = 0; i <= npoints; i++) {
        int count = idx[i];
        if (count == 0 && i < npoints) {
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
    for (i = 0; i < nelm; i++) {
        for (j = 0; j < 2; j++) {
            int pnt = pElm[3 * i + j];
            pelem[idx[pnt]] = i;
            idx[pnt]++;
        }
    }

    /*----------------------------------------------------------------------------
    | Restore 'idx' field:
    ----------------------------------------------------------------------------*/
    start = 0;
    for (i = 0; i < npoints; i++) {
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

    for (i = 0; i < nElm; i++) {
        flag[i] = false;
    }

    for (i = 0; i < nElm; i++) {
        if (flag[i]) {
            continue;
        }

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

        while (true) {
            beg = pntidx[ptidxs];
            end = pntidx[ptidxs + 1];

            eidxp = ortElm[ortidx - 1];

            for (j = beg; j < end; j++) {
                eidx = pntelm[j];

                if (eidx == eidxp) {
                    continue;
                } else if (eidx == fstidx) {
                    newloop = true;
                    break;
                }

                ortElm[ortidx] = eidx;
                flag[eidx] = true;

                for (k = 0; k < 2; k++) {
                    int idx = pElm[3 * eidx + k];
                    if (idx != ptidxs) {
                        ptidxe = idx;
                        pElmN[2 * ortidx + 0] = ptidxs;
                        pElmN[2 * ortidx + 1] = ptidxe;

                        ptidxs = ptidxe;
                        ortidx++;

                        break;
                    }
                }
            }

            if (newloop) {
                break;
            }

        } // while true
    }

    loopidx[cntl] = nElm;

    spdlog::info("loop number: {}\n", cntl);

    // orient each loop according to its area
    double max = 0.0;
    int maxi;
    for (i = 0; i < cntl; i++) {
        double t = 0.0, x1, y1, x2, y2;
        int idx1, idx2;
        beg = loopidx[i], end = loopidx[i + 1];
        for (j = beg; j < end; j++) {
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
        if (fabsf(t) > max) {
            max = fabsf(t);
            maxi = i;
        }
    }

    // orient
    for (i = 0; i < cntl; i++) {
        bool borient = false;
        if (i == maxi) {
            if (looparea[i] < 0.0) {
                borient = true;
            }
        } else {
            if (looparea[i] > 0.0) {
                borient = true;
            }
        }

        if (borient) {
            beg = loopidx[i], end = loopidx[i + 1];
            int tmp;
            for (j = beg; j < end; j++) {
                tmp = pElmN[2 * j + 0];
                pElmN[2 * j + 0] = pElmN[2 * j + 1];
                pElmN[2 * j + 1] = tmp;
            }
        }
    }
}

#define FILENAMESIZE  120
#define DEFAULTCFNAME "TRAN_CONFIG"

char config[FILENAMESIZE];

#define ERRORCONFIG     -1
#define DEFAULTCONFIG   1
#define INPUTCONFIG     0
#define NODEFAULTCONFIG -2

#define N0_LICENSE_FILE        -1
#define INCORRECT_LICENSE_FILE 0
#define CORRECT_LICENSE_FILE   1
#define EXPIRED_LICENSE        2
#define INCORRECT_SYSTEM_TIME  3
#define _CONFIG_FILE
struct ConfigArgc {
    char filenam[1024]; // the project name
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

double *pt = NULL;
int *elm = NULL, npt, nelm, *idx, *pelem;

std::tuple<
    std::string,
    std::vector<std::array<double, 3>>,
    std::vector<std::array<int, 4>>,
    std::vector<int>,
    std::vector<double>> blpre(blpreConfig blcf,
                               std::string &f,
                               std::vector<std::array<double, 3>> points)
{
    int n = blcf.n;
    double len = blcf.len;
    double Ro = blcf.Ro;
    vector<int> box = blcf.box;
    vector<int> wall = blcf.wall;
    vector<int> symm = blcf.symm;
    vector<int> match = blcf.match;
    vector<int> per_face;
    vector<int> adjacent_face = blcf.adjacent;

    for (auto i : blcf.per) {
        per_face.push_back(i);
        symm.push_back(i);
    }
    int faceCount = box.size() + wall.size() + symm.size() + match.size() + adjacent_face.size();

    spdlog::info("begin to covert...");

    stringstream file(f);
    stringstream fout;
    string line;
    string trash;
    getline(file, line);
    stringstream ss(line);

    int nelm, npt;
    ss >> nelm >> npt >> trash >> trash >> trash >> trash;
    std::vector<std::array<double, 3>> pt(npt);
    std::vector<std::array<int, 4>> elm(nelm);
    std::vector<std::array<double, 2>> bc(nelm);
    std::vector<int> wbc(nelm);

    // store symm
    if (symm.size()) {
        cout << "symm face";
        for (auto i : symm) {
            cout << i << " ";
        }
        cout << endl;
    }

    double x, y, z;
    for (int i = 0; i < npt; i++) {
        if (points.empty()) {
            throw std::runtime_error("blpre points is empty ");
        } else {
            pt[i][0] = points[i][0];
            pt[i][1] = points[i][1];
            pt[i][2] = points[i][2];
        }
    }

    std::set<int> nmpts;
    std::set<int>::iterator nmptit;
    std::multimap<int, int> mapFaceElms;

    int p1, p2, p3, fi, ci;
    for (int i = 0; i < nelm; i++) {
        file >> trash >> p1 >> p2 >> p3 >> fi;

        elm[i][0] = p1;
        elm[i][1] = p3;
        elm[i][2] = p2;
        elm[i][3] = fi;

        int nnmlt = 0;

        auto iterBox = std::find(box.begin(), box.end(), fi);
        auto iterSym = std::find(symm.begin(), symm.end(), fi);
        auto iterMatch = std::find(match.begin(), match.end(), fi);
        auto iterPer = std::find(per_face.begin(), per_face.end(), fi);
        auto iterAdjacent = std::find(adjacent_face.begin(), adjacent_face.end(), fi);

        if (iterBox != box.end()) {
            wbc[i] = 0;
        } 
        else if (iterSym != symm.end()) {
            wbc[i] = 2;

            mapFaceElms.insert({fi, i});

            if (iterPer != per_face.end()) {
                wbc[i] = 5;
            }
        } 
        else if (iterMatch != match.end()) {
            wbc[i] = 3;
        } 
        else if (iterAdjacent != adjacent_face.end()) {
            wbc[i] = 5;
        } 
        else {
            wbc[i] = 1;
        }
    }

    // 建立mesh，检查网格单元
    MeshInfo mshInfo(MeshType::MESH_3D);
    mshInfo.Initialize(npt, nelm, pt, elm);

    auto bdryPt = mshInfo.getBoundaryPoints();
    auto bdryEdge = mshInfo.getBoundaryEdges();

    // determines the symmetry planes and coordinate of each symmetry plane
    int nTtlSymBdry = 0;

    // 对称面信息 需要和blmesh中setboundary一起修改，都是无用工
    std::vector<int> nSymBdry(symm.size(), 0);
    std::vector<std::vector<int>> pSymBdry(symm.size()); // 每个对称面存成扁平数组：[p0,p1,p0,p1,...]
    std::vector<int> pSymFidx(symm.size(), -1);
    std::vector<int> symaxis(symm.size(), 3);
    std::vector<double> symcoord(symm.size(), 0.0);

    int fcnt = 0;
    for (int fidx : symm) {
        double xcoord = 0.0, ycoord = 0.0, zcoord = 0.0;
        double eps = 1.0e-5;

        double xttl = 0.0, yttl = 0.0, zttl = 0.0;
        double xavg = 0.0, yavg = 0.0, zavg = 0.0;
        double xeps = 0.0, yeps = 0.0, zeps = 0.0;

        int icnt = 0;

        pSymFidx[fcnt] = fidx;

        auto mapbeg = mapFaceElms.lower_bound(fidx);
        auto mapend = mapFaceElms.upper_bound(fidx);

        while (mapbeg != mapend) {
            int eidx = mapbeg->second;
            int pidx = elm[eidx][0];

            xcoord = pt[pidx][0];
            ycoord = pt[pidx][1];
            zcoord = pt[pidx][2];

            xttl += xcoord;
            yttl += ycoord;
            zttl += zcoord;

            xavg = xttl / (icnt + 1);
            yavg = yttl / (icnt + 1);
            zavg = zttl / (icnt + 1);

            if (icnt != 0) {
                xeps = std::max(xeps, std::fabs(xavg - xcoord));
                yeps = std::max(yeps, std::fabs(yavg - ycoord));
                zeps = std::max(zeps, std::fabs(zavg - zcoord));
            }

            ++icnt;
            ++mapbeg;
        }

        double mineps = min(xeps, min(yeps, zeps));

        if (std::fabs(xeps - mineps) < eps) {
            symaxis[fcnt] = 0;
            symcoord[fcnt] = xcoord;
        } else if (std::fabs(yeps - mineps) < eps) {
            symaxis[fcnt] = 1;
            symcoord[fcnt] = ycoord;
        } else if (std::fabs(zeps - mineps) < eps) {
            symaxis[fcnt] = 2;
            symcoord[fcnt] = zcoord;
        } else {
            symaxis[fcnt] = 3;
            symcoord[fcnt] = xcoord;
        }

        // 假设返回 vector<pair<int,int>> 或 vector<array<int,2>>
        auto edges = mshInfo.getBoundaryEdges(fidx);

        nSymBdry[fcnt] = static_cast<int>(edges.size());
        pSymBdry[fcnt].reserve(2 * nSymBdry[fcnt]);

        for (const auto &e : edges) {
            // 如果是 pair<int,int>，就改成 e.first / e.second
            pSymBdry[fcnt].push_back(e[0]);
            pSymBdry[fcnt].push_back(e[1]);
        }

        nTtlSymBdry += nSymBdry[fcnt];
        ++fcnt;
    }

    std::vector<double> sym_coord;
    sym_coord.reserve(symm.size());

    fout << nelm << " " << npt << " " << symm.size() << " " << nTtlSymBdry << " " << bdryPt.size() << endl;

    nmptit = nmpts.begin();
    int i = 0;
    while (nmptit != nmpts.end()) {
        ++i;
        int pidx = *nmptit;
        fout << i << " " << pidx << " " << 2 << endl;
        ++nmptit;
    }

    for (int i = 0; i < symm.size(); i++) {
        sym_coord.push_back(symcoord[i]);

        fout << nSymBdry[i] << " " << symaxis[i] << " " << pSymFidx[i] << endl;
        for (int j = 0; j < nSymBdry[i]; j++) {
            fout << j + 1 << " " << pSymBdry[i][j * 2 + 0] + 1 << " " << pSymBdry[i][j * 2 + 1] + 1
                 << endl;
        }
    }

    for (i = 0; i < bdryPt.size(); i++) {
        fout << bdryPt[i] + 1 << " ";
        if ((i + 1) % 10 == 0) {
            fout << endl;
        }
    }

    fout << endl;

    spdlog::info("finished converting files!");

    return std::make_tuple(fout.str(), pt, elm, wbc, sym_coord);
}


std::tuple<
    std::string,
    double*,
    int*,
    int*,
    std::vector<double>> temptransform(std::tuple<
                                       std::string,
                                       std::vector<std::array<double, 3>>,
                                       std::vector<std::array<int, 4>>,
                                       std::vector<int>,
                                       std::vector<double>>& in)
{
    const auto& outstr    = std::get<0>(in);
    const auto& pt        = std::get<1>(in);
    const auto& elm       = std::get<2>(in);
    const auto& wbc       = std::get<3>(in);
    const auto& sym_coord = std::get<4>(in);

    double* pt_raw = nullptr;
    int* elm_raw = nullptr;
    int* wbc_raw = nullptr;

    try {
        // pt: vector<array<double,3>> -> double*
        pt_raw = new double[pt.size() * 3];
        for (size_t i = 0; i < pt.size(); ++i) {
            pt_raw[i * 3 + 0] = pt[i][0];
            pt_raw[i * 3 + 1] = pt[i][1];
            pt_raw[i * 3 + 2] = pt[i][2];
        }

        // elm: vector<array<int,4>> -> int*
        elm_raw = new int[elm.size() * 4];
        for (size_t i = 0; i < elm.size(); ++i) {
            elm_raw[i * 4 + 0] = elm[i][0];
            elm_raw[i * 4 + 1] = elm[i][1];
            elm_raw[i * 4 + 2] = elm[i][2];
            elm_raw[i * 4 + 3] = elm[i][3];
        }

        // wbc: vector<int> -> int*
        wbc_raw = new int[wbc.size()];
        for (size_t i = 0; i < wbc.size(); ++i) {
            wbc_raw[i] = wbc[i];
        }
    }
    catch (...) {
        delete[] pt_raw;
        delete[] elm_raw;
        delete[] wbc_raw;
        throw;
    }

    return std::make_tuple(outstr, pt_raw, elm_raw, wbc_raw, sym_coord);
}
}

