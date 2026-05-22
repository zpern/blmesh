
#include "BLMesh.h"
#include "BLFront3D.h"
#include "BLLineMesh.h"
#include "BLMesh_define.h"
#include "BLNode.h"
#include "genmesh2d.h"
#include "mesh_orient.h"
#include "mshgen3d_def.h"
#include "spdlog/spdlog.h"
#include <set>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifndef __APPLE__
#include <omp.h>
#endif
#include "../postprocess/postprocess.h"
#include "dataio.h"
#include <algorithm>
#include <atomic>
#include <time.h>

#include "../include/VoronoiNormalSmoothStrategy.h"
#include "../include/WangzhiNormalSmoothStrategy.h"
#include "../include/distancecalculation.h"

#include "../2dremesh/SurfaceRemesh.h"
#include "../2dremesh/Types.h"

#include "igl/list_to_matrix.h"
#include "igl/remove_unreferenced.h"

#include <algorithm>
#include <set>

#ifdef infect
#undef infect
#endif
#ifdef triface
#undef triface
#endif
#include "../third/tetgen/tetgen.h"
#include "tiger_volmtetrahedralizer.h"

#include "../common/singleton_terminate.h"
#include <Eigen/Dense>
#include <igl/bfs_orient.h>
#include <igl/list_to_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/triangle/triangulate.h>
#include <igl/volume.h>
// #define _FAST_MULTIPOLE

#define _GEOM_NORMAL
#define _VECTOR_FIELD
// #define OLD
// for test(need to be deleted later)
int printflag;
int tempcalculate = 0;
// a callback function of threads
unsigned bccal(void *para)
{
    BCArg *bcArg = (BCArg *)para;
    return 0;
}

void BLMesh::DeleteOctree()
{

    if (m_ocAgent) {
        delete m_ocAgent;
        m_ocAgent = nullptr;
    }
    if (m_ocTree) {
        delete m_ocTree;
        m_ocTree = nullptr;
    }
}
void BLMesh::DeleteOctreeSymm()
{

    // if (m_ocAgent_symm)
    //{
    //	delete m_ocAgent_symm;
    //	m_ocAgent_symm = nullptr;
    // }
    // if (m_ocTree_symm)
    //{
    //	delete m_ocTree_symm;
    //	m_ocTree_symm = nullptr;
    // }
}
BLMesh::~BLMesh(void)
{

    DeleteOctree();
    DeleteOctreeSymm();
    FreeMemoryInFrontAndNode();

    if (m_pNodes) {

        delete[] (m_pNodes);
        m_pNodes = nullptr;
    }
    if (m_pElems) {
        free(m_pElems);
        m_pElems = nullptr;
    }
    m_mapSymline.clear();
    distance_ratio.clear();
    m_vBdyFront.clear();
    m_vecData.clear();
    if (m_symVals) {
        delete[] m_symVals;
        m_symVals = nullptr;
    }
    if (m_symFidx) {
        delete[] m_symFidx;
        m_symFidx = nullptr;
    }
    if (m_pSymBdrys) {
        free(m_pSymBdrys);
        m_pSymBdrys = nullptr;
    }
    if (m_pPntIdx) {
        delete[] m_pPntIdx;
        m_pPntIdx = nullptr;
    }
    if (m_pPntElm) {
        delete[] m_pPntElm;
        m_pPntElm = nullptr;
    }
    if (m_pBElm) {
        delete[] m_pBElm;
        m_pBElm = nullptr;
    }
    if (m_pU) {
        delete[] m_pU;
        m_pU = nullptr;
    }
    if (m_pNorm) {
        delete[] m_pNorm;
        m_pNorm = nullptr;
    }

    if (m_pBdryPnt) {
        free(m_pBdryPnt);
        m_pBdryPnt = nullptr;
    }
    if (outbdry) {
        delete[] outbdry;
        outbdry = nullptr;
    }

    // delete m_pPotentialBEM;
}
auto intersect_two = [](const std::vector<int> &a, const std::vector<int> &b) {
    std::vector<int> result;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
    return result;
};
/**

  *  @author yhf

  *  @brief  output the problem version and time stamp

  *  @param[in]  inArgName input argument description.type name:purpose(Range)

  *  @param[out] outArgName output argument description.

  *  @Return:    :description of the return value

  *  @note   matters needing attention

  *  @see    other functions

  */
void BLMesh::InitBLMesh()
{

    m_blFrontList = std::shared_ptr<BLFrontList>(new BLFrontList);
    m_blFrontListAll = std::vector<std::shared_ptr<BLFrontList>>(m_nTotalLayer);
    GEOM_FUNC::exactinit();
    time_t tt = time(NULL); // return a time stamp
    tm *t = localtime(&tt);
}

double BLMesh::GetAvergeLayer() { return m_nPrism * 1.0 / first_layer_FrtNods; }

#define _NEW_SYMM

/**

  *  @author yhf

  *  @brief  read and handle input file

  *  @param[in]  inArgName input argument description.type name:purpose(Range)

  *  @param[out] outArgName output argument description.

  *  @Return:  whether the code run successfully

  *  @note   matters needing attention

  *  @see    ReadBoundary(string filename)

  */

int BLMesh::SetBoundary(INPUTFORMAT file)
{

    //	spdlog::info("boufilename=%s", filename);
    // FILE *fin = nullptr;
    int i, j, npt, nelm, nboupt, nnmlt, tmp, nconn, nbdrypt, nbdry, dim, nbdout;
    int nsymloop, nsymbdry, pidx;
    BLVector box_max(-1e200, -1e200, -1e200);
    BLVector box_min(1e200, 1e200, 1e200);
    // 	int nbc, *ielm = nullptr;
    // 	double *pt = nullptr, *bcv = nullptr, *dnorm = nullptr, *u = nullptr;

    /*fin = fopen(filename, "r");
    if (!fin)
    {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }*/

    stringstream *fin = new stringstream(std::get<0>(file));
    nboupt = 0;
    nnmlt = 0;
    nbdrypt = 0;
    nbdry = 0;
    nsymloop = 0;

    // fscanf(fin, "%d %d %d %d %d\n", &nelm, &npt, &nsymloop, &nsymbdry, &nbdrypt);
    *fin >> nelm >> npt >> nsymloop >> nsymbdry >> nbdrypt;
    // nbdrypt = 0;

    m_nSurfNodes = npt;
    m_nSurfElems = nelm;

    m_vBdyFront.reserve(int(m_nSurfElems * 3)); // 预分配内存，加快速度
    m_nBdryPnt = nbdrypt;

    m_nSymLoop = nsymloop;

    m_nNodes = npt;
    m_nElems = nelm;

    // 	bcv = new double[2*nelm];
    // 	ielm = new int[2*nelm];
    // 	pt = new double[2*npt];
    // 	dnorm = new double[nelm];
    // 	u = new double[nelm];
    PreparePotential();

    m_pNodes = new MBLNode[npt];
    if (m_pNodes == nullptr) {
        spdlog::info("Can not allocate memory for Nodes\n");
    }
    m_nAllocNodes = m_nNodes;

    m_pElems = (Elem *)malloc(sizeof(Elem) * nelm);
    if (m_pElems == nullptr) {
        spdlog::info("Can not allocate memory for Elements\n");
    }
    m_nAllocElems = m_nElems;

    m_symVals = new double[m_nSymLoop];
    m_symFidx = new int[m_nSymLoop];

    if (nbdrypt > 0) {
        m_pBdryPnt = (int *)malloc(sizeof(int) * nbdrypt);
        if (m_pBdryPnt == nullptr) {
            spdlog::info("Can not allocate memory for Boundary points\n");
        }
    }

#ifndef _NEW_SYMM
    m_pSymBdrys = (Elem *)malloc(sizeof(Elem) * nbdry);
    if (m_pSymBdrys == nullptr) {
        spdlog::info("Can not allocate memory for Boundary points\n");
        return OUT_OF_MEMORY;
    }
    m_nAllocSymBdrys = nbdry;
#endif

#ifdef _CHECK_INTERSECTION
    // intersection data
    spdlog::info("number of surface element = {}"); // << nelm << endl;
    m_TriElm.Reserve(10 * nelm);

    // m_nTriElm = nelm;
    //  end of intersection data
#endif

    // stores the boundary of the region needed to be meshed with unstructured mesh
    outbdry = new int[m_nSurfElems * 12];
    noutbdry = 0;
    nbdout = 0;
    dim = 3;

    // 预分配内存
#ifdef USE_MEMORY_POOL
    BLFront::preAllocate(nelm * 3);
    BLNode::preAllocate(npt * 3);
#endif

    Eigen::MatrixXd V;
    V.resize(npt, 3);
    for (i = 0; i < npt; i++) {

        m_pNodes[i].coord[0] = std::get<1>(file)[3 * i + 0];
        m_pNodes[i].coord[1] = std::get<1>(file)[3 * i + 1];
        m_pNodes[i].coord[2] = std::get<1>(file)[3 * i + 2];

        for (int k = 0; k < 3; k++) {
            V(i, k) = m_pNodes[i].coord[k];
        }

        // create new front node
        BLNode *blNode = new BLNode(i);
        blNode->SetDecentID(i);
        blNode->SetLayerNum(0);

        m_pNodes[i].pointer = (void *)blNode;
        m_pNodes[i].bsysm = false;
        m_pNodes[i].bfarfield = false;
        m_pNodes[i].bwall = false;
        m_pNodes[i].reserved = -1;
    }

    double ave_front_size = 0;
    int front_count = 0;
    int cnt = 0;
    int *ele = get<2>(file);
    vector<int> l2g(nelm, -1);

    std::vector<double> point_size(npt + 1, 0);
    std::vector<int> point_manifold_size(npt + 1, 0);
    std::map<int, std::vector<std::array<int, 3>>> symmface;

    std::vector<int> per_point(npt, 0);
    for (i = 0; i < nelm; i++) {
        BLFront *blFront = nullptr;
        int bct;
        int idx0, idx1, idx2, bc1, ifc;
        double bc2;

        m_pElems[i].nconn = 3;
        nconn = 3;

        //*fin >> tmp >> idx0 >> idx1 >> idx2 >> ifc >> bct >> bc1 >> bc2;

        m_pElems[i].conn[0] = ele[4 * i + 0];
        m_pElems[i].conn[1] = ele[4 * i + 1];
        m_pElems[i].conn[2] = ele[4 * i + 2];
        m_pElems[i].topo = BLEntityTopology::TRIANGLE;

        m_pElems[i].igom = ele[4 * i + 3];

        // m_pElems[i].neig[0] = m_pElems[i].neig[1] = m_pElems[i].neig[2] = -1;
        bct = get<3>(file)[i];

        if (bct == BoundaryType::per) {
            for (j = 0; j < nconn; j++) {
                int idx = m_pElems[i].conn[j];
                if (per_point[idx] / 10 == 0) {
                    per_point[idx] += 10;
                }
            }
        }
#ifdef _CHECK_INTERSECTION
        // intersection
        if (true) {
            // m_pTriElm[intercnt * 3 + 0] = ele[4 * i + 0];
            // m_pTriElm[intercnt * 3 + 1] = ele[4 * i + 1];
            // m_pTriElm[intercnt * 3 + 2] = ele[4 * i + 2];
            std::array<int, 3> element;
            for (int k = 0; k < 3; k++) {
                element[k] = ele[4 * i + k];
            }

            for (int k = 0; k < 3; k++) {
                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 0]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 0]].coord[k]);

                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 1]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 1]].coord[k]);

                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 2]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 2]].coord[k]);
            }

            // int elmSize = intercnt;
            // if (m_pTriElm[elmSize * 3 + 0] > m_pTriElm[elmSize * 3 + 1]) {
            //	swap(m_pTriElm[elmSize * 3 + 0], m_pTriElm[elmSize * 3 + 1]);
            // }

            // if (m_pTriElm[elmSize * 3 + 1] > m_pTriElm[elmSize * 3 + 2]) {
            //	swap(m_pTriElm[elmSize * 3 + 1], m_pTriElm[elmSize * 3 + 2]);
            // }

            // if (m_pTriElm[elmSize * 3 + 0] > m_pTriElm[elmSize * 3 + 1]) {
            //	swap(m_pTriElm[elmSize * 3 + 0], m_pTriElm[elmSize * 3 + 1]);

            //}

            sort(element.begin(), element.end());
            int id = 0;
            if (bct != BoundaryType::symmetry && bct != BoundaryType::per) {
                id = m_TriElm.AddElem(element);
                m_vecData.push_back(id);
                l2g[id] = i;
            }
            // else{
            //	id = m_TriElm.AddElem(element);
            //	m_vecData_symm.push_back(id);
            //	l2g[id] = i;
            // }

            // intercnt++;
        }
#endif

        if (bct == BoundaryType::wall || bct == BoundaryType::match || bct == BoundaryType::adjacent) {
            if (bct == BoundaryType::adjacent) {
                for (j = 0; j < dim; j++) {
                    m_pNodes[m_pElems[i].conn[j]].badjacent = true;
                    BLNode *blNodTmp = (BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer;
                    blNodTmp->adjacent = true;
                }
            }
            blFront = new BLFront();
            // m_blFrontList->AddFront(blFront);
            int id = m_TriElm.GetSize() - 1;
            // set the neighbors of front node
            for (j = 0; j < nconn; j++) {
                int idx = m_pElems[i].conn[j];
                if (per_point[idx] % 10 == 0) {
                    per_point[idx] += 1;
                }
                BLNode *blNodTmp = (BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer;
                blNodTmp->AddNeigFronts(blFront);
                assert(blNodTmp->GetNodIdx() >= 0);
                blFront->AddBLFrontNods(j, blNodTmp);
            }

            blFront->SetLayerNum(0);
            CalFrontSize(blFront);

            ave_front_size += blFront->GetFrontSize();
            for (j = 0; j < nconn; j++) {
                point_manifold_size[m_pElems[i].conn[j]]++;
                point_size[m_pElems[i].conn[j]] += blFront->GetFrontSize();
            }
            front_count++;
            blFront->CalNormal(m_pNodes);

            for (j = 0; j < dim; j++) {
                if (bct == BoundaryType::wall) {
                    m_pNodes[m_pElems[i].conn[j]].bwall = true;
                }
            }
            if (bct == BoundaryType::match) {
                for (j = 0; j < nconn; j++) {
                    if (!m_pNodes[m_pElems[i].conn[j]].bwall) {
                        ((BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer)->SetStopFlag();
                    }
                }

            } else {
                for (j = 0; j < nconn; j++) {
                    ((BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer)->SetStopFlag(false);
                }
            }
#ifdef _CHECK_INTERSECTION
            blFront->SetTriIdx(m_TriElm.GetSize() - 1);

            blFront->SetSurfaceElmIdx(m_TriElm.GetSize() - 1);

            // end of intersection
#endif
        } else if (bct == BoundaryType::farfield) {
#if 1
            outbdry[3 * noutbdry + 0] = m_pElems[i].conn[0];
            outbdry[3 * noutbdry + 1] = m_pElems[i].conn[2];
            outbdry[3 * noutbdry + 2] = m_pElems[i].conn[1];

            if (m_pElems[i].conn[0] < 0 || m_pElems[i].conn[1] < 0 || m_pElems[i].conn[2] < 0) {
                throw std::runtime_error("try to input error");
            }
            noutbdry++;
#else
            SetElmDelete(i);
            nbdout++;
#endif

            for (j = 0; j < dim; j++) {
                m_pNodes[m_pElems[i].conn[j]].bfarfield = true;
                //	m_pNodes[m_pElems[i].conn[j]].pointer = nullptr;
            }
        } else if (bct == BoundaryType::symmetry || bct == BoundaryType::per) {
            std::array<int, 3> element;
            for (int k = 0; k < 3; k++) {
                element[k] = ele[4 * i + k];
            }

            symmface[m_pElems[i].igom].push_back(element);

            for (j = 0; j < dim; j++) {
                m_pNodes[m_pElems[i].conn[j]].bsysm = true;

                auto &vec = m_pNodes[m_pElems[i].conn[j]].isymfc;
                if (auto it = std::find(vec.begin(), vec.end(), m_pElems[i].igom); it == vec.end()) {
                    vec.push_back(m_pElems[i].igom);
                }

#ifndef _NEW_SYMM
                if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[0];
                } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[1];
                } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[2];
                }

                cnt++;
#endif
            }

            SetElmDelete(i);
        }

        m_pElems[i].pointer = (void *)blFront;
    }

    std::vector<int> per_idxs;
    std::vector<Eigen::Vector3d> per_coords;
    for (int i = 0; i < npt; i++) {
        if (per_point[i] == 11) // peroid boundary point
        {
            per_idxs.push_back(i);
            per_coords.push_back(Eigen::Vector3d(m_pNodes[i].coord));
        }
    }

    // perform one-to-one pairing among the periodic boundary points;
    std::vector<std::pair<int, int>> perpairs;
    for (int i = 0; i < per_coords.size(); i++) {
        Eigen::Vector3d np = cf.rotate_matrix * (per_coords[i]);
        np += cf.shift_vec;
        double mindis = 1e20;
        int minidx = -1;

        for (int j = 0; j < per_coords.size(); j++) {
            if (j != i) {

                double distance = (np - per_coords[j]).norm();
                if (mindis > distance && distance < 1e-2 * (np - per_coords[i]).norm()) {
                    mindis = distance;
                    minidx = j;
                }
            }
        }

        if (minidx >= 0) {

            int j = minidx;
            BLNode *f((BLNode *)m_pNodes[per_idxs[i]].pointer);
            BLNode *e((BLNode *)m_pNodes[per_idxs[j]].pointer);
            f->setPerNode(e, BLNode::PerType::Forward);
            e->setPerNode(f, BLNode::PerType::Reverse);
            int fid = f->GetNodIdx();
            int eid = e->GetNodIdx();
        }
        //	spdlog::info(f->GetNodIdx() << " " << e->GetNodIdx());
    }

    for (auto i : symmface) {
        int faceid = i.first;
        Eigen::MatrixXi F;
        igl::list_to_matrix(i.second, F);
        Eigen::MatrixXi NF, J;
        Eigen::MatrixXd NV;
        igl::remove_unreferenced(V, F, NV, NF, J);
        faceid2sp[faceid].init(NV, NF);
    }

    if (cf.iscompresslen) {
        // do step length modification
        DistanceCalculator dc;
        dc.ReadInput(file, cf);
        auto projection = dc.getHeightProjection();

        if (!projection.empty()) {
            for (int i = 0; i < npt; i++) {

                double full_length = (cf.step_len - cf.step_len * pow(cf.ratio, cf.layer_num - 1)) / (1 - cf.ratio);
                if (projection[i] >= 0 && projection[i] < full_length * 0.6) {
                    if (m_pNodes[i].bsysm > 0 || m_pNodes[i].bfarfield) {
                        continue;
                    }
                    auto node = ((BLNode *)m_pNodes[i].pointer);
                    if (node && !node->iso_stop) {

                        node->recommend_height = (-0.4 * (projection[i] - 1) * (projection[i] - 1) + 1) * cf.step_len;
                    }
                }
            }
        }
    }

#ifndef _NEW_SYMM
    if (cnt > 0) {
        m_sysValue /= cnt;
    }
#endif

#ifdef _CHECK_INTERSECTION
    // intersection
    double longest = 0;
    for (int k = 0; k < 3; k++) {
        longest = std::max(longest, box_max[k] - box_min[k]);
    }

    if (front_count) {
        ave_front_size /= front_count;
    } else {
        ave_front_size = longest / 100;
    }

    max_depth_ += int(log(m_TriElm.GetSize()) / log(10) / 1.2);
    // cout << max_depth_;
    m_ocAgent = new OctreeAgent(m_TriElm, m_pNodes);
    m_ocTree = new OCT::Octree(m_ocAgent, max_depth_);

    // cout << max_depth_;
    //   m_ocAgent_symm = new OctreeAgent(m_TriElm, m_pNodes);
    // m_ocTree_symm = new OCT::Octree(m_ocAgent_symm, max_depth_);

    for (int k = 0; k < 3; k++) {
        box_max[k] += ave_front_size * 30;
        box_min[k] -= ave_front_size * 30;
    }
    m_cbCube.upper.x = box_max.x;
    m_cbCube.upper.y = box_max.y;
    m_cbCube.upper.z = box_max.z;
    m_cbCube.lower.x = box_min.x;
    m_cbCube.lower.y = box_min.y;
    m_cbCube.lower.z = box_min.z;

    cf.min_volumn_eps = ave_front_size * cf.step_len * cf.step_len * 1.0 / 6.0 * 0.1;

    m_ocTree->buildOcTree(m_cbCube, m_vecData, max_obj_);
    // m_ocTree_symm->buildOcTree(m_cbCube, m_vecData_symm, max_obj_);

    model_centor = (m_cbCube.upper + m_cbCube.lower) / 2;
    // m_ocTree->saveOCTreeVectorVTK("octree.vtk");
    // BFS
    stack<TreeNode *> q;
    q.push(m_ocTree->getRootNode());
    vector<bool> unique(nelm, false);
    int order[8] = {0, 1, 3, 2, 6, 7, 5, 4};
    while (!q.empty()) {
        TreeNode *node = q.top();
        q.pop();
        if (node->getChild(0)) {
            for (int i = 0; i < 8; i++) {
                q.push(node->getChild(order[i]));
            }
        } else {
            for (auto i : node->order_data) {
                if ((!unique[l2g[i.id]]) && m_pElems[l2g[i.id]].pointer) {
                    m_blFrontList->AddFront((BLFront *)(m_pElems[l2g[i.id]].pointer));
                    unique[l2g[i.id]] = true;
                }
            }
        }
    }

    // build m_ocTree_symm

    // end of intersection
    m_blFrontList->RestoreFront();

    // while (m_blFrontList->HasNextFront()) {
    //	BLFront* blFront = m_blFrontList->GetNextFront();
    //	int idx=blFront->GetTriIdx();
    //	m_ocTree->insertPreProcess(idx);
    // }
#endif

#ifndef _NEW_SYMM
    UpdateSymnode();
#endif

    m_blFrontList->RestoreFront();
    // spdlog::info("noutbdry: {}\n", /*nbdout*/noutbdry);

#ifndef _NEW_SYMM
    int mtype;
    double vrtio;
    BLNode *bNod;
    for (i = 0; i < nboupt; i++) {
        fscanf(fin, "%d %d %lf\n", &tmp, &pidx, &vrtio);
        bNod = (BLNode *)m_pNodes[pidx - 1].pointer;
        bNod->SetHightRatio(vrtio - 0.2);
    }

    for (i = 0; i < nnmlt; i++) {
        fscanf(fin, "%d %d %d\n", &tmp, &pidx, &mtype);
        bNod = (BLNode *)m_pNodes[pidx - 1].pointer;
        bNod->SetNormalMethod(mtype);
    }

    for (i = 0; i < nbdrypt; i++) {
        fscanf(fin, "%d ", &m_pBdryPnt[i]);
        m_pBdryPnt[i] -= 1;
    }
#else
    // nbdry = 176;
    m_pSymBdrys = (Elem *)malloc(sizeof(Elem) * nsymbdry);
    if (m_pSymBdrys == nullptr) {
        spdlog::info("Can not allocate memory for Boundary points\n");
    }
    m_nAllocSymBdrys = nsymbdry;

#endif

    // read symmetry boundary point(2D)/segment(3D)

    CUBECOORD centor = (m_cbCube.lower + m_cbCube.upper) / 2;
    m_nSymBdrys = 0;
#ifdef _NEW_SYMM
    for (i = 0; i < nsymloop; i++) {
        int nbdry_i, iaxis;
        int fidx = -1;
        double fsymval;
        // fscanf(fin, "%d %d %d %lf\n", &nbdry_i, &iaxis, &fidx, &fsymval);
        *fin >> nbdry_i >> iaxis >> fidx;
        fsymval = std::get<4>(file)[i];
        m_sysPlane = (SymmetryPlane)(iaxis + 1);
        m_symVals[i] = fsymval;
        m_symFidx[i] = fidx;

        SymmetryPlaneV spv;
        double coord = centor[iaxis];
        spv.st = (SymmType)(iaxis + 1);
        spv.value = fsymval;

        if (coord < spv.value) {
            spv.st = (SymmType)(iaxis + 1 + 3);
        }

        sj.addPlane(spv);

        // m_sysValue = fsymval;

        for (j = 0; j < nbdry_i; j++) {
            int id1, id2;
            // fscanf(fin, "%d %d %d\n", &tmp, &id1, &id2);

            *fin >> tmp >> id1 >> id2;
            if (m_pNodes[id1 - 1].bsysm && m_pNodes[id2 - 1].bsysm) {

                // m_pNodes[id1 - 1].isymfc.push_back(fidx);		//temporally setting

                BLNode *blNod1 = (BLNode *)m_pNodes[id1 - 1].pointer;
                blNod1->SetBSys(true, iaxis);

                //	m_pNodes[id2 - 1].isymfc.push_back(fidx); // temporally setting

                BLNode *blNod2 = (BLNode *)m_pNodes[id2 - 1].pointer;
                blNod2->SetBSys(true, iaxis);

                m_pSymBdrys[m_nSymBdrys].nconn = 2;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].conn[1] = id2 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall && m_pNodes[id2 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].nconn = -2; // need to be deleted later;
                }

                m_nSymBdrys++;
            }
        }
    }

    for (i = 0; i < nbdrypt; i++) {
        // fscanf(fin, "%d ", &pidx);
        *fin >> pidx;
        m_pBdryPnt[i] = pidx;

        // set boundary point flag
        BLNode *pNod = (BLNode *)m_pNodes[pidx].pointer;
        // assert(pNod != nullptr);
        if (pNod) {
            pNod->SetBdryPt(true);
        }
    }
#else
    for (i = 0; i < nbdry; i++) {
        int id1, id2;
        if (m_blType == BLMType::blm2d) {
            fscanf(fin, "%d %d\n", &tmp, &id1);
            if (m_pNodes[id1 - 1].bsysm) {
                m_pSymBdrys[m_nSymBdrys].nconn = 1;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].reserved = 1; // need to be deleted later;
                } else {
                    m_pSymBdrys[m_nSymBdrys].reserved = 0;
                }

                m_nSymBdrys++;
            }
        } else if (m_blType == BLMType::blm3d) {
            fscanf(fin, "%d %d %d\n", &tmp, &id1, &id2);
            if (m_pNodes[id1 - 1].bsysm && m_pNodes[id2 - 1].bsysm) {
                m_pSymBdrys[m_nSymBdrys].nconn = 2;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].conn[1] = id2 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall && m_pNodes[id2 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].reserved = 1; // need to be deleted later;
                } else {
                    m_pSymBdrys[m_nSymBdrys].reserved = 0;
                }

                m_nSymBdrys++;
            }
        }
    }
#endif
    // fclose(fin);
    // fin = nullptr;
    delete fin;
    UpdateFrontNeig();
    // SetFrontSymm();
    SetBdryFront();

    // 删除包围盒节点的内存
    for (int i = 0; i < m_nNodes; i++) {

        if (m_pNodes[i].bfarfield) {
            if (m_pNodes[i].pointer) {
                // node_to_delete_.push_back((BLNode*)(m_pNodes[i].pointer));
                delete (BLNode *)(m_pNodes[i].pointer);
            }
            m_pNodes[i].pointer = nullptr;
        }
    }

    if (true) {

        if (get<3>(file)) {
            delete[] get<3>(file);
            get<3>(file) = nullptr;
        }

        if (get<1>(file)) {
            delete[] get<1>(file);
            get<1>(file) = nullptr;
        }

        if (get<2>(file)) {
            delete[] get<2>(file);
            get<2>(file) = nullptr;
        }
    }
    return 0;
}

int BLMesh::ReadBoundary(const INPUTFORMAT file, bool clear)
{

    //	spdlog::info("boufilename=%s", filename);
    // FILE *fin = nullptr;
    int i, j, npt, nelm, nboupt, nnmlt, tmp, nconn, nbdrypt, nbdry, dim, nbdout;
    int nsymloop, nsymbdry, pidx;
    BLVector box_max(-1e200, -1e200, -1e200);
    BLVector box_min(1e200, 1e200, 1e200);
    // 	int nbc, *ielm = nullptr;
    // 	double *pt = nullptr, *bcv = nullptr, *dnorm = nullptr, *u = nullptr;

    /*fin = fopen(filename, "r");
    if (!fin)
    {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }*/

    stringstream *fin = new stringstream(std::get<0>(file));
    nboupt = 0;
    nnmlt = 0;
    nbdrypt = 0;
    nbdry = 0;
    nsymloop = 0;

    // fscanf(fin, "%d %d %d %d %d\n", &nelm, &npt, &nsymloop, &nsymbdry, &nbdrypt);
    *fin >> nelm >> npt >> nsymloop >> nsymbdry >> nbdrypt;
    // nbdrypt = 0;

    m_nSurfNodes = npt;
    m_nSurfElems = nelm;

    m_vBdyFront.reserve(int(m_nSurfElems * 3)); // 预分配内存，加快速度
    m_nBdryPnt = nbdrypt;

    m_nSymLoop = nsymloop;

    m_nNodes = npt;
    m_nElems = nelm;

    // 	bcv = new double[2*nelm];
    // 	ielm = new int[2*nelm];
    // 	pt = new double[2*npt];
    // 	dnorm = new double[nelm];
    // 	u = new double[nelm];
    PreparePotential();

    m_pNodes = new MBLNode[npt];
    if (m_pNodes == nullptr) {
        spdlog::info("Can not allocate memory for Nodes\n");
        return OUT_OF_MEMORY;
    }
    m_nAllocNodes = m_nNodes;

    m_pElems = (Elem *)malloc(sizeof(Elem) * nelm);
    if (m_pElems == nullptr) {
        spdlog::info("Can not allocate memory for Elements\n");
        return OUT_OF_MEMORY;
    }
    m_nAllocElems = m_nElems;

    m_symVals = new double[m_nSymLoop];
    m_symFidx = new int[m_nSymLoop];

    if (nbdrypt > 0) {
        m_pBdryPnt = (int *)malloc(sizeof(int) * nbdrypt);
        if (m_pBdryPnt == nullptr) {
            spdlog::info("Can not allocate memory for Boundary points\n");
            return OUT_OF_MEMORY;
        }
    }

#ifndef _NEW_SYMM
    m_pSymBdrys = (Elem *)malloc(sizeof(Elem) * nbdry);
    if (m_pSymBdrys == nullptr) {
        spdlog::info("Can not allocate memory for Boundary points\n");
        return OUT_OF_MEMORY;
    }
    m_nAllocSymBdrys = nbdry;
#endif

#ifdef _CHECK_INTERSECTION
    // intersection data
    cout << "number of surface element = " << nelm << endl;
    m_TriElm.Reserve(10 * nelm);

    // m_nTriElm = nelm;
    m_nAllocTriElems = nelm;
    // end of intersection data
#endif

    // stores the boundary of the region needed to be meshed with unstructured mesh
    outbdry = new int[m_nSurfElems * 12];
    noutbdry = 0;
    nbdout = 0;
    dim = DIM3;

    // 预分配内存
#ifdef USE_MEMORY_POOL
    BLFront::preAllocate(nelm * 3);
    BLNode::preAllocate(npt * 3);
#endif

    for (i = 0; i < npt; i++) {

        m_pNodes[i].coord[0] = std::get<1>(file)[3 * i + 0];
        m_pNodes[i].coord[1] = std::get<1>(file)[3 * i + 1];
        m_pNodes[i].coord[2] = std::get<1>(file)[3 * i + 2];

        // create new front node
        BLNode *blNode = new BLNode(i);
        blNode->SetDecentID(i);
        blNode->SetLayerNum(0);

        m_pNodes[i].pointer = (void *)blNode;
        m_pNodes[i].bsysm = false;
        m_pNodes[i].bfarfield = false;
        m_pNodes[i].bwall = false;
        m_pNodes[i].reserved = -1;

#ifdef _CHECK_INTERSECTION
        // intersection

        // end of intersection data
#endif
    }

#ifdef _CHECK_INTERSECTION

    // end of intersection
#endif
    double ave_front_size = 0;
    int front_count = 0;
    int cnt = 0;
    int *ele = get<2>(file);
    vector<int> l2g(nelm, -1);

    std::vector<double> point_size(npt + 1, 0);
    std::vector<int> point_manifold_size(npt + 1, 0);

    for (i = 0; i < nelm; i++) {
        BLFront *blFront = nullptr;
        int bct;
        int idx0, idx1, idx2, bc1, ifc;
        double bc2;

        m_pElems[i].nconn = 3;
        nconn = 3;

        //*fin >> tmp >> idx0 >> idx1 >> idx2 >> ifc >> bct >> bc1 >> bc2;

        m_pElems[i].conn[0] = ele[4 * i + 0];
        m_pElems[i].conn[1] = ele[4 * i + 1];
        m_pElems[i].conn[2] = ele[4 * i + 2];

        m_pElems[i].topo = BLEntityTopology::TRIANGLE;

        m_pElems[i].igom = ele[4 * i + 3];
        ;

        // m_pElems[i].neig[0] = m_pElems[i].neig[1] = m_pElems[i].neig[2] = -1;
        bct = get<3>(file)[i];

#ifdef _CHECK_INTERSECTION
        // intersection
        if (true) {
            // m_pTriElm[intercnt * 3 + 0] = ele[4 * i + 0];
            // m_pTriElm[intercnt * 3 + 1] = ele[4 * i + 1];
            // m_pTriElm[intercnt * 3 + 2] = ele[4 * i + 2];
            std::array<int, 3> element;
            for (int k = 0; k < 3; k++) {
                element[k] = ele[4 * i + k];
            }

            for (int k = 0; k < 3; k++) {
                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 0]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 0]].coord[k]);

                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 1]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 1]].coord[k]);

                box_max[k] = std::max(box_max[k], m_pNodes[ele[4 * i + 2]].coord[k]);
                box_min[k] = std::min(box_min[k], m_pNodes[ele[4 * i + 2]].coord[k]);
            }

            // int elmSize = intercnt;
            // if (m_pTriElm[elmSize * 3 + 0] > m_pTriElm[elmSize * 3 + 1]) {
            //	swap(m_pTriElm[elmSize * 3 + 0], m_pTriElm[elmSize * 3 + 1]);
            // }

            // if (m_pTriElm[elmSize * 3 + 1] > m_pTriElm[elmSize * 3 + 2]) {
            //	swap(m_pTriElm[elmSize * 3 + 1], m_pTriElm[elmSize * 3 + 2]);
            // }

            // if (m_pTriElm[elmSize * 3 + 0] > m_pTriElm[elmSize * 3 + 1]) {
            //	swap(m_pTriElm[elmSize * 3 + 0], m_pTriElm[elmSize * 3 + 1]);

            //}
            sort(element.begin(), element.end());

            int id;
            if (bct != BoundaryType::symmetry && bct != BoundaryType::per) {
                id = m_TriElm.AddElem(element);
                m_vecData.push_back(id);

                l2g[id] = i;
            } else {
                id = m_TriElm.AddElem(element);
                m_vecData_symm.push_back(id);
                l2g[id] = i;
            }

            // intercnt++;
        }
#endif

        if (bct == BoundaryType::wall || bct == BoundaryType::match) {
            blFront = new BLFront();
            // m_blFrontList->AddFront(blFront);
            int id = m_TriElm.GetSize() - 1;
            // set the neighbors of front node
            for (j = 0; j < nconn; j++) {
                BLNode *blNodTmp = (BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer;
                blNodTmp->AddNeigFronts(blFront);
                assert(blNodTmp->GetNodIdx() >= 0);
                blFront->AddBLFrontNods(j, blNodTmp);
            }

            blFront->SetLayerNum(0);
            CalFrontSize(blFront);

            ave_front_size += blFront->GetFrontSize();
            for (j = 0; j < nconn; j++) {
                point_manifold_size[m_pElems[i].conn[j]]++;
                point_size[m_pElems[i].conn[j]] += blFront->GetFrontSize();
            }
            front_count++;
            blFront->CalNormal(m_pNodes);

            for (j = 0; j < dim; j++) {
                m_pNodes[m_pElems[i].conn[j]].bwall = true;
            }
            if (bct == BoundaryType::match) {
                for (j = 0; j < nconn; j++) {
                    ((BLNode *)m_pNodes[m_pElems[i].conn[j]].pointer)->SetStopFlag();
                }
            }
#ifdef _CHECK_INTERSECTION
            blFront->SetTriIdx(m_TriElm.GetSize() - 1);

            // end of intersection
#endif
        } else if (bct == BoundaryType::farfield) {
#if 1
            outbdry[3 * noutbdry + 0] = m_pElems[i].conn[0];
            outbdry[3 * noutbdry + 1] = m_pElems[i].conn[2];
            outbdry[3 * noutbdry + 2] = m_pElems[i].conn[1];

            if (m_pElems[i].conn[0] < 0 || m_pElems[i].conn[1] < 0 || m_pElems[i].conn[2] < 0) {
                throw std::runtime_error("try to input error");
            }
            noutbdry++;
#else
            SetElmDelete(i);
            nbdout++;
#endif

            for (j = 0; j < dim; j++) {
                m_pNodes[m_pElems[i].conn[j]].bfarfield = true;
            }
        } else if (bct == BoundaryType::symmetry || bct == BoundaryType::per) {
            for (j = 0; j < dim; j++) {
                m_pNodes[m_pElems[i].conn[j]].bsysm = true;

#ifndef _NEW_SYMM
                if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[0];
                } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[1];
                } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
                    m_sysValue += m_pNodes[m_pElems[i].conn[j]].coord[2];
                }

                cnt++;
#endif
            }

            SetElmDelete(i);
        }

        m_pElems[i].pointer = (void *)blFront;
    }

#ifndef _NEW_SYMM
    if (cnt > 0) {
        m_sysValue /= cnt;
    }
#endif

#ifdef _CHECK_INTERSECTION
    // intersection
    double longest = 0;
    for (int k = 0; k < 3; k++) {
        longest = std::max(longest, box_max[k] - box_min[k]);
    }

    if (front_count) {
        ave_front_size /= front_count;
    } else {
        ave_front_size = longest / 100;
    }

    max_depth_ += int(log(m_TriElm.GetSize()) / log(10) / 1.2);
    // cout << max_depth_;
    m_ocAgent = new OctreeAgent(m_TriElm, m_pNodes);
    m_ocTree = new OCT::Octree(m_ocAgent, max_depth_);

    // m_ocAgent_symm = new OctreeAgent(m_TriElm, m_pNodes);
    // m_ocTree_symm = new OCT::Octree(m_ocAgent, max_depth_);
    for (int k = 0; k < 3; k++) {
        box_max[k] += ave_front_size * 30;
        box_min[k] -= ave_front_size * 30;
    }
    m_cbCube.upper.x = box_max.x;
    m_cbCube.upper.y = box_max.y;
    m_cbCube.upper.z = box_max.z;
    m_cbCube.lower.x = box_min.x;
    m_cbCube.lower.y = box_min.y;
    m_cbCube.lower.z = box_min.z;

    cf.min_volumn_eps = ave_front_size * cf.step_len * cf.step_len * 1.0 / 6.0 * 0.1;

    m_ocTree->buildOcTree(m_cbCube, m_vecData, max_obj_);
    // m_ocTree_symm->buildOcTree(m_cbCube, m_vecData_symm, max_obj_);
    // m_ocTree->saveOCTreeVectorVTK("octree.vtk");
    // BFS
    stack<TreeNode *> q;
    q.push(m_ocTree->getRootNode());
    vector<bool> unique(nelm, false);
    int order[8] = {0, 1, 3, 2, 6, 7, 5, 4};
    while (!q.empty()) {
        TreeNode *node = q.top();
        q.pop();
        if (node->getChild(0)) {
            for (int i = 0; i < 8; i++) {
                q.push(node->getChild(order[i]));
            }
        } else {
            for (auto i : node->order_data) {
                if ((!unique[l2g[i.id]]) && m_pElems[l2g[i.id]].pointer) {
                    m_blFrontList->AddFront((BLFront *)(m_pElems[l2g[i.id]].pointer));
                    unique[l2g[i.id]] = true;
                }
            }
        }
    }

    // end of intersection
    m_blFrontList->RestoreFront();

    // while (m_blFrontList->HasNextFront()) {
    //	BLFront* blFront = m_blFrontList->GetNextFront();
    //	int idx=blFront->GetTriIdx();
    //	m_ocTree->insertPreProcess(idx);
    // }
#endif

#ifndef _NEW_SYMM
    UpdateSymnode();
#endif

    m_blFrontList->RestoreFront();
    // spdlog::info("noutbdry: {}\n", /*nbdout*/noutbdry);

#ifndef _NEW_SYMM
    int mtype;
    double vrtio;
    BLNode *bNod;
    for (i = 0; i < nboupt; i++) {
        fscanf(fin, "%d %d %lf\n", &tmp, &pidx, &vrtio);
        bNod = (BLNode *)m_pNodes[pidx - 1].pointer;
        bNod->SetHightRatio(vrtio - 0.2);
    }

    for (i = 0; i < nnmlt; i++) {
        fscanf(fin, "%d %d %d\n", &tmp, &pidx, &mtype);
        bNod = (BLNode *)m_pNodes[pidx - 1].pointer;
        bNod->SetNormalMethod(mtype);
    }

    for (i = 0; i < nbdrypt; i++) {
        fscanf(fin, "%d ", &m_pBdryPnt[i]);
        m_pBdryPnt[i] -= 1;
    }
#else
    // nbdry = 176;
    m_pSymBdrys = (Elem *)malloc(sizeof(Elem) * nsymbdry);
    if (m_pSymBdrys == nullptr) {
        spdlog::info("Can not allocate memory for Boundary points\n");
        return OUT_OF_MEMORY;
    }
    m_nAllocSymBdrys = nsymbdry;

#endif

    // read symmetry boundary point(2D)/segment(3D)

    m_nSymBdrys = 0;
#ifdef _NEW_SYMM
    for (i = 0; i < nsymloop; i++) {
        int nbdry_i, iaxis, fidx;
        double fsymval;
        // fscanf(fin, "%d %d %d %lf\n", &nbdry_i, &iaxis, &fidx, &fsymval);
        *fin >> nbdry_i >> iaxis >> fidx;
        fsymval = std::get<4>(file)[i];
        m_sysPlane = (SymmetryPlane)(iaxis + 1);
        m_symVals[i] = fsymval;
        m_symFidx[i] = fidx;
        // m_sysValue = fsymval;

        for (j = 0; j < nbdry_i; j++) {
            int id1, id2;
            // fscanf(fin, "%d %d %d\n", &tmp, &id1, &id2);

            *fin >> tmp >> id1 >> id2;
            if (m_pNodes[id1 - 1].bsysm && m_pNodes[id2 - 1].bsysm) {
                //  m_pNodes[id1 - 1].isymfc.push_back(fidx);
                // temporally setting

                BLNode *blNod1 = (BLNode *)m_pNodes[id1 - 1].pointer;
                blNod1->SetBSys(true, iaxis);

                //	m_pNodes[id2 - 1].isymfc.push_back(fidx); // temporally setting

                BLNode *blNod2 = (BLNode *)m_pNodes[id2 - 1].pointer;
                blNod2->SetBSys(true, iaxis);

                m_pSymBdrys[m_nSymBdrys].nconn = 2;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].conn[1] = id2 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall && m_pNodes[id2 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].nconn = -2; // need to be deleted later;
                }

                m_nSymBdrys++;
            }
        }
    }

    for (i = 0; i < nbdrypt; i++) {
        // fscanf(fin, "%d ", &pidx);
        *fin >> pidx;
        m_pBdryPnt[i] = pidx - 1;

        // set boundary point flag
        BLNode *pNod = (BLNode *)m_pNodes[pidx - 1].pointer;
        // assert(pNod != nullptr);
        if (pNod) {
            pNod->SetBdryPt(true);
        }
    }
#else
    for (i = 0; i < nbdry; i++) {
        int id1, id2;
        if (m_blType == BLMType::blm2d) {
            fscanf(fin, "%d %d\n", &tmp, &id1);
            if (m_pNodes[id1 - 1].bsysm) {
                m_pSymBdrys[m_nSymBdrys].nconn = 1;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].reserved = 1; // need to be deleted later;
                } else {
                    m_pSymBdrys[m_nSymBdrys].reserved = 0;
                }

                m_nSymBdrys++;
            }
        } else if (m_blType == BLMType::blm3d) {
            fscanf(fin, "%d %d %d\n", &tmp, &id1, &id2);
            if (m_pNodes[id1 - 1].bsysm && m_pNodes[id2 - 1].bsysm) {
                m_pSymBdrys[m_nSymBdrys].nconn = 2;
                m_pSymBdrys[m_nSymBdrys].conn[0] = id1 - 1;
                m_pSymBdrys[m_nSymBdrys].conn[1] = id2 - 1;
                m_pSymBdrys[m_nSymBdrys].pointer = nullptr;

                if (m_pNodes[id1 - 1].bwall && m_pNodes[id2 - 1].bwall) {
                    m_pSymBdrys[m_nSymBdrys].reserved = 1; // need to be deleted later;
                } else {
                    m_pSymBdrys[m_nSymBdrys].reserved = 0;
                }

                m_nSymBdrys++;
            }
        }
    }
#endif
    // fclose(fin);
    // fin = nullptr;
    delete fin;
    UpdateFrontNeig();
    // SetFrontSymm();
    SetBdryFront();

    // 删除包围盒节点的内存
    for (int i = 0; i < m_nNodes; i++) {

        if (m_pNodes[i].bfarfield) {
            if (m_pNodes[i].pointer) {
                // node_to_delete_.push_back((BLNode*)(m_pNodes[i].pointer));
                delete (BLNode *)(m_pNodes[i].pointer);
            }
            m_pNodes[i].pointer = nullptr;
        }
    }

    return BLM_SUCCESS;
}
void BLMesh::destroyNode() {}
void BLMesh::destroyFront()
{

    for (auto bdy_front_ptr : m_vBdyFront) {
        if (bdy_front_ptr->GetLowerFront()) {
            delete bdy_front_ptr->GetLowerFront();
        }
        delete bdy_front_ptr;
    }
    for (auto bdy_front_ptr : front_to_delete_) {
        bdy_front_ptr->RmAllNode();
        if (bdy_front_ptr->GetLowerFront()) {
            delete bdy_front_ptr->GetLowerFront();
        }
        delete bdy_front_ptr;
    }
    for (auto bdy_node_ptr : node_to_delete_) {
        delete bdy_node_ptr;
    }

    m_vBdyFront.clear();
    m_vBdyFront.shrink_to_fit();
    node_to_delete_.clear();
    front_to_delete_.clear();
    front_to_delete_.shrink_to_fit();
}
void BLMesh::CalZeroNorm()
{
    int nflt, nnod;
    BLFront *blFrnts[DIM3];
    BLNode *blnod[MAX_NCONN];

    set<BLFront *> virtualFront;
    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        if (blFront->GetNormal().magnitude() < 0.5) {
            virtualFront.insert(blFront);

            blFront->GetNodes(&nnod, blnod);
            // double dis[3] = {0.0};
            for (int j = 0; j < 3; j++) {
                double distance = 0;
                for (int k = 0; k < 3; k++) {
                    distance += pow(m_pNodes[blnod[j]->GetNodIdx()].coord[k] - m_pNodes[blnod[(j + 1) % 3]->GetNodIdx()].coord[k], 2);
                }
                if (distance < 1e-8) {
                    blnod[j]->SetVirtualFlag(true);
                    blnod[(j + 1) % 3]->SetVirtualFlag(true);

                    // cout << "123";
                }
            }
        }
    }
    m_blFrontList->RestoreFront();
    while (!virtualFront.empty()) {
        vector<BLFront *> toDel;
        for (auto &i : virtualFront) {
            i->GetNeigbourFronts(&nflt, blFrnts);
            BLVector n;
            for (int j = 0; j < nflt; j++) {
                n += blFrnts[j]->GetNormal();
            }
            n.normalize();
            if (i->GetNormal() * n > 0.999) {
                i->SetNormal(n);
                toDel.push_back(i);
            } else {
                i->SetNormal(n);
            }
        }
        for (auto i : toDel) {
            virtualFront.erase(i);
        }
    }
}

void BLMesh::setpntelm3d(int npt, int nElm, Elem *pElem, int **pntidx, int **pntelm)
{
    int i, j, start, *idx, *pelem;
    int npoints, nelm, dim;

    dim = DIM3;

    nelm = nElm;
    npoints = npt;

    idx = new int[npoints + 1];
    pelem = new int[nelm * dim];
    *pntidx = idx;
    *pntelm = pelem;

    for (i = 0; i <= npoints; i++) {
        idx[i] = 0;
    }

    /*----------------------------------------------------------------------------
    | Count the number of elements for each point:
    ----------------------------------------------------------------------------*/
    for (i = 0; i < nelm; i++) {
        for (j = 0; j < dim; j++) {
            idx[pElem[i].conn[j]]++;
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
        for (j = 0; j < dim; j++) {
            int pnt = pElem[i].conn[j];
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

void BLMesh::updateneig(int nElm, Elem *pElem, int *pntidx, int *pntelm, vector<std::array<int, 3>> &neigh)
{
    int i, j, k, nsurf, beg, end;
    int pidx1, pidx2, eidx;
    int dim = DIM3;

    /*遍历每个曲面建立网格单元邻接关系*/
    for (i = 0; i < nElm; i++) {
        for (k = 0; k < dim; k++) {
            if (neigh[i][k] == -1) {
                /*网格点(k+1)%dim的共享单元中找出k对应的邻接单元*/
                pidx1 = pElem[i].conn[(k + 1) % dim];
                // if(m_blType == BLMType::blm3d)
                pidx2 = pElem[i].conn[(k + 2) % dim];
                int pbeg = pntidx[pidx1];
                int pend = pntidx[pidx1 + 1];

                /*遍历共享网格点pidx1的网格单元*/
                for (int j1 = pbeg; j1 < pend; j1++) {
                    int peidx = pntelm[j1], lidx;

                    if (peidx != i && iselmedge(peidx, pidx1, pidx2, &lidx, pElem)) {

                        neigh[i][k] = peidx;
                        neigh[peidx][lidx] = i;
                        break;
                    }
                }
            }
        }
    }
}

bool BLMesh::iselmedge(int eidx, int idx1, int idx2, int *pidx, Elem *pElm)
{
    int i, p1, p2, dim = DIM3;

    for (i = 0; i < dim; i++) {
        p1 = pElm[eidx].conn[i];
        p2 = pElm[eidx].conn[(i + 1) % 3];
        *pidx = (i + 2) % dim; // 第三个点的局部编号

        if ((p1 == idx1 && p2 == idx2) || (p1 == idx2 && p2 == idx1)) {
            return true;
        }
    }
    return false;
}

void BLMesh::UpdateFrontNeig()
{
    int *pntidx, *pntelm, i, j, eidx;
    int dim = DIM3;
    vector<std::array<int, 3>> neigh(m_nSurfElems, std::array<int, 3>{-1, -1, -1});
    setpntelm3d(m_nSurfNodes, m_nSurfElems, m_pElems, &pntidx, &pntelm);
    updateneig(m_nSurfElems, m_pElems, pntidx, pntelm, neigh);

    m_pPntIdx = pntidx;
    m_pPntElm = pntelm;

    for (i = 0; i < m_nSurfElems; i++) {
        BLFront *blFlt = (BLFront *)m_pElems[i].pointer;

        if (blFlt != nullptr) {
            for (j = 0; j < dim; j++) {
                eidx = neigh[i][j]; // m_pElems[i].neig[j];
                if (eidx >= 0) {
                    BLFront *blTmp = (BLFront *)m_pElems[eidx].pointer;
                    if (blTmp != nullptr) {
                        blFlt->AddNeigbourFronts(blTmp);
                    }
                }
            }
        }
    }
}

void BLMesh::SetFrontSymm()
{
    int i, nNods, pidx, cnt, idx;
    BLNode *blNods[MAX_POINT_IN_ELEMENT];
    BLFront *blFront;

    while (m_blFrontList->HasNextFront()) {
        blFront = m_blFrontList->GetNextFront();
        blFront->GetNodes(&nNods, blNods);

        cnt = 0;
        for (i = 0; i < nNods; i++) {
            idx = blNods[i]->GetNodIdx();
            if (m_pNodes[idx].bsysm) {
                cnt++;
            } else {
                pidx = i;
            }
        }

        if (cnt == 2) {
            blFront->SetSymm(pidx);
        }
    }
    m_blFrontList->RestoreFront();
}

void BLMesh::AddSymmSegment(BLFront *blFront)
{
    int idx, conn[DIM3], pid1, pid2, nNods;
    BLNode *blNods[MAX_POINT_IN_ELEMENT];
    idx = blFront->GetSymm();

    blFront->GetNodes(&nNods, blNods);

    pid1 = blNods[(idx + 1) % DIM3]->GetNodIdx();
    pid2 = blNods[(idx + 2) % DIM3]->GetNodIdx();
    conn[0] = pid1;
    conn[1] = pid2;
    AddSymBdry(2, conn, BLEntityTopology::LINE);

    if (blNods[(idx + 1) % DIM3]->HasUpperNode()) {
        pid1 = blNods[(idx + 1) % DIM3]->GetNodIdx();
        pid2 = blNods[(idx + 1) % DIM3]->GetUpperNode()->GetNodIdx();
        conn[0] = pid1;
        conn[1] = pid2;
        AddSymBdry(2, conn, BLEntityTopology::LINE);
    }

    if (blNods[(idx + 2) % DIM3]->HasUpperNode()) {
        pid1 = blNods[(idx + 2) % DIM3]->GetNodIdx();
        pid2 = blNods[(idx + 2) % DIM3]->GetUpperNode()->GetNodIdx();
        conn[0] = pid1;
        conn[1] = pid2;
        AddSymBdry(2, conn, BLEntityTopology::LINE);
    }
}

void BLMesh::SetBdryFront()
{
    int i, j, idx1, idx2, idx3, beg, end, eidx;
    bool belm = false;

    for (i = 0; i < m_nSymBdrys; i++) {
        if (m_pSymBdrys[i].nconn < 0) {
            belm = false;
            idx1 = m_pSymBdrys[i].conn[0];

            idx2 = m_pSymBdrys[i].conn[1];

            beg = m_pPntIdx[idx1];
            end = m_pPntIdx[idx1 + 1];
            for (j = beg; j < end; j++) {
                eidx = m_pPntElm[j];
                for (int k = 0; k < DIM3; k++) {
                    idx3 = (m_pElems[eidx].conn[0] + m_pElems[eidx].conn[1] + m_pElems[eidx].conn[2]) - (idx1 + idx2);
                    if ((idx2 == m_pElems[eidx].conn[k]) && m_pNodes[idx3].bwall && (m_pElems[eidx].pointer != nullptr)) {
                        belm = true;
                        break;
                    }
                }
                if (belm) {
                    break;
                }
            }

            if (!belm) {
                throw(std::string("No matched element!"));
            }

            BLFront *blFront = (BLFront *)m_pElems[eidx].pointer;
            m_pSymBdrys[i].pointer = blFront;
            // 				if(i == 15)
            // 					spdlog::info("front idx: {} ({} {} {})\n", eidx, m_pElems[eidx].conn[0], m_pElems[eidx].conn[1],
            // m_pElems[eidx].conn[2]);
        }
    }
}
/*将一条直线加入到待选的边界中*/
void BLMesh::insertSymFace(pair<int, int> p) { m_mapSymline[p.first].insert(p.second); }

/* remove some sym line if tow node in symm face inner the symm face */
/*  1-------2----7---------
 *    \     /
 *	  \   /
 *       4
 *        \
 *          6
 *           \
 *
 *    remove (2,4)
 */
#ifdef _DEBUG
void BLMesh::RemoveSymLineByManifoldCretiria()
{
    std::map<int, set<int>> sym_lines;
    for (auto i : m_mapSymline) {
        for (auto j : i.second) {
            sym_lines[i.first].insert(j);
            sym_lines[j].insert(i.first);
        }
    }
    for (auto i : sym_lines) {
        if (i.second.size() > 2) {
            for (auto j : i.second) {
                if (sym_lines[j].size() > 2) {
                    m_mapSymline[min(i.first, j)].erase(max(i.first, j));
                }
            }
        }
    }
}
#endif
#ifdef _DEBUG
bool BLMesh::closeToPoint(int node_index, BLVector coord)
{
    double p = 0;
    for (int i = 0; i < 3; i++) {
        p += pow(m_pNodes[node_index].coord[i] - coord[i], 2);
    }
    return p < 1e-3;
}

#endif

// #define _COLOR_DEF
#ifdef APIFUNC
void BLMesh::SaveBLMeshAndFreeMemory(VM &v)
{
    MeshTopoStatics();
    int i, j, npt, nelm, sidx;
    npt = m_nNodes;
    v.ppdMNC = new double[3 * npt];
    for (auto i = 0; i < npt; i++) {
        for (int j = 0; j < 3; j++) {
            v.ppdMNC[3 * i + j] = m_pNodes[i].coord[j];
        }
    }
    if (m_pNodes) {
        delete[] (m_pNodes);
        m_pNodes = nullptr;
    }
    v.pnMN = npt;

    int sum = 0;
    int index = 0;
    m_nPrism = 0;
    m_nPyramid = 0;
    m_nTet = 0;
    for (int i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                m_nPrism++;
            }
            if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                if (generate_pyramid) {
                    m_nPyramid++;
                }
            }
            if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                m_nTet++;
            }
        }
    }
    v.pnME = m_nPrism + m_nPyramid + m_nTet;
    v.ppnMETp = new int[v.pnME];
    v.ppnMEFm = new int[6 * m_nPrism + 5 * m_nPyramid + 4 * m_nTet];
    index = 0;
    int index_p = 0;
    for (int i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                v.ppnMETp[index_p++] = m_pElems[i].topo;
                v.ppnMEFm[index++] = m_pElems[i].conn[0];
                v.ppnMEFm[index++] = m_pElems[i].conn[2];
                v.ppnMEFm[index++] = m_pElems[i].conn[1];
                v.ppnMEFm[index++] = m_pElems[i].conn[3];
                v.ppnMEFm[index++] = m_pElems[i].conn[5];
                v.ppnMEFm[index++] = m_pElems[i].conn[4];
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                if (generate_pyramid) {
                    v.ppnMETp[index_p++] = m_pElems[i].topo;
                    v.ppnMEFm[index++] = m_pElems[i].conn[0];
                    v.ppnMEFm[index++] = m_pElems[i].conn[2];
                    v.ppnMEFm[index++] = m_pElems[i].conn[1];
                    v.ppnMEFm[index++] = m_pElems[i].conn[3];
                    v.ppnMEFm[index++] = m_pElems[i].conn[4];
                }
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                v.ppnMETp[index_p++] = m_pElems[i].topo;
                v.ppnMEFm[index++] = m_pElems[i].conn[0];
                v.ppnMEFm[index++] = m_pElems[i].conn[1];
                v.ppnMEFm[index++] = m_pElems[i].conn[2];
                v.ppnMEFm[index++] = m_pElems[i].conn[3];
            }
        }
    }
    if (m_pElems) {
        free(m_pElems);
        m_pElems = nullptr;
    }

    return;
}
#endif
int BLMesh::SaveTetMesh(char *filename)
{
    int i, j, npt, nelm, sidx;
    FILE *fout = nullptr;

    npt = m_nNodes;

    MeshTopoStatics();

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }

    // vtk
    fprintf(fout, "# vtk DataFile Version 2.0\n");
    fprintf(fout, "boundary layer mesh\n");
    fprintf(fout, "ASCII\n");
    fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout, "POINTS %d double\n", npt);

    for (i = 0; i < npt; i++) {
        fprintf(fout, "%.10lf %.10lf %.10lf\n", m_pNodes[i].coord[0], m_pNodes[i].coord[1], m_pNodes[i].coord[2]);
    }

    // only output volume mesh
// #define TESTPYRA
#ifdef TESTPYRA
    nelm = m_nPyramid;
    fprintf(fout, "CELLS %d %d\n", nelm, m_nPyramid * 6);
#else
    nelm = m_nTet;
    fprintf(fout, "CELLS %d %d\n", nelm, m_nTet * 5);
#endif

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {

            if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fprintf(fout, "4 %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[1], m_pElems[i].conn[2], m_pElems[i].conn[3]);
            }
        }
    }

    fprintf(fout, "CELL_TYPES %d\n", nelm);

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
#ifdef TESTPYRA
            if (m_pElems[i].topo == EntityTopology::PYRAMID) {
                fprintf(fout, "%d\n", 14);
            }
#else

            if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fprintf(fout, "%d\n", 10);
            }
#endif
        }
    }

    fclose(fout);
    fout = nullptr;

    return BLM_SUCCESS;
}
#ifdef _DEBUG
void BLMesh::SaveCurrentLayerBLMesh(std::string filename, int layer)
{
    int i, j, npt, nelm, sidx;

    npt = m_nNodes;

    MeshTopoStatics();
    map<int, int> m;
    map<int, int> rm;
    set<int> valid;
    std::ofstream fout(filename);
    int count = 0;
    int c1 = 0, c2 = 0, c3 = 0;
    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i) && layer == m_pElems[i].iLayer) {
            valid.insert(i);
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                c1++;
                for (int j = 0; j < 6; j++) {
                    if (m.find(m_pElems[i].conn[j]) == m.end()) {
                        m[m_pElems[i].conn[j]] = count;
                        rm[count] = m_pElems[i].conn[j];
                        count++;
                    }
                }
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                c2++;
                for (int j = 0; j < 5; j++) {
                    if (m.find(m_pElems[i].conn[j]) == m.end()) {
                        m[m_pElems[i].conn[j]] = count;
                        rm[count] = m_pElems[i].conn[j];
                        count++;
                    }
                }
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                c3++;
                for (int j = 0; j < 4; j++) {
                    if (m.find(m_pElems[i].conn[j]) == m.end()) {
                        m[m_pElems[i].conn[j]] = count;
                        rm[count] = m_pElems[i].conn[j];
                        count++;
                    }
                }
            }
        }
    }

    // vtk
    fout << "# vtk DataFile Version 2.0" << endl;
    fout << "boundary layer mesh" << endl;
    fout << "ASCII" << endl;
    fout << "DATASET UNSTRUCTURED_GRID" << endl;
    fout << "POINTS " << rm.size() << " double" << endl;

    for (auto j : rm) {
        int i = j.second;
        fout << m_pNodes[i].coord[0] << " " << m_pNodes[i].coord[1] << " " << m_pNodes[i].coord[2] << endl;
    }

    // only output volume mesh
    // #define TESTPYRA
    nelm = valid.size();
    fout << "CELLS " << nelm << " " << c1 * 7 + c3 * 5 + c2 * 6 << endl;

    for (auto i : valid) {
        if (!IsElmDelete(i)) {

#ifdef TESTPYRA
            if (m_pElems[i].topo == EntityTopology::PYRAMID) {
                fprintf(fout, "5 %d %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[2], m_pElems[i].conn[1], m_pElems[i].conn[3],
                        m_pElems[i].conn[4]);
            }
#else
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                fout << "6 " << m[m_pElems[i].conn[0]] << " " << m[m_pElems[i].conn[2]] << " " << m[m_pElems[i].conn[1]] << " "
                     << m[m_pElems[i].conn[3]] << " " << m[m_pElems[i].conn[5]] << " " << m[m_pElems[i].conn[4]] << endl;
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                fout << "5 " << m[m_pElems[i].conn[0]] << " " << m[m_pElems[i].conn[2]] << " " << m[m_pElems[i].conn[1]] << " "
                     << m[m_pElems[i].conn[3]] << " " << m[m_pElems[i].conn[4]] << endl;
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fout << "4 " << m[m_pElems[i].conn[0]] << " " << m[m_pElems[i].conn[1]] << " " << m[m_pElems[i].conn[2]] << " "
                     << m[m_pElems[i].conn[3]] << endl;
            }
#endif;
        }
    }

    fout << "CELL_TYPES " << nelm << endl;

    for (auto i : valid) {
        if (!IsElmDelete(i)) {
#ifdef TESTPYRA
            if (m_pElems[i].topo == EntityTopology::PYRAMID) {
                fprintf(fout, "%d\n", 14);
            }
#else

            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                fout << "13" << endl;
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                fout << "14" << endl;
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fout << "10" << endl;
            }
#endif
        }
    }

    fout.close();

    return;
}
#endif
int BLMesh::SaveBLMesh(char *filename)
{

    int i, j, npt, nelm, sidx;
    FILE *fout = nullptr;

    npt = m_nNodes;

    MeshTopoStatics();

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }

    // vtk
    fprintf(fout, "# vtk DataFile Version 2.0\n");
    fprintf(fout, "boundary layer mesh\n");
    fprintf(fout, "ASCII\n");
    fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout, "POINTS %d double\n", npt);

    for (i = 0; i < npt; i++) {
        fprintf(fout, "%.10lf %.10lf %.10lf\n", m_pNodes[i].coord[0], m_pNodes[i].coord[1], m_pNodes[i].coord[2]);
    }

    // only output volume mesh
// #define TESTPYRA
#ifdef TESTPYRA
    nelm = m_nPyramid;
    fprintf(fout, "CELLS %d %d\n", nelm, m_nPyramid * 6);
#else
    nelm = m_nPrism + m_nPyramid + m_nTet;
    if (!generate_pyramid) {
        nelm = m_nPrism + m_nTet;
    }
    if (generate_pyramid) {
        fprintf(fout, "CELLS %d %d\n", nelm, m_nPrism * 7 + m_nTet * 5 + m_nPyramid * 6);
    } else {
        fprintf(fout, "CELLS %d %d\n", nelm, m_nPrism * 7 + m_nTet * 5);
    }
#endif

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {

#ifdef TESTPYRA
            if (m_pElems[i].topo == EntityTopology::PYRAMID) {
                fprintf(fout, "5 %d %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[2], m_pElems[i].conn[1], m_pElems[i].conn[3],
                        m_pElems[i].conn[4]);
            }
#else
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                fprintf(fout, "6 %d %d %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[2], m_pElems[i].conn[1], m_pElems[i].conn[3],
                        m_pElems[i].conn[5], m_pElems[i].conn[4]);
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                if (generate_pyramid) {
                    fprintf(fout, "5 %d %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[2], m_pElems[i].conn[1], m_pElems[i].conn[3],
                            m_pElems[i].conn[4]);
                }
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fprintf(fout, "4 %d %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[1], m_pElems[i].conn[2], m_pElems[i].conn[3]);
            }
#endif;
        }
    }

    fprintf(fout, "CELL_TYPES %d\n", nelm);

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
#ifdef TESTPYRA
            if (m_pElems[i].topo == EntityTopology::PYRAMID) {
                fprintf(fout, "%d\n", 14);
            }
#else

            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                fprintf(fout, "%d\n", 13);
            } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
                if (generate_pyramid) {
                    fprintf(fout, "%d\n", 14);
                }
            } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
                fprintf(fout, "%d\n", 10);
            }
#endif
        }
    }

    fclose(fout);
    fout = nullptr;

    return BLM_SUCCESS;
}

int BLMesh::SaveBdry(VM &mesh)
{
    int i, j, nelm, cnt = 0;
    FILE *fout = nullptr;

    set<int> sym_face;
    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {

            if (m_pElems[i].topo == BLEntityTopology::TRIANGLE) {
                cnt++;
            } else if (m_pElems[i].topo == BLEntityTopology::QUADRILATERAL) {
                cnt++;
            }
        }
    }
    mesh.num_boundary_face = cnt;
    mesh.boundary_mesh = new int[4 * cnt];
    mesh.boundary_face = new int[cnt];
    cnt = 0;
    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {

            if (m_pElems[i].topo == BLEntityTopology::TRIANGLE) {
                for (int k = 0; k < 3; k++) {
                    mesh.boundary_mesh[4 * cnt + k] = m_pElems[i].conn[k];
                    for (int l = 0; l < 3; l++) {}
                }
                mesh.boundary_mesh[4 * cnt + 3] = -1;
                mesh.boundary_face[cnt] = m_pElems[i].igom;
                cnt++;
            } else if (m_pElems[i].topo == BLEntityTopology::QUADRILATERAL) {
                for (int k = 0; k < 4; k++) {
                    mesh.boundary_mesh[4 * cnt + k] = m_pElems[i].conn[k];
                }
                sym_face.insert(m_pElems[i].igom);
                mesh.boundary_face[cnt] = m_pElems[i].igom;
                cnt++;
            }
        }
    }
    // Try to uniform unit winding in_symm
    int wind_count = 0;
    for (auto sym_id : sym_face) {
        for (int i = 0; i < cnt; i++) {
            if (mesh.boundary_face[i] == sym_id) {
                BLVector node[3];
                for (int k = 0; k < 3; k++) {
                    for (int j = 0; j < 3; j++) {
                        node[k][j] = m_pNodes[mesh.boundary_mesh[4 * i + k]].coord[j];
                    }
                }
                BLVector normal = (node[1] - node[0]) ^ (node[2] - node[0]);
                BLVector directon = model_centor - node[0];
                //	cout << model_centor<< normal<< directon << normal*directon << endl;
                if (normal * directon > 0) {
                    swap(mesh.boundary_mesh[4 * i + 0], mesh.boundary_mesh[4 * i + 2]);
                    wind_count++;
                } else {
                }
            }
        }
    }

    spdlog::info("wind_count={}", wind_count);
    return BLM_SUCCESS;
}

int BLMesh::SaveBdry(char *filename)
{
    int i, j, nelm, cnt;
    FILE *fout = nullptr;

    MeshTopoStatics();

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }

    cnt = 0;
    // only output volume mesh
    nelm = m_nTri + m_nQuad;

    fprintf(fout, "%d %d %d\n", nelm, m_nTri, m_nQuad);

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::TRIANGLE) {
                fprintf(fout, "%d 3 %d %d %d %d\n", cnt + 1, m_pElems[i].conn[0] + 1, m_pElems[i].conn[1] + 1, m_pElems[i].conn[2] + 1,
                        m_pElems[i].igom);

                cnt++;
            } else if (m_pElems[i].topo == BLEntityTopology::QUADRILATERAL) {
                fprintf(fout, "%d 4 %d %d %d %d %d\n", cnt + 1, m_pElems[i].conn[0] + 1, m_pElems[i].conn[1] + 1, m_pElems[i].conn[2] + 1,
                        m_pElems[i].conn[3] + 1, m_pElems[i].igom);

                cnt++;
            }
        }
    }

    fclose(fout);
    fout = nullptr;

    return BLM_SUCCESS;
}

int BLMesh::SaveSurfGrid(char *filename)
{
    int i, j, npt, nelm, sidx;
    FILE *fout = nullptr;

    npt = m_nNodes;
    // nelm = m_nElems;
    nelm = m_nSurfElems;

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }

    // vtk
    fprintf(fout, "# vtk DataFile Version 2.0\n");
    fprintf(fout, "boundary layer mesh\n");
    fprintf(fout, "ASCII\n");
    fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout, "POINTS %d double\n", npt);

    for (i = 0; i < npt; i++) {
        fprintf(fout, "%f %f %f\n", m_pNodes[i].coord[0], m_pNodes[i].coord[1], m_pNodes[i].coord[2]);
    }

    // sidx = m_nSurfElems;
    sidx = 0;

    // sidx = 2914;
    // nelm = sidx + 100;
    fprintf(fout, "CELLS %d %d\n", nelm - sidx, (nelm - sidx) * 4);
    for (i = 0; i < nelm - sidx; i++)
    // for (i=sidx; i<nelm; i++)
    {
        fprintf(fout, "3 %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[1], m_pElems[i].conn[2]);
    }

    fprintf(fout, "CELL_TYPES %d\n", nelm - sidx);

    for (i = 0; i < nelm - sidx; i++) {
        fprintf(fout, "%d\n", 5);
    }

#if 1
    // for testing
    fprintf(fout, "POINT_DATA %d\n", npt);
    // 		fprintf(fout, "SCALARS fixed double\n");
    // 		fprintf(fout, "LOOKUP_TABLE default\n");
    // 		for (i=0; i<npt; i++)
    // 		{
    // 			fprintf(fout, "%lf\n", m_pNodes[i].uvalue);
    // 		}
    fprintf(fout, "NORMALS node_normals double\n");
    for (i = 0; i < npt; i++) {
        BLNode *blNod = (BLNode *)m_pNodes[i].pointer;
        fprintf(fout, "%lf %lf %lf\n", blNod->GetNormal().x, blNod->GetNormal().y, blNod->GetNormal().z);
    }
#endif
#if 0
		//for testing
		fprintf(fout, "CELL_DATA %d\n", nelm);
		// 		fprintf(fout, "SCALARS fixed double\n");
		// 		fprintf(fout, "LOOKUP_TABLE default\n");
		// 		for (i=0; i<npt; i++)
		// 		{
		// 			fprintf(fout, "%lf\n", m_pNodes[i].uvalue);
		// 		}
		fprintf(fout, "NORMALS face_normals double\n");
		for (i = 0; i < nelm; i++)
		{
			BLVector norm = m_pElems[i].norm;
			fprintf(fout, "%lf %lf %lf\n", norm.x, norm.y, norm.z);
		}
#endif

    fclose(fout);
    fout = nullptr;

    return BLM_SUCCESS;
}

// require from zjj
int BLMesh::SaveSurfGrid2(char *filename)
{
    // return 0;
    int i, j, npt, nelm, sidx;
    FILE *fout = nullptr;

    npt = m_nSurfNodes;
    // nelm = m_nElems;
    nelm = m_nSurfElems;

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", filename);
        return OPEN_FILE_FAILD;
    }

    // vtk
    fprintf(fout, "# vtk DataFile Version 3.0\n");
    fprintf(fout, "Really cool data\n");
    fprintf(fout, "ASCII\n");
    fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout, "POINTS %d double\n", npt);

    for (i = 0; i < npt; i++) {
        fprintf(fout, "%f %f %f\n", m_pNodes[i].coord[0], m_pNodes[i].coord[1], m_pNodes[i].coord[2]);
    }

    sidx = m_nSurfElems;
    // sidx = 0;
    // sidx = 2914;
    // nelm = sidx + 100;
    fprintf(fout, "CELLS %d %d\n", sidx, sidx * 4);
    for (i = 0; i < sidx; i++)
    // for (i=sidx; i<nelm; i++)
    {
        fprintf(fout, "3 %d %d %d\n", m_pElems[i].conn[0], m_pElems[i].conn[1], m_pElems[i].conn[2]);
    }

    fprintf(fout, "CELL_TYPES %d\n", sidx);

    for (i = 0; i < sidx; i++) {
        fprintf(fout, "%d\n", 5);
    }

#if 1
    // for testing
    fprintf(fout, "POINT_DATA %d\n", npt);
    // fprintf(fout, "SCALARS pointtype double\n");
    // fprintf(fout, "LOOKUP_TABLE default\n");
    // for (i = 0; i < npt; i++)
    //{

    //	fprintf(fout, "%d\n", m_pNodes[i].pointType < 0 ? 0 : m_pNodes[i].pointType);
    //}
    // fprintf(fout, "POINT_DATA %d\n", npt);
    // fprintf(fout, "SCALARS internum double\n");
    // fprintf(fout, "LOOKUP_TABLE default\n");
    // for (i = 0; i < npt; i++)
    //{
    //	fprintf(fout, "%d\n", m_pNodes[i].numOfInter < 0 ? 0 : m_pNodes[i].numOfInter);
    //}
    fprintf(fout, "SCALARS beita double\n");
    fprintf(fout, "LOOKUP_TABLE default\n");
    for (i = 0; i < npt; i++) {
        // fprintf(fout, "%lf\n", m_pNodes[i].distance);
    }
#endif
#if 0
	fprintf(fout, "NORMALS vectors double\n");
	for (i = 0; i < npt; i++)
	{
		BLNode *blNod = (BLNode *)m_pNodes[i].pointer;
		BLVector norm;
		if (blNod)
			norm = blNod->GetNormal();
		//norm.normalize();
		if (norm.magnitude() < 0.5)
			fprintf(fout, "%lf %lf %lf\n", 1.0f, 0.0f, 0.0f);
		else
			fprintf(fout, "%lf %lf %lf\n", norm.x, norm.y, norm.z);
	}
	fprintf(fout, "NORMALS upper_vector double\n");
	for (i = 0; i < npt; i++)
	{
		BLNode *blNod = (BLNode *)m_pNodes[i].pointer;
		BLVector norm;
		if (blNod)
			if (blNod->GetUpperNode())
				norm = blNod->GetUpperNode()->GetNormal();
		//norm.normalize();
		if (norm.magnitude() < 0.5)
			fprintf(fout, "%lf %lf %lf\n", 1.0f, 0.0f, 0.0f);
		else
			fprintf(fout, "%lf %lf %lf\n", norm.x, norm.y, norm.z);
	}

#endif

    fclose(fout);
    fout = nullptr;

    return BLM_SUCCESS;
}

void BLMesh::MeshTopoStatics()
{
    int i = 0, nTri = 0, nQuad = 0, nTet = 0, nPrism = 0, nPyramid = 0;
    BLEntityTopology mshtopo;

    for (i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            mshtopo = m_pElems[i].topo;
            switch (mshtopo) {
                case BLEntityTopology::TRIANGLE:
                    nTri++;
                    break;

                case BLEntityTopology::QUADRILATERAL:
                    nQuad++;
                    break;

                case BLEntityTopology::TETRAHEDRON:
                    nTet++;
                    break;

                case BLEntityTopology::PRISM:
                    nPrism++;
                    break;

                case BLEntityTopology::PYRAMID:
                    nPyramid++;
                    break;

                default:
                    break;
            }
        }
    }

    m_nTri = nTri;
    m_nQuad = nQuad;
    m_nTet = nTet;
    m_nPrism = nPrism;
    m_nPyramid = nPyramid;
}

int BLMesh::VTK2Pls(char *filename)
{
    int i, j, npt, nelm, idx;
    FILE *fout = nullptr, *fin = nullptr;
    char filein[256], fileout[256], strtmp[256];
    double *pt;

    memset(filein, 0, strlen(filein));
    memset(fileout, 0, strlen(fileout));

    memcpy(filein, filename, strlen(filename));
    memcpy(fileout, filename, strlen(filename));

    strcat(filein, ".vtk");
    strcat(fileout, ".pl3");

    fin = fopen(filein, "r");
    if (fin == nullptr) {
        spdlog::info("Can not open file %s\n", filein);
        return OPEN_FILE_FAILD;
    }

    fout = fopen(fileout, "w");
    if (fout == nullptr) {
        spdlog::info("Can not open file %s\n", fileout);
        return OPEN_FILE_FAILD;
    }

    // vtk
    fscanf(fin, "%s %s %s %s %s\n", strtmp, strtmp, strtmp, strtmp, strtmp);
    fscanf(fin, "%s %s %s\n", strtmp, strtmp, strtmp);
    fscanf(fin, "%s\n", strtmp);
    fscanf(fin, "%s %s\n", strtmp, strtmp);
    fscanf(fin, "%s %d %s", strtmp, &npt, strtmp);

    pt = new double[3 * npt];

    for (i = 0; i < npt; i++) {
        fscanf(fin, "%f %f %f\n", &pt[i * 3 + 0], &pt[i * 3 + 1], &pt[i * 3 + 2]);
    }

    fscanf(fin, "%s %d %s\n", strtmp, &nelm, strtmp);

    fprintf(fout, "%d %d 0", nelm, npt);
    for (i = 0; i < npt; i++) {
        fprintf(fout, "%d %f %f %f\n", i + 1, pt[i * 3 + 0], pt[i * 3 + 1], pt[i * 3 + 2]);
    }

    for (i = 0; i < nelm; i++) {
        int p1, p2, p3, p4, p5, p6;
        fscanf(fin, "%d %d %d %d %d %d %d\n", &idx, &p1, &p2, &p3, &p4, &p5, &p6);

        fprintf(fout, "%d %d %d %d %d %d %d\n", i + 1, p1, p2, p3, p4, p5, p6);
    }

    fclose(fin);
    fclose(fout);
    fout = nullptr;
    fin = nullptr;

    return BLM_SUCCESS;
}

void BLMesh::PreparePotential()
{
    // m_pPt = new double[3 * m_nSurfNodes];
    // m_pPt = new double[3*m_nNodes];
#ifdef _COLOR_DEF
    m_pBc = new double[4 * m_nSurfElems];
#else
#endif
    //	m_pBElm = new int[MAX_NCONN*m_nSurfElems];
    //	m_pU = new double[m_nSurfElems];
    //	m_pNorm = new double[MAX_NORMAL_COMPONENT*m_nSurfElems];
}

int BLMesh::AddNode(BLVector pnt, double space)
{
    MBLNode *pNewNodes = nullptr;
    int i, nAlloc = 0, nodeSize;

    double d = pow(647.001 - pnt[0], 2) + pow(-921.032 - pnt[1], 2) + pow(-6.91823 - pnt[2], 2);
    if (d < 1e-3) {
        cout << "find the point!" << endl;
    }
    nodeSize = m_nNodes;

    // check memory for nodes
    if (nodeSize >= m_nAllocNodes) {
#ifdef _OUTPUT_LOG
        spdlog::info("try to reallocate memory for nodes\n");
#endif
        nAlloc = m_nAllocNodes;
        do {
            nAlloc += m_nSurfNodes + 0.5 * nAlloc;
        } while (nAlloc <= nodeSize);

        pNewNodes = new MBLNode[nAlloc];
        std::copy(m_pNodes, m_pNodes + m_nNodes, pNewNodes);
        delete[] m_pNodes;
        if (!pNewNodes) {
            spdlog::info("cannot reallocate memory!\n");
            return OUT_OF_MEMORY;
        }

        m_pNodes = pNewNodes;
        m_nAllocNodes = nAlloc;

#ifdef _CHECK_INTERSECTION
        // intersection
        if (m_ocAgent) {
            m_ocAgent->setNod(m_pNodes);
        }
        // m_ocAgent_symm->setNod(m_pNodes);
        // end of intersection
#endif
    }

    m_pNodes[nodeSize].coord[0] = pnt.x;
    m_pNodes[nodeSize].coord[1] = pnt.y;
    m_pNodes[nodeSize].coord[2] = pnt.z;
    // if (abs(m_pNodes[nodeSize].coord[0] - 369.453) < 5e-3&&abs(m_pNodes[nodeSize].coord[2] - 105.134) < 5e-3) {
    //	spdlog::info("=================nodeid=" << nodeSize + 1);
    //	spdlog::info(m_pNodes[nodeSize].coord[2]);
    // }
    m_pNodes[nodeSize].reserved = -1;
    m_nNodes++;

    return nodeSize;
}

int inline BLMesh::AddElem(int nconn, int *conn, BLEntityTopology topu, int igeom)
{
    Elem *pNewElems = nullptr;
    int i, nAlloc = 0, elmSize;

    elmSize = m_nElems;

    // check memory for elements
    if (elmSize >= m_nAllocElems) {
#ifdef _DEBUG
        spdlog::info("try to reallocate memory for elements\n");
#endif
        nAlloc = m_nAllocElems;
        do {
            nAlloc += 2 * m_nSurfElems + 0.3 * nAlloc;
        } while (nAlloc <= elmSize);

        pNewElems = (Elem *)realloc(m_pElems, sizeof(Elem) * nAlloc);

        if (nAlloc > m_nAllocElems) {
            memset(&pNewElems[m_nAllocElems], 0, sizeof(Elem) * (nAlloc - m_nAllocElems));
        }

        m_pElems = pNewElems;
        m_nAllocElems = nAlloc;
    }
    for (i = 0; i < nconn; i++) {
        m_pElems[elmSize].conn[i] = conn[i];
        assert(conn[i] > -1);
    }
    if (m_pElems[elmSize].conn[0] < 0) {
        spdlog::info(m_pElems[i].conn[0]);
    }
    m_pElems[elmSize].nconn = nconn;
    m_pElems[elmSize].topo = topu;
    m_pElems[elmSize].igom = igeom;
#ifdef _DEBUG
    m_pElems[elmSize].iLayer = m_nCurrLayer;
#endif
    m_nElems++;

    return elmSize;
}

int BLMesh::AddSymBdry(int nconn, int *conn, int topu)
{
    Elem *pNewElems = nullptr;
    int i, nAlloc = 0, elmSize;

    elmSize = m_nSymBdrys;

    // check memory for elements
    if (elmSize >= m_nAllocSymBdrys) {
#ifdef _OUTPUT_LOG
        spdlog::info("try to reallocate memory for symmetry boundaries\n");
#endif
        nAlloc = m_nAllocSymBdrys;
        do {
            nAlloc += m_nSymBdrys / 10 + 100;
        } while (nAlloc <= elmSize);

        pNewElems = (Elem *)realloc(m_pSymBdrys, sizeof(Elem) * nAlloc);
        if (!pNewElems) {
            spdlog::info("cannot reallocate memory!\n");
            return OUT_OF_MEMORY;
        }

        if (nAlloc > m_nAllocSymBdrys) {
            memset(&pNewElems[m_nAllocSymBdrys], 0, sizeof(Elem) * (nAlloc - m_nAllocSymBdrys));
        }

        m_pSymBdrys = pNewElems;
        m_nAllocSymBdrys = nAlloc;
    }
    for (i = 0; i < nconn; i++) {
        m_pSymBdrys[elmSize].conn[i] = conn[i];
    }

    m_pSymBdrys[elmSize].nconn = nconn;

    m_nSymBdrys++;

    return elmSize;
}

INDEX_TYPE BLMesh::AddTriElem(std::array<int, 3> &conn)
{
    if (conn[0] > conn[1]) {
        swap(conn[0], conn[1]);
    }
    if (conn[1] > conn[2]) {
        swap(conn[1], conn[2]);
        if (conn[0] > conn[1]) {
            swap(conn[0], conn[1]);
        }
    }
    return m_TriElm.AddElem(conn);
}
INDEX_TYPE BLMesh::AddTriElem(int nconn, int *conn)
{
    static std::array<int, 3> element;

    for (int i = 0; i < 3; i++) {

        element[i] = conn[i];
    }
    if (element[0] > element[1]) {
        swap(element[0], element[1]);
    }
    if (element[1] > element[2]) {
        swap(element[1], element[2]);
        if (element[0] > element[1]) {
            swap(element[0], element[1]);
        }
    }
    // cout << m_TriElm.GetSize() << endl;;
    return m_TriElm.AddElem(element);
}

bool BLMesh::CheckTri(int nconn, int *conn, double coord[3][3])
{
    int i, j, k, cnt = 0, idx, nds = 0;
    double eps = 1.0e-6;
    for (i = 0; i < nconn; i++) {
        idx = conn[i];
        for (j = 0; j < 3; j++) {
            cnt = 0;
            for (k = 0; k < 3; k++) {
                if (fabs(m_pNodes[idx].coord[k] - coord[j][k]) < eps) {
                    cnt++;
                }
            }
            if (cnt == 3) {
                nds++;
                break;
            }
        }
    }

    if (nds == 3) {
        return true;
    }
    return false;
}

int BLMesh::RmvSymBdry(int p1, int p2)
{
    int i;
    for (i = 0; i < m_nSymBdrys; i++) {
        if ((m_pSymBdrys[i].conn[0] == p1 && m_pSymBdrys[i].conn[1] == p2) ||
            (m_pSymBdrys[i].conn[0] == p2 && m_pSymBdrys[i].conn[1] == p1)) {
            SetSymBdryDelete(i);
        }
    }
    return 0;
}

void BLMesh::SetSymBdryDelete(int i) { m_pSymBdrys[i].nconn = -1; }

bool BLMesh::IsSymBdryDelete(int i) { return m_pSymBdrys[i].nconn < 0; }

/**
 * @brief Key function in generating boundary layer mesh
 * @author yhf
 * @note Big function
 */

void BLMesh::GenerateBLMesh()
{
    int iLayer = 0, i, nNods, iNod, iNodNew, cnt = 0, nFrtNods;
    double tmls, tmle, tmts, tmte;
    BLVector newpos, normal, normn;
    BLNode *blNods[MAX_FRONT_NODES], *blNod = nullptr, **blFrtNods = nullptr, **blFrtNodsNew = nullptr;
    bool bCreateFront;
    std::unordered_set<int> setblNods;
    double *angle = nullptr;

    static vector<int> p;
#define NORMALTYPE 3
    double min_height_first_layer = std::numeric_limits<double>::max();
    double max_height_first_layer = std::numeric_limits<double>::lowest();
    double ave_height_first_layer = 0;
    tmts = clock();

    while (iLayer < m_nTotalLayer) {
        if (checkterminate()) {
            FreeMemoryInFrontAndNode();
            resetterminate();
            return;
        }
        string tmp = to_string(iLayer);

        if (iLayer < 10) {
            tmp = '0' + tmp;
        }
        // SaveBLMesh(const_cast<char*>((tmp+".vtk").c_str()));
        m_blNxtFList = std::shared_ptr<BLFrontList>(new BLFrontList);
        cnt = 0;
        m_averUb = 0.0;

        tmls = clock();

        // 1. Get all nodes need to be propagated
        auto num_nodes_lower_layer = setblNods.size();
        setblNods.clear();
        setblNods.reserve(num_nodes_lower_layer);
        ave_front_size = 0;
        int nfront = 0;

#ifdef _DEBUG
        int num_front = 0;
        m_blFrontList->RestoreFront();
        while (m_blFrontList->HasNextFront()) {
            BLFront *blFront = m_blFrontList->GetNextFront();
            num_front++;
        }
        p.push_back(num_front);
#endif

        if (iLayer == 0) {
            m_blFrontList->RestoreFront();
            while (m_blFrontList->HasNextFront()) {
                // 计算front宽度
                BLFront *blFront = m_blFrontList->GetNextFront();
                CalFrontSize(blFront);
                ave_front_size += blFront->GetFrontSize();
                nfront++;

                blFront->GetNodes(&nNods, blNods);
                for (i = 0; i < nNods; i++) {
                    blNod = blNods[i];

                    if (m_pNodes[blNod->GetNodIdx()].bfarfield) {
                        blNod->SetStopFlag(true);
                    }
                    // else
                    // int iiii = blNod->GetNodIdx();

                    setblNods.insert(blNod->GetNodIdx());
                }
            }
            nFrtNods = setblNods.size();
            if (nFrtNods <= 0) {
                break;
            }
            if (nfront) {
                ave_front_size /= nfront;
            }

            if (blFrtNods) {
                delete[] blFrtNods;
                blFrtNods = nullptr;
            }
            blFrtNods = new BLNode *[nFrtNods];

            i = 0;
            auto setit = setblNods.begin();
            while (setit != setblNods.end()) {

                BLNode *node = (BLNode *)m_pNodes[*setit].pointer;
                if (node) {
                    blFrtNods[i++] = node;
                }
                ++setit;
            }
            nFrtNods = i;

            // 处理最小层数
            resetZeroHeight(nFrtNods, blFrtNods);
        } else { // not the first layer
            int count = 0;
            for (int i = 0; i < nFrtNods; i++) {
                if (blFrtNods[i]->GetUpperNode()) {
                    blFrtNods[count] = blFrtNods[i]->GetUpperNode();
                    count++;
                }
            }
            nFrtNods = count;
            // std::ofstream fout("1.txt");
            // for (int i = 0; i < count; i += 1)
            //	fout << blFrtNods[i]->GetCoord(m_pNodes)[0]<<" "<< blFrtNods[i]->GetCoord(m_pNodes)[1] << "	";
            // fout.close();
        }

        if (nFrtNods == 0) {
            break;
        }

#ifdef _DEBUG
        if (m_nCurrLayer == 0) {
            m_ocTree->num_inter = 0;
        }
        // p.push_back(m_ocTree->num_inter);
        cout << m_ocTree->num_inter << " ";
#endif

        spdlog::info("generating layer {} ...", iLayer + 1); // propagating each front in the list

                                                             // 2. Calculate normal vectors for all front nodes
#if 1

#ifndef _DEBUG

#endif

        // omp_set_num_threads(10);
        int nConcave = 0;
        int nConvex = 0;

// #ifndef DEBUG
#ifdef USE_OPENMP
        BLNode **node_array = new BLNode *[nFrtNods];
        for (int i = 0; i < nFrtNods; i++) {
            node_array[i] = new BLNode();
        }

#pragma omp parallel for private(blNod, iNod, iNodNew, newpos, normal, normn) if (nFrtNods > 100000)
#endif // USE_OPENMP

        // #endif

        for (int i = 0; i < nFrtNods; i++) {
            double pt[DIM3];

            blNod = blFrtNods[i];
            iNod = blNod->GetNodIdx();

#ifdef _GEOM_NORMAL
            BLNode *blNods[MAX_FRONT_NODES], *neigNods[2]; // , * blNodNew;
            normn = blNod->GetNormal(m_pNodes, NORMALTYPE);

            int min_layer = 100000000;
            double min_length = 1e20;
            double min_ratio = 1e20;
            for (auto f : blNod->GetNeigFronts()) {
                int eid = f->GetSurfaceElmIdx();
                int id = m_pElems[eid].igom;
                if (id < cf.step_len_vec.size()) {
                    int target_layer = cf.layer_num_vec[id];
                    double target_length = cf.step_len_vec[id];
                    min_layer = std::min(min_layer, target_layer);
                    min_length = std::min(min_length, target_length);
                    min_ratio = std::min(min_ratio, cf.ratio_vec[id]);
                }
            }
            if (min_layer < 1000000) {
                blNod->respect_layer = min_layer;
                blNod->respect_height = min_length;
                blNod->respect_ratio = min_ratio;
            } else {
                blNod->respect_layer = cf.layer_num;
                blNod->respect_height = cf.step_len;
                blNod->respect_ratio = cf.ratio;
            }

            if (cf.iscompresslen) {
                blNod->respect_height = blNod->recommend_height;
            }
            int a = blNod->GetDecentID();
            const double eps = 0.1;
            const double target = 3000.57;

            if (std::abs(m_pNodes[a].coord[0] - target) <= eps && std::abs(m_pNodes[a].coord[1] - 56.7135) <= eps &&
                std::abs(m_pNodes[a].coord[2] - 365.083) <= eps) {
                tempcalculate++;
            }
            if (blNod->GetBSys()) {
                auto ans = blNod->GetHeight();
                int decent_id = blNod->GetDecentID();
                int nodeid = iNod;
                Eigen::RowVector3d start_point(m_pNodes[nodeid].coord[0], m_pNodes[nodeid].coord[1], m_pNodes[nodeid].coord[2]);
                Eigen::RowVector3d normal(ans[0], ans[1], ans[2]);
                std::vector<int> faceid = m_pNodes[nodeid].isymfc;
                if (faceid.size() == 1) {
                    faceid2sp[faceid[0]].adjustNormal(start_point, normal);
                } else {
                    for (int j = 0; j < 1; j++) {
                        faceid2sp[faceid[0]].adjustNormal(start_point, normal, true);
                        faceid2sp[faceid[1]].adjustNormal(start_point, normal, true);
                    }
                }
                BLVector n(normal(0), normal(1), normal(2));
                blNod->SetNormal(n.normalized());
            }

#ifdef CHECK_SKEWNWSS
            if (blNod->GetBeitaVisu() < 0.01 / 180 * PI) {
#else
            if (blNod->GetBeitaVisu() < 0) {
#endif
                blNod->SetStopFlag();
                // cout << blNod->GetBeitaVisu();
                continue;
            }

            blNod->SetHightRatio(0);

            normn = blNod->GetNormal();
#endif
            normal = blNod->GetHeight();

            newpos = normal;
            newpos.x += m_pNodes[iNod].coord[0];
            newpos.y += m_pNodes[iNod].coord[1];
            newpos.z += m_pNodes[iNod].coord[2];

#ifdef USE_OPENMP
            BLNode *blNodNew = node_array[i];
#else
            BLNode *blNodNew = new BLNode();
#endif
            blNodNew->SetBSys(blNod->GetBSys(), blNod->GetSymAxis());

            // blNodNew->SetNodeType(blNod->GetNodeType());
            blNodNew->SetNormal(blNod->GetNormal());
            blNodNew->SetBeitaVisu(0.8);
            blNodNew->beitaVisu = blNod->beitaVisu;
            blNodNew->m_h0 = blNod->m_h0 * 0.95 + cf.step_len * 0.05;
            blNodNew->respect_height = blNod->respect_height;
            blNodNew->respect_layer = blNod->respect_layer;
            blNodNew->respect_ratio = blNod->respect_ratio;
            // blNodNew
            //
            pt[0] = newpos.x;
            pt[1] = newpos.y;
            pt[2] = newpos.z;

#ifdef _VECTOR_FIELD
#ifdef _GEOM_NORMAL
            // geometric method
            // normn = blNod->GetNormal(m_pNodes);
            blNodNew->SetNormal(normn);

            // blNodNew->SetDistanceRatio(blNod->GetDistanceRatio());
            blNodNew->SetVirtualFlag(blNod->GetVirtualFlag());
            blNodNew->SetDecentID(blNod->GetDecentID());

#else

#if 1
#ifndef _FAST_MULTIPOLE
            // regular BEM
            normn = blNodNew->GetNormal(m_pPotentialBEM, pt, m_pPt, m_nSurfNodes, m_pBElm, m_nSurfElems, pBc, m_pNorm, pU,
                                        m_pNodes[iNod].bsysm, iNod);

            if (blNod->GetNormMethod() == NormalGenMethod::PHYSICAL_COPY || blNod->GetNormMethod() == NormalGenMethod::GEOMETRICAL_COPY) {
                normn = blNod->GetNormal();
                blNodNew->SetNormal(normn);
            }
#else
            // fast multipole method

            normn = blNodNew->GetNormal(m_pDomain, pt, m_pNodes, iLayer);
#endif
#else
            // for test
            BLVector nor1, nor2;
            // fast multipole method

            nor1 = blNodNew->GetNormal(m_pDomain, pt, iLayer);

            // regular BEM
            normn = blNodNew->GetNormal(m_pPotentialBEM, pt, m_pPt, m_nSurfNodes, m_pBElm, m_nSurfElems, pBc, m_pNorm, pU,
                                        m_pNodes[iNod].bsysm, iNod);

            if (blNod->GetNormMethod() == NormalGenMethod::PHYSICAL_COPY || blNod->GetNormMethod() == NormalGenMethod::GEOMETRICAL_COPY) {
                normn = blNod->GetNormal();
                blNodNew->SetNormal(normn);
            }

            angle[i] = AngleNorm(nor1, normn);
#endif

#endif
#else
            normn = blNodNew->GetNormal(m_pPotentialBEM, pt, m_pPt, m_nSurfNodes, m_pBElm, m_nSurfElems, m_pBc, m_pNorm, m_pU);
#endif
            // 			//for test
            // 			if(iNod == 2038)
            // 				blNod->SetStopFlag();
            iNod = blNod->GetDecentID();

            if (CheckStop(blNod, blNodNew, iLayer) && !blNod->GetBSys()) {
                // cout << "stop" << endl;
                blNod->SetStopFlag();
                delete blNodNew;
                continue;
            } else {
                blNodNew->AddLowerNode(blNod);
                blNod->AddUpperNode(blNodNew);

                if (blNod->getPerNode()) {
                    if (blNod->getPerNode()->GetUpperNode()) {
                        blNodNew->setPerNode(blNod->getPerNode()->GetUpperNode(), blNod->getPerNodeType());
                        blNod->getPerNode()->GetUpperNode()->setPerNode(blNodNew, blNod->getPerNode()->getPerNodeType());
                    }
                }
            }
        } // 2. Calculate normal vectors of all front nodes
#endif

        // 2.1 Smooth normal vectors of all front nodes

        double time_smooth_start, time_smooth_stop;
        time_smooth_start = clock();

        // cout << "preprocess cost" << (time_smooth_start - tmls) / CLOCKS_PER_SEC << endl;
        int ipass = iLayer;

        NormalSmoothStrategy *ns_strategy = (new SimpleNormalSmoothStrategy(blFrtNods, m_pNodes, nFrtNods));

        // NormalSmoothStrategy *ns_strategy = (new VORONOISMOOTHING::VoronoiNormalSmoothStategy(blFrtNods, m_pNodes, nFrtNods));

        ns_strategy->SetFaceidSP(faceid2sp);
        //	ns_strategy->SetSmoothTimes(ipass);
        double angle = 70;
#ifdef VIRTUAL_MESH
        for (i = 0; i < nFrtNods; i++) {
            blNod = blFrtNods[i];
            if (blNod->GetVirtualFlag()) {
                int nblnod;
                BLNode *vir_nod[MAX_NCONN];
                blNod->GetVirtualNeigNods(vir_nod, &nblnod, m_pNodes);
                vir_nod[nblnod++] = blNod;
                // cout << nblnod << endl;
                BLVector ave_normal(0, 0, 0);
                for (int j = 0; j < nblnod; j++) {
                    ave_normal += vir_nod[j]->GetNormal();
                }
                ave_normal.normalize();
                for (int j = 0; j < nblnod; j++) {
                    double d = vir_nod[j]->GetNormal() * ave_normal;
                    if (d > cos(angle * PI / 180) && d < 0.9999) {
                        BLVector new_normal = vir_nod[j]->GetNormal();
                        BLVector third = ave_normal ^ new_normal;
                        BLVector orth = third ^ ave_normal;
                        orth.normalize();
                        new_normal = orth * tan(angle * PI / 360) + ave_normal;
                        vir_nod[j]->SetNormal(new_normal.normalized());
                    }
                }
            }
        }
#endif
        ns_strategy->SmoothNormal();

        //	delete ns_strategy;
        //	ns_strategy = (new VORONOISMOOTHING::VoronoiNormalSmoothStategy(blFrtNods, m_pNodes, nFrtNods));

        //	ns_strategy = (new WZSMOOTHING::WangZhiNormalSmoothStategy(blFrtNods, m_pNodes, nFrtNods));
        //	ns_strategy->SetSmoothTimes(ipass);

        if (checkterminate()) {
            FreeMemoryInFrontAndNode();
            resetterminate();
            return;
        }

        /// good
        delete ns_strategy;
#ifdef VIRTUAL_MESH
        ipass = 1;
        while (ipass > 0) {
            for (i = 0; i < nFrtNods; i++) {
                blNod = blFrtNods[i];
                SmoothVirtualFrontHeight(blNod, m_pNodes);
            }
            ipass--;
        }
#endif
        ipass = 1;
        // caculate font height fix
        if (iLayer == 0) {
            for (int i = 0; i < nFrtNods; i++) {
                blNod = blFrtNods[i];
                FixHightRatio(blNod, m_pNodes);
            }
        }

        while (ipass > 0) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif // USE_OPENMP
            for (int i = 0; i < nFrtNods; i++) {
                blNod = blFrtNods[i];
                SmoothHeightRatio(blNod, m_pNodes);
            }
            ipass--;
        }

        /*if (iLayer >= 1)
            for (i = 0; i < nFrtNods; i++)
            {
                blNod = blFrtNods[i];
                SmoothHorNodeNorm(blNod, m_pNodes);
            }*/

        time_smooth_stop = clock();
#ifdef _DEBUG
        cout << "smooth process cost" << (time_smooth_stop - time_smooth_start) / CLOCKS_PER_SEC << endl;
#endif
        smooth_time += time_smooth_stop - time_smooth_start;

#ifdef CHANGE_STEP_BY_DISTANCE
        if (iLayer == 0) {
            double distance_start = clock();
#ifndef DEBUG
#pragma omp parallel for
#endif

            for (int i = 0; i < nFrtNods; i++) {
                blNod = blFrtNods[i];
                m_pNodes[blNod->GetNodIdx()].distance = GetDistance(blNod);
                blNod->SetDistanceRatio(max(0.20, m_pNodes[blNod->GetNodIdx()].distance));
            }

            // double ipass = 40;
            int nNeigNods;
            std::vector<BLNode *> blNodes;
            queue<BLNode *> idx;
            for (int i = 0; i < nFrtNods; i++) {
                idx.push(blFrtNods[i]);
            }
            const double distance_step_control = 0.13;
            while (!idx.empty()) {
                blNod = idx.front();
                idx.pop();
                blNodes = blNod->GetNeigNods();
                nNeigNods = blNodes.size();
                double ave = 0;
                for (int j = 0; j < nNeigNods; j++) {
                    ave += blNodes[j]->GetDistanceRatio();
                }
                if (nNeigNods) {
                    ave /= nNeigNods;
                } else {
                    ave = 1;
                }
                if (blNod->GetDistanceRatio() - ave > distance_step_control) {
                    for (int j = 0; j < nNeigNods; j++) {
                        idx.push(blNodes[j]);
                    }
                    idx.push(blNod);
                    blNod->SetDistanceRatio(min(ave, blNod->GetDistanceRatio()));
                    m_pNodes[blNod->GetNodIdx()].distance = GetDistance(blNod);
                }
            }
            double distance_stop = clock();
            spdlog::info("caculating distance cost" << (distance_stop - distance_start) / CLOCKS_PER_SEC);
        }
#endif
        double genera_prism_start = clock();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif // USE_OPENMP

        for (int i = 0; i < nFrtNods; i++) {

            BLFront *blNeigFrts[MAX_NCONN * 2];
            BLNode *blNeigNodes[MAX_NCONN * 2];
            int nblNeigFrts;
            int nblNodes;
            if (m_bIsotropicStop > 1e-5) {
                //	cout << m_bIsotropicStop << " ";

                const std::vector<BLFront *> &blNeigFrts = blFrtNods[i]->GetNeigFronts();
                nblNeigFrts = blNeigFrts.size();

                bool iso_stop = true;
                for (int j = 0; j < nblNeigFrts; j++) {
                    CalFrontSize(blNeigFrts[j]);
                    if (max(0.1 * ave_front_size + 0.9 * blNeigFrts[j]->GetFrontSize(),
                            0.30 * ave_front_size + 0.7 * blNeigFrts[j]->GetFrontSize()) >
                        0.95 * blFrtNods[j]->GetHeightLength() / m_bIsotropicStop) {
                        iso_stop = false;
                    }
                    if (blFrtNods[i]->GetHeightLength() / m_bIsotropicStop > blNeigFrts[j]->GetSqrtFrontSize() * 1.3) {
                        iso_stop = true;
                        break;
                    }
                    if (blFrtNods[i]->GetHeightLength() / m_bIsotropicStop > blNeigFrts[j]->GetMinFrontSize() * 1.8) {
                        iso_stop = true;
                        break;
                    }
                }
                blFrtNods[i]->iso_stop = iso_stop;
            }
        }

        //    		++iLayer;
        //    		continue;
#if 0
		//for test
		OutputAngle("normalangle", angle, nFrtNods, iLayer);
		//exit(0);
		//for test
#endif

        // 3. Propagate all front nodes

        // 这段不能多线程并行，所以要单独运行
        int count = 0;
        for (int i = 0; i < nFrtNods; i++) {
            if (blFrtNods[i]->respect_layer <= iLayer) {
                blFrtNods[i]->SetStopFlag(true);
            }

            if (CheckIsotroStop(blFrtNods[i]) && !blFrtNods[i]->GetBSys()) {
#ifdef _DEBUG
                // cout << "=====================";
                // cout << "isostop";
                // cout << " " << blFrtNods[i]->GetDecentID() << endl;
#endif
                blFrtNods[i]->SetStopFlag();
                if (blFrtNods[i]->GetUpperNode()) {
                    delete blFrtNods[i]->GetUpperNode();
                    blFrtNods[i]->RmvUpperNod(blFrtNods[i]->GetUpperNode());
                }
                continue;
            }
        }

        for (int i = 0; i < nFrtNods; i++) {
            blNod = blFrtNods[i];
            iNod = blNod->GetDecentID();

            if (blNod->getPerNode() && blNod->GetUpperNode() && !blNod->GetStopFlag()) {
                if (!blNod->getPerNode()->GetUpperNode() || blNod->getPerNode()->GetStopFlag()) {
                    blNod->SetStopFlag(true);
                    blNod->GetUpperNode()->SetLowerNode(nullptr);
                    blNod->SetUpperNode(nullptr);
                }
            }
        }
#ifdef USE_OPENMP
#pragma omp parallel for private(blNod, normal)
#endif // USE_OPENMP

        for (int i = 0; i < nFrtNods; i++) {

            blNod = blFrtNods[i];

            // 			if (iLayer == 0)
            // 				normal = blNod->GetNormal();
            // 			else
            normal = blNod->GetHeight();
            iNod = blNod->GetNodIdx();
            if (iLayer == 0 && !blNod->GetStopFlag()) {
                min_height_first_layer = std::min(min_height_first_layer, normal.magnitude());
                max_height_first_layer = std::max(max_height_first_layer, normal.magnitude());
                ave_height_first_layer += normal.magnitude();
                // cout << normal.magnitude() << endl;
                count++;
            }
            PropagateNode(blNod, normal, iLayer);
        }
        if (iLayer == 0) {
            ave_height_first_layer /= count;
        }

        // 4. Form current layer elements
        double propegate_start = clock();
        /*if(m_nCurrLayer>25)
        m_ocTree->printElement("before.pls");*/

        Propagate();
        if (checkterminate()) {
            FreeMemoryInFrontAndNode();
            resetterminate();
            return;
        }

        for (int i = 0; i < nFrtNods; i++) {

            if (blFrtNods[i]->getPerNode()) {
                int pernodeidx = blFrtNods[i]->getPerNode()->GetDecentID();

                if ((!blFrtNods[i]->GetUpperNode()) || blFrtNods[i]->GetUpperNode()->GetNeigFronts().empty()) { // deleted
                    if (blFrtNods[i]->getPerNode()->GetUpperNode()) {
                        // blFrtNods[i]->getPerNode()->RmvUpperNod(blFrtNods[i]->getPerNode()->GetUpperNode());
                        // blFrtNods[i]->getPerNode()->GetNeigFronts();
                        RmvUperNeigFronts(blFrtNods[i]);
                        RmvUperNeigFronts(blFrtNods[i]->getPerNode());
                        if (blFrtNods[i]->GetUpperNode()) {
                            blFrtNods[i]->RmvUpperNod(blFrtNods[i]->GetUpperNode());
                        }
                        if (blFrtNods[i]->getPerNode()->GetUpperNode()) {
                            BLNode *pernode = blFrtNods[i]->getPerNode();
                            pernode->SetStopFlag();
                            auto fs = pernode->GetNeigFronts();
                            for (auto f : fs) {
                                if (f->GetUpperFront()) {
                                    f->SetUpperFront(nullptr);
                                }
                                if (f->GetElmIdx() >= 0) {
                                    SetElmDelete(f->GetElmIdx());
                                }
                            }
                            pernode->RmvUpperNod(pernode->GetUpperNode());
                        }
                    }
                }
            }
        }

        /*if (m_nCurrLayer > 25)
        m_ocTree->printElement("after.pls");*/

        double genera_prism_end = clock();
        check_prism_quality_time += genera_prism_end - genera_prism_start;
#ifdef _DEBUG
        double propegate_stop = clock();
        cout << "propagate process cost" << (propegate_stop - propegate_start) / CLOCKS_PER_SEC << endl;
#endif
#ifdef _DEBUG
        cout << "numberof tri elements=" << m_TriElm.GetSize() << " memory=" << m_TriElm.GetCapacity() << endl;
#endif

        if (iLayer == 0) {
            first_layer_FrtNods = cnt;
        }
        // 4.1 smooth normal vectors of fronts

        // 5. update node's uppernode status by checking all its neighbor fronts

        int nBLFronts;
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif // USE_OPENMP
        for (i = 0; i < nFrtNods; i++) {

            auto blNod = blFrtNods[i];
            std::vector<BLFront *> &blFronts = blNod->GetNeigFronts();
            int nBLFronts = blFronts.size();

            bool bHasUpper = false;
            for (int j = 0; j < nBLFronts; j++) {
                BLFront *tmp = blFronts[j]->GetUpperFront();
                if (tmp != nullptr) {
                    bHasUpper = true;
                    break;
                }
            }
            BLNode *bup = blNod->GetUpperNode();
            if (!bHasUpper && bup) {
                blNod->RmvUpperNod(bup);
                blNod->SetStopFlag();
                delete bup;
            }
        }

        // Postpropcess
        m_blFrontList->RestoreFront();
        while (m_blFrontList->HasNextFront()) {
            BLFront *blFront = m_blFrontList->GetNextFront();
            if (blFront->GetLowerFront()) {
                for (int k = 0; k < 3; k++) {
#ifdef USE_DYNAMIC_ARRAY
                    for (int j = 0; j < 2; j++) {
                        int idx = blFront->GetLowerFront()->GetSTriIdx(k, j);
                        if (idx > 0) {
                            m_TriElm.DeleteElem(idx);
                        }
                    }
#endif
                }

                BLFront *b[4];
                int num_front;
                blFront->GetLowerFront()->GetNeigbourFronts(&num_front, b);
                if (num_front != 3) {
                    // front_to_delete_.push_back(blFront->GetLowerFront());
                } else {
                    m_TriElm.DeleteElem(blFront->GetLowerFront()->GetTriIdx());
                    delete blFront->GetLowerFront();
                    blFront->SetLowerFront(nullptr);
                }
            }
        }

        m_blFrontList->RestoreFront();
        m_blFrontListAll[iLayer] = m_blFrontList;
        m_blFrontList = m_blNxtFList;
        int a = m_blFrontList->Count();

#if 0
		if (m_blType == BLMType::blm3d)
		{
			m_blFrontList->RestoreFront();
			while (m_blFrontList->HasNextFront())
			{
				BLFront* blFront = m_blFrontList->GetNextFront();

				PostChckIntersect(blFront);
			}
			m_blFrontList->RestoreFront();
		}
#else
        double generate_pyramid_start = clock();
        //	if(generate_pyramid)

        CreateTransitionElements();
        UpdateSymmetry();

        // if (this->m_nCurrLayer > 16)
        //	m_ocTree->printElement("afterpyramid.pls");
        double generate_pyramid_end = clock();
        generate_pyramid_time += generate_pyramid_end - generate_pyramid_start;
#endif

        // m_ocTree->printElement("after.pls");
        // m_ocTree->printElement("after.pls");
        ++iLayer;
        ++m_nCurrLayer;

        tmle = clock();
        // if (setblNods.empty()) {
        //	spdlog::info("no element need to generate boundary layer mesh\n");
        //	break;
        // }
        // cout << "post process cost" << (tmle - time_smooth_stop) / CLOCKS_PER_SEC << endl;
        // m_ocTree->printElement(string("outele.pls"));
        // m_ocTree->printNodeSize();
        // this->SaveBLMesh(const_cast<char*>(("part/"+std::to_string(m_nCurrLayer)+".vtk").c_str()));
        spdlog::info("finished layer {}, wall-clock time used: {}", iLayer, (tmle - tmls) / CLOCKS_PER_SEC);

        fflush(stdout);
    }

    // SaveBLMesh("boundary_final.vtk");
    // m_ocTree->printElement("after.pls");
    /*free some memory*/
#ifdef _DEBUG
    cout << "number of front" << endl;
    for (auto i : p) {
        cout << i << " ";
    }
    cout << endl << "============================";
#endif
    if (blFrtNods) {
        delete[] blFrtNods;
        blFrtNods = nullptr;
    }
#ifdef WIN32
#ifdef MSC VER
    _ASSERTE(_CrtCheckMemory());
#endif
#endif
    // m_ocTree->printElement("after.pls");
    RemoveOutbdry(boundary_to_delete_.size() / 3, boundary_to_delete_.data());

    tmte = clock();
    printf("Total CPU time, wall-clock time used: %fs\n", (tmte - tmts) / CLOCKS_PER_SEC);
    fflush(stdout);
    m_blFrontList->RestoreFront();
    RemoveOverlapFace(outbdry, noutbdry);
    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        blFront->GetNodes(&nNods, blNods);
        if (blNods[0]->GetNodIdx() < 0 || blNods[1]->GetNodIdx() < 0 || blNods[2]->GetNodIdx() < 0) {
            continue;
        }
        m_vBdyFront.push_back(blFront); // 将最后一层剩余front 加入待选的front
    }

    // RemoveNonManifoldFront();

    for (i = 0; i < m_vBdyFront.size(); i++) {
        BLFront *ft = m_vBdyFront[i];
        ft->GetNodes(&nNods, blNods);
        if (blNods[0]->GetNodIdx() < 0 || blNods[1]->GetNodIdx() < 0 || blNods[2]->GetNodIdx() < 0) {
            continue;
        }

        outbdry[noutbdry * 3 + 0] = blNods[0]->GetNodIdx();
        outbdry[noutbdry * 3 + 1] = blNods[2]->GetNodIdx();
        outbdry[noutbdry * 3 + 2] = blNods[1]->GetNodIdx();
        if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
            throw std::runtime_error("try to input error");
        }

        ++noutbdry;
    }
    PstprecsMergedElm();
    m_ocTree->printTransferCost();
    DeleteOctree();
    DeleteOctreeSymm();
    std::cout.setf(ios::left, ios::adjustfield);
    std::cout.fill(' ');
    // spdlog::info(std::setiosflags(ios::left) << std::setfill(' ') );
    spdlog::info("             Mesh Quality    ");
    spdlog::info("+================================================+");
    spdlog::info("|             Item                |number of cell|");
    spdlog::info("-------------------------------------------------");
    spdlog::info("|         min height(first layer) |{}", min_height_first_layer); // << '|' );
    spdlog::info("|         max height(first layer) |{}", max_height_first_layer); // << '|' );
    spdlog::info("|     average height(first layer) |{}", ave_height_first_layer); // << '|' );
    spdlog::info("|    expected height(first layer) |{}", cf.step_len);
    spdlog::info("+================================================+");

    // free memory (due to fix)
}

int BLMesh::ElmBdryPtCnt(int eidx)
{
    int i, j, dim = DIM3, pidx[DIM3], cnt;

    if (eidx < 0) {
        return 0;
    }

    for (i = 0; i < dim; i++) {
        pidx[i] = m_pElems[eidx].conn[i];
    }

    cnt = 0;
    for (i = 0; i < m_nBdryPnt; i++) {
        for (j = 0; j < dim; j++) {
            if (m_pBdryPnt[i] == pidx[j]) {
                cnt++;
            }
        }
    }

    return cnt;
}

void BLMesh::UpdateBdryNorm()
{
    int i, j, k, pidx, eidx, beg, end, cnt;
    BLVector norms[MAX_NCONN * 2];
    bool flag[MAX_NCONN * 2];
    double angle_eps = rad(15);
    BLVector norm;

    for (i = 0; i < m_nBdryPnt; i++) {
        std::vector<std::vector<BLVector>> normgroup;
        cnt = 0;
        pidx = m_pBdryPnt[i];

        if (m_pNodes[pidx].bsysm) {
            continue;
        }

        beg = m_pPntIdx[pidx];
        end = m_pPntIdx[pidx + 1];

        for (j = beg; j < end; j++) {
            eidx = m_pPntElm[j];

            norms[cnt] = norm;
            flag[cnt] = 0;
            cnt++;

            if (cnt > MAX_NCONN * 2) {
                throw(std::string("error cut"));
            }
        }

        // classify norms
        for (j = 0; j < cnt; j++) {
            for (k = 0; k < normgroup.size(); k++) {
                std::vector<BLVector> vtmp = normgroup[k];
                BLVector ntmp = vtmp[0];

                if (ntmp * norms[j] > cos(angle_eps)) {
                    normgroup[k].push_back(norms[j]);
                    flag[j] = 1;
                }
            }

            if (!flag[j]) {
                std::vector<BLVector> vtmp;
                vtmp.push_back(norms[j]);
                normgroup.push_back(vtmp);
                flag[j] = 1;
            }
        }

        // update norms
        BLVector fnorm(0, 0, 0);
        for (j = 0; j < normgroup.size(); j++) {
            BLVector ntmp(0, 0, 0);
            std::vector<BLVector> vtmp = normgroup[j];
            for (k = 0; k < vtmp.size(); k++) {
                ntmp += vtmp[k];
            }
            ntmp.x /= vtmp.size();
            ntmp.y /= vtmp.size();
            ntmp.z /= vtmp.size();

            ntmp.normalize();

            fnorm += ntmp;
        }
        fnorm.normalize();

        for (j = beg; j < end; j++) {
            eidx = m_pPntElm[j];

#if 0 //before \
	  //omit if the end points of this element are all boundary points
			if (CheckElmPtsBdry(eidx))
				continue;

			m_pElems[eidx].norm = fnorm;
#else // update

#endif
        }
    }
}

// void BLMesh
void BLMesh::SmoothNodeNormAndRatio(BLNode *blNod, MBLNode *pNodes)
{

#define POWER 1
    double thredhold = 30;
    std::vector<BLNode *> blNodes;
    std::set<BLNode *> stBlnods;
    int nNeigFrts, nNeigNods;

    BLVector norms[MAX_NCONN * 2];
    if (blNod->GetBeitaVisu(blNod->GetNormal()) * 180 / PI < thredhold) {
        return;
    }

    if (blNod->IsComplexNode()) {
        return;
    }

    // if (blNod->GetBSys())
    if (blNod->GetBdryPt()) {
        blNodes = blNod->GetNeigNods();
        nNeigNods = blNodes.size();
    } else {
        blNodes = blNod->GetNeigNods();
        nNeigNods = blNodes.size();
    }

    nNeigFrts = blNod->NeighFrontCount();

    if (nNeigNods > MAX_NCONN * 2) {
        spdlog::info("Error: the neighborhood is not right for this node ({}-{})!\n", nNeigNods, nNeigFrts);
        // throw(std::string("Error: the neighborhood is not right"));
    }

    BLVector fnorm(0, 0, 0), norm0(0, 0, 0);
    bool nearCliff = false;

    double height = 0;
    double lengthSum = 0;
#ifdef AVERC
    BLVector averC(0.0, 0.0, 0.0);
    BLFront *blNeigFrts[MAX_NCONN * 2];
    int nblNeigFrts;
    blNod->GetNeigFronts(blNeigFrts, &nblNeigFrts);
    double pk = 1e10; // 附近front的平均宽高比
    for (int i = 0; i < nblNeigFrts; i++) {
        pk = min(pk, (blNeigFrts[i]->GetFrontSize() / blNod->GetHeight().magnitude()));
    }
    // cout << pk << endl;
    if (pk > 1e9) {
        pk = 1;
    }

#endif
    for (int i = 0; i < nNeigNods; i++) {
        BLNode *bnod = blNodes[i];
        // averC.x+=bnod.
        int id = bnod->GetNodIdx();
        if (bnod->GetBeitaVisu(blNod->GetNormal()) * 180 / PI < thredhold) {
            nearCliff = true;
        }
        BLVector ntmp = bnod->GetNormal();
        fnorm += ntmp;
#ifdef AVERC
        BLVector h = bnod->GetHeight();
        averC += h + BLVector(pNodes[id].coord[0], pNodes[id].coord[1], pNodes[id].coord[2]);
#endif
    }

#ifdef POWER
    fnorm = BLVector(0, 0, 0);
#endif
    blNodes = blNod->GetNeigNods();
    nNeigNods = blNodes.size();
    for (int i = 0; i < nNeigNods; i++) {
        BLNode *bnod = blNodes[i];
        // averC.x+=bnod.
        int id = bnod->GetNodIdx();
#ifdef POWER
        double length = 1.0 / BLVector(pNodes[id].coord[0] - pNodes[blNod->GetNodIdx()].coord[0],
                                       pNodes[id].coord[1] - pNodes[blNod->GetNodIdx()].coord[1],
                                       pNodes[id].coord[2] - pNodes[blNod->GetNodIdx()].coord[2])
                                  .magnitude();
        height += bnod->GetHightRatio() * length;
        lengthSum += length;
        BLVector ntmp = bnod->GetNormal();
        ntmp.normalize();
        fnorm += ntmp * length * length;
#else
        height += bnod->GetHightRatio();

#endif
    }
#ifdef AVERC
    if (nNeigNods > 0) {
        averC = averC / nNeigNods -
                BLVector(pNodes[blNod->GetNodIdx()].coord[0], pNodes[blNod->GetNodIdx()].coord[1], pNodes[blNod->GetNodIdx()].coord[2]);
    }
    averC.normalize();
#endif
#ifdef POWER
    if (lengthSum > 0) {
        height /= lengthSum;
    }
#else
    height /= nNeigNods;
#endif
    fnorm.normalize();

    norm0 = blNod->GetNormal();

#ifdef AVERC

    if (blNod->GetBdryPt()) {
        norm0 = 0.8 * norm0 + 0.2 * averC;
    } else {
        norm0 = 0.8 * norm0 + 0.2 * averC;
    }
#else
    if (blNod->GetBdryPt()) {
        norm0 = 0.6 * norm0 + 0.4 * fnorm;
    } else {
        norm0 = 0.3 * norm0 + 0.7 * fnorm;
    }

#endif
    norm0.normalize();
    if (!(blNod->GetBdryPt())) {
        double heightratio = blNod->GetHightRatio();
        if (!nearCliff) {
            heightratio = 0.5 * height + 0.5 * heightratio;
        }
        blNod->SetHightRatio(heightratio);
    }
    blNod->SetNormal(norm0);
}

void BLMesh::SmoothVirtualFrontHeight(BLNode *blNod, MBLNode *pNodes)
{
    BLNode *blNodes[MAX_NCONN * 2];
    int nNeigFrts, nNeigNods;

    BLVector norms[MAX_NCONN * 2];

    if (blNod->GetBdryPt()) {
        return;
    }

    blNod->GetVirtualNeigNods(blNodes, &nNeigNods, m_pNodes);
    double height = blNod->GetHightRatio();
    for (int i = 0; i < nNeigNods; i++) {
        height += blNodes[i]->GetHightRatio();
    }
    if (!nNeigNods) {
        return;
    }
    blNod->SetVirtualFlag(true);
    height /= nNeigNods + 1;
    blNod->SetHightRatio(height);
}

void BLMesh::SmoothHorNodeNorm(BLNode *blNod, MBLNode *pNodes)
{
    if (!blNod->GetLowerNode()) {
        return;
    }
    BLVector lowerNormal = blNod->GetLowerNode()->GetNormal();
    BLVector newnormal = blNod->GetNormal();
    ;
    double angle;
    if (blNod->GetVirtualFlag()) {
        angle = 1 + m_nCurrLayer / 10;
    } else {
        angle = 3 + 2 * m_nCurrLayer / 3;
    }
    if (lowerNormal * newnormal < acos(angle / 180 * PI)) { // acos(3.6/180*PI)
        BLVector v1 = lowerNormal ^ newnormal;
        BLVector v2 = -1 * lowerNormal ^ v1;
        v2.normalize();
        newnormal = tan(angle / 180 * PI) * v2 + lowerNormal;
    }

    blNod->SetNormal(newnormal.normalized());
}

/// #define _MULTI_THREAD_CAL

double BLMesh::GetDistance(BLNode *nod)
{
    if (nod->GetVirtualFlag()) {
        return 1;
    }
    BLVector normal;
    std::vector<BLNode *> front_arr;
    int nfront;
    front_arr = nod->GetNeigNods();
    nfront = front_arr.size();
    double size = 0;
    double q = cf.ratio2;
    double a0 = cf.step_len;
    double ans = 1.0;
    BLVector start(m_pNodes[nod->GetNodIdx()].coord[0], m_pNodes[nod->GetNodIdx()].coord[1], m_pNodes[nod->GetNodIdx()].coord[2]);
    for (int i = 0; i < nfront; i++) {
        BLVector neop(m_pNodes[front_arr[i]->GetNodIdx()].coord[0], m_pNodes[front_arr[i]->GetNodIdx()].coord[1],
                      m_pNodes[front_arr[i]->GetNodIdx()].coord[2]);
        double size = (neop - start).magnitude();
        size = min(size, cf.step_len * pow(cf.ratio2, cf.layer_num));
        double sum = (size * q - a0) * 3 / (q - 1);
        normal = front_arr[i]->GetNormal();
        ans = min(ans, m_ocTree->chckIntersectWithLine(start + normal * 1e-3, start + normal * sum));
    }
    return ans;

    // return m_ocTree->chckIntersectWithLine();
}

void BLMesh::PropagateNode(BLNode *blNod, BLVector normal, int iLayer)
{
    BLVector newpos;
    int iNod, iNodNew;

    iNod = blNod->GetNodIdx();
    /////////////////// 可删除
    // if ((!m_pNodes[blNod->GetNodIdx()].bsysm)&& normal.y + m_pNodes[iNod].coord[1] < -0.01)
    //	blNod->SetStopFlag();
    // if (m_pNodes[blNod->GetNodIdx()].bsysm) {
    //	normal.y = 0;
    // }
    /////////////////// 可删除

    if (!blNod->GetStopFlag()) {
        newpos = normal;
        newpos.x += m_pNodes[iNod].coord[0];
        newpos.y += m_pNodes[iNod].coord[1];
        newpos.z += m_pNodes[iNod].coord[2];

        double pt[DIM3];
        pt[0] = newpos.x;
        pt[1] = newpos.y;
        pt[2] = newpos.z;

        if (!m_pNodes[blNod->GetNodIdx()].bsysm) {
            // if (!sj.judge(pt)) {
            //	blNod->SetStopFlag();
            //	return;
            // }
        }

        // BLNode *blNodNew = new BLNode(m_blType);

        BLNode *blNodNew = blNod->GetUpperNode();
        // if (!blNodNew) {
        //	blNod->SetStopFlag(true);
        //	return;
        // }
        if (blNod->getPerNode()) {

            if (blNod->getPerNode()->GetUpperNode() && blNod->getPerNode()->GetUpperNode()->GetLayerNum() == iLayer + 1) {
                BLNode *per_node = blNod->getPerNode()->GetUpperNode();
                blNodNew->setPerNode(per_node, blNod->getPerNodeType());
                per_node->setPerNode(blNodNew, blNod->getPerNode()->getPerNodeType());
                Eigen::Vector3d c_coord(newpos.x, newpos.y, newpos.z);
                Eigen::Vector3d per_coord(per_node->GetCoord(m_pNodes).x, per_node->GetCoord(m_pNodes).y, per_node->GetCoord(m_pNodes).z);
                Eigen::Vector3d nnode;
                if (blNod->getPerNodeType() == BLNode::PerType::Forward) {
                    nnode = per_coord - cf.shift_vec;
                    nnode = cf.reverse_matrix * nnode;
                } else {
                    nnode = cf.rotate_matrix * per_coord;
                    nnode = nnode + cf.shift_vec;
                }
                if ((nnode - c_coord).norm() < 1) {

                    newpos.x = nnode(0);
                    newpos.y = nnode(1);
                    newpos.z = nnode(2);
                } else {
                }
            }
        }

        // set uvalue and gvalue
        iNodNew = AddNode(newpos, 0.0);

        m_pNodes[iNodNew].pointer = (void *)blNodNew;

        // need to fix next version
        if (m_pNodes[iNod].bsysm) {
            m_pNodes[iNodNew].isymfc = m_pNodes[iNod].isymfc;
            m_pNodes[iNodNew].bsysm = m_pNodes[iNod].bsysm;

            //	blNodNew->SetBSys(true, m_pNodes[iNod].symaxis);
        }

        // cout << m_pNodes[iNodNew].height << endl;;

        blNodNew->SetNodIdx(iNodNew);
        blNodNew->SetLayerNum(iLayer + 1);
        blNodNew->SetHightRatio(blNod->GetHightRatio());
        blNodNew->SetFixedHightRatio(blNod->GetFixedHightRatio());
        blNodNew->SetBdryPt(blNod->GetBdryPt());

        if (!blNod->GetDescentNod()) {
            blNodNew->SetDescentNode(blNod);
        } else {
            blNodNew->SetDescentNode(blNod->GetDescentNod());
        }
    } else { /*
               BLNode *blNodNew = blNod->GetUpperNode();
               if (blNodNew)
               {
                   blNodNew->RmvLowerNod(blNod);
                   blNod->RmvUpperNod(blNodNew);
               }*/
        StopPropagateNode(blNod);
    }
}

void BLMesh::StopPropagateNode(BLNode *blNod)
{
    if (blNod == nullptr) {
        spdlog::info("Warning: trying to stop propagating a nullptr BLNode !\n");
        return;
    }
    BLNode *pTmp = blNod->GetUpperNode();
    if (pTmp) {
        pTmp->RmvLowerNod(blNod);
        blNod->RmvUpperNod(pTmp);
        delete pTmp;
    }
    blNod->SetStopFlag();
    if (blNod->getPerNode() && blNod->getPerNode()->GetUpperNode()) {
        StopPropagateNode(blNod->getPerNode());
    }
}

void printprism(MBLNode *node, BLFront *blFront)
{
    BLEntityTopology topu;
    int i, j, iLayer, nNods, iNodNew;
    BLNode *blNods[MAX_NCONN], *blNod;
    int nconn, conn[MAX_NCONN], neigs;
    BLFront *neigFrts[DIM3];
    bool bflag[DIM3];
    int itrix;

    BLVector front_normal = blFront->GetNormal();

    iLayer = blFront->GetLayerNum();
    blFront->GetNodes(&nNods, blNods);
    nconn = nNods;
    bool isvalid = true;
    // check validity
    int pidx;
    for (i = 0; i < nNods; i++) {
        blNod = blNods[i];
        conn[i] = blNod->GetNodIdx();
        conn[nconn + i] = blNod->GetUpperNode()->GetNodIdx();
    }

    int idx0, idx1;
    FILE *fout = nullptr;

    fout = fopen("prism.pls", "w");

    fprintf(fout, "%d %d 0 0 0 0\n", 2, 6);
    for (i = 0; i < 6; i++) {
        pidx = conn[i];
        fprintf(fout, "%d %20.15lf %20.15lf %20.15lf \n", i + 1, node[pidx].coord[0], node[pidx].coord[1], node[pidx].coord[2]);
    }

    fprintf(fout, "1 %d %d %d %d\n", 1, 2, 3, 1);
    fprintf(fout, "2 %d %d %d %d\n", 4, 5, 6, 1);
    fprintf(fout, "2 %d %d %d %d\n", 1, 2, 4, 1);
    fprintf(fout, "2 %d %d %d %d\n", 2, 4, 5, 1);
    fprintf(fout, "2 %d %d %d %d\n", 2, 3, 5, 1);
    fprintf(fout, "2 %d %d %d %d\n", 3, 5, 6, 1);
    fprintf(fout, "2 %d %d %d %d\n", 3, 1, 6, 1);
    fprintf(fout, "2 %d %d %d %d\n", 1, 6, 4, 1);

    fclose(fout);
}

void BLMesh::CheckInsertSideSuface(BLFront *blFront)
{
    BLFront *upper_front = blFront->GetUpperFront();
    if (!upper_front) {
        return;
    }

    int id = blFront->GetTriIdx();

    // return false;
    int i, j, neigs, ntri = 0, tris[3];
    int itri[13]; /// 新增的面
    int idx1, idx2, nNods;
    BLFront *neigFrts[DIM3];
    // bool flag[DIM3];
    bool is_inserect = false;
    bool used_by_neigh_front[3] = {false};
    int *conn;
    BLNode *blNods[DIM3];

    blFront->GetNodes(&nNods, blNods);

    blFront->GetNeigbourFronts(&neigs, neigFrts);

    /*获取底层front周围的三个front的多出来的侧面的三角形ID，存在ithird中*/
    for (i = 0; i < neigs; i++) {
        BLFront *ft = neigFrts[i];
        if (ft->GetUpperFront()) {
            for (j = 0; j < DIM3; j++) {
                if (!ft->IncludeNode(blNods[(j) % DIM3])) {
                    used_by_neigh_front[j] = true; // 标记被使用过
                    break;
                }
            }
        }
    }
    for (i = 0; i < DIM3; i++) {
        if (!used_by_neigh_front[i]) {
            itri[ntri++] = blFront->GetSTriIdx(i, 0);
            itri[ntri++] = blFront->GetSTriIdx(i, 1);
        }
    }

    itri[ntri++] = upper_front->GetTriIdx();

    // m_ocTree->setNodeBefore(upper_front->GetOuterNode());
    // itri[ntri++] = upper_front->GetTriIdx();
    for (int k = 0; k < ntri; k++) {
        if (m_ocTree->chckIntersectPreProcess(itri[k])) {
            is_inserect = true;
            break;
        }
        // if (!blNods[0]->GetBSys() && !blNods[1]->GetBSys() && !blNods[2]->GetBSys()) {
        //	if (m_ocTree_symm->chckIntersectPreProcess(itri[k]))
        //	{
        //		is_inserect = true;
        //		break;
        //	}
        // }
    }

    if (is_inserect) {
#ifdef _DEBUG
        cout << "side inter surface id = " << blNods[0]->GetDecentID() << endl;
#endif
        scheck.clear();
        for (i = 0; i < DIM3; i++) {

            RmvUperNeigFronts(blNods[i]);
            StopPropagateNode(blNods[i]);
        }

        for (auto k : scheck) {
            inser_queue.push_back(k);
        }
    } else {
    }
}
void BLMesh::CheckInsertSuface(BLFront *blFront)
{
    BLFront *upper_front = blFront->GetUpperFront();
    if (!upper_front) {
        return;
    }

    int id = blFront->GetTriIdx();

    // return false;
    int i, j, neigs, ntri = 0, tris[3];
    int itri[13]; /// 新增的面
    int idx1, idx2, nNods;
    BLFront *neigFrts[DIM3];
    // bool flag[DIM3];
    bool is_inserect = false;
    bool used_by_neigh_front[3] = {false};
    BLNode *blNods[DIM3];

    blFront->GetNodes(&nNods, blNods);

    blFront->GetNeigbourFronts(&neigs, neigFrts);

    for (int k = 0; k < 3; k++) {
        if (m_ocTree->check_intersection_in_set(blNods[k]->GetUpperNode()->GetNodIdx())) {
            is_inserect = true;
#ifdef _DEBUG
            cout << "side inter wall or far" << endl;
#endif
            for (i = 0; i < neigs; i++) {
                inser_queue.push_back(neigFrts[i]);
            }
            inser_queue.push_back(blFront);
            return;
        }
        //		if (!blNods[0]->GetBSys() && !blNods[1]->GetBSys() && !blNods[2]->GetBSys()) {
        //			if (m_ocTree_symm->check_intersection_in_set(blNods[k]->GetUpperNode()->GetNodIdx()))
        //			{
        //				is_inserect = true;
        // #ifdef _DEBUG
        //				cout << "=====================" << endl;
        //				cout << "side inter symm";
        // #endif
        //				for (i = 0; i < neigs; i++) {
        //					inser_queue.push_back(neigFrts[i]);
        //				}
        //				inser_queue.push_back(blFront);
        //				return;
        //			}
        //		}
    }
}
void BLMesh::insertAndRmTriInOctree(BLFront *blFront)
{
    int id = blFront->GetTriIdx();

    // return false;
    int i, j, neigs, ntri = 0, tris[3], istri[12], nstri = 0;
    int itri[13]; /// 新增的面
    int idx1, idx2, nNods;
    BLFront *neigFrts[DIM3];
    // bool flag[DIM3];
    bool is_inserect = false;
    bool used_by_neigh_front[3] = {false};
    BLNode *blNods[DIM3];
    blFront->GetNodes(&nNods, blNods);

    blFront->GetNeigbourFronts(&neigs, neigFrts);

    /*获取底层front周围的三个front的多出来的侧面的三角形ID，存在ithird中*/
    for (i = 0; i < neigs; i++) {
        BLFront *ft = neigFrts[i];
        if (ft->GetUpperFront()) {
            for (j = 0; j < DIM3; j++) {
                if (!ft->IncludeNode(blNods[j % DIM3])) {
                    used_by_neigh_front[j] = true; // 标记被使用过
                    break;
                }
            }
        }
    }

    m_ocTree->rmDataPreProcess(id);

    // m_ocTree->setNodeBefore(blFront->GetOuterNode());
    for (i = 0; i < DIM3; i++) {
        if (!used_by_neigh_front[i]) {

            m_ocTree->insertPreProcess(blFront->GetSTriIdx(i, 0));
            m_ocTree->insertPreProcess(blFront->GetSTriIdx(i, 1));
        }
    }

    m_ocTree->insertPreProcess(blFront->GetUpperFront()->GetTriIdx());

    // blFront->GetUpperFront()->SetOuterNode(m_ocTree->getNodeBefore());
    m_ocTree->setNodeBefore(blFront->GetOuterNode());
}
void BLMesh::Propagate()
{

    m_ocTree->record_stone_position();

    m_blFrontList->RestoreFront();
    // outer boundary
    std::set<int> blNodStop;

    insert_container.clear();
    rm_container.clear();
    int cnt = 0, i;
    int nNods;
    BLNode *blNods[3];
    bool bCreateFront = false;
    std::vector<BLFront *> vec_create_fronts;
    vec_create_fronts.reserve(m_nSurfElems);
    while (m_blFrontList->HasNextFront()) {
        cnt++;

        bCreateFront = true;
        BLFront *blFront = m_blFrontList->GetNextFront();
        blFront->GetNodes(&nNods, blNods);

        for (i = 0; i < 3; i++) {
            int iidn = blNods[i]->GetNodIdx();

            if (blNods[i]->GetStopFlag() || !blNods[i]->GetUpperNode()) {
                bCreateFront = false;
                blNods[i]->SetStopFlag();
                // store stop node
                blNodStop.insert(blNods[i]->GetNodIdx());
            }
        }

        if (bCreateFront) {
            vec_create_fronts.push_back(blFront);
        } else {

            m_vBdyFront.push_back(blFront);
        }
    }
    int size = vec_create_fronts.size();
#ifdef _DEBUG
    cout << "number of active fornt=" << size << endl;
#endif

#ifdef USE_OPENMP
    omp_set_num_threads(omp_get_num_procs());
#pragma omp parallel for if (size > 3000) schedule(guided, 100)
#endif
    for (int i = 0; i < size; i++) {
        PreCheckPrismValid(vec_create_fronts[i]);
    }
    for (int i = 0; i < size; i++) {
        prePropagate(vec_create_fronts[i]);
    }

    inser_queue = deque<BLFront *>();

    m_ocTree->pOctreeAgent->mset.clear();

    if (size > m_nSurfElems * 0.2) {
        m_ocTree->sort_telerant = true;
    }

    double check_intersection_start = clock();
    m_blFrontList->RestoreFront();

    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        if (blFront->GetUpperFront()) {
            insertAndRmTriInOctree(blFront);
        }
    }

    vector<BLFront *> &bdy_front = vec_create_fronts;

#ifdef USE_OPENMP
#pragma omp parallel for if (size > 3000) schedule(guided, 100)
#endif
    if (size > m_nSurfElems * 0.2) {
        m_ocTree->find_intersected_triangles();
        m_ocTree->sort_telerant = false;
        for (int i = 0; i < size; i++) {
            CheckInsertSuface(bdy_front[i]);
        }
    } else {

        for (int i = 0; i < size; i++) {
            CheckInsertSideSuface(bdy_front[i]);
        }
    }

    while (!inser_queue.empty()) {
        BLFront *blFront = inser_queue.front();
        inser_queue.pop_front();
        CheckInsertSideSuface(blFront);
    }
#ifdef BRUTE_FORCE_CHECK
    for (int i = 0; i < size; i++) {
        if (bdy_front[i]->interact_with_other && bdy_front[i]->GetUpperFront()) {
            throw std::logic_error("not deleted!");
        }
    }

#endif

    m_ocTree->sorted = true;
    double check_intersection_end = clock();
    this->check_prism_time += check_intersection_end - check_intersection_start;
}
void BLMesh::getConnection(BLFront *blFront)
{
    BLEntityTopology topu;
    int i, j, iLayer, nNods, iNodNew;
    BLNode *blNods[DIM3], *blNod;
    int nconn, neigs;
    BLFront *neigFrts[DIM3];
    bool bflag[DIM3];
    int itrix;

    BLVector front_normal = blFront->GetNormal();

    iLayer = blFront->GetLayerNum();
    blFront->GetNodes(&nNods, blNods);
    nconn = nNods;
    bool isvalid = true;
    // check validity
    int pidx;
    for (i = 0; i < nNods; i++) {
        blNod = blNods[i];
        blFront->conn[i] = blNod->GetNodIdx();
        blFront->conn[nconn + i] = blNod->GetUpperNode()->GetNodIdx();
    }
    int tmp;

    tmp = blFront->conn[1];
    blFront->conn[1] = blFront->conn[2];
    blFront->conn[2] = tmp;
    tmp = blFront->conn[4];
    blFront->conn[4] = blFront->conn[5];
    blFront->conn[5] = tmp;
}
void BLMesh::PreCheckPrismValid(BLFront *blFront)
{

    getConnection(blFront);
    BLEntityTopology topu;
    int i, j, iLayer, nNods, iNodNew;
    BLNode *blNods[DIM3], *blNod;
    int nconn, *conn, neigs;
    BLFront *neigFrts[DIM3];
    bool bflag[DIM3];
    int itrix;

    BLVector front_normal = blFront->GetNormal();
    blFront->GetNodes(&nNods, blNods);
    nconn = nNods;
    // check validity
    conn = blFront->conn.data();

    BLVector e1 = BLVector(m_pNodes[conn[5]].coord[0] - m_pNodes[conn[3]].coord[0], m_pNodes[conn[5]].coord[1] - m_pNodes[conn[3]].coord[1],
                           m_pNodes[conn[5]].coord[2] - m_pNodes[conn[3]].coord[2])
                      .normalized();
    BLVector e2 = BLVector(m_pNodes[conn[4]].coord[0] - m_pNodes[conn[3]].coord[0], m_pNodes[conn[4]].coord[1] - m_pNodes[conn[3]].coord[1],
                           m_pNodes[conn[4]].coord[2] - m_pNodes[conn[3]].coord[2])
                      .normalized();
    BLVector e3 = (e1 - e2).normalized();

    BLVector up_front_normal = -1 * e1 ^ e2;

#ifdef CHECK_SKEWNWSS
    double throushold = sin((0.002 * (m_nCurrLayer + 1)) * PI / 180);
#else
    double throushold = 0.0;
#endif
    if (m_nCurrLayer > 10) {
        throushold = sin((0.02 + 0.03 * (m_nCurrLayer - 10)) * PI / 180);
    }
    blFront->is_prism_valid = 1;

    // 上层法向与节点法向反向
    for (int i = 0; i < nNods; i++) {
        blNod = blNods[i];
        if (up_front_normal * blNods[i]->GetNormal() < 1e-6) {
            blFront->is_prism_valid = 0;
#ifdef _DEBUG
            cout << "upper face normal opposite ";
            cout << "id= " << blNod->GetDecentID() << endl;
#endif
            return;
        }
    }

    // 上层面几何退化
    if (up_front_normal.magnitude2() < throushold * throushold || (e2 ^ e3).magnitude2() < throushold * throushold ||
        (e3 ^ e1).magnitude2() < throushold * throushold) {

        blFront->is_prism_valid = 0;
#ifdef _DEBUG
        cout << "upper face degenerate ";
        cout << "id = " << blNods[0]->GetDecentID() << endl;
#endif
        return;
    }

    // 节点法向与上下两层面法向夹角过大
    up_front_normal.normalize();
    for (i = 0; i < nNods; i++) {
        blNod = blNods[i];
        if (blNod->GetNormal() * blFront->GetNormal() < throushold || blNod->GetNormal() * up_front_normal < throushold) {
            blFront->is_prism_valid = 0;
#ifdef _DEBUG
            cout << "node normal deviates form face normal ";
            cout << "id = " << blNods[i]->GetDecentID() << endl;
#endif
            return;
        }
    }

    if (!CheckPrismVolumn(2 * nconn, conn) || !CheckPrismSkewness(2 * nconn, conn) || !CheckPrismOrth(blFront, false)) {
        blFront->is_prism_valid = 0;
#ifdef _DEBUG
        if (!CheckPrismVolumn(2 * nconn, conn)) {
            cout << "volume id=";
        }
        if (!CheckPrismSkewness(2 * nconn, conn)) {
            cout << "skewness id=";
        }
        if (!CheckPrismOrth(blFront, false)) {
            cout << "orth id=";
        }

        cout << blNods[0]->GetDecentID() << endl;

#endif
    }

    // 前沿探测留空
    if (cf.clearance > 0) {
        up_front_normal = up_front_normal.normalized() * cf.clearance + up_front_normal;
    } else {
        up_front_normal = up_front_normal * (blFront->m_pBLNods[0]->GetHeightLength());
    }

    // 取三角形三个点（这里仍按你原来的 conn[i+3] 作为三点）
    BLVector triPos[3];
    BLVector startpos(0, 0, 0);

    for (int i = 0; i < 3; ++i) {
        triPos[i] = BLVector(m_pNodes[conn[i + 3]].coord[0], m_pNodes[conn[i + 3]].coord[1], m_pNodes[conn[i + 3]].coord[2]);
        startpos += triPos[i];
    }
    startpos = startpos / 3.0;

    // 封装一次探测：从 p + n*a 到 p + n*b
    auto probeBlocked = [&](const BLVector &p, double a, double b) -> bool {
        return m_ocTree->chckIntersectWithLine(p + up_front_normal * a, p + up_front_normal * b) <
               1 /*|| m_ocTree_symm->chckIntersectWithLine(p + up_front_normal * a, p + up_front_normal * (b-0.3)) < 1*/;
    };

    bool blocked = false;

    // 1) 先做你原来的“中点”探测（保持到 1.3）
    blocked = probeBlocked(startpos, 0.05, 1.3);

    // 2) 再对三角形三个顶点分别做 0.7 的探测
    if (!blocked) {
        for (int i = 0; i < 3; ++i) {
            if (probeBlocked(triPos[i], 0.05, 0.7)) {
                blocked = true;
#ifdef _DEBUG
                cout << "forward point id = ";
                cout << blNods[i]->GetDecentID() << endl;
#endif
                break;
            }
        }
    }

    if (blocked) {
        blFront->is_prism_valid = 0;
#ifdef _DEBUG
        cout << "forward surface id = ";
        cout << blNods[0]->GetDecentID() << endl;
#endif
    }

    return;
}
void BLMesh::prePropagate(BLFront *blFront)
{
    BLEntityTopology topu;
    int i, nNods, iNodNew;
    BLNode *blNods[DIM3], *blNod;
    int nconn, *conn, neigs;
    BLFront *neigFrts[DIM3];
    bool bflag[DIM3];
    int itrix;
    blFront->GetNodes(&nNods, blNods);

    bool deleted = false;
    for (int i = 0; i < nNods; i++) {
        if (!blNods[i]->GetUpperNode()) {
            deleted = true;
        }
    }
    if (deleted) {
        m_vBdyFront.push_back(blFront);
    } else if (blFront->is_prism_valid != 1) {
        for (int j = 0; j < DIM3; j++) {
            RmvUperNeigFronts(blNods[j]);
            StopPropagateNode(blNods[j]);
        }
        // add to the exposed fronts and symmplane segments
        m_vBdyFront.push_back(blFront);
    } else {

        auto conn = blFront->conn.data();
        BLNode *blNod;
        int neigs;
        BLFront *neigFrts[DIM3];
        int itrix = 0;
        CreateTriangles(blFront, itrix);

        // create a new front
        BLFront *blFrontNew = nullptr;

        blFrontNew = new BLFront();
        blFrontNew->SetTriIdx(itrix);
        blFrontNew->SetSurfaceElmIdx(blFront->GetSurfaceElmIdx());

        // set the relationship between the new propagated nodes and the new front
        for (int j = 0; j < nNods; j++) {
            blNod = blNods[j];

            blFrontNew->AddBLFrontNods(j, blNod->GetUpperNode());
            blNod->GetUpperNode()->AddNeigFronts(blFrontNew);
        }

        blFront->SetUpperFront(blFrontNew);
        blFrontNew->SetLowerFront(blFront);
        blFrontNew->SetLayerNum(blFront->GetLayerNum() + 1);
        blFrontNew->SetSymm(blFront->GetSymm());

        m_blNxtFList->AddFront(blFrontNew);

        // set neighbors
        blFront->GetNeigbourFronts(&neigs, neigFrts);
        for (int j = 0; j < neigs; j++) {
            BLFront *ft = nullptr, *upft = nullptr;
            ft = neigFrts[j];
            upft = ft->GetUpperFront();
            if (upft) {
                blFrontNew->AddNeigbourFronts(upft);
                upft->AddNeigbourFronts(blFrontNew);
            }
        }

        CalFrontSize(blFrontNew);
        BLVector up_front_normal =
            -1 * BLVector(m_pNodes[conn[5]].coord[0] - m_pNodes[conn[3]].coord[0], m_pNodes[conn[5]].coord[1] - m_pNodes[conn[3]].coord[1],
                          m_pNodes[conn[5]].coord[2] - m_pNodes[conn[3]].coord[2]) ^
            BLVector(m_pNodes[conn[4]].coord[0] - m_pNodes[conn[3]].coord[0], m_pNodes[conn[4]].coord[1] - m_pNodes[conn[3]].coord[1],
                     m_pNodes[conn[4]].coord[2] - m_pNodes[conn[3]].coord[2]);
        blFrontNew->SetNormal(up_front_normal.normalized());
        // create a new mesh element
        int k = AddElem(6, conn, BLEntityTopology::PRISM);
        m_nPrism++;

        blFront->SetElmIdx(k);
        m_pElems[k].pointer = blFront;
    }
}
bool BLMesh::CheckPrismOrth(BLFront *blFront, bool include_top_face /*=false*/
)
{
    if (!blFront) {
        return false;
    }

    // 确保 conn 已经填好（你 PreCheckPrismValid 里也会调用，但这里再调用一次更稳）
    getConnection(blFront);
    int *conn = blFront->conn.data();

    auto vpos = [&](int nid) -> BLVector {
        if (nid < 0) {
            double nanv = std::numeric_limits<double>::quiet_NaN();
            return BLVector(nanv, nanv, nanv);
        }
        return BLVector(m_pNodes[nid].coord[0], m_pNodes[nid].coord[1], m_pNodes[nid].coord[2]);
    };

    auto isNan = [&](const BLVector &v) -> bool { return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z); };

    // 点：底(0,1,2) 顶(3,4,5)
    BLVector p0 = vpos(conn[0]), p1 = vpos(conn[1]), p2 = vpos(conn[2]);
    BLVector q0 = vpos(conn[3]), q1 = vpos(conn[4]), q2 = vpos(conn[5]);
    if (isNan(p0) || isNan(p1) || isNan(p2) || isNan(q0) || isNan(q1) || isNan(q2)) {
        return false;
    }

    // prism centroid
    BLVector C = (p0 + p1 + p2 + q0 + q1 + q2) / 6.0;

    auto triCentroid = [&](const BLVector &a, const BLVector &b, const BLVector &c) -> BLVector { return (a + b + c) / 3.0; };
    auto triAreaNormal = [&](const BLVector &a, const BLVector &b, const BLVector &c) -> BLVector { return (b - a) ^ (c - a); };
    auto quadCentroid = [&](const BLVector &a, const BLVector &b, const BLVector &c, const BLVector &d) -> BLVector {
        return (a + b + c + d) / 4.0;
    };
    auto quadAreaNormal = [&](const BLVector &a, const BLVector &b, const BLVector &c, const BLVector &d) -> BLVector {
        return ((b - a) ^ (c - a)) + ((c - a) ^ (d - a));
    };

    auto safeCos = [&](const BLVector &A, const BLVector &d) -> double {
        double AA2 = A.magnitude2();
        double dd2 = d.magnitude2();
        if (AA2 <= 1e-30 || dd2 <= 1e-30) {
            return 0.0;
        }
        double v = std::fabs(A * d) / (std::sqrt(AA2) * std::sqrt(dd2));
        if (v < 0.0) {
            v = 0.0;
        }
        if (v > 1.0) {
            v = 1.0;
        }
        return v;
    };

    // 由某个 front 的 conn 直接算它对应棱柱中心
    auto prismCentroidFromConn = [&](BLFront *f, BLVector &Cc) -> bool {
        if (!f || !f->GetUpperFront()) {
            return false;
        }
        getConnection(f);
        int *c = f->conn.data();
        BLVector a0 = vpos(c[0]), a1 = vpos(c[1]), a2 = vpos(c[2]);
        BLVector b0 = vpos(c[3]), b1 = vpos(c[4]), b2 = vpos(c[5]);
        if (isNan(a0) || isNan(a1) || isNan(a2) || isNan(b0) || isNan(b1) || isNan(b2)) {
            return false;
        }
        Cc = (a0 + a1 + a2 + b0 + b1 + b2) / 6.0;
        return true;
    };

    // 内部面 dual-cos：min( cos(A·(Cnb-C)), cos(A·(F-C)) )
    auto orthFaceDual = [&](const BLVector &A, const BLVector &Fc, const BLVector *Cnb) -> double {
        double cos_cf = safeCos(A, Fc - C);
        if (!Cnb) {
            return cos_cf;
        }
        double cos_cc = safeCos(A, (*Cnb) - C);
        return std::min(cos_cc, cos_cf);
    };

    // 找共享边(ida,idb)的邻 front（在 3 个邻居里找；可能非流形则取最差）
    auto worstOrthAcrossEdge = [&](int ida, int idb, const BLVector &Af, const BLVector &Fc) -> double {
        int neigs = 0;
        BLFront *neigFrts[DIM3] = {nullptr, nullptr, nullptr};
        blFront->GetNeigbourFronts(&neigs, neigFrts);

        double best = 1.0;
        bool hasInternal = false;

        for (int i = 0; i < neigs; i++) {
            BLFront *nb = neigFrts[i];
            if (!nb || !nb->GetUpperFront()) {
                continue;
            }

            // 判断 nb 是否包含 ida 和 idb
            int nn = 0;
            BLNode *ns[3] = {nullptr, nullptr, nullptr};
            nb->GetNodes(&nn, ns);
            if (nn != 3) {
                continue;
            }
            int a = ns[0]->GetNodIdx(), b = ns[1]->GetNodIdx(), c = ns[2]->GetNodIdx();
            bool hasA = (a == ida) || (b == ida) || (c == ida);
            bool hasB = (a == idb) || (b == idb) || (c == idb);
            if (!hasA || !hasB) {
                continue;
            }

            BLVector Cnb;
            if (!prismCentroidFromConn(nb, Cnb)) {
                continue;
            }

            best = std::min(best, orthFaceDual(Af, Fc, &Cnb));
            hasInternal = true;
        }

        if (!hasInternal) {
            best = orthFaceDual(Af, Fc, nullptr); // boundary
        }
        return best;
    };

    double minOrth = 1.0;

    // ===== 1) 底面（三角）=====
    {
        BLVector Fc = triCentroid(p0, p1, p2);
        BLVector Af = triAreaNormal(p0, p1, p2);

        BLFront *lower = blFront->GetLowerFront();
        BLVector Cnb;
        if (lower && prismCentroidFromConn(lower, Cnb)) {
            minOrth = std::min(minOrth, orthFaceDual(Af, Fc, &Cnb));
        } else {
            minOrth = std::min(minOrth, orthFaceDual(Af, Fc, nullptr));
        }
    }

    // ===== 2) 三个侧面（四边形）=====
    // 侧面按 wedge 标准连法： (0,1,4,3) (1,2,5,4) (2,0,3,5)
    {
        // face (0,1,4,3) 共享边(0,1)
        BLVector Fc = quadCentroid(p0, p1, q1, q0);
        BLVector Af = quadAreaNormal(p0, p1, q1, q0);
        minOrth = std::min(minOrth, worstOrthAcrossEdge(conn[0], conn[1], Af, Fc));
    }
    {
        // face (1,2,5,4) 共享边(1,2)
        BLVector Fc = quadCentroid(p1, p2, q2, q1);
        BLVector Af = quadAreaNormal(p1, p2, q2, q1);
        minOrth = std::min(minOrth, worstOrthAcrossEdge(conn[1], conn[2], Af, Fc));
    }
    {
        // face (2,0,3,5) 共享边(2,0)
        BLVector Fc = quadCentroid(p2, p0, q0, q2);
        BLVector Af = quadAreaNormal(p2, p0, q0, q2);
        minOrth = std::min(minOrth, worstOrthAcrossEdge(conn[2], conn[0], Af, Fc));
    }

    // ===== 3) 顶面（三角）可选（通常先不算/按边界近似）=====
    if (include_top_face) {
        BLVector Fc = triCentroid(q0, q1, q2);
        BLVector Af = triAreaNormal(q0, q1, q2);
        minOrth = std::min(minOrth, orthFaceDual(Af, Fc, nullptr));
    }

    return (minOrth >= 1 - cf.max_orth[0]);
}
BLVector BLMesh::TransNorm(BLVector norm, double *angl)
{
    BLVector normx(1.0, 0.0, 0.0), vtmp;
    double dval, thd, cosa;

    vtmp = normx ^ norm;

    if (vtmp.z >= 0) {
        cosa = normx * norm / (normx.magnitude() * norm.magnitude());
        thd = acos(cosa);
        thd = fmod(thd * 4.0, 2 * PI);
    } else {
        cosa = normx * norm / (normx.magnitude() * norm.magnitude());
        thd = 2 * PI - acos(cosa);
        thd = fmod(thd * 4.0, 2 * PI);
    }
    *angl = thd / PI * 180;

    vtmp.z = 0;
    vtmp.x = normx.x * cos(thd) - normx.y * sin(thd);
    vtmp.y = normx.x * sin(thd) + normx.y * cos(thd);

    vtmp.normalize();
    return vtmp;
}

double BLMesh::NormAngle(BLVector norm)
{
    BLVector normx(1.0, 0.0, 0.0), vtmp;
    double dval, thd, cosa, angl;

    vtmp = normx ^ norm;

    if (vtmp.z >= 0) {
        cosa = normx * norm / (normx.magnitude() * norm.magnitude());
        thd = acos(cosa);
    } else {
        cosa = normx * norm / (normx.magnitude() * norm.magnitude());
        thd = 2 * PI - acos(cosa);
    }
    angl = thd / PI * 180;

    return angl;
}

BLVector BLMesh::CalNorm(double coord0[3], double coord1[3], double coord2[3])
{
    BLVector vect, vecn, vectmp1, vectmp2;
    vectmp1.x = coord1[0] - coord0[0];
    vectmp1.y = coord1[1] - coord0[1];
    vectmp1.z = coord1[2] - coord0[2];

    vectmp2.x = coord2[0] - coord0[0];
    vectmp2.y = coord2[1] - coord0[1];
    vectmp2.z = coord2[2] - coord0[2];

    // outward normal
    vecn = vectmp2 ^ vectmp1;

    // inward normal
    // vecn = vectmp1^vectmp2;
    vecn.normalize();

    return vecn;
}

double BLMesh::AngleNorm(BLVector norm1, BLVector norm2)
{
    double cos, angle;
    cos = norm1 * norm2 / (norm1.magnitude() * norm2.magnitude());
    angle = acos(cos);
    return (180 * angle) / (PI);
}

void BLMesh::OutputAngle(char *filename, double *angle, int nFrtNods, int idx)
{
    char file1[256];
    double fw = 0.01, min = FLT_MAX, max = FLT_MIN;
    int i, cnt, *qcnt = nullptr;

    memset(file1, 0, sizeof(file1));
    sprintf(file1, "%s_%d.csv", filename, idx);

    cnt = (double)1.0 / fw;
    qcnt = new int[cnt + 1];

    for (i = 0; i < cnt + 1; i++) {
        qcnt[i] = 0;
    }

    for (i = 0; i < nFrtNods; i++) {
        cnt = angle[i] / fw;
        qcnt[cnt]++;

        if (angle[i] > max) {
            max = angle[i];
        }

        if (angle[i] < min) {
            min = angle[i];
        }
    }

    FILE *fout = fopen(file1, "w");
    if (!fout) {
        spdlog::info("error: cann't open file!\n");
        return;
    }
    cnt = (double)1.0 / fw;
    for (i = 0; i < cnt + 1; i++) {
        double raito = (double)qcnt[i] / nFrtNods;
        fprintf(fout, "%f, %d, %f\n", fw * i + 0.001, qcnt[i], raito * 100);
    }
    fprintf(fout, "min, %f, max, %f\n", min, max);

    fclose(fout);
    fout = nullptr;
}
// 3dversion

// 2dversion
bool BLMesh::IsConvexNode(double coord0[3], double coord1[3], double coord2[3])
{
    BLVector vecn, vectmp1, vectmp2;

    vectmp1.x = coord0[0] - coord2[0];
    vectmp1.y = coord0[1] - coord2[1];
    vectmp1.z = 0.0;

    vectmp2.x = coord1[0] - coord2[0];
    vectmp2.y = coord1[1] - coord2[1];
    vectmp2.z = 0.0;

    // outward normal
    vecn = vectmp2 ^ vectmp1;

    return vecn.z > 0;
}

void BLMesh::test()
{
    double A[3] = {1.266530, -0.256233, -0.040692}, B[3] = {1.339154, -0.255740, -0.016062}, C[3] = {1.280364, -0.253356, 0.051173},
           O[3] = {1.289233, -0.255227, 0.030649}, P[3] = {1.303460, -0.272334, 0.038875}, Q[3] = {1.295262, -0.254155, 0.024839};
    int k = tri_tri_inter(A, B, C, Q, P, O);
    spdlog::info("return={}", k);
    exit(0);
}

void BLMesh::CalFrontSize(BLFront *blFront)
{
    if (blFront->GetFrontSize() > 0) {
        return;
    }
    int i, pidx[MAX_NCONN], nNods;
    BLNode *blNods[MAX_NCONN];

    blFront->GetNodes(&nNods, blNods);
    BLVector vlen1(0, 0, 0), vlen2(0, 0, 0), vlen3(0, 0, 0);
    for (i = 0; i < nNods; i++) {
        pidx[i] = blNods[i]->GetNodIdx();
    }

    vlen1.x = m_pNodes[pidx[0]].coord[0] - m_pNodes[pidx[1]].coord[0];
    vlen1.y = m_pNodes[pidx[0]].coord[1] - m_pNodes[pidx[1]].coord[1];
    vlen1.z = m_pNodes[pidx[0]].coord[2] - m_pNodes[pidx[1]].coord[2];

    vlen2.x = m_pNodes[pidx[0]].coord[0] - m_pNodes[pidx[2]].coord[0];
    vlen2.y = m_pNodes[pidx[0]].coord[1] - m_pNodes[pidx[2]].coord[1];
    vlen2.z = m_pNodes[pidx[0]].coord[2] - m_pNodes[pidx[2]].coord[2];

    vlen3.x = m_pNodes[pidx[2]].coord[0] - m_pNodes[pidx[1]].coord[0];
    vlen3.y = m_pNodes[pidx[2]].coord[1] - m_pNodes[pidx[1]].coord[1];
    vlen3.z = m_pNodes[pidx[2]].coord[2] - m_pNodes[pidx[1]].coord[2];

    double vsize = (vlen1.magnitude() + vlen2.magnitude() + vlen3.magnitude()) / 3.0;
    double minvsize = std::min(std::min(vlen1.magnitude(), vlen2.magnitude()), vlen3.magnitude());
    double sqrtsize = pow(vlen1.magnitude() * vlen2.magnitude() * vlen3.magnitude(), 1.0 / 3);
    blFront->SetMinFrontSize(minvsize);
    blFront->SetSqrtFrontSize(sqrtsize);
    blFront->SetFrontSize(vsize);
}

// check whether is it need to generate pyramid
#if 0
void BLMesh::SetFrontStopflag(BLFront* blFront)
{
	int nNods, i, j;
	BLNode* blNods[MAX_FRONT_NODES], *blNod;
	int nflt, nnods, inei, conn[MAX_NCONN];
	BLFront* blFrnts[DIM3], *blLowerFrt = nullptr;
	BLNode* blnds[MAX_FRONT_NODES];
	bool bstop = false, bcrepramid = false;

	blFront->GetNodes(&nNods, blNods);

	//check neighbor fronts stopping flag
	blLowerFrt = blFront->GetLowerFront();
	if (blLowerFrt != nullptr)
	{
		int idx1, idx2, idx;
		blLowerFrt->GetNeigbourFronts(&nflt, blFrnts);
		for (i = 0; i < nflt; i++)
		{
			bcrepramid = false;
			bstop = false;

			BLFront* ft = blFrnts[i];

			NeighIdx(blLowerFrt, ft, &inei);
			idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
			idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();

			ft->GetNodes(&nnods, blnds);

			for (int j = 0; j < nnods; j++)
			{
				blNod = blnds[j];
				if (blNod->GetNodIdx() != idx1 && blNod->GetNodIdx() != idx2)
					idx = blNod->GetNodIdx();

				if (blNod->GetStopFlag())
				{
					bstop = true;
					bcrepramid = true;
					//break;
				}
			}

			if (bstop)
			{
				for (int j = 1; j < nNods; j++)
				{
					blNod = blNods[(inei + j) % DIM3];
					//blNod->SetStopFlag();
					StopPropagateNode(blNod);
				}
			}
#if 1
			if (bcrepramid)
			{
				//create a new mesh element
				conn[0] = idx1;
				conn[1] = blNods[(inei + 2) % DIM3]->GetNodIdx();
				conn[2] = blNods[(inei + 1) % DIM3]->GetNodIdx();
				conn[3] = idx2;
				conn[4] = idx;
				AddElem(5, conn, EntityTopology::PYRAMID);
				m_nPyramid++;

				std::vector<BLFront*>::iterator vit;
				vit = find(m_vBdyFront.begin(), m_vBdyFront.end(), ft);
				if (vit != m_vBdyFront.end())
					m_vBdyFront.erase(vit);

				outbdry[noutbdry * 3 + 0] = conn[1];
				outbdry[noutbdry * 3 + 1] = idx;
				outbdry[noutbdry * 3 + 2] = conn[2];
				noutbdry++;

				if (!(m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm))
				{
					outbdry[noutbdry * 3 + 0] = conn[1];
					outbdry[noutbdry * 3 + 1] = idx2;
					outbdry[noutbdry * 3 + 2] = idx;
					noutbdry++;
				}

				if (m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm)
				{

					if (idx < conn[1])
						m_mapSymline.insert(std::make_pair(idx, conn[1]));
					else
						m_mapSymline.insert(std::make_pair(conn[1], idx));

					/* 2019/03/26
					int con[DIM2];
					con[0] = idx;
					con[1] = conn[1];
					AddSymBdry(2, con, EntityTopology::LINE);

					if(m_pNodes[idx2].bsysm)
					{
						RmvSymBdry(idx, idx2);
						RmvSymBdry(idx2, conn[1]);
					}
					*/
				}

				if (!(m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm))
				{
					outbdry[noutbdry * 3 + 0] = conn[2];
					outbdry[noutbdry * 3 + 1] = idx;
					outbdry[noutbdry * 3 + 2] = idx1;
					noutbdry++;
				}

#ifdef _CHECK_INTERSECTION
				//intersection
				int cons[3], ii;
				cons[0] = conn[1];
				cons[1] = idx;
				cons[2] = conn[2];
				ii = AddTriElem(3, cons);
				m_ocTree->insert(ii, m_ocTree->getRootNode(), m_cbCube);

				cons[0] = conn[1];
				cons[1] = idx2;
				cons[2] = idx;
				ii = AddTriElem(3, cons);
				m_ocTree->insert(ii, m_ocTree->getRootNode(), m_cbCube);

				cons[0] = conn[2];
				cons[1] = idx;
				cons[2] = idx1;
				ii = AddTriElem(3, cons);
				m_ocTree->insert(ii, m_ocTree->getRootNode(), m_cbCube);

				m_ocTree->rmData(ft->GetTriIdx(), m_ocTree->getRootNode(), m_cbCube);
				ii = blLowerFrt->GetSTriIdx(inei, 0);
				m_ocTree->rmData(ii, m_ocTree->getRootNode(), m_cbCube);
				ii = blLowerFrt->GetSTriIdx(inei, 1);
				m_ocTree->rmData(ii, m_ocTree->getRootNode(), m_cbCube);
				//end of intersection
#endif

				if (m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm)
				{

					if (idx < conn[2])
						m_mapSymline.insert(std::make_pair(idx, conn[2]));
					else
						m_mapSymline.insert(std::make_pair(conn[2], idx));

					/* 2019/03/26
					int con[DIM2];
					con[0] = idx;
					con[1] = conn[2];
					AddSymBdry(2, con, EntityTopology::LINE);

					if(m_pNodes[idx1].bsysm)
					{
						RmvSymBdry(idx, idx1);
						RmvSymBdry(idx1, conn[2]);
					}
					*/
				}
			}
#endif
		}
		/*
		if (bstop)
		{
			for (i=0; i<nNods; i++)
			{
				blNod = blNods[i];
				blNod->SetStopFlag();
			}
		}*/
	}
}
#else
// check whether is it need to create pyramid elements in the neighbor of blFront

#endif

#if 0
void BLMesh::PostChckIntersect(BLFront* blFront)
{
	int nNods, i, j, tridx, ntri, itri[9], nlNods;
	BLNode* blNods[MAX_FRONT_NODES], *blNod;
	int nflt, nnods, inei, conn[MAX_NCONN];
	BLFront* blFrnts[DIM3], *blLowerFrt = nullptr;
	BLNode* blnds[MAX_FRONT_NODES], *blLoNods[MAX_FRONT_NODES];
	bool bstop = false, bcrepramid, ret = false;

	blFront->GetNodes(&nNods, blNods);

	ntri = 0;
	//check neighbor fronts stopping flag
	blLowerFrt = blFront->GetLowerFront();
	if (blLowerFrt != nullptr)
	{
		int idx1, idx2, idx;
		blLowerFrt->GetNeigbourFronts(&nflt, blFrnts);
		blLowerFrt->GetNodes(&nlNods, blLoNods);
		for (i = 0; i < nflt; i++)
		{
			bcrepramid = false;
			bstop = false;

			BLFront* ft = blFrnts[i];

			NeighIdx(blLowerFrt, ft, &inei);
			idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
			idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();

			ft->GetNodes(&nnods, blnds);

			for (int j = 0; j < nnods; j++)
			{
				blNod = blnds[j];
				if (blNod->GetNodIdx() != idx1 && blNod->GetNodIdx() != idx2)
					idx = blNod->GetNodIdx();

				if (blNod->GetStopFlag())
				{
					bstop = true;
					bcrepramid = true;
					//break;
				}
			}

			if (bcrepramid)
			{
				//create a new mesh element
				conn[0] = idx1;
				conn[1] = blNods[(inei + 2) % DIM3]->GetNodIdx();
				conn[2] = blNods[(inei + 1) % DIM3]->GetNodIdx();
				conn[3] = idx2;
				conn[4] = idx;

#ifdef _CHECK_INTERSECTION
				//intersection
				int cons[3], ii;
				cons[0] = conn[1];
				cons[1] = idx;
				cons[2] = conn[2];
				tridx = itri[ntri++] = AddTriElem(3, cons);
				m_ocTree->insertPreProcess(tridx);

				cons[0] = conn[1];
				cons[1] = idx2;
				cons[2] = idx;
				tridx = itri[ntri++] = AddTriElem(3, cons);
				m_ocTree->insertPreProcess(tridx);

				cons[0] = conn[2];
				cons[1] = idx;
				cons[2] = idx1;
				tridx = itri[ntri++] = AddTriElem(3, cons);
				m_ocTree->insertPreProcess(tridx);

				
				//m_ocTree->rmData(ft->GetTriIdx(), m_ocTree->getRootNode(), m_cbCube);
				//end of intersection
#endif
/*
				ret = false;
				for (int j=0; j<3; j++)
					ret |= m_ocTree->chckIntersect(itri[j], m_cbCube);

				if(ret)
				{
					for (int k=0; k<DIM3; k++)
					{
						blLoNods[k]->SetStopFlag();
						RmvUperNeigFronts(blLoNods[k], false);
					}

					//add to the exposed fronts segments
					//m_vBdyFront.push_back(blLowerFrt);

					for (int k=0; k<DIM3; k++)
						m_ocTree->rmData(itri[k], m_ocTree->getRootNode(), m_cbCube);

					break;
				}
				else
				{
					//remove
				}*/
			}
		}

		for (int j = 0; j < ntri; j++)
			ret |= m_ocTree->chckIntersectPreProcess(itri[j]);

		if (ret)
		{
#if 0
			for (int k = 0; k < DIM3; k++)
			{
				RmvUperNeigFronts(blLoNods[k], false);
				/*
				blLoNods[k]->SetStopFlag();

				BLNode* tmp = blLoNods[k]->GetUpperNode();
				if(tmp)
				{
					tmp->RmvLowerNod(blLoNods[k]);
					blLoNods[k]->RmvUpperNod(tmp);
				}*/
				StopPropagateNode(blLoNods[k]);
			}
#else

			for (int k = 1; k < DIM3; k++)
			{
				RmvUperNeigFronts(blLoNods[(inei + k) % DIM3], false);
#ifdef DEBUG
				cout << "Debug info post check" << blNod->GetDescentNod()->GetNodIdx() << endl;
#endif
				StopPropagateNode(blLoNods[(inei + k) % DIM3]);
			}
#endif

			for (int k = 0; k < ntri; k++) {
				m_ocTree->rmDataPreProcess(itri[k]);
			}
		}
	}
}
#endif
double BLMesh::CheckPyramidVolumn(double coordinates[][3])
{
    double volumn1 = 0, volumn2 = 0;
    BLVector side1, side2, side3;

    /* Update: New calculation (average two ways of division of bottom face) */

    // divide the pyramid into 2 tets (0124 + 0234) and calculate each
    // volume of the 0124 tet
    side1.set(coordinates[1][0] - coordinates[0][0], coordinates[1][1] - coordinates[0][1], coordinates[1][2] - coordinates[0][2]);

    side2.set(coordinates[2][0] - coordinates[0][0], coordinates[2][1] - coordinates[0][1], coordinates[2][2] - coordinates[0][2]);

    side3.set(coordinates[4][0] - coordinates[0][0], coordinates[4][1] - coordinates[0][1], coordinates[4][2] - coordinates[0][2]);

    volumn1 = (side3 * (side1 ^ side2)) / 6.0;

    // volume of the 0234 tet
    side1.set(coordinates[2][0] - coordinates[0][0], coordinates[2][1] - coordinates[0][1], coordinates[2][2] - coordinates[0][2]);

    side2.set(coordinates[3][0] - coordinates[0][0], coordinates[3][1] - coordinates[0][1], coordinates[3][2] - coordinates[0][2]);

    side3.set(coordinates[4][0] - coordinates[0][0], coordinates[4][1] - coordinates[0][1], coordinates[4][2] - coordinates[0][2]);

    volumn1 += (side3 * (side1 ^ side2)) / 6.0;

    // divide the pyramid into 2 tets (0134 + 1234) and calculate each
    // volume of the 0134 tet
    side1.set(coordinates[1][0] - coordinates[0][0], coordinates[1][1] - coordinates[0][1], coordinates[1][2] - coordinates[0][2]);

    side2.set(coordinates[3][0] - coordinates[0][0], coordinates[3][1] - coordinates[0][1], coordinates[3][2] - coordinates[0][2]);

    side3.set(coordinates[4][0] - coordinates[0][0], coordinates[4][1] - coordinates[0][1], coordinates[4][2] - coordinates[0][2]);

    volumn2 = (side3 * (side1 ^ side2)) / 6.0;

    // volume of the 1234 tet
    side1.set(coordinates[2][0] - coordinates[1][0], coordinates[2][1] - coordinates[1][1], coordinates[2][2] - coordinates[1][2]);

    side2.set(coordinates[3][0] - coordinates[1][0], coordinates[3][1] - coordinates[1][1], coordinates[3][2] - coordinates[1][2]);

    side3.set(coordinates[4][0] - coordinates[1][0], coordinates[4][1] - coordinates[1][1], coordinates[4][2] - coordinates[1][2]);

    volumn2 += (side3 * (side1 ^ side2)) / 6.0;

    if (volumn1 < volumn2) {
        return volumn1 / 2;
    }

    return volumn2 / 2;
}
bool BLMesh::CheckPyramidSkewness(double coordinates[][3])
{
    using Eigen::Vector3d;
    static constexpr const double pi = 3.14159265358979323846;
    Vector3d a(coordinates[0][0], coordinates[0][1], coordinates[0][2]);
    Vector3d b(coordinates[1][0], coordinates[1][1], coordinates[1][2]);
    Vector3d c(coordinates[2][0], coordinates[2][1], coordinates[2][2]);
    Vector3d d(coordinates[3][0], coordinates[3][1], coordinates[3][2]);
    Vector3d e(coordinates[4][0], coordinates[4][1], coordinates[4][2]);

    // 向量夹角（弧度）
    auto angleBetween = [](const Vector3d &u, const Vector3d &v) -> double {
        double nu = u.norm();
        double nv = v.norm();
        if (nu <= std::numeric_limits<double>::epsilon() || nv <= std::numeric_limits<double>::epsilon()) {
            return 0.0; // 退化角，当 0 处理
        }
        double c = u.dot(v) / (nu * nv);
        c = std::max(std::min(c, 1.0), -1.0);
        return std::acos(c);
    };

    // 给一组角度 + 理想角，算等角偏斜（单位：弧度）
    auto equiangleSkewFromAngles = [](const std::vector<double> &angles, double idealAngle) -> double {
        if (angles.empty()) {
            return 0.0;
        }

        double theta_min = angles[0];
        double theta_max = angles[0];
        for (size_t i = 1; i < angles.size(); ++i) {
            if (angles[i] < theta_min) {
                theta_min = angles[i];
            }
            if (angles[i] > theta_max) {
                theta_max = angles[i];
            }
        }
        double s1 = (theta_max - idealAngle) / (pi - idealAngle);
        double s2 = (idealAngle - theta_min) / idealAngle;

        double skew = s1 > s2 ? s1 : s2;
        if (skew < 0.0) {
            skew = 0.0;
        }
        if (skew > 1.0) {
            skew = 1.0;
        }
        return skew;
    };

    // 底面四边形 a-b-c-d 的等角偏斜
    std::vector<double> quadAngles;
    quadAngles.reserve(4);
    // 顶点 a：角(b-a, d-a)
    quadAngles.push_back(angleBetween(b - a, d - a));
    // 顶点 b：角(c-b, a-b)
    quadAngles.push_back(angleBetween(c - b, a - b));
    // 顶点 c：角(d-c, b-c)
    quadAngles.push_back(angleBetween(d - c, b - c));
    // 顶点 d：角(a-d, c-d)
    quadAngles.push_back(angleBetween(a - d, c - d));

    double quadIdeal = pi / 2.0; // 90°
    double baseSkew = equiangleSkewFromAngles(quadAngles, quadIdeal);

    // 四个侧面三角形的等角偏斜
    auto triangleSkew = [&](const Vector3d &p0, const Vector3d &p1, const Vector3d &p2) -> double {
        std::vector<double> angs;
        angs.reserve(3);

        angs.push_back(angleBetween(p1 - p0, p2 - p0)); // at p0
        angs.push_back(angleBetween(p0 - p1, p2 - p1)); // at p1
        angs.push_back(angleBetween(p0 - p2, p1 - p2)); // at p2

        double triIdeal = pi / 3.0;                     // 60°
        return equiangleSkewFromAngles(angs, triIdeal);
    };

    double tri1 = triangleSkew(a, b, e); // 面 (0,1,4)
    double tri2 = triangleSkew(b, c, e); // 面 (1,2,4)
    double tri3 = triangleSkew(c, d, e); // 面 (2,3,4)
    double tri4 = triangleSkew(d, a, e); // 面 (3,0,4)

    double equal_angle_skewness = std::max({baseSkew, tri1, tri2, tri3, tri4});

    if (cf.max_skewness[1] < 0.1) {
        throw(std::logic_error("maximum equal skewnwass is too small!"));
    }
    if (equal_angle_skewness > cf.max_skewness[1]) {
#ifdef _DEBUG
        cout << "pyramid equal angle skewness failed: " << equal_angle_skewness << endl;
#endif
        return false;
    }

    return true;
}
bool BLMesh::CheckPyramidOrth(double coordinates[][3], double neighcenter[3])
{
    using Eigen::Vector3d;

    const Vector3d p0(coordinates[0][0], coordinates[0][1], coordinates[0][2]);
    const Vector3d p1(coordinates[1][0], coordinates[1][1], coordinates[1][2]);
    const Vector3d p2(coordinates[2][0], coordinates[2][1], coordinates[2][2]);
    const Vector3d p3(coordinates[3][0], coordinates[3][1], coordinates[3][2]);
    const Vector3d p4(coordinates[4][0], coordinates[4][1], coordinates[4][2]);

    const Vector3d Cnei(neighcenter[0], neighcenter[1], neighcenter[2]);

    auto safeNormalize = [](const Vector3d &v) -> Vector3d {
        double n = v.norm();
        if (n <= std::numeric_limits<double>::epsilon()) {
            return Vector3d::Zero();
        }
        return v / n;
    };

    // cell centroid（简单平均足够用于正交性）
    const Vector3d C = (p0 + p1 + p2 + p3 + p4) / 5.0;

    // shared face: base quad (0,1,2,3)
    const Vector3d F = (p0 + p1 + p2 + p3) / 4.0;

    // base normal（用两个三角的法向相加，较稳）
    Vector3d n = (p1 - p0).cross(p2 - p0) + (p2 - p0).cross(p3 - p0);
    Vector3d n_u = safeNormalize(n);
    if (n_u.isZero(0)) {
#ifdef _DEBUG
        std::cout << "pyramid orth failed: base normal degenerate\n";
#endif
        return false;
    }

    Vector3d dCC_u = safeNormalize(Cnei - C); // centroid -> neighbor centroid
    if (dCC_u.isZero(0)) {
#ifdef _DEBUG
        std::cout << "pyramid orth failed: neighbor center coincides with cell center\n";
#endif
        return false;
    }

    Vector3d dCF_u = safeNormalize(F - C); // centroid -> face centroid
    if (dCF_u.isZero(0)) {
#ifdef _DEBUG
        std::cout << "pyramid orth failed: face centroid coincides with cell center\n";
#endif
        return false;
    }

    double cos_cc = std::abs(n_u.dot(dCC_u));
    double cos_cf = std::abs(n_u.dot(dCF_u));

    // orthogonality for this shared face
    double orth = std::min(cos_cc, cos_cf);

    // 可选：数值夹紧
    orth = std::max(0.0, std::min(1.0, orth));

    // 阈值检查（你自己把名字对上）
    if (cf.max_orth[1] < 1e-6) {
        throw std::logic_error("min_pyramid_orth is too small / not set!");
    }

    if (orth < 1 - cf.max_orth[1]) {
#ifdef _DEBUG
        std::cout << "pyramid orth failed: orth=" << 1 - orth << "\n";
#endif
        return false;
    }

    return true;
}

bool BLMesh::ChckIntersectforTransit(BLFront *blFront)
{
    int nNods, ntri = 0;
    int nlowNods;
    int nnods; // 邻接面点数
    int inei;  // 当前正在使用的点i
    int tridx, itri[9];
    int conn[MAX_NCONN];
    BLNode *blNods[MAX_FRONT_NODES], *blLowNods[MAX_FRONT_NODES], *blLownNods[MAX_FRONT_NODES];
    BLNode *blNod;
    BLFront *blLowFrt = nullptr;
    bool bstop, bcrepramid;
    bool ret = false;

    // Get nods and front
    blFront->GetNodes(&nNods, blNods);
    blLowFrt = blFront->GetLowerFront(); // 下一层的front

    if (blLowFrt != nullptr) {
        int idx1, idx2, idx3, idx = 0;
        auto &blLowNeigFrnts = (blLowFrt)->m_arrNeigFronts;
        blLowFrt->GetNodes(&nlowNods, blLowNods);
        bool used[3] = {false};

        for (int i = 0; i < blLowFrt->m_nNeiFront; i++) {
            bcrepramid = false;
            bstop = false;
            BLFront *neig_front = blLowNeigFrnts[i];

            // 寻找neig不包含的那个点
            NeighIdx(blLowFrt, neig_front, &inei);
            used[inei] = true;

            idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
            idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();
            idx3 = blNods[(inei + 3) % DIM3]->GetLowerNode()->GetNodIdx();

            neig_front->GetNodes(&nnods, blLownNods);
            for (int j = 0; j < nnods; j++) {
                blNod = blLownNods[j];
                if (blNod->GetNodIdx() != idx1 && blNod->GetNodIdx() != idx2) {
                    idx = blNod->GetNodIdx();
                }

                if (blNod->GetStopFlag()) {
                    bstop = true;
                    bcrepramid = true;
                    // break;
                }
            }

            if (bcrepramid) {
                // create a new mesh element
                conn[0] = idx;                                    // 下层邻接面不相交点
                conn[1] = idx1;                                   // 邻接1号点
                conn[2] = idx2;                                   // 邻接2号点
                conn[3] = idx3;                                   // 下层本面不相交点
                conn[4] = blNods[(inei + 1) % DIM3]->GetNodIdx(); // 上层本面1号点
                conn[5] = blNods[(inei + 2) % DIM3]->GetNodIdx(); // 上层本面2号点
                conn[6] = blNods[(inei + 3) % DIM3]->GetNodIdx(); // 上层本面不相交点
#ifdef CHECK_VOLUMN
                // volumn
                double coordp[5][3];
                double neighCenter[3] = {0};
                for (int k = 0; k < 3; k++) {
                    coordp[0][k] = m_pNodes[conn[1]].coord[k];
                    coordp[1][k] = m_pNodes[conn[4]].coord[k];
                    coordp[2][k] = m_pNodes[conn[5]].coord[k];
                    coordp[3][k] = m_pNodes[conn[2]].coord[k];
                    coordp[4][k] = m_pNodes[conn[0]].coord[k];
                }
                for (int k = 0; k < 3; k++) {
                    neighCenter[k] += m_pNodes[conn[1]].coord[k];
                    neighCenter[k] += m_pNodes[conn[2]].coord[k];
                    neighCenter[k] += m_pNodes[conn[3]].coord[k];
                    neighCenter[k] += m_pNodes[conn[4]].coord[k];
                    neighCenter[k] += m_pNodes[conn[5]].coord[k];
                    neighCenter[k] += m_pNodes[conn[6]].coord[k];
                    neighCenter[k] /= 6;
                }

                double volumn = CheckPyramidVolumn(coordp);
                if (volumn < 1e-25) {
                    ret = true;
#ifdef _DEBUG
                    cout << "pyramid volumn=" << volumn << endl;
#endif // DEBUG
                }

                if (!CheckPyramidSkewness(coordp) || !CheckPyramidOrth(coordp, neighCenter)) {
                    ret = true;
                }
#endif

                int cons[3];
                cons[0] = conn[0];
                cons[1] = conn[1];
                cons[2] = conn[4];
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                // m_ocTree_symm->insertPreProcess(tridx);

                cons[0] = conn[0];
                cons[1] = conn[2];
                cons[2] = conn[5];
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                // m_ocTree_symm->insertPreProcess(tridx);

                cons[0] = conn[0];
                cons[1] = conn[4];
                cons[2] = conn[5];
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                // m_ocTree_symm->insertPreProcess(tridx);
            }
        }

        if (cf.max_layer_diff != 1) {
            int i;
            for (i = 0; i < 3; i++) {
                if (blNods[i]->GetBSys()) {
                    break;
                }
            }
            if (i == 3) {
                for (i = 0; i < 3; i++) {
                    if (!used[i]) {
                        bcrepramid = true;
                        bstop = true;

                        inei = i;
                        used[inei] = true;
                        idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
                        idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();
                        idx3 = blNods[(inei + 3) % DIM3]->GetLowerNode()->GetNodIdx();
                        idx = blLowFrt->ps.node_idx[inei];
                        int layer = blLowFrt->ps.exposed_layer[inei];
                        int sum = 0;

                        if (bcrepramid) {
                            // create a new mesh element
                            conn[0] = idx;
                            conn[1] = idx1;
                            conn[2] = idx2;
                            conn[3] = idx3;
                            conn[4] = blNods[(inei + 1) % DIM3]->GetNodIdx();
                            conn[5] = blNods[(inei + 2) % DIM3]->GetNodIdx();
                            conn[6] = blNods[(inei + 3) % DIM3]->GetNodIdx();
#ifdef CHECK_VOLUMN
                            // volumn
                            double coordp[5][3];
                            double neighCenter[3] = {0};
                            for (int k = 0; k < 3; k++) {
                                coordp[0][k] = m_pNodes[conn[1]].coord[k];
                                coordp[1][k] = m_pNodes[conn[4]].coord[k];
                                coordp[2][k] = m_pNodes[conn[5]].coord[k];
                                coordp[3][k] = m_pNodes[conn[2]].coord[k];
                                coordp[4][k] = m_pNodes[conn[0]].coord[k];
                            }
                            for (int k = 0; k < 3; k++) {
                                neighCenter[k] += m_pNodes[conn[1]].coord[k];
                                neighCenter[k] += m_pNodes[conn[2]].coord[k];
                                neighCenter[k] += m_pNodes[conn[3]].coord[k];
                                neighCenter[k] += m_pNodes[conn[4]].coord[k];
                                neighCenter[k] += m_pNodes[conn[5]].coord[k];
                                neighCenter[k] += m_pNodes[conn[6]].coord[k];
                                neighCenter[k] /= 6;
                            }
                            double volumn = CheckPyramidVolumn(coordp);

                            if (volumn < 1e-25) {
                                ret = true;
#ifdef _DEBUG
                                cout << "pyramid volumn=" << volumn << endl;
#endif // DEBUG
                            }
                            if (!CheckPyramidSkewness(coordp) || !CheckPyramidOrth(coordp, neighCenter)) {
                                ret = true;
                            }

#endif
                            // 通过两个策略控制 1. 最大层差 2.层差和横向的比值
                            int k = 1;
                            BLNode *ptr = blFront->m_pBLNods[inei];
                            const double growth_ratio = ptr->respect_ratio > 0.0 ? ptr->respect_ratio : cf.ratio2;
                            double thick = ptr->GetHeightLength();
                            const int layer_diff = ptr->GetLayerNum() - layer;
                            if (layer_diff > 0) {
                                thick /= pow(growth_ratio, layer_diff);
                            }
                            double thick_total = thick;
                            while (thick_total < blFront->GetMinFrontSize() * cf.max_ratio_diff) {
                                thick *= growth_ratio;
                                k++;
                                thick_total += thick;
                            }
                            // cout << min(k, cf.max_layer_diff)<<" "<<k<<":";
                            if (blFront->GetLayerNum() > layer + min(k, cf.max_layer_diff)) {
                                ret = true;
                            }

                            int cons[3];
                            cons[0] = conn[0];
                            cons[1] = conn[1];
                            cons[2] = conn[4];
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);
                            // m_ocTree_symm->insertPreProcess(tridx);

                            cons[0] = conn[0];
                            cons[1] = conn[4];
                            cons[2] = conn[5];
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);
                            // m_ocTree_symm->insertPreProcess(tridx);

                            cons[0] = conn[0];
                            cons[1] = conn[5];
                            cons[2] = conn[2];
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);
                            // m_ocTree_symm->insertPreProcess(tridx);
                        }
                    }
                }
            }
        }

#ifdef _DEBUG
        m_ocTree->num_inter += ntri;
#endif
        for (int j = 0; j < ntri; j++) {
            if (!ret) {
                ret = m_ocTree->chckIntersectPreProcess(itri[j]) /*|| m_ocTree_symm->chckIntersectPreProcess(itri[j])*/;
            }
        }

        if (ret) // intersections happen
        {
            for (int k = 1; k < DIM3; k++) {
                RmvUperNeigFrontsAndFreeNode(blLowNods[(inei + k) % DIM3]);
                StopPropagateNode(blLowNods[(inei + k) % DIM3]);
            }

            for (int k = 0; k < ntri; k++) {
                m_ocTree->rmDataPreProcess(itri[k]);
                // m_ocTree_symm->rmDataPreProcess(itri[k]);
            }

            // m_ocTree->insert(ft->GetTriIdx(), m_ocTree->getRootNode(), m_cbCube);
        } else {
            if (bcrepramid) {
                blFront->SetPyramidFlag(true);
            }
        }
    }

    return ret;
}

void BLMesh::CreatePyramid(BLFront *blFront)
{
    int nNods, i, j;
    BLNode *blNods[MAX_FRONT_NODES], *blNod;
    int nflt, nnods, inei, conn[MAX_NCONN];
    BLFront *blLowerFrt = nullptr;
    BLNode *blnds[MAX_FRONT_NODES];
    bool bstop = false, bcrepramid;

    blFront->GetNodes(&nNods, blNods);

    // check neighbor fronts stopping flag
    blLowerFrt = blFront->GetLowerFront();
    bool used[3] = {false};
    if (blLowerFrt != nullptr) {
        int idx1, idx2, idx = 0;
        auto &blFrnts = (blLowerFrt)->m_arrNeigFronts;
        blFront->ps = blLowerFrt->ps;
        for (i = 0; i < blLowerFrt->m_nNeiFront; i++) {
            bcrepramid = false;
            bstop = false;

            BLFront *lower_neigbour = blFrnts[i]; // lower neigbour

            NeighIdx(blLowerFrt, lower_neigbour, &inei);
            used[inei] = true;
            idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
            idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();

            lower_neigbour->GetNodes(&nnods, blnds);
            int sum = 0;
            for (int j = 0; j < 3; j++) {
                blNod = blnds[j];
                if (blNod->GetNodIdx() != idx1 && blNod->GetNodIdx() != idx2) {
                    idx = blNod->GetNodIdx();
                }
                /// idx 是下层的第三点

                if (blNod->GetStopFlag()) {
                    bstop = true;
                    bcrepramid = true;
                    // break;
                }
            }
            if (bstop) {
                bool in_sys = false;
                for (int j = 1; j < nNods; j++) {
                    blNod = blNods[(inei + j) % DIM3];
                    // cout << "Debug info create pyramid" << blNod->GetDescentNod()->GetNodIdx() << endl;
                    if (cf.max_layer_diff == 1) {
                        StopPropagateNode(blNod);
                    }
                    if (blNod->GetBSys()) {
                        in_sys = true;
                    }
                }
                if (in_sys || blNods[inei]->GetBSys()) {
                    for (int j = 1; j < nNods; j++) {
                        blNod = blNods[(inei + j) % DIM3];
                        StopPropagateNode(blNod);
                    }
                }
            }

            if (bcrepramid) {

                // create a new mesh element
                conn[0] = idx1;
                conn[1] = blNods[(inei + 2) % DIM3]->GetNodIdx();
                conn[2] = blNods[(inei + 1) % DIM3]->GetNodIdx();
                conn[3] = idx2;
                conn[4] = idx;
                AddElem(5, conn, BLEntityTopology::PYRAMID);
                m_nPyramid++;
                blFront->ps.node_idx[inei] = idx;
                blFront->ps.exposed_layer[inei] = blNods[(inei + 2) % DIM3]->GetLayerNum() - 1;

                lower_neigbour->is_boundary_ = false;

                outbdry[noutbdry * 3 + 0] = conn[1];
                outbdry[noutbdry * 3 + 1] = idx;
                outbdry[noutbdry * 3 + 2] = conn[2];

                if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
                    throw std::runtime_error("try to input error");
                }

                noutbdry++;

                if (!(m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm)) {
                    outbdry[noutbdry * 3 + 0] = conn[1];
                    outbdry[noutbdry * 3 + 1] = idx2;
                    outbdry[noutbdry * 3 + 2] = idx;
                    if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
                        throw std::runtime_error("try to input error");
                    }

                    noutbdry++;
                }

                if (m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm) {
                    auto smaller = std::min(idx, conn[1]);
                    auto biger = std::max(idx, conn[1]);
                    // m_mapSymline.insert(std::make_pair(smaller, biger));
                    m_pyramid_symline.insert(std::array<int, 2>{smaller, biger});
                    // insertSymFace(std::make_pair(smaller, biger));
                    int conn_sym[3];
                    conn_sym[0] = idx;
                    conn_sym[1] = conn[1];
                    conn_sym[2] = idx2;

                    auto face = intersect_two(intersect_two(m_pNodes[conn_sym[0]].isymfc, m_pNodes[conn_sym[1]].isymfc),
                                              m_pNodes[conn_sym[2]].isymfc);
                    if (face.size() == 0) {
                        outbdry[noutbdry * 3 + 0] = conn[1];
                        outbdry[noutbdry * 3 + 1] = idx2;
                        outbdry[noutbdry * 3 + 2] = idx;
                        if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
                            throw std::runtime_error("try to input error");
                        }
                        noutbdry++;
                    } else {
                        AddElem(3, conn_sym, BLEntityTopology::TRIANGLE, face.size() ? face[0] : 0);
                    }
                }

                if (!(m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm)) {
                    outbdry[noutbdry * 3 + 0] = conn[2];
                    outbdry[noutbdry * 3 + 1] = idx;
                    outbdry[noutbdry * 3 + 2] = idx1;
                    if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
                        throw std::runtime_error("try to input error");
                    }
                    noutbdry++;
                }

                if (m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm) {

                    auto smaller = std::min(idx, conn[2]);
                    auto biger = std::max(idx, conn[2]);
                    // insertSymFace(std::make_pair(smaller, biger));
                    m_pyramid_symline.insert(std::array<int, 2>{smaller, biger});

                    int conn_sym[3];
                    conn_sym[0] = idx;
                    conn_sym[1] = conn[2];
                    conn_sym[2] = idx1;

                    auto face = intersect_two(intersect_two(m_pNodes[conn_sym[0]].isymfc, m_pNodes[conn_sym[1]].isymfc),
                                              m_pNodes[conn_sym[2]].isymfc);
                    if (face.size() == 0) {
                        outbdry[noutbdry * 3 + 0] = conn[2];
                        outbdry[noutbdry * 3 + 1] = idx;
                        outbdry[noutbdry * 3 + 2] = idx1;
                        if (outbdry[noutbdry * 3 + 0] < 0 || outbdry[noutbdry * 3 + 1] < 0 || outbdry[noutbdry * 3 + 2] < 0) {
                            throw std::runtime_error("try to input error");
                        }
                        noutbdry++;
                    } else {
                        AddElem(3, conn_sym, BLEntityTopology::TRIANGLE, face.size() ? face[0] : 0);
                    }
                }

#ifdef _CHECK_INTERSECTION
                int tridx;
                int itri[13];
                int ntri = 0;
                int cons[3], ii;
                /*cons[0] = conn[1];
                cons[1] = idx;
                cons[2] = conn[2];
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);

                if (m_ocTree->chckIntersectPreProcess(tridx)) {
                    auto b = (ChckIntersectforTransit(blFront));
                    cout << blFront->GetElmIdx() << blFront->GetLayerNum() << endl;
                    cout << b << endl;
                    exit(-77);
                }
                cons[0] = conn[1];
                cons[1] = idx2;
                cons[2] = idx;
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                if (m_ocTree->chckIntersectPreProcess(tridx)) {
                    ChckIntersectforTransit(blFront);
                }
                cons[0] = conn[2];
                cons[1] = idx;
                cons[2] = idx1;
                tridx = itri[ntri++] = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                if (m_ocTree->chckIntersectPreProcess(tridx)) {
                    ChckIntersectforTransit(blFront);
                }*/

                int itrii;
                m_ocTree->rmDataPreProcess(lower_neigbour->GetTriIdx());

                itrii = blLowerFrt->GetSTriIdx(inei, 0);
                m_ocTree->rmDataPreProcess(itrii);
                itrii = blLowerFrt->GetSTriIdx(inei, 1);
                m_ocTree->rmDataPreProcess(itrii);

#endif
            }
        }
        if (cf.max_layer_diff != 1) {

            for (i = 0; i < 3; i++) {
                if (blNods[i]->GetBSys()) {
                    break;
                }
            }
            if (i == 3) {
                for (i = 0; i < 3; i++) {
                    if (!used[i]) {
                        bcrepramid = true;
                        bstop = true;

                        inei = i;
                        used[inei] = true;
                        idx1 = blNods[(inei + 1) % DIM3]->GetLowerNode()->GetNodIdx();
                        idx2 = blNods[(inei + 2) % DIM3]->GetLowerNode()->GetNodIdx();

                        idx = blLowerFrt->ps.node_idx[inei];
                        int sum = 0;
                        if (bstop) {
                            for (int j = 1; j < nNods; j++) {
                                blNod = blNods[(inei + j) % DIM3];
                                // cout << "Debug info create pyramid" << blNod->GetDescentNod()->GetNodIdx() << endl;
                            }
                        }

                        if (bcrepramid) {

                            // create a new mesh element
                            conn[0] = idx1;
                            conn[1] = blNods[(inei + 2) % DIM3]->GetNodIdx();
                            conn[2] = blNods[(inei + 1) % DIM3]->GetNodIdx();
                            conn[3] = idx2;
                            conn[4] = idx;
                            AddElem(5, conn, BLEntityTopology::PYRAMID);
                            m_nPyramid++;
                            blFront->ps.node_idx[inei] = idx;

                            outbdry[noutbdry * 3 + 0] = conn[1];
                            outbdry[noutbdry * 3 + 1] = idx;
                            outbdry[noutbdry * 3 + 2] = conn[2];
                            noutbdry++;

                            boundary_to_delete_.push_back(idx2);
                            boundary_to_delete_.push_back(idx);
                            boundary_to_delete_.push_back(idx1);

                            if (!(m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm)) {
                                outbdry[noutbdry * 3 + 0] = conn[1];
                                outbdry[noutbdry * 3 + 1] = idx2;
                                outbdry[noutbdry * 3 + 2] = idx;
                                noutbdry++;
                            }

                            if (m_pNodes[idx].bsysm && m_pNodes[conn[1]].bsysm) {
                                auto smaller = std::min(idx, conn[1]);
                                auto biger = std::max(idx, conn[1]);
                                // m_mapSymline.insert(std::make_pair(smaller, biger));
                                m_pyramid_symline.insert(std::array<int, 2>{smaller, biger});
                                m_pyramid_symline.erase(std::array<int, 2>{std::min(idx, idx2), std::max(idx, idx2)});
                                // insertSymFace(std::make_pair(smaller, biger));
                                int conn_sym[3];
                                conn_sym[0] = idx;
                                conn_sym[1] = conn[1];
                                conn_sym[2] = idx2;

                                auto face = intersect_two(intersect_two(m_pNodes[conn_sym[0]].isymfc, m_pNodes[conn_sym[1]].isymfc),
                                                          m_pNodes[conn_sym[2]].isymfc);
                                if (face.size() == 0) {
                                    for (int k = 0; k < 3; k++) {
                                        for (auto l : m_pNodes[conn_sym[k]].isymfc) {
                                            face.push_back(l);
                                        }
                                    }
                                }
                                AddElem(3, conn_sym, BLEntityTopology::TRIANGLE, face.size() ? face[0] : 0);
                            }

                            if (!(m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm)) {
                                outbdry[noutbdry * 3 + 0] = conn[2];
                                outbdry[noutbdry * 3 + 1] = idx;
                                outbdry[noutbdry * 3 + 2] = idx1;
                                noutbdry++;
                            }

                            if (m_pNodes[idx].bsysm && m_pNodes[conn[2]].bsysm) {

                                auto smaller = std::min(idx, conn[2]);
                                auto biger = std::max(idx, conn[2]);
                                // insertSymFace(std::make_pair(smaller, biger));
                                m_pyramid_symline.insert(std::array<int, 2>{smaller, biger});
                                m_pyramid_symline.erase(std::array<int, 2>{std::min(idx, idx1), std::max(idx, idx1)});

                                int conn_sym[3];
                                conn_sym[0] = idx;
                                conn_sym[1] = conn[2];
                                conn_sym[2] = idx1;
                                auto face = intersect_two(intersect_two(m_pNodes[conn_sym[0]].isymfc, m_pNodes[conn_sym[1]].isymfc),
                                                          m_pNodes[conn_sym[2]].isymfc);
                                if (face.size() == 0) {
                                    for (int k = 0; k < 3; k++) {
                                        for (auto l : m_pNodes[conn_sym[k]].isymfc) {
                                            face.push_back(l);
                                        }
                                    }
                                }
                                AddElem(3, conn_sym, BLEntityTopology::TRIANGLE, face.size() ? face[0] : 0);
                            }

#ifdef _CHECK_INTERSECTION
                            int tridx;
                            int itri[13];
                            int ntri = 0;
                            int cons[3], ii;
                            /*cons[0] = conn[1];
                            cons[1] = idx;
                            cons[2] = conn[2];
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);

                            if (m_ocTree->chckIntersectPreProcess(tridx)) {
                                auto b = (ChckIntersectforTransit(blFront));
                                cout << blFront->GetElmIdx() << blFront->GetLayerNum() << endl;
                                cout << b << endl;
                                exit(-77);
                            }
                            cons[0] = conn[1];
                            cons[1] = idx2;
                            cons[2] = idx;
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);
                            if (m_ocTree->chckIntersectPreProcess(tridx)) {
                                ChckIntersectforTransit(blFront);
                            }
                            cons[0] = conn[2];
                            cons[1] = idx;
                            cons[2] = idx1;
                            tridx = itri[ntri++] = AddTriElem(3, cons);
                            m_ocTree->insertPreProcess(tridx);
                            if (m_ocTree->chckIntersectPreProcess(tridx)) {
                                ChckIntersectforTransit(blFront);
                            }*/

                            int itrii;

                            itrii = blLowerFrt->GetSTriIdx(inei, 0);
                            m_ocTree->rmDataPreProcess(itrii);
                            itrii = blLowerFrt->GetSTriIdx(inei, 1);
                            m_ocTree->rmDataPreProcess(itrii);

#endif
                        }
                    }
                }
            }
        }
    }
    // if (bcrepramid != blFront->IsPyramid())
    //{
    //	cout << bcrepramid << " " << blFront->IsPyramid();
    //	ChckIntersectforTransit(blFront);
    //	CreatePyramid(blFront);
    //	exit(-98);
    // }
}

void BLMesh::UpdateSymmetry()
{

    m_blFrontList->RestoreFront();
    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        BLNode *nei[3];
        int num_nei;
        blFront->GetNodes(&num_nei, nei);

        for (int i = 0; i < num_nei; i++) {
            int idx1 = nei[i]->GetNodIdx();
            int idx2 = nei[(i + 1) % 3]->GetNodIdx();
            int idx3 = nei[(i + 2) % 3]->GetNodIdx();
            //			const auto& isymfc1 = m_pNodes[idx1].isymfc;
            // const auto& isymfc2 = m_pNodes[idx2].isymfc;

            // 获取 isymfc 的两个节点的 set
            auto &isymfc1 = m_pNodes[idx1].isymfc;
            auto &isymfc2 = m_pNodes[idx2].isymfc;

            std::sort(isymfc1.begin(), isymfc1.end());
            std::sort(isymfc2.begin(), isymfc2.end());

            // 使用 set_intersection 检查是否有公共单元
            std::set<int> intersection;
            std::set_intersection(isymfc1.begin(), isymfc1.end(), isymfc2.begin(), isymfc2.end(),
                                  std::inserter(intersection, intersection.begin()));

            if (m_pNodes[idx1].bsysm && m_pNodes[idx2].bsysm && !intersection.empty()) {
                int id[4] = {idx2, idx1, nei[i]->GetLowerNode()->GetNodIdx(), nei[(i + 1) % 3]->GetLowerNode()->GetNodIdx()};

                // 合并相交面
                auto face = intersect_two(intersect_two(isymfc1, isymfc2), m_pNodes[idx3].isymfc);

                // if (face.empty()) {
                //	for (int k = 0; k < 3; k++) {
                //		auto addNodeFace = [&](int nodeIdx) {
                //			for (auto l : m_pNodes[nodeIdx].isymfc) {
                //				face.push_back(l);
                //			}
                //			};

                //		addNodeFace(nei[(i + k) % 3]->GetNodIdx());
                //		addNodeFace(nei[(i + k) % 3]->GetLowerNode()->GetNodIdx());
                //	}
                //}

                AddElem(4, id, BLEntityTopology::QUADRILATERAL, *intersection.begin());
            }
        }
    }
}
void BLMesh::CreateTransitionElements()

{
    BLFront *blFront;
    bool bInter = false;

#ifdef GEN_PYRAMID

    // take out blFront
    m_blFrontList->RestoreFront();
    queue<BLFront *> createPyramid_queue;
    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        createPyramid_queue.push(blFront);
    }

    scheck.clear();
    while (!createPyramid_queue.empty()) {

        auto front = createPyramid_queue.front();
        if (ChckIntersectforTransit(front)) {
            for (auto k : scheck) {
                if (k->GetUpperFront()) {
                    createPyramid_queue.push(k->GetUpperFront());
                }
            }
            scheck.clear();
        }
        createPyramid_queue.pop();
    }

    std::map<std::vector<int>, BLFront *> pair_map;

    while (true) {
        auto front_list = m_blFrontListAll[m_nCurrLayer];
        bool is_break = true;
        front_list->RestoreFront();
        while (front_list->HasNextFront()) {
            BLFront *blFront = front_list->GetNextFront();
            int nnode;
            BLNode *nodes[3];
            blFront->GetNodes(&nnode, nodes);
            std::vector<int> count;
            for (int i = 0; i < nnode; i++) {
                if (nodes[i]->getPerNode()) {
                    count.push_back(nodes[i]->GetNodIdx());
                }
            }
            if (count.size() >= 2) {
                sort(count.begin(), count.end());
                pair_map[count] = blFront;
            }
        }

        front_list->RestoreFront();
        while (front_list->HasNextFront()) {
            BLFront *blFront = front_list->GetNextFront();

            int nnode;
            BLNode *nodes[3];
            blFront->GetNodes(&nnode, nodes);
            std::set<BLFront *> common;
            BLFront *per_front = nullptr;
            std::vector<int> percount;
            for (int i = 0; i < nnode; i++) {
                if (nodes[i]->getPerNode()) {
                    percount.push_back(nodes[i]->getPerNode()->GetNodIdx());
                }
            }

            if (percount.size() >= 2) {
                sort(percount.begin(), percount.end());
                if (pair_map.find(percount) != pair_map.end()) {
                    per_front = pair_map[percount];
                } else {
                    is_break = false;
                    pair_map[percount] = nullptr;
                    for (int i = 0; i < nnode; i++) {
                        RmvUperNeigFrontsAndFreeNode(nodes[i]);
                        StopPropagateNode(nodes[i]);
                    }
                }
            }

            if (per_front) {

                if (blFront->GetUpperFront() && !IsElmDelete(blFront->GetElmIdx())) {
                    if (!per_front->GetUpperFront() || IsElmDelete(per_front->GetElmIdx())) {
                        is_break = false;
                        for (int i = 0; i < nnode; i++) {
                            RmvUperNeigFrontsAndFreeNode(nodes[i]);
                            StopPropagateNode(nodes[i]);
                        }
                    }
                }

                if (per_front->GetUpperFront() && !IsElmDelete(per_front->GetElmIdx())) {
                    if (!blFront->GetUpperFront() || IsElmDelete(blFront->GetElmIdx())) {
                        is_break = false;
                        int pernnode;
                        BLNode *pernodes[3];
                        per_front->GetNodes(&pernnode, pernodes);
                        for (int i = 0; i < pernnode; i++) {
                            RmvUperNeigFrontsAndFreeNode(pernodes[i]);
                            StopPropagateNode(pernodes[i]);
                        }
                    }
                }
            }
        }
        if (is_break) {
            break;
        }
    }

    // do
    //{
    //	cout << "123#";
    //	bInter = false;
    //	m_blFrontList->RestoreFront();
    //	for (auto i : rm_flag)
    //	{
    //		m_ocTree->rmDataPreProcess(i);
    //	}
    //	rm_flag.clear();
    //	while (m_blFrontList->HasNextFront())
    //	{
    //		BLFront *blFront = m_blFrontList->GetNextFront();

    //		//if(!blFront->IsPyramid())
    //		{
    //			//if (blFront->GetTriIdx() == 1043686 && blFront->GetLayerNum() == 9)
    //			//	cout << "##########" << ChckIntersectforTransit(blFront) << endl;
    //			bInter |=
    //		}
    //	}
    //	if (bInter)
    //	{
    //	}
    //} while (bInter);
#endif

    // generate pyramid
    m_blFrontList->RestoreFront();
    while (m_blFrontList->HasNextFront()) {
        BLFront *blFront = m_blFrontList->GetNextFront();
        CreatePyramid(blFront);
    }
    int start = 0;
    for (int i = 0; i < m_vBdyFront.size(); i++) {
        m_vBdyFront[start] = m_vBdyFront[i];
        if (m_vBdyFront[i]->is_boundary_) {
            start++;
        } else {
            front_to_delete_.push_back(m_vBdyFront[i]);
        }
    }
    m_vBdyFront.resize(start);

    m_blFrontList->RestoreFront();
}

bool BLMesh::CheckIsotroStop(BLNode *blNod)
{
    // TODO：需要完善
    BLVector vecC, vecN, vecD;
    BLNode *blDescentNod;
    BLFront *blNeigFrts[MAX_NCONN * 2];
    BLNode *blNeigNodes[MAX_NCONN * 2];
    int nblNeigFrts;
    int nblNodes;
    double flag, cosa, angle;

    if (m_bIsotropicStop > 1e-5) {

        auto Neig_vec = blNod->GetNeigNods();
        bool iso_stop = true;
        for (auto &i : Neig_vec) {
            if (!i->iso_stop) {
                iso_stop = false;
            }
        }
        return iso_stop;
    }
    return false;
}
bool BLMesh::CheckStop(BLNode *blNod, BLNode *blNodNew, int iLayer)
{
    BLVector vecC, vecN, vecD;
    BLNode *blDescentNod;
    BLFront *blNeigFrts[MAX_NCONN * 2];
    BLNode *blNeigNodes[MAX_NCONN * 2];
    int nblNeigFrts;
    int nblNodes;
    double flag, cosa, angle;

    vecC = blNod->GetNormal();
    vecN = blNodNew->GetNormal();

    flag = vecC * vecN;

    // stop propagating a node when its normal inverts
    if (flag < 0.0) {
        return true;
    }

    // stop propagating a node when its normal deviate its descentnod normal a lot
    // blDescentNod = blNod->GetDescentNod();
    // if (blDescentNod != nullptr)
    //{
    //	// 		if(blDescentNod->GetUvalue() < 0.03 && iLayer > 5)
    //	// 			return true;

    //	vecD = blDescentNod->GetNormal();

    //	angle = 90/*45*/;
    //	cosa = vecN * vecD / (vecN.magnitude()*vecD.magnitude());
    //	if (cosa < cos(rad(angle))) {
    //		//cout << "3";
    //		return true;
    //	}
    //}

#if 1

#endif

    return false;
}
bool BLMesh::CheckPrismEveryVolumn(int nconn, int *conn)
{
#ifdef CHECK_VOLUMN
    const double threshold = 1e-15 * cf.step_len;
    BLVector gprism[6];
    for (int i = 0; i < 6; i++) {
        gprism[i].x = m_pNodes[conn[i]].coord[0];
        gprism[i].y = m_pNodes[conn[i]].coord[1];
        gprism[i].z = m_pNodes[conn[i]].coord[2];
    }

    // 所有四面体体积
    double vols[] = {
        (gprism[3] - gprism[0]) * ((gprism[1] - gprism[0]) ^ (gprism[2] - gprism[0])) / 6, // volumn0123
        (gprism[4] - gprism[0]) * ((gprism[1] - gprism[0]) ^ (gprism[2] - gprism[0])) / 6, // volumn0124
        (gprism[5] - gprism[0]) * ((gprism[1] - gprism[0]) ^ (gprism[2] - gprism[0])) / 6, // volumn0125

        (gprism[0] - gprism[3]) * ((gprism[5] - gprism[3]) ^ (gprism[4] - gprism[3])) / 6, // volumn3450
        (gprism[1] - gprism[3]) * ((gprism[5] - gprism[3]) ^ (gprism[4] - gprism[3])) / 6, // volumn3451
        (gprism[2] - gprism[3]) * ((gprism[5] - gprism[3]) ^ (gprism[4] - gprism[3])) / 6, // volumn3452

        (gprism[5] - gprism[0]) * ((gprism[3] - gprism[0]) ^ (gprism[1] - gprism[0])) / 6, // volumn0135
        (gprism[5] - gprism[0]) * ((gprism[4] - gprism[0]) ^ (gprism[1] - gprism[0])) / 6, // volumn0145

        (gprism[3] - gprism[1]) * ((gprism[4] - gprism[1]) ^ (gprism[5] - gprism[1])) / 6, // volumn1243
        (gprism[3] - gprism[1]) * ((gprism[5] - gprism[1]) ^ (gprism[2] - gprism[1])) / 6, // volumn1253

        (gprism[4] - gprism[2]) * ((gprism[5] - gprism[2]) ^ (gprism[0] - gprism[2])) / 6, // volumn0254
        (gprism[4] - gprism[2]) * ((gprism[3] - gprism[2]) ^ (gprism[0] - gprism[2])) / 6  // volumn0234
    };

    for (int i = 0; i < sizeof(vols) / sizeof(double); ++i) {
        if (fabs(vols[i]) < threshold) {

            std::cout << "Negative or tiny prism sub-volume detected: " << vols[i] << std::endl;
            return false;
        }
    }
    return true;
#else
    return true;
#endif
}

bool BLMesh::CheckPrismVolumn(int nconn, int *conn)
{
#ifdef CHECK_VOLUMN
    BLVector gprism[6];
    for (int i = 0; i < 6; i++) {
        gprism[i].x = m_pNodes[conn[i]].coord[0];
        gprism[i].y = m_pNodes[conn[i]].coord[1];
        gprism[i].z = m_pNodes[conn[i]].coord[2];
    }
    double volumn1 = (gprism[3] - gprism[0]) * ((gprism[1] - gprism[0]) ^ (gprism[2] - gprism[0])) / 6 +
                     (gprism[3] - gprism[1]) * ((gprism[4] - gprism[1]) ^ (gprism[5] - gprism[1])) / 6 +
                     (gprism[3] - gprism[1]) * ((gprism[5] - gprism[1]) ^ (gprism[2] - gprism[1])) / 6;

    double volumn2 = (gprism[4] - gprism[1]) * ((gprism[2] - gprism[1]) ^ (gprism[0] - gprism[1])) / 6 +
                     (gprism[4] - gprism[2]) * ((gprism[5] - gprism[2]) ^ (gprism[3] - gprism[2])) / 6 +
                     (gprism[4] - gprism[2]) * ((gprism[3] - gprism[2]) ^ (gprism[0] - gprism[2])) / 6;
    // cout << volumn << endl;
    if (min(volumn1, volumn2) < 1e-15 * cf.step_len) {
        cout << "negetive prism volumn" << endl;
    }
    return min(volumn1, volumn2) > 1e-15 * cf.step_len;
#else
    return true;
#endif
}
bool BLMesh::CheckPrismValidity(int nconn, int *conn, int *pidx)
{
    int i, idx;
    double fprismqual = FLT_MAX, ftmp, Jackbin, mjkb1, mjkb2, mjkb3;
    BLVector gprism[6];
    constexpr int eps_eta_kxi[6][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1}
    }; // 六个顶点
    double para = 2.0 * sqrtf(3.0) /*pow( 3,0.5 )*/, jack[6], eps = 1e-6 /*1e-6*/;

    for (i = 0; i < nconn; i++) {
        gprism[i].x = m_pNodes[conn[i]].coord[0];
        gprism[i].y = m_pNodes[conn[i]].coord[1];
        gprism[i].z = m_pNodes[conn[i]].coord[2];
    }

    BLVector t2(gprism[0] - gprism[2]);
    BLVector t3(gprism[1] - gprism[0]);

    BLVector s1(gprism[3] - gprism[0]);
    BLVector s2(gprism[4] - gprism[1] - s1);
    BLVector s3(gprism[5] - gprism[2] - s1);

    for (int i = 0; i < 6; ++i) {
        BLVector jackb1(t3 + s2 * eps_eta_kxi[i][2]);
        BLVector jackb2(t2 * (-1) + s3 * eps_eta_kxi[i][2]);
        // 		jackb1 = t2 + ( s2 - s1 ) * eps_eta_kxi[i][2];
        // 		jackb2 = t3 + ( s3 - s1 ) * eps_eta_kxi[i][2];
        BLVector jackb3(s1 + s2 * eps_eta_kxi[i][0] + s3 * eps_eta_kxi[i][1]);

        Jackbin = jackb1.x * jackb2.y * jackb3.z + jackb2.x * jackb3.y * jackb1.z + jackb3.x * jackb1.y * jackb2.z -
                  jackb2.x * jackb1.y * jackb3.z - jackb1.x * jackb3.y * jackb2.z - jackb3.x * jackb2.y * jackb1.z;

        Jackbin *= para;

        mjkb1 = jackb1.magnitude2();
        mjkb2 = jackb2.magnitude2();
        mjkb3 = jackb3.magnitude();
        double mjkb12 = (jackb1 - jackb2).magnitude2();

        /*
         */
        if (!(mjkb3 >= 0. && (mjkb1 + mjkb2 + mjkb12) >= 0.)) {
            spdlog::info("prism error!  \n");

            spdlog::info("%lf %lf %lf %lf\n", mjkb1, mjkb2, mjkb3, mjkb12);
            return false;
            // exit(0);
        }

        ftmp = Jackbin / (mjkb3 * (mjkb1 + mjkb2 + mjkb12));
        jack[i] = ftmp;
        // if( ftmp > 1 || ftmp < -1 )
        if (ftmp < -1) {
            spdlog::info("ERROR:prism_qulity(VPRISM gprism)\n");
            return false;
        }
        if (ftmp < fprismqual) {
            fprismqual = ftmp;
            idx = i;
        }
    }

    idx = idx > 2 ? (idx - 3) : idx;
    if (idx == 1) {
        idx = 2;
    } else if (idx == 2) {
        idx = 1;
    }

    *pidx = idx;

    if (fprismqual < eps) {
        spdlog::info("qual: {} {}\n", fprismqual, idx);
        spdlog::info("triangle: ({} {} {})\n", conn[0], conn[1], conn[2]);
        return false;
    }

    return true;
}

bool BLMesh::CheckPrismSkewness(int nconn, int *conn)
{
    using Eigen::Vector3d;
    static constexpr const double pi = 3.14159265358979323846;
    Vector3d a(m_pNodes[conn[0]].coord[0], m_pNodes[conn[0]].coord[1], m_pNodes[conn[0]].coord[2]);
    Vector3d b(m_pNodes[conn[1]].coord[0], m_pNodes[conn[1]].coord[1], m_pNodes[conn[1]].coord[2]);
    Vector3d c(m_pNodes[conn[2]].coord[0], m_pNodes[conn[2]].coord[1], m_pNodes[conn[2]].coord[2]);
    Vector3d d(m_pNodes[conn[3]].coord[0], m_pNodes[conn[3]].coord[1], m_pNodes[conn[3]].coord[2]);
    Vector3d e(m_pNodes[conn[4]].coord[0], m_pNodes[conn[4]].coord[1], m_pNodes[conn[4]].coord[2]);
    Vector3d f(m_pNodes[conn[5]].coord[0], m_pNodes[conn[5]].coord[1], m_pNodes[conn[5]].coord[2]);

    // 向量夹角（弧度）
    auto angleBetween = [](const Vector3d &u, const Vector3d &v) -> double {
        double nu = u.norm();
        double nv = v.norm();
        if (nu <= std::numeric_limits<double>::epsilon() || nv <= std::numeric_limits<double>::epsilon()) {
            return 0.0;
        }
        double c = u.dot(v) / (nu * nv);
        if (c > 1.0) {
            c = 1.0;
        }
        if (c < -1.0) {
            c = -1.0;
        }
        return std::acos(c);
    };

    // 给一组角度 + 理想角（弧度），算等角偏斜
    auto equiangleSkewFromAngles = [](const std::vector<double> &angles, double idealAngle) -> double {
        if (angles.empty()) {
            return 0.0;
        }

        double theta_min = angles[0];
        double theta_max = angles[0];
        for (size_t i = 1; i < angles.size(); ++i) {
            theta_min = std::min(theta_min, angles[i]);
            theta_max = std::max(theta_max, angles[i]);
        }

        double s1 = (theta_max - idealAngle) / (pi - idealAngle);
        double s2 = (idealAngle - theta_min) / idealAngle;

        double skew = std::max(s1, s2);
        skew = std::max(std::min(skew, 1.0), 0.0);
        return skew;
    };

    // 三角形等角偏斜
    auto triangleSkew = [&](const Vector3d &p0, const Vector3d &p1, const Vector3d &p2) -> double {
        std::vector<double> angs;
        angs.reserve(3);

        angs.push_back(angleBetween(p1 - p0, p2 - p0)); // at p0
        angs.push_back(angleBetween(p0 - p1, p2 - p1)); // at p1
        angs.push_back(angleBetween(p0 - p2, p1 - p2)); // at p2

        double ideal = pi / 3.0;                        // 60°
        return equiangleSkewFromAngles(angs, ideal);
    };

    // 四边形等角偏斜：按环绕顺序 (p0,p1,p2,p3)
    auto quadSkew = [&](const Vector3d &p0, const Vector3d &p1, const Vector3d &p2, const Vector3d &p3) -> double {
        std::vector<double> angs;
        angs.reserve(4);

        // 顶点 p0：角(p1-p0, p3-p0)
        angs.push_back(angleBetween(p1 - p0, p3 - p0));
        // 顶点 p1：角(p2-p1, p0-p1)
        angs.push_back(angleBetween(p2 - p1, p0 - p1));
        // 顶点 p2：角(p3-p2, p1-p2)
        angs.push_back(angleBetween(p3 - p2, p1 - p2));
        // 顶点 p3：角(p0-p3, p2-p3)
        angs.push_back(angleBetween(p0 - p3, p2 - p3));

        double ideal = pi / 2.0; // 90°
        return equiangleSkewFromAngles(angs, ideal);
    };

    // 底三角: (a,b,c)
    // 顶三角: (d,e,f)
    // 侧四边形: (a,b,e,d), (b,c,f,e), (c,a,d,f)
    double s_tri0 = triangleSkew(a, b, c);
    double s_tri1 = triangleSkew(d, e, f);
    double s_q0 = quadSkew(a, b, e, d);
    double s_q1 = quadSkew(b, c, f, e);
    double s_q2 = quadSkew(c, a, d, f);

    double equal_angle_skewness = std::max({s_tri0, s_tri1, s_q0, s_q1, s_q2});

    // 判定 1 (带动态阈值)
    if (cf.max_skewness[0] < 0.1) {
        throw(std::logic_error("maximum equal skewnwass is too small!"));
    }

    // 原逻辑：根据层数放宽阈值
    // cf.max_equal_skewness = 1.0 - (1.00 - cf.max_equal_skewness) * ((m_nCurrLayer > 5 ? 5 : m_nCurrLayer) * 0.2);

    if (equal_angle_skewness > cf.max_skewness[0]) {
        return false;
    }

    return true;
}

bool BLMesh::CheckIntersection(BLFront *blFront, bool *used_by_neigh_front, int *tridx)
{

    int id = blFront->GetTriIdx();

    // return false;
    int i, j, neigs, ntri, tris[3], istri[12], nstri = 0;
    int itri[13]; /// 新增的面
    int idx1, idx2, ithird[DIM3], nNods, nconn, conn[DIM3 * 2];
    // bool flag[DIM3];
    atomic_bool is_inserect;
    BLNode *blNods[DIM3];

    nconn = 3;
    blFront->GetNodes(&nNods, blNods);

    for (i = 0; i < DIM3; i++) {
        // blNods[i] = (BLNode*)m_pNodes[conn[i]].pointer;
        conn[i] = blNods[i]->GetNodIdx();
        conn[nconn + i] = blNods[i]->GetUpperNode()->GetNodIdx();
        used_by_neigh_front[i] = false;
    }

    auto &neigFrts = (blFront)->m_arrNeigFronts;
    neigs = blFront->m_nNeiFront;

    /*获取底层front周围的三个front的多出来的侧面的三角形ID，存在ithird中*/
    for (i = 0; i < neigs; i++) {
        BLFront *ft = neigFrts[i];
        if (ft->GetUpperFront()) {
            for (j = 0; j < DIM3; j++) {
                if (ft->IncludeNode(blNods[(j + 1) % DIM3]) && ft->IncludeNode(blNods[(j + 2) % DIM3])) {
                    used_by_neigh_front[j] = true; // 标记被使用过
                    ithird[j] = ft->GetThirdNodIdx(blNods[(j + 1) % DIM3], blNods[(j + 2) % DIM3]);
                    istri[nstri++] = ft->GetSTriIdx(ithird[j], 0);
                    istri[nstri++] = ft->GetSTriIdx(ithird[j], 1);
                }
            }
        }
    }

    ntri = 0;
    for (i = 0; i < DIM3; i++) {
        if (!used_by_neigh_front[i]) {
            // add triangles
            //  type 1
            tris[0] = conn[(i + 1) % DIM3];
            tris[1] = conn[(i + 2) % DIM3];
            tris[2] = conn[(i + 1) % DIM3 + 3];
            idx1 = itri[ntri++] = AddTriElem(3, tris);

            tris[0] = conn[(i + 2) % DIM3];
            tris[1] = conn[(i + 1) % DIM3 + 3];
            tris[2] = conn[(i + 2) % DIM3 + 3];
            idx2 = itri[ntri++] = AddTriElem(3, tris);

            blFront->SetSTriIdx(i, 0, idx1);
            blFront->SetSTriIdx(i, 1, idx2);
        }
    }
    tris[0] = conn[3];
    tris[1] = conn[4];
    tris[2] = conn[5];
    itri[ntri++] = *tridx = AddTriElem(3, tris);

    // is_inserect = false;

#ifndef _DEBUG

#endif

    /*三棱柱内部求相交*/

    // for (int k = 0; k < ntri - 1; k++) {
    //	for (int j = k + 1; j < ntri - 1; j++) {
    //		if (k % 2 == 0 && j == k + 1)//这里相邻的话一定是不相交的
    //			continue;
    //		if (!is_inserect&&Octree::tri_overlap_test(itri[k], itri[j], m_ocTree->pOctreeAgent->getNod(),
    //m_ocTree->pOctreeAgent->getElm())) { 			is_inserect = true; 			break;
    //			//throw string("self interact in prism") + std::to_string(k) + std::to_string(j);
    //		}
    //	}
    // }

    return false;
}

void BLMesh::CreateTriangles(BLFront *blFront, int &itrix)
{
    int id = blFront->GetTriIdx();

    // return false;
    int i, j, neigs, ntri, nstri = 0;
    static std::array<int, 3> tris;
    bool used_by_neigh_front[3] = {0}; // 标记被使用过
    int idx1, idx2, ithird[DIM3], nNods, nconn = 3, conn[DIM3 * 2];
    // bool flag[DIM3];
    atomic_bool is_inserect;
    BLNode *blNods[DIM3];
    blFront->GetNodes(&nNods, blNods);

    for (i = 0; i < DIM3; i++) {
        // blNods[i] = (BLNode*)m_pNodes[conn[i]].pointer;
        conn[i] = blNods[i]->GetNodIdx();
        conn[nconn + i] = blNods[i]->GetUpperNode()->GetNodIdx();
    }
    auto &neigFrts = (blFront)->m_arrNeigFronts;
    neigs = blFront->m_nNeiFront;

    /*获取底层front周围的三个front的多出来的侧面的三角形ID，存在ithird中*/
    for (i = 0; i < neigs; i++) {
        BLFront *ft = neigFrts[i];
        if (ft->GetUpperFront()) {
            for (j = 0; j < DIM3; j++) {
                if (!ft->IncludeNode(blNods[j])) {
                    used_by_neigh_front[j] = i; // 标记被使用过
                }
            }
        }
    }
    for (i = 0; i < DIM3; i++) {
        if (!used_by_neigh_front[i]) {
            // add triangles
            //  type 1
            tris[0] = conn[(i + 1) % DIM3];
            tris[1] = conn[(i + 2) % DIM3];
            tris[2] = conn[(i + 1) % DIM3 + 3];
            idx1 = AddTriElem(tris);

            tris[0] = conn[(i + 2) % DIM3];
            tris[1] = conn[(i + 1) % DIM3 + 3];
            tris[2] = conn[(i + 2) % DIM3 + 3];
            idx2 = AddTriElem(tris);

            blFront->SetSTriIdx(i, 0, idx1);
            blFront->SetSTriIdx(i, 1, idx2);
#ifdef _DEBUG
            m_ocTree->num_inter += 2;
#endif
        } else {
            int third;

            auto nei = neigFrts[used_by_neigh_front[i]];
            NeighIdx(nei, blFront, &third);
            blFront->SetSTriIdx(i, 0, nei->GetSTriIdx(third, 0));
            blFront->SetSTriIdx(i, 1, nei->GetSTriIdx(third, 1));
        }
    }
    tris[0] = conn[3];
    tris[1] = conn[4];
    tris[2] = conn[5];
    itrix = AddTriElem(tris);
#ifdef _DEBUG
    m_ocTree->num_inter += 1;
#endif

    return;
}

void BLMesh::FixHightRatio(BLNode *blNod, MBLNode *pNodes)
{
    const auto &node_array = blNod->GetNeigNods();
    if (node_array.empty()) {
        return;
    }
    double length = 0;
    BLVector COORD(pNodes[blNod->GetNodIdx()].coord[0], pNodes[blNod->GetNodIdx()].coord[1], pNodes[blNod->GetNodIdx()].coord[2]);
    for (auto i : node_array) {
        length += ((BLVector(pNodes[i->GetNodIdx()].coord[0], pNodes[i->GetNodIdx()].coord[1], pNodes[i->GetNodIdx()].coord[2])) - COORD)
                      .magnitude();
    }
    length /= node_array.size();

    if (length * 0.8 < cf.step_len) {
        blNod->SetFixedHightRatio(length * 0.3 / cf.step_len - 1);
    }
}
#ifdef OLD
void BLMesh::SmoothHeightRatio(BLNode *blNod, MBLNode *pNodes)
{
    const auto &node_array = blNod->GetNeigNods();
    if (node_array.empty()) {
        return;
    }
    double length = 0;
    BLVector COORD(pNodes[blNod->GetNodIdx()].coord[0], pNodes[blNod->GetNodIdx()].coord[1], pNodes[blNod->GetNodIdx()].coord[2]);
    for (auto i : node_array) {
        length += ((BLVector(pNodes[i->GetNodIdx()].coord[0], pNodes[i->GetNodIdx()].coord[1], pNodes[i->GetNodIdx()].coord[2]) +
                    i->GetHeight()) -
                   COORD) *
                  blNod->GetNormal();
    }
    double suppose_height = blNod->GetHeightLength() / (1 + blNod->GetHightRatio()) / (1 - blNod->GetFixedHightRatio());
    double min_dos = 1.0;

    for (auto i : node_array) {
        min_dos = min(min_dos, blNod->GetNormal() * i->GetNormal());
    }
    min_dos = 1 - min_dos;
    length /= node_array.size();
    length = (length - suppose_height) / suppose_height;

    double alength = 0;
    // if(length>0)
    //	alength = atan(length * (1 + min_dos / 2)) * 2 / PI;
    // else {
    //	alength = atan(length * (1 + min_dos / 2)) * (0.5) * 2 / PI;
    // }

    alength = 1.0 / (1.0 + std::exp(-0.01 * length)) - 0.5;

    blNod->SetHightRatio(alength);

    // smooth fixed height ratio

    queue<BLNode *> q;
    q.push(blNod);
    while (!q.empty()) {

        BLNode *n = q.front();
        q.pop();
        auto nei = n->GetNeigNods();
        for (auto i : nei) {
            if (i->GetFixedHightRatio() > n->GetFixedHightRatio() + 0.1) {
                i->SetFixedHightRatio(n->GetFixedHightRatio() + 0.1);
                q.push(i);
            }
        }
    }
}
#else
void BLMesh::SmoothHeightRatio(BLNode *blNod, MBLNode *pNodes)
{
    const auto &node_array = blNod->GetNeigNods();
    if (node_array.empty()) {
        return;
    }
    double length = 0;
    BLVector COORD(pNodes[blNod->GetNodIdx()].coord[0], pNodes[blNod->GetNodIdx()].coord[1], pNodes[blNod->GetNodIdx()].coord[2]);
    for (auto i : node_array) {
        length += ((BLVector(pNodes[i->GetNodIdx()].coord[0], pNodes[i->GetNodIdx()].coord[1], pNodes[i->GetNodIdx()].coord[2]) +
                    i->GetHeight()) -
                   COORD) *
                  blNod->GetNormal();
    }
    double suppose_height = blNod->GetHeightLength() / (1 + blNod->GetHightRatio()) / (1 + blNod->GetFixedHightRatio());
    double min_dos = 1.0;

    for (auto i : node_array) {
        min_dos = min(min_dos, blNod->GetNormal() * i->GetNormal());
    }
    min_dos = 1 - min_dos;
    length /= node_array.size();
    length = (length - suppose_height) / suppose_height;

    double alength = 0;
    // if(length>0)
    //	alength = atan(length * (1 + min_dos / 2)) * 2 / PI;
    // else {
    //	alength = atan(length * (1 + min_dos / 2)) * (0.5) * 2 / PI;
    // }

    alength = 1.0 / (1.0 + std::exp(-0.5 * length)) - 0.5;

    blNod->SetHightRatio(alength);

    // smooth fixed height ratio

    queue<BLNode *> q;
    q.push(blNod);
    while (!q.empty()) {

        BLNode *n = q.front();
        q.pop();
        auto nei = n->GetNeigNods();
        for (auto i : nei) {
            if (i->GetFixedHightRatio() > n->GetFixedHightRatio() + 0.1) {
                i->SetFixedHightRatio(n->GetFixedHightRatio() + 0.1);
                q.push(i);
            }
        }
    }
}
#endif
int BLMesh::GetOuterBoundary(int *npt, int *nlem, double **pt, int **elm, int **l_to_g, bool add_symm)
{
    int i, j, *nodmap = nullptr, idx;
    double *bpnt = nullptr;
    int nbpt, nbelm, *belm = nullptr;
    int nbdry, *pbdry = nullptr;
    int nbdryelm = 0, *pbdryelm = nullptr;

    std::set<int> setbdry;
    std::set<int>::iterator iter;

    // get elements on symmetry

    int *nbdryi = nullptr, **pbdryi = nullptr, nTtlInitSymBdrys;
    int **pbdryelmi = nullptr;

    nbdryi = new int[m_nSymLoop];
    pbdryi = new int *[m_nSymLoop];
    std::vector<int> nbdryelmi(m_nSymLoop, 0);
    pbdryelmi = new int *[m_nSymLoop];
    for (i = 0; i < m_nSymLoop; i++) {
        pbdryi[i] = nullptr;
        pbdryelmi[i] = nullptr;
    }

    if (m_nSymBdrys > 0) {
#ifdef _NEW_SYMM
        m_nTtlInitSymBdrys = m_nSymBdrys;
        std::set<int> perid;
        for (auto id : cf.perfc) {
            perid.insert(id);
        }
        for (i = 0; i < m_nSymLoop; i++) {
            if (!perid.size()) {
                if (perid.find(m_symFidx[i]) == perid.end()) {
                    CalSymplnBdry(&nbdryi[i], &(pbdryi[i]), m_symFidx[i], add_symm);
#ifdef _DEBUG
                    OutputSymplnBdry("testbdry.vtk", nbdryi[i], pbdryi[i], m_symFidx[i]);
#endif
                    if (nbdryi[i] && add_symm) {
                        // OutputSymplnBdry("testbdry.vtk", nbdryi[i], pbdryi[i], m_symFidx[i]);
                        CalSymplnMsh(nbdryi[i], pbdryi[i], &nbdryelmi[i], &pbdryelmi[i], m_symVals[i], m_symFidx[i]);

                        nbdryelm += nbdryelmi[i];
                    }
                }
            }
        }

        RemvDeletedSymBdry();

        pbdryelm = new int[nbdryelm * 3];
        int cnt = 0;
        for (i = 0; i < m_nSymLoop; i++) {
            for (j = 0; j < nbdryelmi[i]; j++) {
                pbdryelm[cnt * 3 + 0] = pbdryelmi[i][j * 3 + 0];
                pbdryelm[cnt * 3 + 1] = pbdryelmi[i][j * 3 + 1];
                pbdryelm[cnt * 3 + 2] = pbdryelmi[i][j * 3 + 2];

                cnt++;
            }
        }
#else
        CalSymplnBdry(&nbdry, &pbdry);
        OutputSymplnBdry("testbdry.fr2", nbdry, pbdry);
        CalSymplnMsh(nbdry, pbdry, &nbdryelm, &pbdryelm);
#endif

        for (i = 0; i < nbdryelm; i++) {

            setbdry.insert(pbdryelm[i * 3 + 0]);
            setbdry.insert(pbdryelm[i * 3 + 1]);
            setbdry.insert(pbdryelm[i * 3 + 2]);
        }
    }

    // get exposed elements
    for (i = 0; i < noutbdry; i++) {
        setbdry.insert(outbdry[3 * i + 0]);
        setbdry.insert(outbdry[3 * i + 1]);
        setbdry.insert(outbdry[3 * i + 2]);
    }

    // set local to global mapping and global to local mapping
    *l_to_g = new int[setbdry.size()];

    nodmap = new int[m_nNodes];
    for (i = 0; i < m_nNodes; i++) {
        nodmap[i] = -1;
    }

    idx = 0;
    iter = setbdry.begin();
    while (iter != setbdry.end()) {
        nodmap[*iter] = idx;
        (*l_to_g)[idx] = *iter;

        ++idx;
        ++iter;
    }

    //
    *npt = nbpt = setbdry.size();
    *nlem = nbelm = noutbdry + nbdryelm;

    bpnt = new double[3 * nbpt];
    belm = new int[3 * nbelm];

    *pt = bpnt;
    *elm = belm;

    i = 0;
    iter = setbdry.begin();
    while (iter != setbdry.end()) {
        int ipt0 = *iter;

        bpnt[i * 3 + 0] = m_pNodes[ipt0].coord[0];
        bpnt[i * 3 + 1] = m_pNodes[ipt0].coord[1];
        bpnt[i * 3 + 2] = m_pNodes[ipt0].coord[2];

        ++iter;
        ++i;
    }

    // set boundary elements for outer unstructured mesh generation
    if (m_nSymBdrys > 0) {
        int ipt0, ipt1, ipt2;
        for (i = 0; i < nbdryelm; i++) {
            ipt0 = pbdryelm[i * 3 + 0];
            ipt1 = pbdryelm[i * 3 + 1];
            ipt2 = pbdryelm[i * 3 + 2];

            belm[i * 3 + 0] = nodmap[ipt0];
            belm[i * 3 + 1] = nodmap[ipt1];
            belm[i * 3 + 2] = nodmap[ipt2];
        }
    }

    for (i = nbdryelm, j = 0; i < nbdryelm + noutbdry; i++, j++) {
        int ipt0, ipt1, ipt2;

        ipt0 = outbdry[3 * j + 0];
        ipt1 = outbdry[3 * j + 1];
        ipt2 = outbdry[3 * j + 2];

        belm[i * 3 + 0] = nodmap[ipt0];
        belm[i * 3 + 1] = nodmap[ipt1];
        belm[i * 3 + 2] = nodmap[ipt2];
    }

    // free
    if (nodmap) {
        delete[] nodmap;
    }
    if (pbdry) {
        delete[] pbdry;
    }
    if (pbdryelm) {
        delete[] pbdryelm;
    }

    for (i = 0; i < m_nSymLoop; i++) {
        if (pbdryi[i]) {
            delete[] (pbdryi[i]);
        }
        if (pbdryelmi[i]) {
            delete[] (pbdryelmi[i]);
        }
    }

    if (pbdryi) {
        delete[] pbdryi;
    }
    if (pbdryelmi) {
        delete[] pbdryelmi;
    }
    if (nbdryi) {
        delete[] nbdryi;
    }

    return 0;
}

void BLMesh::OutputFr2(string file, int npt, int nelm, double *pt, int *elm)
{
    int i;
    FILE *fout = nullptr;
    const char *filename = file.data();
    fout = fopen(filename, "w");
    if (fout == nullptr) {
        throw(std::string("no authority to write fr2 file!"));
    }
    fprintf(fout, "# vtk DataFile Version 2.0\n");
    fprintf(fout, "test0.vtk\n");
    fprintf(fout, "ASCII\n");
    fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fout, "POINTS %d double\n", npt);

    // fprintf(fout, "%d %d 0 0 0 0 0\n", npt, nelm);
    for (i = 0; i < npt; i++) {
        fprintf(fout, "%lf %lf 0.0\n", pt[2 * i + 0], pt[2 * i + 1]);
    }
    fprintf(fout, "CELLS %d %d\n", nelm, 3 * nelm);
    for (i = 0; i < nelm; i++) {
        fprintf(fout, "2 %d %d\n", elm[2 * i + 0], elm[2 * i + 1]);
    }
    fprintf(fout, "CELL_TYPES %d\n", nelm);

    for (i = 0; i < nelm; i++) {
        fprintf(fout, "3\n");
    }

    fclose(fout);
    fout = nullptr;
}

void BLMesh::OutputPls(string file, int npt, int nelm, double *pt, int *elm)
{
    int i;
    const char *filename = file.data();
    FILE *fout = nullptr;
    fout = fopen(filename, "w");
    if (fout == nullptr) {
        throw(std::string("Error: no authority to wirte pls"));
    }

    fprintf(fout, "%d %d 0 0 0 0\n", nelm, npt);
    for (i = 0; i < npt; i++) {
        fprintf(fout, "%d %.10lf %.10lf %.10lf\n", i + 1, pt[3 * i + 0], pt[3 * i + 1], pt[3 * i + 2]);
    }

    for (i = 0; i < nelm; i++) {
        fprintf(fout, "%d %d %d %d 1\n", i + 1, elm[3 * i + 0] + 1, elm[3 * i + 1] + 1, elm[3 * i + 2] + 1);
    }
    fclose(fout);
    fout = nullptr;
}
int BLMesh::OutputOuterBoundary(string filename, int npt, int nelm, double *pt, int *elm)
{
    std::string fs(filename);
    BLMesh::OutputPls(fs, npt, nelm, pt, elm);

    return 0;
}

void BLMesh::CheckElemNum(int itmck)
{
    if (itmck == 0 && m_nElems > 10000) {
        throw(std::string("too much element"));
    }
}

void BLMesh::PstprecsMergedElm()
{
    int i, j, sidx, npt, nelm, ipt1, ipt2;
    int nconn, nsm, cnt;

    npt = m_nNodes;
    nelm = m_nElems;
    sidx = m_nSurfElems;
    // handle merged points from collapse operation
    for (i = sidx; i < nelm; i++) {
        nconn = m_pElems[i].nconn;
        nsm = 0;
        for (j = 0; j < nconn; j++) {
            ipt1 = m_pElems[i].conn[j];
            if (m_pNodes[ipt1].reserved > 0) {
                m_pElems[i].conn[j] = m_pNodes[ipt1].reserved;
                nsm++;
            }
        }

        if (nsm > 0) {
            int kk = 0;
            bool bsm = false;
            for (j = 0; j < nconn; j++) {
                bsm = false;
                ipt1 = m_pElems[i].conn[j];
                for (int k = 0; k < j; k++) {
                    ipt2 = m_pElems[i].conn[k];
                    if (ipt1 == ipt2) {
                        bsm = true;
                        break;
                    }
                }
                if (!bsm) {
                    m_pElems[i].conn[kk] = m_pElems[i].conn[j];
                    kk++;
                }
            }
            m_pElems[i].nconn = kk;
        }
    }
    int dim = DIM3;

    for (i = 0; i < noutbdry; i++) {
        for (j = 0; j < dim; j++) {
            if (m_pNodes[outbdry[2 * i + j]].reserved > 0) {
                outbdry[2 * i + 0] = m_pNodes[outbdry[2 * i + j]].reserved;
            }
        }
    }
}

void BLMesh::FreeMemory()
{

    if (outbdry) {
        noutbdry = 0;
        delete[] outbdry;
        outbdry = nullptr;
    }

    m_vBdyFront.clear();

    if (m_pPntIdx) {
        delete[] m_pPntIdx;
        m_pPntIdx = nullptr;
    }

    if (m_pPntElm) {
        delete[] m_pPntElm;
        m_pPntElm = nullptr;
    }
}

std::vector<int> BLMesh::RemvNodElm()
{
    int i, j, sidx, npt, nelm, ipt1, ipt2;
    int nconn, nsm, cnt, nd;

    npt = m_nNodes;
    nelm = m_nElems;
    sidx = m_nSurfElems;

    std::vector<int> o_to_n(npt, -1);
    // remove deleted element
    cnt = 0 /*sidx*/;
    nd = 0;
    for (i = 0 /*sidx*/; i < nelm; i++) {
        if (!IsElmDelete(i)) {
            m_pElems[cnt] = m_pElems[i];
            cnt++;
        } else {
            nd++;
        }
    }
    spdlog::info("Num of removed elements: {} ({} {})\n", nd, m_nElems, cnt);
    nelm = cnt;
    m_nElems = cnt;

    // label unreferenced nodes
    //  	for (i=0; i<m_nSurfNodes; i++)
    //  		m_pNodes[i].reserved = 1;

    for (i = 0 /*m_nSurfNodes*/; i < npt; i++) {
        m_pNodes[i].reserved = -1;
    }
#ifndef _DEBUG

#endif
    for (int i = 0 /*sidx*/; i < nelm; i++) {
        for (int j = 0; j < m_pElems[i].nconn; j++) {
            m_pNodes[m_pElems[i].conn[j]].reserved += 1;
        }
    }

    // remove unreferenced nodes
    cnt = 0;
    nd = 0;

    for (i = 0; i < npt; i++) {
        if (m_pNodes[i].reserved >= 0) {
            m_pNodes[cnt].bsysm = m_pNodes[i].bsysm;

            m_pNodes[cnt].coord[0] = m_pNodes[i].coord[0];
            m_pNodes[cnt].coord[1] = m_pNodes[i].coord[1];
            m_pNodes[cnt].coord[2] = m_pNodes[i].coord[2];
            m_pNodes[cnt].pointer = m_pNodes[i].pointer;
            m_pNodes[i].reserved = cnt;
            o_to_n[i] = cnt;
            cnt++;
        } else {
            nd++;
        }
    }
    spdlog::info("Num of removed points: {}\n", nd);
    m_nNodes = cnt;
#ifndef _DEBUG

#endif
    for (int i = 0 /*sidx*/; i < nelm; i++) {
        for (int j = 0; j < m_pElems[i].nconn; j++) {
            m_pElems[i].conn[j] = m_pNodes[m_pElems[i].conn[j]].reserved;
        }
    }

    for (i = 0; i < noutbdry; i++) {
        for (j = 0; j < DIM3; j++) {
            if (m_pNodes[outbdry[DIM3 * i + j]].reserved >= 0) {
                outbdry[DIM3 * i + j] = m_pNodes[outbdry[DIM3 * i + j]].reserved;
            }
            // else
            //	spdlog::info("error in RemvNodElm()\n");
        }
    }
    return o_to_n;
}
/*
 * 该函数建立点和边的map关系，为什么不用std::map<int,std::array<int,2>>呢，简单明了,绕来绕去我吐了
 */
void BLMesh::setpntelm(int npt, int nElm, int *pElm, int **pntidx, int **pntelm)
{
    int i, j, start, *idx, *pelem;
    int npoints, nelm;

    nelm = nElm;
    npoints = npt;
    map<int, vector<int>> mp;
    idx = new int[npoints + 1];
    pelem = new int[nelm * 2];
    *pntidx = idx;
    *pntelm = pelem;

    for (i = 0; i <= npoints; i++) {
        idx[i] = 0;
    }

    /*----------------------------------------------------------------------------
    | Count the number of elements for each point:
    ----------------------------------------------------------------------------*/
    for (i = 0; i < nelm; i++) {
        for (j = 0; j < 2; j++) {
            idx[pElm[2 * i + j]]++;
            mp[pElm[2 * i + j]].push_back(pElm[2 * i + 1 - j]);
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

        if (count > 2) {
            spdlog::info("Node {} is connected more than 2 element!\n", i);
            for (auto j : mp[i]) {
                cout << j << " ";
            }
            RemvNodElm();
            SaveBLMesh("TEMP_BOUNDARY.vtk");
            throw(std::string("non manifold"));
            // return;
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
            int pnt = pElm[2 * i + j];
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

void BLMesh::OrientLoop(int nPt, int nElm, double *pPt, int *element_array, int **elmN, int **ortEl, vector<std::array<int, 2>> &subholes)
{
    int *oriented_element = nullptr, *element_new_array = nullptr, i, j, point_index, ptidxe;
    int begin_num, end_num, edge_index_point, edge_index, k, start_edge, cntl = 0;

    /*标识是否被其他环使用了*/
    vector<bool> isused(nElm, false);
    bool newloop = false;
    // int *point_array = nullptr, *point_element_map = nullptr;
    std::map<int, int> loop_array;
    std::map<int, vector<int>> point_edge_map;
    double looparea[10];

    /*排好序的单元存放与此*/
    oriented_element = new int[nElm];

    element_new_array = new int[2 * nElm];

    for (int i = 0; i < nElm; i++) {
        point_edge_map[element_array[2 * i + 0]].push_back(i);
        point_edge_map[element_array[2 * i + 1]].push_back(i);
    }
    for (auto i : point_edge_map) {
        if (i.second.size() != 2) {
            cout << "point index=" << i.first << " size=" << i.second.size() << endl;
            throw std::logic_error("non manifold point in symmetry plane!");
        }
    }
    *ortEl = oriented_element;
    *elmN = element_new_array;

    // setpntelm(nPt, nElm, element_array, &point_array, &point_element_map);

    /*新的数组的index*/
    int orient_loop_index = 0;

    for (i = 0; i < nElm; i++) {
        if (isused[i]) {
            continue;
        }

        /*新的loop*/
        loop_array[cntl++] = orient_loop_index;

        /*这个loop起始单元是什么*/
        oriented_element[orient_loop_index] = i;

        start_edge = i;

        /*第i个单元的第二个点的id*/
        point_index = element_array[i * 2 + 1];

        /**/
        edge_index_point = i;

        /*将老的elementarray中的元素存入新的*/
        element_new_array[orient_loop_index * 2 + 0] = element_array[i * 2 + 0];
        element_new_array[orient_loop_index * 2 + 1] = element_array[i * 2 + 1];

        newloop = false;
        isused[i] = true;
        orient_loop_index++;
        int next_edge = start_edge;
        /*循环找环，通过点-边 映射关系找*/
        while (true) {

            next_edge = point_edge_map[point_index][0] + point_edge_map[point_index][1] - next_edge;
            if (next_edge == start_edge) {
                break;
            }
            /*找到当前loop的index*/
            edge_index_point = oriented_element[orient_loop_index - 1];
            /*这个边被使用过了哦*/
            isused[next_edge] = true;

            element_new_array[2 * orient_loop_index + 0] = point_index;
            int other_index = element_array[2 * next_edge + 0] + element_array[2 * next_edge + 1] - point_index;
            if (other_index == element_array[2 * next_edge + 0]) {
                swap(element_array[2 * next_edge + 0], element_array[2 * next_edge + 1]);
            }
            element_new_array[2 * orient_loop_index + 1] = other_index;

            /*往后找，链表指向下一个*/
            point_index = other_index;

            /*index++，存下一个单位*/
            orient_loop_index++;
        } // while true
    }

    loop_array[cntl] = orient_loop_index;

    std::vector<int> element_copy(2 * orient_loop_index);
    for (int i = 0; i < element_copy.size(); i++) {
        element_copy[i] = element_new_array[i];
    }

    // reset orientation by geometry
    int before = 0;
    for (int i = 0; i < cntl; i++) {
        int elementstart = before;
        int elementend = loop_array[i];
        before = elementend;

        double min_geometry = std::numeric_limits<double>::max();
        int index;

        for (int j = elementstart; j < elementend; j++) {
            double d = pPt[2 * element_new_array[2 * j + 0] + 0] + pPt[2 * element_new_array[2 * j + 0] + 1];
            if (min_geometry > d) {
                d = min_geometry;
                index = j;
            }
        }

        for (int j = elementstart; j < elementend; j++) {
            int pindex = index + j - elementstart;
            if (pindex >= elementend) {
                pindex -= elementend;
            }
            element_new_array[2 * j + 0] = element_copy[2 * pindex + 0];
            element_new_array[2 * j + 1] = element_copy[2 * pindex + 1];
        }
    }

    spdlog::info("loop number: {}\n", cntl);

    // orient each loop according to its area
    double max_value = 0.0;
    int maxi;
    for (i = 0; i < cntl; i++) {
        double t = 0.0, x1, y1, x2, y2;
        int idx1, idx2;
        begin_num = loop_array[i], end_num = loop_array[i + 1];
        for (j = begin_num; j < end_num; j++) {
            idx1 = element_new_array[2 * j + 0];
            idx2 = element_new_array[2 * j + 1];

            x1 = pPt[2 * idx1 + 0];
            y1 = pPt[2 * idx1 + 1];
            x2 = pPt[2 * idx2 + 0];
            y2 = pPt[2 * idx2 + 1];

            t += x1 * y2 - x2 * y1;
        }

        looparea[i] = t;
        spdlog::info("i{}: %lf\n", i, t);
        if (abs(t) > max_value) {
            max_value = abs(t);
            maxi = i;
        }
    }

    // orient
    for (i = 0; i < cntl; i++) {
        bool borient = false;
        if (i != maxi) {
            subholes.push_back(std::array<int, 2>{loop_array[i], loop_array[i + 1]});
        }
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
            begin_num = loop_array[i], end_num = loop_array[i + 1];
            int s1 = 2 * begin_num, s2 = 2 * end_num - 1;
            while (s2 > s1) {
                swap(element_new_array[s2], element_new_array[s1]);
                s2--;
                s1++;
            }
        }
    }
    // output new array;
}
void BLMesh::RemoveBox(int npt, int &nbdry, double *pt, int *elm) { MeshOrient::removeBox(npt, pt, nbdry, elm); }
void BLMesh::OrientFace(int npt, int nbdry, double *pt, int *elm, int **elmN, int **ortEl, vector<double> &hole)
{
    int *sign = nullptr;
    sign = new int[nbdry];

    MeshOrient::orientTriangularSurface(npt, pt, nbdry, elm, sign, nullptr, hole);

    if (sign) {
        delete[] sign;
        sign = nullptr;
    }
}

void BLMesh::ReadAddBdry(char *filename, bool bfirst)
{
    int i, npt, nelm, tmp, pidx, *ptmap = nullptr;
    FILE *fin = nullptr;
    BLVector pnt;

    fin = fopen(filename, "r");
    if (!fin) {
        spdlog::info("Can not open file %s\n", filename);
    }

    fscanf(fin, "%d %d 0 0 0 0 \n", &nelm, &npt);

    if (bfirst) {
        m_pAddBdry = new int[nelm * 30];
    }

    ptmap = new int[npt];

    for (i = 0; i < npt; i++) {
        double coord[3];
        fscanf(fin, "%d %lf %lf %lf\n", &tmp, &coord[0], &coord[1], &coord[2]);
        pnt.x = coord[0];
        pnt.y = coord[1];
        pnt.z = coord[2];

        pidx = AddNode(pnt, 0.0);
        if (i == 0) {
            spdlog::info("Begin idx: {}\n", pidx);
        }

        ptmap[i] = pidx;
    }

    int cnt = 0;
    for (i = 0; i < nelm; i++) {
        int conn[3], nconn = 3;
        fscanf(fin, "%d %d %d %d %d\n", &tmp, &conn[0], &conn[1], &conn[2], &tmp);

        conn[0] = ptmap[conn[0] - 1];
        conn[1] = ptmap[conn[1] - 1];
        conn[2] = ptmap[conn[2] - 1];

        AddElem(nconn, conn, BLEntityTopology::TRIANGLE);

        if (bfirst) {
            outbdry[3 * noutbdry + 0] = conn[0];
            outbdry[3 * noutbdry + 1] = conn[1];
            outbdry[3 * noutbdry + 2] = conn[2];
            noutbdry++;
        }

        m_pAddBdry[3 * m_nAddBdry + 0] = conn[0];
        m_pAddBdry[3 * m_nAddBdry + 1] = conn[1];
        m_pAddBdry[3 * m_nAddBdry + 2] = conn[2];
        m_nAddBdry++;
    }
}

int BLMesh::GetAddBoundary(int *npt, int *nlem, double **pt, int **elm, int **l_to_g)
{
    int i, j, *nodmap = nullptr, idx;
    double *bpnt = nullptr;
    int nbpt, nbelm, *belm = nullptr;
    int nbdry, *pbdry = nullptr;
    int nbdryelm = 0, *pbdryelm = nullptr;

    std::set<int> setbdry;
    std::set<int>::iterator iter;

    // get exposed elements
    for (i = 0; i < m_nAddBdry; i++) {
        setbdry.insert(m_pAddBdry[3 * i + 0]);
        setbdry.insert(m_pAddBdry[3 * i + 1]);
        setbdry.insert(m_pAddBdry[3 * i + 2]);
    }

    // set local to global mapping and global to local mapping
    *l_to_g = new int[setbdry.size()];

    nodmap = new int[m_nNodes];
    for (i = 0; i < m_nNodes; i++) {
        nodmap[i] = -1;
    }

    idx = 0;
    iter = setbdry.begin();
    while (iter != setbdry.end()) {
        nodmap[*iter] = idx;
        (*l_to_g)[idx] = *iter;

        idx++;
        ++iter;
    }

    //
    *npt = nbpt = setbdry.size();
    *nlem = nbelm = m_nAddBdry;

    bpnt = new double[3 * nbpt];
    belm = new int[3 * nbelm];

    *pt = bpnt;
    *elm = belm;

    i = 0;
    iter = setbdry.begin();
    while (iter != setbdry.end()) {
        int ipt0 = *iter;

        bpnt[i * 3 + 0] = m_pNodes[ipt0].coord[0];
        bpnt[i * 3 + 1] = m_pNodes[ipt0].coord[1];
        bpnt[i * 3 + 2] = m_pNodes[ipt0].coord[2];

        ++iter;
        ++i;
    }

    // set boundary elements for outer unstructured mesh generation
    for (i = 0, j = 0; i < m_nAddBdry; i++, j++) {
        int ipt0, ipt1, ipt2;

        ipt0 = m_pAddBdry[3 * j + 0];
        ipt1 = m_pAddBdry[3 * j + 1];
        ipt2 = m_pAddBdry[3 * j + 2];

        belm[i * 3 + 0] = nodmap[ipt0];
        belm[i * 3 + 1] = nodmap[ipt1];
        belm[i * 3 + 2] = nodmap[ipt2];
    }

    // free
    if (nodmap) {
        delete[] nodmap;
    }
    if (pbdry) {
        delete[] pbdry;
    }
    if (pbdryelm) {
        delete[] pbdryelm;
    }

    return 0;
}
void BLMesh::GenTopMesh(VM &v, bool add_symm)
{
    int nbpt = 0, nbelm = 0, *belm = nullptr, *belmN = nullptr, *orientElm = nullptr;
    int npt = 0, nelm = 0, *elm = nullptr, *l_to_g = nullptr;
    int nbdry, *pbdry = nullptr;
    double *bpt = nullptr, *pt = nullptr;

#if 1
    // ReadAddBdry("add.pls", true);
#endif

    double *pt_size;
    GetOuterBoundary(&nbpt, &nbelm, &bpt, &belm, &l_to_g, add_symm);
    CalOriginSize(nbpt, nbelm, bpt, belm, pt_size, l_to_g);
    v.sizing = pt_size;
    v.l2g = l_to_g;
    vector<double> hole;
    OrientFace(nbpt, nbelm, bpt, belm, &belmN, &orientElm, hole);

    /*add top face*/
    v.pnSNO = nbelm;
    v.ppnSFFmO = new int[3 * nbelm];
    v.ppnSFTpO = new int[nbelm];
    v.nSN0 = nbpt;
    v.ppdSNC0 = new double[3 * nbpt];
    for (int i = 0; i < nelm; i++) {
        v.ppnSFTpO[i] = TRIANGLE;
    }
    for (int i = 0; i < 3 * nbpt; i++) {
        v.ppdSNC0[i] = bpt[i];
    }
    // cout << "nbelm=" << nbelm << endl;
    for (int i = 0; i < 3 * nbelm; i++) {
        v.ppnSFFmO[i] = belm[i];
    }
    // try {
    //	v.convert2Manifold();
    //	v.checkManifold();
    // }
    // catch (std::logic_error error) {
    //	spdlog::info(error.what());
    // }
}
int disfind(vector<int> &arr, int index)
{
    if (arr[index] < 0) {
        return index;
    }
    return arr[index] = disfind(arr, arr[index]);
}
void unino(vector<int> &arr, int l1, int l2)
{
    if (disfind(arr, l1) != disfind(arr, l2)) {
        arr[disfind(arr, l2)] += arr[disfind(arr, l1)];
        arr[disfind(arr, l1)] = disfind(arr, l2);
    }
}
void BLMesh::RemoveOverlapNodeAndElement()
{
    dSet = vector<int>(m_nNodes, -1);
    int conn[10];
    int deToTetra = 0;
    int deToPyra = 0;
    int deTriangle = 0;

    for (int i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::TRIANGLE) {

                map<vector<double>, int> points;
                for (int j = 0; j < 3; j++) {
                    int index = m_pElems[i].conn[j];
                    vector<double> point;
                    for (int k = 0; k < 3; k++) {
                        point.push_back(m_pNodes[index].coord[k]);
                    }

                    if (points.find(point) != points.end()) {
                        unino(dSet, points[point], index);
                    } else {
                        points[point] = index;
                    }
                }
            }
        }
    }

    /*remove zero triangles*/
    for (int i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::TRIANGLE) {
                for (int j = 0; j < 3; j++) {
                    m_pElems[i].conn[j] = disfind(dSet, m_pElems[i].conn[j]);
                }
                if (m_pElems[i].conn[0] == m_pElems[i].conn[1] || m_pElems[i].conn[0] == m_pElems[i].conn[2] ||
                    m_pElems[i].conn[1] == m_pElems[i].conn[2]) {
                    SetElmDelete(i);
                    deTriangle++;
                    m_nTri--;
                }
            }
        }
    }

    /*remove zero triangles in outbdy*/
    int nbdry = noutbdry;
    int *bndry = outbdry;
    int i, cnt, i1, i2, i3;
    vector<int> flag(nbdry);
    for (i = 0; i < nbdry; i++) {
        flag[i] = 1;
    }

    // flag deleted element
    cnt = 0;
    for (i = 0; i < nbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        if (disfind(dSet, i1) == disfind(dSet, i3) || disfind(dSet, i1) == disfind(dSet, i2) || disfind(dSet, i2) == disfind(dSet, i3)) {
        } else {
            bndry[cnt * 3 + 0] = disfind(dSet, bndry[i * 3 + 0]);
            bndry[cnt * 3 + 1] = disfind(dSet, bndry[i * 3 + 1]);
            bndry[cnt * 3 + 2] = disfind(dSet, bndry[i * 3 + 2]);
            cnt++;
        }
    }

    // delete element
    noutbdry = cnt;

    /*remove zero prism*/

    for (int i = 0; i < m_nElems; i++) {
        if (!IsElmDelete(i)) {
            if (m_pElems[i].topo == BLEntityTopology::PRISM) {
                for (int j = 0; j < 3; j++) {
                    m_pElems[i].conn[j] = disfind(dSet, m_pElems[i].conn[j]);
                }
                if (m_pElems[i].conn[0] == m_pElems[i].conn[1] && m_pElems[i].conn[0] == m_pElems[i].conn[2]) {
                    SetElmDelete(i);
                    deToTetra++;

                    conn[0] = m_pElems[i].conn[3];
                    conn[1] = m_pElems[i].conn[5];
                    conn[2] = m_pElems[i].conn[4];
                    conn[3] = m_pElems[i].conn[0];
                    AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                    m_nPrism--;
                    m_nTet++;
                } else if (m_pElems[i].conn[0] == m_pElems[i].conn[1] || m_pElems[i].conn[0] == m_pElems[i].conn[2] ||
                           m_pElems[i].conn[1] == m_pElems[i].conn[2]) {
                    SetElmDelete(i);
                    BLVector co;
                    for (int j = 0; j < 6; j++) {
                        co += BLVector(m_pNodes[m_pElems[i].conn[j]].coord);
                    }
                    co = co / 6;
                    int iNodNew = AddNode(co, 0.0);
                    if (m_pElems[i].conn[0] == m_pElems[i].conn[1]) {
                        conn[0] = m_pElems[i].conn[1];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = m_pElems[i].conn[2];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[2];
                        conn[1] = m_pElems[i].conn[3];
                        conn[2] = m_pElems[i].conn[5];
                        conn[3] = m_pElems[i].conn[0];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[0];
                        conn[1] = m_pElems[i].conn[3];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[3];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                    } else if (m_pElems[i].conn[0] == m_pElems[i].conn[2]) {
                        conn[0] = m_pElems[i].conn[1];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = m_pElems[i].conn[2];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[0];
                        conn[1] = m_pElems[i].conn[4];
                        conn[2] = m_pElems[i].conn[3];
                        conn[3] = m_pElems[i].conn[1];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[2];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[3];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[3];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                    } else if (m_pElems[i].conn[1] == m_pElems[i].conn[2]) {
                        conn[0] = m_pElems[i].conn[2];
                        conn[1] = m_pElems[i].conn[3];
                        conn[2] = m_pElems[i].conn[5];
                        conn[3] = m_pElems[i].conn[0];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[0];
                        conn[1] = m_pElems[i].conn[4];
                        conn[2] = m_pElems[i].conn[3];
                        conn[3] = m_pElems[i].conn[1];
                        conn[4] = iNodNew;
                        AddElem(5, conn, BLEntityTopology::PYRAMID, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[2];
                        conn[1] = m_pElems[i].conn[4];
                        conn[2] = m_pElems[i].conn[5];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                        conn[0] = m_pElems[i].conn[3];
                        conn[1] = m_pElems[i].conn[5];
                        conn[2] = m_pElems[i].conn[4];
                        conn[3] = iNodNew;
                        AddElem(4, conn, BLEntityTopology::TETRAHEDRON, m_pElems[i].igom);
                    }
                    deToPyra += 2;
                    m_nPrism--;
                    m_nTet += 2;
                    m_nPyramid += 2;
                    deToTetra += 2;
                }
            }
        }
    }
    spdlog::info("number of element degenerate to point {}", deTriangle);
    spdlog::info("number of element degenerate to tetra {}", deToTetra);
    spdlog::info("number of element degenerate to Pyramid {}", deToPyra);
}
void BLMesh::FreeMemoryInFrontAndNode()
{
    destroyFront();
    m_TriElm.FreeMemory();
#ifdef USE_MEMORY_POOL
    BLFront::freeMemoryInPool();
    BLNode::freeMemoryInPool();
#endif
}

double createSizingFunc(double x, double y, double z)
{

    if (BLMesh::GetSizingFunction()) {
        return BLMesh::GetSizingFunction()(std::array<double, 3>{x, y, z});
    }
    // return std::numeric_limits<double>::max();
}
void BLMesh::NarrowPointperturbation() { auto points = GetNarrowConstrainedPointidx(); }
std::set<int> BLMesh::GetNarrowConstrainedPointidx()
{
    std::set<int> ans;
    VM v;
    GenTopMesh(v);
    constexpr double threadhold = 179.5;
    /*
    v.pnSNO = nbelm;
    v.ppnSFFmO = new int[3 * nbelm];
    v.ppnSFTpO = new int[nbelm];
    v.nSN0 = nbpt;
    v.ppdSNC0 = new double[3 * nbpt];
    */

    TETGEN_BLMESH_::tetgenbehavior behavior;
    behavior.quality = 1;
    behavior.nobisect = 1;
    behavior.vtkview = 1;

    // behavior.refine = 1;
    behavior.plc = 1;
    behavior.quality = 1;
    // behavior.nobisect = 1;
    behavior.minratio = 1.2;
    behavior.docheck = 1;

    TETGEN_BLMESH_::tetgenio surface;
    TETGEN_BLMESH_::tetgenio volume;

    vector<int> BndFcts, PtMap;
    vector<double> BndPts;
    int &bndPtNum = v.nSN0;
    auto &coord = v.ppdSNC0;
    auto &bndFctNum = v.pnSNO;
    auto &bndFcts = v.ppnSFFmO;

    surface.pointlist = coord;
    surface.numberofpoints = bndPtNum;
    surface.numberoffacets = bndFctNum;
    surface.facetlist = new TETGEN_BLMESH_::tetgenio::facet[bndFctNum];
    for (int i = 0; i < bndFctNum; i++) {
        auto f = &surface.facetlist[i];
        TETGEN_BLMESH_::tetgenio::init(f);
        f->numberofpolygons = 1;
        f->polygonlist = new TETGEN_BLMESH_::tetgenio::polygon[1];
        auto p = &f->polygonlist[0];
        TETGEN_BLMESH_::tetgenio::init(p);
        // Read the number of vertices, it should be greater than 0.
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        for (int k = 0; k < 3; k++) {
            p->vertexlist[k] = bndFcts[3 * i + k];
        }
    }

    TETGEN_BLMESH_::tetrahedralize(&behavior, &surface, &volume);

    int npt = volume.numberofpoints;
    int nelm = volume.numberoftetrahedra;

    double max_angle = 0;
    for (int i = 0; i < nelm; i++) {
        BLVector coord[4];
        volume.numberofcorners;
        for (int k = 0; k < 4; k++) {
            for (int j = 0; j < 3; j++) {
                coord[k][j] = volume.pointlist[3 * volume.tetrahedronlist[4 * i + k] + j];
            }
        }
        BLVector normal[4];
        for (int k = 0; k < 4; k++) {
            normal[k] = ((coord[(k + 1) % 4] - coord[k]) ^ (coord[(k + 2) % 4] - coord[k])).normalized();
        }
        double max_diff = 0;
        for (int k = 0; k < 4; k++) {
            max_diff = max(max_diff, abs(normal[k] * normal[(k + 1) % 4]));
        }
        max_angle = max(max_angle, 180 - acos(max_diff) * 180 / PI);
        if (threadhold < 180 - acos(max_diff) * 180 / PI) {
            int dd = 0;
            for (int k = 0; k < 4; k++) {
                if (surface.numberofpoints > volume.tetrahedronlist[4 * i + k]) {
                    cout << volume.tetrahedronlist[4 * i + k] << endl;
                    cout << coord[k] << " ";
                    ans.insert(v.l2g[volume.tetrahedronlist[4 * i + k]]);
                    dd++;
                }
            }
            cout << "number of point in boundary=" << dd << endl;
        }
    }

    v.saveVTK("test2.vtk");
    cout << "max angle=" << max_angle << endl;
    /*
    double *coordX = new double[npt];
    double *coordY = new double[npt];
    double *coordZ = new double[npt];

    for (int i = 0; i < npt; i++)
    {
        coordX[i] = volume.pointlist[3 * i];
        coordY[i] = volume.pointlist[3 * i + 1];
        coordZ[i] = volume.pointlist[3 * i + 2];
    }

    auto elm = nullptr;
    std::swap(elm, volume.tetrahedronlist);
    volume.numberoftetrahedra = 0;
    */
    return ans;
}
/**
 * @brief 生成各向同性的全四面体网格
 *
 * @param smooth_attempt
 */
void BLMesh::GenOuterMesh(int smooth_attempt)
{

    int nbpt = 0, nbelm = 0, *belm = nullptr, *belmN = nullptr, *orientElm = nullptr;
    int npt = 0, nelm = 0, *elm = nullptr, *l_to_g = nullptr;
    int nbdry, *pbdry = nullptr;
    double *bpt = nullptr, *pt = nullptr;
    double *pt_size = nullptr;

#if 1
    // ReadAddBdry("add.pls", true);
#endif
    /*生成外网格*/
    GetOuterBoundary(&nbpt, &nbelm, &bpt, &belm, &l_to_g);

    CalOriginSize(nbpt, nbelm, bpt, belm, pt_size, l_to_g);

    vector<double> hole;
    OrientFace(nbpt, nbelm, bpt, belm, &belmN, &orientElm, hole);

    FreeMemoryInFrontAndNode();

#ifdef _DEBUG
    OutputOuterBoundary("test.pls", nbpt, nbelm, bpt, belm);
    // exit(-4);
#endif
#ifdef USE_TETGEN
    TETGEN_BLMESH_::tetgenbehavior behavior;
    behavior.quality = 1;
    behavior.nobisect = 1;
    behavior.vtkview = 1;

    // behavior.refine = 1;
    behavior.plc = 1;
    behavior.quality = 1;
    // behavior.nobisect = 1;
    behavior.minratio = 1.2;

    TETGEN_BLMESH_::tetgenio surface;
    TETGEN_BLMESH_::tetgenio volume;

    surface.numberofholes = hole.size() / 3;
    surface.holelist = new double[hole.size()];

    for (int i = 0; i < hole.size(); i++) {
        surface.holelist[i] = hole[i];
    }

    vector<int> BndFcts, PtMap;
    vector<double> BndPts;
    int &bndPtNum = nbpt;
    auto &bndPts = bpt;
    auto &bndFctNum = nbelm;
    auto &bndFcts = belm;
    double *coord = new double[bndPtNum * 3];
    for (int i = 0; i < 3 * bndPtNum; i++) {
        coord[i] = bndPts[i];
    }
    surface.pointlist = coord;
    surface.numberofpoints = bndPtNum;
    surface.numberoffacets = bndFctNum;
    surface.facetlist = new TETGEN_BLMESH_::tetgenio::facet[bndFctNum];
    for (int i = 0; i < bndFctNum; i++) {
        auto f = &surface.facetlist[i];
        TETGEN_BLMESH_::tetgenio::init(f);
        f->numberofpolygons = 1;
        f->polygonlist = new TETGEN_BLMESH_::tetgenio::polygon[1];
        auto p = &f->polygonlist[0];
        TETGEN_BLMESH_::tetgenio::init(p);
        // Read the number of vertices, it should be greater than 0.
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        for (int k = 0; k < 3; k++) {
            p->vertexlist[k] = bndFcts[3 * i + k];
        }
    }
    TETGEN_BLMESH_::tetrahedralize(&behavior, &surface, &volume);

    npt = volume.numberofpoints;
    nelm = volume.numberoftetrahedra;
    double *coordX = new double[npt];
    double *coordY = new double[npt];
    double *coordZ = new double[npt];

    for (int i = 0; i < npt; i++) {
        coordX[i] = volume.pointlist[3 * i];
        coordY[i] = volume.pointlist[3 * i + 1];
        coordZ[i] = volume.pointlist[3 * i + 2];
    }

    elm = nullptr;
    std::swap(elm, volume.tetrahedronlist);
    volume.numberoftetrahedra = 0;

#else
    DTIso3DConfig config;
    config.nSmoothAttempts = smooth_attempt;

    if (expan_ratio < 1.0) {
        throw std::invalid_argument("expan_ratio should be greater than 1!");
    }

    int re_obj_handler;
    std::array<int, 3> intCstNum;
    intCstNum[0] = 0;
    intCstNum[1] = 0;
    intCstNum[2] = 0;
    if (checkterminate()) {
        FreeMemoryInFrontAndNode();
        return;
    }
    try {
        API_Tetrahedralize_Single_Domain(nbpt, bpt, nbelm, belm, intCstNum.data(), nullptr, smooth_attempt, nullptr, &re_obj_handler);
    } catch (std::exception e) {
        // 发生了错误

        cout << "Catch an exception! Trying to free memory" << endl;
        ;
        API_DelTetrahedraObj(re_obj_handler);

        if (bpt) {
            delete[] bpt;
            bpt = nullptr;
        }

        if (belm) {
            delete[] belm;
            belm = nullptr;
        }
        if (l_to_g) {
            delete[] l_to_g;
            l_to_g = nullptr;
        }
        if (belmN) {
            delete[] belmN;
            belmN = nullptr;
        }
        if (orientElm) {
            delete[] orientElm;
            orientElm = nullptr;
        }
        if (elm) {
            delete[] elm;
            elm = nullptr;
        }
        if (pbdry) {
            delete[] pbdry;
            pbdry = nullptr;
        }
        if (pt_size) {
            delete[] pt_size;
            pt_size = nullptr;
        }

        throw e;
    }
    // delete SizingFuncHelper<1>::getSingleton();
    // SizingFuncHelper<1>::setSingleton(nullptr);
    // API_Tetrahedralize_Multi_Domain(
    //	nbpt, bpt, nbelm, belm, &re_obj_handler
    //);
    // exit(-1);

    npt = API_GetTetrahedraPntNum(re_obj_handler);
    nelm = API_GetTetrahedraElemNum(re_obj_handler);
    double *pts = new double[3 * npt];
    double *coordX = new double[npt];
    double *coordY = new double[npt];
    double *coordZ = new double[npt];
    elm = new int[nelm * 4];
    for (int i = 0; i < npt; i++) {
        API_GetTetrahedraPointCoord(re_obj_handler, i + 1, pts[3 * i + 0], pts[3 * i + 1], pts[3 * i + 2]);

        if (API_GetTetrahedraPointCoord(re_obj_handler, i + 1, coordX[i], coordY[i], coordZ[i])) {
            throw("get coord error");
        }
    }
    for (int i = 0; i < nelm; i++) {
        if (API_GetTetrahedraElemPntIdx(re_obj_handler, i + 1, elm[4 * i + 0], elm[4 * i + 1], elm[4 * i + 2], elm[4 * i + 3])) {
            throw("get coord error");
        }
    }
    API_DelTetrahedraObj(re_obj_handler);
#endif

    // GenMesh3D(smooth_attempt, nbpt, nbelm, bpt, belm, &npt, &pt, &nelm, &elm);
    cout << "number of Tet mesh in boundary layer mesh=" << m_nTet << endl;
    m_nTet += nelm;
    cout << "number of Tet mesh=" << m_nTet << endl;

    // TETGEN_BLMESH_::tetgenio volume;

    // volume.numberofpoints = npt;
    // volume.numberoftetrahedra = nelm;
    // volume.pointlist = pts;
    // volume.tetrahedronlist = elm;

    // double threadhold = 179;
    // double max_angle = 0;
    // for (int i = 0; i < nelm; i++) {
    //	BLVector coord[4];
    //	volume.numberofcorners;
    //	for (int k = 0; k < 4; k++)
    //		for (int j = 0; j < 3; j++) {
    //			coord[k][j] = volume.pointlist[3 * volume.tetrahedronlist[4 * i + k] + j];
    //		}
    //	BLVector normal[4];
    //	for (int k = 0; k < 4; k++)
    //		normal[k] = ((coord[(k + 1) % 4] - coord[k]) ^ (coord[(k + 2) % 4] - coord[k])).normalized();
    //	double max_diff = 0;
    //	for (int k = 0; k < 4; k++) {
    //		max_diff = max(max_diff, abs(normal[k] * normal[(k + 1) % 4]));
    //	}
    //	max_angle = max(max_angle, 180 - acos(max_diff) * 180 / PI);
    //	if (threadhold < 180 - acos(max_diff) * 180 / PI) {
    //		int dd = 0;
    //		for (int k = 0; k < 4; k++) {
    //			if (nbpt > volume.tetrahedronlist[4 * i + k]) {
    //				cout << volume.tetrahedronlist[4 * i + k] ;
    //				cout << coord[k]<<" "<< 180 - acos(max_diff) * 180 / PI << " " << endl;
    //				//ans.insert(v.l2g[volume.tetrahedronlist[4 * i + k]]);
    //				dd++;

    //			}

    //		}
    //		cout << "number of point in boundary=" << dd << endl;
    //	}
    //}

    UpdateDomainGrid(npt, nbpt, nelm, coordX, coordY, coordZ, elm, &l_to_g);
    // UpdateDomainGrid(npt, nbpt, nelm, pt, elm, &l_to_g, 0);
    if (pt_size) {
        delete[] pt_size;
        pt_size = nullptr;
    }
    if (coordX) {
        delete[] coordX;
        coordX = nullptr;
    }

    if (coordY) {
        delete[] coordY;
        coordY = nullptr;
    }

    if (coordZ) {
        delete[] coordZ;
        coordZ = nullptr;
    }

    if (bpt) {
        delete[] bpt;
        bpt = nullptr;
    }

    if (belm) {
        delete[] belm;
        belm = nullptr;
    }
    if (l_to_g) {
        delete[] l_to_g;
        l_to_g = nullptr;
    }
    if (belmN) {
        delete[] belmN;
        belmN = nullptr;
    }
    if (orientElm) {
        delete[] orientElm;
        orientElm = nullptr;
    }
    if (elm) {
        delete[] elm;
        elm = nullptr;
    }
    if (pbdry) {
        delete[] pbdry;
        pbdry = nullptr;
    }
    /*volume.pointlist = nullptr;
    volume.tetrahedronlist = nullptr;
    volume.numberofpoints = 0;
    volume.numberoftetrahedra = 0;*/

    cout << "Well done!  Generating tet mesh is completed." << endl;
}

void BLMesh::CalOriginSize(int npt, int nlem, double *pt, int *elm, double *&pt_size, int *l_to_g)
{
    pt_size = new double[npt];
    std::vector<std::vector<double>> size_on_point(npt);

    for (int i = 0; i < nlem; i++) {

        double size[3] = {0};

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                size[j] += pow(pt[3 * elm[3 * i + j] + k] - pt[3 * elm[3 * i + (j + 1) % 3] + k], 2);
            }
        }
        for (int k = 0; k < 3; k++) {
            size[k] = sqrt(size[k]);
            size_on_point[elm[3 * i + k]].push_back(size[k]);
        }
    }

    for (int i = 0; i < size_on_point.size(); i++) {
        double ave_size = 1.0;
        double min_size = std::numeric_limits<double>::max();
        double max_size = -1e100;

        for (int j = 0; j < size_on_point[i].size(); j++) {
            ave_size *= size_on_point[i][j];
            min_size = std::min(size_on_point[i][j], min_size);
            max_size = std::max(size_on_point[i][j], max_size);
        }
        ave_size = std::pow(ave_size, 1.0 / size_on_point[i].size());
        // ave_size = sqrt(min_size*max_size);
        // double target_size = ave_size;
        // BLNode *bln =
        //                  reinterpret_cast<BLNode*>(m_pNodes[l_to_g[i]].pointer);
        // if(bln&&bln->GetLowerNode()){
        //                      int lowerid =
        //                          bln->GetLowerNode()->GetNodIdx();
        //				int upperid=l_to_g[i];
        // target_size= sqrt(	pow(m_pNodes[lowerid].coord[0]-m_pNodes[upperid].coord[0],2)+
        //                                                  pow(m_pNodes[lowerid]
        //                                                              .coord[1] -
        //                                                          m_pNodes[upperid]
        //                                                              .coord[1],
        //                                                      2) +
        //                                                  pow(m_pNodes[lowerid]
        //                                                              .coord[2] -
        //                                                          m_pNodes[upperid]
        //                                                              .coord[2],
        //                                                      2));
        //}
        // else{
        //}
        // pt_size[i] = std::min(ave_size,target_size);
        pt_size[i] = ave_size;
    }
}

void BLMesh::UpdateDomainGrid(int ngp, int nbp, int nel, double *g_coordx, double *g_coordy, double *g_coordz, int *gtopu, int **l_to_g)
{
    int i, j, pidx, eidx, conn[4];
    BLVector pnt;
    BLEntityTopology topu;

    if (ngp > nbp) {
        int *newpt = (int *)realloc(*l_to_g, sizeof(int) * ngp);

        (*l_to_g) = newpt;
        memset((*l_to_g) + nbp, -1, sizeof(int) * (ngp - nbp));

        int nAlloc = m_nAllocNodes + ngp - nbp;
        MBLNode *pNewNodes = new MBLNode[nAlloc];
        std::copy(m_pNodes, m_pNodes + m_nNodes, pNewNodes);
        delete[] m_pNodes;

        m_pNodes = pNewNodes;
        m_nAllocNodes = nAlloc;
#ifndef _DEBUG

#endif
        for (int i = nbp; i < ngp; i++) {
            int idx = m_nNodes + i - nbp;
            m_pNodes[idx].coord[0] = g_coordx[i];
            m_pNodes[idx].coord[1] = g_coordy[i];
            m_pNodes[idx].coord[2] = g_coordz[i];
            m_pNodes[idx].reserved = -1;
            (*l_to_g)[i] = m_nNodes + i - nbp;
        }
        m_nNodes += ngp - nbp;

        // for (i = nbp; i < ngp; i++)
        //{
        //	pnt.x = g_coordx[i];
        //	pnt.y = g_coordy[i];
        //	pnt.z = g_coordz[i];

        //	pidx = AddNode(pnt, 0.0);

        //	(*l_to_g)[i] = pidx;
        //}
    }
    int elmSize = m_nElems;
    int nAlloc = m_nAllocElems;
    nAlloc += nel;

    Elem *pNewElems = (Elem *)realloc(m_pElems, sizeof(Elem) * nAlloc);

    if (nAlloc > m_nAllocElems) {
        memset(&pNewElems[m_nAllocElems], 0, sizeof(Elem) * (nAlloc - m_nAllocElems));
    }

    m_pElems = pNewElems;
    m_nAllocElems = nAlloc;

    for (int i = 0; i < nel; i++) {
        for (int k = 0; k < 4; k++) {
            m_pElems[m_nElems + i].conn[k] = (*l_to_g)[gtopu[4 * i + k]];
        }
        m_pElems[m_nElems + i].topo = TETRAHEDRON;
        m_pElems[m_nElems + i].nconn = 4;
    }

    m_nElems += nel;
    for (int i = 0; i < m_nElems; i++) {
        if (m_pElems[i].nconn == 4) {
            if (GEOM_FUNC::orient3d(m_pNodes[m_pElems[i].conn[0]].coord, m_pNodes[m_pElems[i].conn[1]].coord,
                                    m_pNodes[m_pElems[i].conn[2]].coord, m_pNodes[m_pElems[i].conn[3]].coord) > 0) {
                //	cout<<GEOM_FUNC::orient3d(m_pNodes[m_pElems[i].conn[0]].coord, m_pNodes[m_pElems[i].conn[1]].coord,
                //m_pNodes[m_pElems[i].conn[2]].coord, m_pNodes[m_pElems[i].conn[3]].coord); 	cout << " " <<  i << endl;
                // throw(std::exception());
            }
        }
    }
    m_nOutTri += nel;
}

void BLMesh::UpdateDomainGrid(int ngp, int nbp, int nel, double *g_coord, int *gtopu, int **l_to_g, int index)
{
    int i, j, pidx, eidx, conn[4];
    BLVector pnt;
    BLEntityTopology topu;

    if (ngp > nbp) {
        int *newpt = (int *)realloc(*l_to_g, sizeof(int) * ngp);
        (*l_to_g) = newpt;
        memset((*l_to_g) + nbp, -1, sizeof(int) * (ngp - nbp));

        for (i = nbp; i < ngp; i++) {
            pnt.x = g_coord[i * 3 + 0];
            pnt.y = g_coord[i * 3 + 1];
            pnt.z = g_coord[i * 3 + 2];

            pidx = AddNode(pnt, 0.0);

            (*l_to_g)[i] = pidx;
        }
    }

    for (i = 0; i < nel; i++) {
        j = i * 4;

        conn[0] = (*l_to_g)[gtopu[j + 0]];
        conn[1] = (*l_to_g)[gtopu[j + 1]];
        conn[2] = (*l_to_g)[gtopu[j + 2]];
        conn[3] = (*l_to_g)[gtopu[j + 3]];

        topu = BLEntityTopology::TETRAHEDRON;
        eidx = AddElem(4, conn, topu);

        m_pElems[eidx].igom = index;
    }

    m_nOutTri += nel;
}

void BLMesh::UpdateSymplnGrid(int ngp, int nbp, int nel, double *g_coord, int *gtopu, int **l_to_g, int ifc)
{
    int i, j, pidx, eidx, conn[3];
    BLVector pnt;
    BLEntityTopology topu;

    if (ngp > nbp) {
        int *newpt = (int *)realloc(*l_to_g, sizeof(int) * ngp);
        (*l_to_g) = newpt;
        memset((*l_to_g) + nbp, -1, sizeof(int) * (ngp - nbp));

        for (i = nbp; i < ngp; i++) {
            pnt.x = g_coord[i * 3 + 0];
            pnt.y = g_coord[i * 3 + 1];
            pnt.z = g_coord[i * 3 + 2];
            ;

            pidx = AddNode(pnt, 0.0);

            (*l_to_g)[i] = pidx;
        }
    }

    for (i = 0; i < nel; i++) {
        j = i * 3;

        conn[0] = (*l_to_g)[gtopu[j + 0]];
        conn[1] = (*l_to_g)[gtopu[j + 1]];
        conn[2] = (*l_to_g)[gtopu[j + 2]];

        topu = BLEntityTopology::TRIANGLE;
        eidx = AddElem(3, conn, topu, ifc);
    }
}
void BLMesh::RemoveNonManifoldEdgeInSymFace()
{
    std::map<pair<int, int>, int> edge_hash;
    int i1, i2, i3;
    for (int i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        edge_hash[pair<int, int>(i1, i2)]++;
        edge_hash[pair<int, int>(i2, i3)]++;
        edge_hash[pair<int, int>(i1, i3)]++;
    }
    spdlog::info("begin to remove non-manifold edges in sym face...\n");
    auto it = m_vBdyFront.begin();

    spdlog::info("finished removing non-manifold edges in sym face!\n");
}
void BLMesh::RemoveNonManifoldFrontByForce()
{
    std::map<pair<int, int>, vector<int>> edge_hash;
    int i1, i2, i3;
    for (int i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        edge_hash[pair<int, int>(i1, i2)].push_back(i);
        edge_hash[pair<int, int>(i2, i3)].push_back(i);
        edge_hash[pair<int, int>(i1, i3)].push_back(i);
    }
    int valid = 0;
    for (auto i : edge_hash) {
        if (i.second.size() > 2) {
            valid++;
        }
    }
    if (!valid) {
        return;
    } else {
        map<int, int> no_manifold_edges_count;
        for (auto i : edge_hash) {
            no_manifold_edges_count[i.second.size()]++;
        }

        cout << "still have " << valid << " non manifold edge" << endl;
        cout << "face count ------  nonmanifold edge count" << endl;
        for (auto i : no_manifold_edges_count) {
            cout << i.first << " | " << i.second << endl;
        }
        cout << "============================================" << endl;
        cout << "use brute-force algorithm to try to remove" << endl;
        cout << "warning,this algorithm is not robust, once dtiso3d support nonmanifold edge" << endl;
        cout << "please remove this function" << endl;
    }
    vector<int> pRmBdry;
    int addbdry = 0; // 增加计数
    for (auto i : edge_hash) {
        if (i.second.size() > 2) {
            int start = i.first.first;
            int end = i.first.second;
            auto triangles = i.second;
            int *ptr = outbdry;
            auto third = [start, end, triangles, ptr](int id) {
                for (int l = 0; l < 3; l++) {
                    if (ptr[3 * triangles[id] + l] != start && ptr[3 * triangles[id] + l] != end) {
                        return ptr[3 * triangles[id] + l];
                    }
                    return -1;
                }
            };
            bool succeed = false;
            for (int j = 0; j < triangles.size(); j++) {
                for (int k = j + 1; k < triangles.size(); k++) {
                    int n1 = third(j), n2 = third(k);
                    if (edge_hash.find(pair<int, int>(std::min(n1, n2), std::max(n1, n2))) != edge_hash.end()) {
                        // try to add triangles
                        int tmp1[DIM3], tmp2[DIM3];
                        tmp1[0] = n1;
                        tmp1[1] = n2;
                        tmp1[2] = start;

                        tmp2[0] = n1;
                        tmp2[1] = n2;
                        tmp2[2] = end;

                        auto tridx1 = AddTriElem(3, tmp1);
                        auto tridx2 = AddTriElem(3, tmp2);

                        bool leftin = IsOutbdry(tmp1), rightin = IsOutbdry(tmp2);
                        bool inter = false;
                        if (!leftin) {
                            inter |= m_ocTree->chckIntersectPreProcess(tridx1);
                        }
                        if (!rightin) {
                            inter |= m_ocTree->chckIntersectPreProcess(tridx2);
                        }
                        if (!inter && (!(leftin && rightin))) {
                            if (!leftin) {
                                m_ocTree->insertPreProcess(tridx1);
                                outbdry[noutbdry * 3 + 0] = tmp1[0];
                                outbdry[noutbdry * 3 + 1] = tmp1[1];
                                outbdry[noutbdry * 3 + 2] = tmp1[2];
                                ++noutbdry;
                                addbdry++;
                            } else {
                                pRmBdry.push_back(tmp1[0]);
                                pRmBdry.push_back(tmp1[1]);
                                pRmBdry.push_back(tmp1[2]);
                            }
                            if (!rightin) {
                                m_ocTree->insertPreProcess(tridx2);
                                outbdry[noutbdry * 3 + 0] = tmp2[0];
                                outbdry[noutbdry * 3 + 1] = tmp2[1];
                                outbdry[noutbdry * 3 + 2] = tmp2[2];
                                ++noutbdry;
                                addbdry++;
                            } else {
                                pRmBdry.push_back(tmp2[0]);
                                pRmBdry.push_back(tmp2[1]);
                                pRmBdry.push_back(tmp2[2]);
                            }

                            pRmBdry.push_back(start);
                            pRmBdry.push_back(end);
                            pRmBdry.push_back(n1);

                            pRmBdry.push_back(start);
                            pRmBdry.push_back(end);
                            pRmBdry.push_back(n2);

                            int conn[4];
                            conn[0] = start;
                            conn[2] = end;
                            conn[1] = n1;
                            conn[3] = n2;
                            AddElem(4, conn, BLEntityTopology::TETRAHEDRON, 1);
                            m_nOutTri += addbdry - pRmBdry.size() / 3;
                            succeed = true;
                        }
                    }
                    if (succeed) {
                        break;
                    }
                }
                if (succeed) {
                    break;
                }
            }
            // 失败了强行添加
            if (!succeed) {
                for (int j = 0; j < triangles.size(); j++) {
                    for (int k = j + 1; k < triangles.size(); k++) {
                        int n1 = third(j), n2 = third(k);
                        if (edge_hash.find(pair<int, int>(std::min(n1, n2), std::max(n1, n2))) == edge_hash.end()) {
                            int tmp1[DIM3], tmp2[DIM3];
                            tmp1[0] = n1;
                            tmp1[1] = n2;
                            tmp1[2] = start;

                            tmp2[0] = n1;
                            tmp2[1] = n2;
                            tmp2[2] = end;
                            auto tridx1 = AddTriElem(3, tmp1);
                            auto tridx2 = AddTriElem(3, tmp2);
                            bool inter = false;
                            inter |= m_ocTree->chckIntersectPreProcess(tridx1);
                            inter |= m_ocTree->chckIntersectPreProcess(tridx2);

                            if (!inter) {
                                m_ocTree->insertPreProcess(tridx1);
                                outbdry[noutbdry * 3 + 0] = tmp1[0];
                                outbdry[noutbdry * 3 + 1] = tmp1[1];
                                outbdry[noutbdry * 3 + 2] = tmp1[2];
                                ++noutbdry;
                                addbdry++;
                                m_ocTree->insertPreProcess(tridx2);
                                outbdry[noutbdry * 3 + 0] = tmp2[0];
                                outbdry[noutbdry * 3 + 1] = tmp2[1];
                                outbdry[noutbdry * 3 + 2] = tmp2[2];
                                ++noutbdry;
                                addbdry++;

                                pRmBdry.push_back(start);
                                pRmBdry.push_back(end);
                                pRmBdry.push_back(n1);

                                pRmBdry.push_back(start);
                                pRmBdry.push_back(end);
                                pRmBdry.push_back(n2);

                                int conn[4];
                                conn[0] = start;
                                conn[2] = end;
                                conn[1] = n1;
                                conn[3] = n2;
                                AddElem(4, conn, BLEntityTopology::TETRAHEDRON, 1);
                                m_nOutTri += addbdry - pRmBdry.size() / 3;
                                succeed = true;
                            }
                        }
                        if (succeed) {
                            break;
                        }
                    }
                    if (succeed) {
                        break;
                    }
                }
            }
        }
    }
    RemoveOutbdry(pRmBdry.size() / 3, pRmBdry.data());

    edge_hash.clear();

    for (int i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        edge_hash[pair<int, int>(i1, i2)].push_back(i);
        edge_hash[pair<int, int>(i2, i3)].push_back(i);
        edge_hash[pair<int, int>(i1, i3)].push_back(i);
    }
    valid = 0;
    for (auto i : edge_hash) {
        if (i.second.size() > 2) {
            valid++;
        }
    }
    if (!valid) {
        return;
    } else {
        map<int, int> no_manifold_edges_count;
        for (auto i : edge_hash) {
            no_manifold_edges_count[i.second.size()]++;
        }

        cout << "after operation" << endl;
        cout << "face count ------  nonmanifold edge count" << endl;
        for (auto i : no_manifold_edges_count) {
            cout << i.first << " | " << i.second << endl;
        }
        cout << "============================================" << endl;
    }
}
void BLMesh::RemoveNonManifoldFront()
{
    std::map<pair<int, int>, int> edge_hash;
    int i1, i2, i3;
    for (int i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        edge_hash[pair<int, int>(i1, i2)]++;
        edge_hash[pair<int, int>(i2, i3)]++;
        edge_hash[pair<int, int>(i1, i3)]++;
    }

    int i, nNod, conn[DIM3 + 1];

    int *pRmBdry = nullptr, nRmBdry = 0, tmp1[DIM3], tmp2[DIM3];
    BLFront *blFrt;
    BLNode *blNods[DIM3], *blNodCur;
    std::vector<BLFront *>::iterator it;

    spdlog::info("begin to remove non-manifold edges...\n");

    pRmBdry = new int[noutbdry * 3];

    /*遍历所有前沿front*/
    it = m_vBdyFront.begin();

    while (it != m_vBdyFront.end()) {
        bool succeed = false;
        vector<pair<int, pair<int, int>>> num_neigh_idx;
        blFrt = *it;
        blFrt->GetNodes(&nNod, blNods);
        for (i = 0; i < nNod; i++) {
            if (blNods[i]->GetUpperNode()) {
                int num_nei;
                BLFront *b_fr[MAX_NCONN * 2];
                blNods[i]->GetNeigFronts(b_fr, &num_nei);
                if (m_pNodes[blNods[i]->GetNodIdx()].bsysm) {
                    num_nei += 1;
                }
                num_neigh_idx.push_back(pair<int, pair<int, int>>(num_nei, pair<int, int>(blNods[i]->GetNodIdx(), i)));
            }
        }

        if (blFrt->GetUpperFront()) {
            it++;
            continue;
        }
        BLFront *neigh_front[3];
        int num_nei_front;
        int j;
        blFrt->GetNeigbourFronts(&num_nei_front, neigh_front);
        for (j = 0; j < num_nei_front; j++) {
            if (neigh_front[j]->GetUpperFront()) {
                break;
            }
        }
        if (j != num_nei_front) { // 周围的front有上层front，意味着这个三角形的某个边将会长成四边形，不需要填补
            it++;
            continue;
        }
        if (num_neigh_idx.empty()) { // 无上层节点
            it++;
            continue;
        }
        sort(num_neigh_idx.rbegin(), num_neigh_idx.rend());

        int cidx = num_neigh_idx[0].second.second;
        blNodCur = blNods[cidx];
        int side_idx1 = blNods[(cidx + 1) % 3]->GetNodIdx(); /// 解决多金字塔问题
        int side_idx2 = blNods[(cidx + 2) % 3]->GetNodIdx(); /// 解决多金字塔问题
        BLNode *tmp = blNodCur;
        bool need_to_generate = false;
        while (tmp->GetUpperNode()) {
            int pointl = tmp->GetNodIdx();
            int pointr = tmp->GetUpperNode()->GetNodIdx();
            if (edge_hash[pair<int, int>(std::min(pointl, pointr), std::max(pointl, pointr))] > 2 ||
                (edge_hash[pair<int, int>(std::min(pointl, pointr), std::max(pointl, pointr))] == 2 && m_pNodes[pointl].bsysm)) {
                need_to_generate = true;
            }
            tmp = tmp->GetUpperNode();
        }
        /*只有一个点有上层点*/
        // 第一个三角形
        while (blNodCur->GetUpperNode()) {
            tmp1[0] = blNodCur->GetNodIdx();
            tmp1[1] = blNodCur->GetUpperNode()->GetNodIdx();
            tmp1[2] = side_idx1;
            // 第二个三角形
            tmp2[0] = blNodCur->GetNodIdx();
            tmp2[1] = side_idx2;
            tmp2[2] = blNodCur->GetUpperNode()->GetNodIdx();

            int pointl = blNodCur->GetNodIdx();
            int pointr = blNodCur->GetUpperNode()->GetNodIdx();
            int rmbdry = 0;  // 删除计数
            int addbdry = 0; // 增加计数
            if (need_to_generate) {

                bool leftin = IsOutbdry(tmp1), rightin = IsOutbdry(tmp2);

                if (leftin) {
                    pRmBdry[nRmBdry * 3 + 0] = blNodCur->GetNodIdx();
                    pRmBdry[nRmBdry * 3 + 1] = blNodCur->GetUpperNode()->GetNodIdx();
                    pRmBdry[nRmBdry * 3 + 2] = side_idx1;
                    nRmBdry++;
                    rmbdry++;
                } else {

                    // 这里对称面添加三角形后面有做，所以这里如果是对称面就不需要做了

                    /*consider sym surface*/
                    if (!(m_pNodes[blNodCur->GetNodIdx()].bsysm && m_pNodes[side_idx1].bsysm)) {
                        outbdry[noutbdry * 3 + 0] = blNodCur->GetNodIdx();
                        outbdry[noutbdry * 3 + 1] = blNodCur->GetUpperNode()->GetNodIdx();
                        outbdry[noutbdry * 3 + 2] = side_idx1;
                        ++noutbdry;
                        addbdry++;
                    } else {
                        auto biger = std::max(blNodCur->GetUpperNode()->GetNodIdx(), side_idx1);
                        auto smaller = std::min(blNodCur->GetUpperNode()->GetNodIdx(), side_idx1);
                        // m_mapSymline.insert(std::make_pair(smaller, biger));
                        insertSymFace(std::make_pair(smaller, biger));
                    }
                }
                if (rightin) {
                    pRmBdry[nRmBdry * 3 + 0] = blNodCur->GetNodIdx();
                    pRmBdry[nRmBdry * 3 + 1] = side_idx2;
                    pRmBdry[nRmBdry * 3 + 2] = blNodCur->GetUpperNode()->GetNodIdx();
                    // RemoveOutbdry(conn);
                    nRmBdry++;
                    rmbdry++;
                } else {

                    if (!(m_pNodes[blNodCur->GetNodIdx()].bsysm && m_pNodes[side_idx2].bsysm)) {
                        outbdry[noutbdry * 3 + 0] = blNodCur->GetNodIdx();
                        outbdry[noutbdry * 3 + 1] = side_idx2;
                        outbdry[noutbdry * 3 + 2] = blNodCur->GetUpperNode()->GetNodIdx();
                        ++noutbdry;
                        addbdry++;
                    } else {
                        auto biger = std::max(blNodCur->GetUpperNode()->GetNodIdx(), side_idx2);
                        auto smaller = std::min(blNodCur->GetUpperNode()->GetNodIdx(), side_idx2);
                        // m_mapSymline.insert(std::make_pair(smaller, biger));
                        insertSymFace(std::make_pair(smaller, biger));
                    }

                    // outbdry[noutbdry * 3 + 0] = blNodCur->GetNodIdx();
                    // outbdry[noutbdry * 3 + 1] = blNods[(cidx + 2) % 3]->GetNodIdx();
                    // outbdry[noutbdry * 3 + 2] = blNodCur->GetUpperNode()->GetNodIdx();
                    //++noutbdry; addbdry++;
                }
                outbdry[noutbdry * 3 + 0] = blNodCur->GetUpperNode()->GetNodIdx();
                outbdry[noutbdry * 3 + 1] = side_idx1;
                outbdry[noutbdry * 3 + 2] = side_idx2;
                ++noutbdry;
                addbdry++;

                blFrt->GetNodes(&nNod, blNods);
                pRmBdry[nRmBdry * 3 + 0] = blNodCur->GetNodIdx();
                pRmBdry[nRmBdry * 3 + 1] = side_idx1;
                pRmBdry[nRmBdry * 3 + 2] = side_idx2;
                // RemoveOutbdry(conn);
                nRmBdry++;
                rmbdry++;

                // intersection test
#ifdef _CHECK_INTERSECTION
                // intersection
                int cons[3], tridx, ltriidx, rtriidx;
                cons[0] = blNodCur->GetUpperNode()->GetNodIdx();
                cons[1] = side_idx1;
                cons[2] = side_idx2;
                tridx = AddTriElem(3, cons);
                m_ocTree->insertPreProcess(tridx);
                bool ret = false;
                ret = ret | m_ocTree->chckIntersectPreProcess(tridx);

                if (!leftin) {
                    cons[0] = blNodCur->GetNodIdx();
                    cons[1] = blNodCur->GetUpperNode()->GetNodIdx();
                    cons[2] = side_idx1;
                    ltriidx = AddTriElem(3, cons);
                    m_ocTree->insertPreProcess(ltriidx);
                    ret = ret | m_ocTree->chckIntersectPreProcess(ltriidx);
                }
                if (!rightin) {
                    cons[0] = blNodCur->GetNodIdx();
                    cons[1] = side_idx2;
                    cons[2] = blNodCur->GetUpperNode()->GetNodIdx();
                    rtriidx = AddTriElem(3, cons);
                    m_ocTree->insertPreProcess(rtriidx);
                    ret = ret | m_ocTree->chckIntersectPreProcess(rtriidx);
                }

                if (ret) { // the new created triangle intersect with existing mesh elements
                    cout << ("gave up creating a triangle for potential intersections!") << endl;
                    m_ocTree->rmDataPreProcess(tridx);
                    if (!leftin) {
                        m_ocTree->rmDataPreProcess(ltriidx);
                    }
                    if (!rightin) {
                        m_ocTree->rmDataPreProcess(rtriidx);
                    }
                    noutbdry -= addbdry;
                    nRmBdry -= rmbdry;
                    break;
                } else {
                    // printf("nbdry: %d, alloc: %d\n", noutbdry*3, 6*m_nSurfElems);

                    // create a tetrahedron
                    conn[0] = blNodCur->GetNodIdx();
                    conn[2] = side_idx1;
                    conn[1] = side_idx2;
                    conn[3] = blNodCur->GetUpperNode()->GetNodIdx();
                    int id = AddElem(4, conn, BLEntityTopology::TETRAHEDRON, 1);
                    m_nOutTri += addbdry - rmbdry;
                    succeed = true;
                    cout << "add a tetrahedron mesh id=" << id << " to preseve the manifold outer-most triangle mesh" << endl;
                    // remove blFrt from m_vBdyFront
                }
            }

#endif

            else {

                break;
            }
            blNodCur = blNodCur->GetUpperNode();
        }
        if (succeed) {
            it = m_vBdyFront.erase(it);
        } else {
            it++;
        }
    }
    // RemoveSymLineByManifoldCretiria();
    RemoveOutbdry(nRmBdry, pRmBdry);
    if (pRmBdry) {
        delete[] pRmBdry;
    }
    // RemoveNonManifoldFrontByForce();
    spdlog::info("finished removing non-manifold edges!\n");
}

void BLMesh::RemoveOutbdry(int nRmBdry, int *pRmBdry)
{
    int i, j, cnt;
    int *flag = nullptr;

    //	spdlog::info("begin to remove unused out boundary...\n");

    flag = new int[noutbdry];
    for (i = 0; i < noutbdry; i++) {
        flag[i] = 1;
    }
    map<int, set<pair<int, int>>> remove_tri;

    for (int k = 0; k < nRmBdry; k++) {
        int idx1, idx2, idx3, conn[3];

        conn[0] = pRmBdry[k * 3 + 0];
        conn[1] = pRmBdry[k * 3 + 1];
        conn[2] = pRmBdry[k * 3 + 2];

        idx1 = idx3 = conn[0];

        if (conn[1] < idx1) {
            idx1 = conn[1];
        }
        if (conn[2] < idx1) {
            idx1 = conn[2];
        }

        if (conn[1] > idx3) {
            idx3 = conn[1];
        }
        if (conn[2] > idx3) {
            idx3 = conn[2];
        }

        idx2 = (conn[0] + conn[1] + conn[2]) - (idx1 + idx3);
        remove_tri[idx1].insert(pair<int, int>(idx2, idx3));
        // flag deleted element
    }
    int delete_count = 0;
    for (i = 0; i < noutbdry; i++) {
        int i1, i2, i3;
        GetFacIdx(outbdry, i, &i1, &i2, &i3);

        if (remove_tri[i1].find(pair<int, int>(i2, i3)) != remove_tri[i1].end()) {
            flag[i] = 0;
            delete_count++;
        }
    }
    // delete element
    cnt = 0;
    for (i = 0; i < noutbdry; i++) {
        if (flag[i]) {
            outbdry[cnt * 3 + 0] = outbdry[i * 3 + 0];
            outbdry[cnt * 3 + 1] = outbdry[i * 3 + 1];
            outbdry[cnt * 3 + 2] = outbdry[i * 3 + 2];
            cnt++;
        }
    }
    if (nRmBdry + cnt != noutbdry) {
        cout << endl << "should remove" << nRmBdry << " but remove" << delete_count << " " << noutbdry - cnt << " front" << endl;
    }
    noutbdry = cnt;

    if (flag) {
        delete[] flag;
    }

    //	spdlog::info("finished removing unused out boundary!\n");
}
int BLMesh::GetNumTriangleByEdge(int p1, int p2)
{
    if (p1 > p2) {
        swap(p1, p2);
    }
    int i1, i2, i3;
    int ans = 0;
    for (int i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);
        if (i1 == p1 && i2 == p2) {
            ans++;
        } else if (i2 == p1 && i3 == p2) {
            ans++;
        } else if (i1 == p1 && i3 == p2) {
            ans++;
        }
    }
    return ans;
}
bool BLMesh::IsOutbdry(int conn[3])
{
    int i, j, cnt, idx1, idx2, idx3;
    int i1, i2, i3;

    idx1 = idx3 = conn[0];

    if (conn[1] < idx1) {
        idx1 = conn[1];
    }
    if (conn[2] < idx1) {
        idx1 = conn[2];
    }

    if (conn[1] > idx3) {
        idx3 = conn[1];
    }
    if (conn[2] > idx3) {
        idx3 = conn[2];
    }

    idx2 = (conn[0] + conn[1] + conn[2]) - (idx1 + idx3);

    for (i = 0; i < noutbdry; i++) {
        GetFacIdx(outbdry, i, &i1, &i2, &i3);

        if (i1 == idx1 && i2 == idx2 && i3 == idx3) {
            return true;
        }
    }

    return false;
}

void BLMesh::RemoveOverlapFace(int *bndry, int nbdry)
{
    /*原算法O（N2），超级慢，现换成O（Nlog(N)）去除重复面*/
    int i, j, cnt;
    int i1, i2, i3, j1, j2, j3;
    vector<bool> flag(nbdry, true);
    std::map<std::array<int, 3>, int> overlap_check;

    // flag deleted element
    cnt = 0;
    for (i = 0; i < nbdry; i++) {
        GetFacIdx(bndry, i, &i1, &i2, &i3);
        std::array<int, 3> idx{i1, i2, i3};

        if (overlap_check.find(idx) == overlap_check.end()) {
            overlap_check[idx] = i;
        } else {
            flag[overlap_check[idx]] = false;
            flag[i] = false;
            overlap_check.erase(idx);
        }
    }
    // delete element
    cnt = 0;
    for (i = 0; i < nbdry; i++) {
        if (flag[i]) {
            bndry[cnt * 3 + 0] = bndry[i * 3 + 0];
            bndry[cnt * 3 + 1] = bndry[i * 3 + 1];
            bndry[cnt * 3 + 2] = bndry[i * 3 + 2];
            cnt++;
        }
    }
    noutbdry = cnt;
}

void BLMesh::GetFacIdx(int *bndry, int i, int *i1, int *i2, int *i3)
{
    int idx1, idx2, idx3;

    idx1 = bndry[3 * i + 0];
    idx2 = bndry[3 * i + 1];
    idx3 = bndry[3 * i + 2];

    *i1 = *i3 = idx1;

    if (idx2 < *i1) {
        *i1 = idx2;
    }
    if (idx3 < *i1) {
        *i1 = idx3;
    }

    if (idx2 > *i3) {
        *i3 = idx2;
    }
    if (idx3 > *i3) {
        *i3 = idx3;
    }

    *i2 = (idx1 + idx2 + idx3) - (*i1 + *i3);
}

void BLMesh::RmvUperNeigFrontsAndFreeNode(BLNode *blNod)
{
    int i, j;
    int nblFrts;                                         /*邻近面数量 **/
    BLFront *blFrts[MAX_NCONN * 2] /*邻近面 **/, *bluFt; /*邻近面的上层面 **/
    std::vector<BLNode *> nodes;

    // get neigh fronts
    blNod->GetNeigFronts(blFrts, &nblFrts);

    for (i = 0; i < nblFrts; i++) {
        bool used[3] = {false}; // mark
        bluFt = blFrts[i]->GetUpperFront();

        // delete bluFt
        if (bluFt != nullptr) {
            // 删除周围
            int nNeig, nNod;
            BLFront *neiFts[3];
            BLNode *blNods[3];
            // remove neigbor information
            bluFt->GetNeigbourFronts(&nNeig, neiFts);
            for (j = 0; j < nNeig; j++) {
                neiFts[j]->RmvNeigbourFronts(bluFt);
                if (neiFts[j]->GetLowerFront()) {
                    BLFront *lower = neiFts[j]->GetLowerFront();
                    lower->GetNodes(&nNod, blNods);
                    scheck.insert(neiFts[j]->GetLowerFront());
                }
            }

            // remove node information
            bluFt->GetNodes(&nNod, blNods);
            for (j = 0; j < nNod; j++) {
                blNods[j]->RmvNeigFront(bluFt);
                /// TODO 这里存在内存泄漏
                if (blNods[j]->checkValid()) {
                    nodes.push_back(blNods[j]);
                }
                // if (blNods[j]->checkValid())
                //	node_to_delete_.push_back(blNods[j]);
            }

            blFrts[i]->SetUpperFront(nullptr);
            SetElmDelete(blFrts[i]->GetElmIdx());
            m_nPrism--;

#ifdef _CHECK_INTERSECTION
            // intersection

            if (bluFt->IsPyramid()) {
                bluFt->SetPyramidFlag(false);
            }

            m_ocTree->rmDataPreProcess(bluFt->GetTriIdx());
            m_ocTree->insertPreProcess(blFrts[i]->GetTriIdx());

            blFrts[i]->GetNeigbourFronts(&nNeig, neiFts);
            for (j = 0; j < nNeig; j++) {
                int nidx1, nidx2, sidx1, sidx2;
                NeighIdx(blFrts[i], neiFts[j], &nidx1);
                NeighIdx(neiFts[j], blFrts[i], &nidx2);

                used[nidx1] = true;
                sidx1 = blFrts[i]->GetSTriIdx(nidx1, 0);
                sidx2 = blFrts[i]->GetSTriIdx(nidx1, 1);

                if (neiFts[j]->GetUpperFront()) {
                    ////插入到待处理的边界中
                    // auto it = scheck.begin();
                    // for (; it != scheck.end(); it++)
                    //	if (*it == neiFts[j]) {
                    //		break;
                    //	}
                    // if (it == scheck.end())
                    //	scheck.push_back(neiFts[j]);

                    if (sidx1 > 0 && sidx2 > 0) {
                        m_ocTree->insertPreProcess(sidx1);
                        m_ocTree->insertPreProcess(sidx2);

                        neiFts[j]->SetSTriIdx(nidx2, 0, sidx1);
                        neiFts[j]->SetSTriIdx(nidx2, 1, sidx2);

                        blFrts[i]->SetSTriIdx(nidx1, 0, -1);
                        blFrts[i]->SetSTriIdx(nidx1, 1, -1);
                    } else {
                        sidx1 = neiFts[j]->GetSTriIdx(nidx2, 0);
                        sidx2 = neiFts[j]->GetSTriIdx(nidx2, 1);

                        m_ocTree->insertPreProcess(sidx1);
                        m_ocTree->insertPreProcess(sidx2);
                    }
                } else // 周围前沿无front
                {
                    if (sidx1 > 0 && sidx2 > 0) {
                        m_ocTree->rmDataPreProcess(sidx1);
                        m_ocTree->rmDataPreProcess(sidx2);

                        blFrts[i]->SetSTriIdx(nidx1, 0, -1);
                        blFrts[i]->SetSTriIdx(nidx1, 1, -1);
                    }
                }
            }
            // end of intersection
#endif

            m_blNxtFList->DeleteFront(bluFt);
            m_blFrontList->DeleteFront(bluFt);
            bluFt->SetLowerFront(nullptr);
            bluFt->clearTopo();
            delete bluFt;

            m_vBdyFront.push_back(blFrts[i]);

            // 			if(blFrts[i]->IsSymm())
            // 				AddSymmSegment(blFrts[i]);
            //
            // 			if(bluFt->IsSymm())
            // 			{
            // 				//need to remove segment
            // 			}
            for (int k = 0; k < 3; k++) {
                if (!used[k]) {
                    m_ocTree->rmDataPreProcess(blFrts[i]->GetSTriIdx(k, 0));
                    m_ocTree->rmDataPreProcess(blFrts[i]->GetSTriIdx(k, 1));
                }
            }
        }
    }

    for (i = 0; i < nblFrts; i++) {
        scheck.erase(blFrts[i]);
    }

    for (auto &i : nodes) {
        if (i->GetLowerNode()) {
            i->GetLowerNode()->SetUpperNode(nullptr);
            i->GetLowerNode()->RmvUpperNod(i);
        }

        delete i;
        i = nullptr;
    }
    if (blNod->getPerNode() && blNod->getPerNode()->GetUpperNode()) {
        RmvUperNeigFrontsAndFreeNode(blNod->getPerNode());
    }
}
void BLMesh::RmvUperNeigFronts(BLNode *blNod, bool bnext)
{
    int nblFrts, i, j;

    BLFront *blFrts[MAX_NCONN * 2], *bluFt;
    ;

    blNod->GetNeigFronts(blFrts, &nblFrts);

    for (i = 0; i < nblFrts; i++) {
        int nNeig, nNod;
        BLFront *neiFts[MAX_NCONN];

        blFrts[i]->GetNeigbourFronts(&nNeig, neiFts);
        int min_level = 2;
        for (int j = 0; j < nNeig; j++) {
            min_level = min(min_level, neiFts[j]->is_prism_valid);
        }
        blFrts[i]->is_prism_valid = min(blFrts[i]->is_prism_valid, min_level - 1);
        for (int j = 0; j < nNeig; j++) {
            scheck.insert(neiFts[j]);
        }
    }
    for (i = 0; i < nblFrts; i++) {
        scheck.erase(blFrts[i]);
    }
    for (i = 0; i < nblFrts; i++) {
        bluFt = blFrts[i]->GetUpperFront();

        // delete bluFt
        if (bluFt != nullptr) {
            // 删除周围
            int nNeig, nNod;
            BLFront *neiFts[3];
            BLNode *blNods[3];
            // remove neigbor information
            bluFt->GetNeigbourFronts(&nNeig, neiFts);
            for (j = 0; j < nNeig; j++) {
                neiFts[j]->RmvNeigbourFronts(bluFt);
            }

            // remove node information
            bluFt->GetNodes(&nNod, blNods);
            for (j = 0; j < nNod; j++) {
                blNods[j]->RmvNeigFront(bluFt);
            }

            blFrts[i]->SetUpperFront(nullptr);
            SetElmDelete(blFrts[i]->GetElmIdx());
            m_nPrism--;

#ifdef _CHECK_INTERSECTION
            // intersection

            if (bluFt->IsPyramid()) {
                bluFt->SetPyramidFlag(false);
            }

            m_ocTree->rmDataPreProcess(bluFt->GetTriIdx());
            if (blFrts[i]->GetTriIdx() < 0) {
                cout << "error";
            }
            m_ocTree->insertPreProcess(blFrts[i]->GetTriIdx());

            // if (m_ocTree->chckIntersectPreProcess(blFrts[i]->GetTriIdx())) {
            //	check_insert_again_ = true;
            // }

            blFrts[i]->GetNeigbourFronts(&nNeig, neiFts);
            for (j = 0; j < nNeig; j++) {
                int nidx1, nidx2, sidx1, sidx2;
                NeighIdx(blFrts[i], neiFts[j], &nidx1);
                NeighIdx(neiFts[j], blFrts[i], &nidx2);
                sidx1 = blFrts[i]->GetSTriIdx(nidx1, 0);
                sidx2 = blFrts[i]->GetSTriIdx(nidx1, 1);

                if (neiFts[j]->GetUpperFront()) {
                    ////插入到待处理的边界中
                    // auto it = scheck.begin();
                    // for (; it != scheck.end(); it++)
                    //	if (*it == neiFts[j]) {
                    //		break;
                    //	}
                    // if (it == scheck.end())
                    //	scheck.push_back(neiFts[j]);

                    if (sidx1 > 0 && sidx2 > 0) {
                        m_ocTree->insertPreProcess(sidx1);
                        m_ocTree->insertPreProcess(sidx2);

                        neiFts[j]->SetSTriIdx(nidx2, 0, sidx1);
                        neiFts[j]->SetSTriIdx(nidx2, 1, sidx2);

                        blFrts[i]->SetSTriIdx(nidx1, 0, -1);
                        blFrts[i]->SetSTriIdx(nidx1, 1, -1);
                    } else {
                        sidx1 = neiFts[j]->GetSTriIdx(nidx2, 0);
                        sidx2 = neiFts[j]->GetSTriIdx(nidx2, 1);

                        m_ocTree->insertPreProcess(sidx1);
                        m_ocTree->insertPreProcess(sidx2);
                    }
                    // cout << sidx1 << " " << sidx2 << "###" << endl;
                } else // 周围前沿无front
                {
                    if (sidx1 > 0 && sidx2 > 0) {
                        m_ocTree->rmDataPreProcess(sidx1);
                        m_ocTree->rmDataPreProcess(sidx2);

                        blFrts[i]->SetSTriIdx(nidx1, 0, -1);
                        blFrts[i]->SetSTriIdx(nidx1, 1, -1);
                    }
                }
            }
            // end of intersection
#endif

            m_blNxtFList->DeleteFront(bluFt);
            m_blFrontList->DeleteFront(bluFt);
            bluFt->SetLowerFront(nullptr);
            bluFt->clearTopo();
            delete bluFt;

            m_vBdyFront.push_back(blFrts[i]);

            // 			if(blFrts[i]->IsSymm())
            // 				AddSymmSegment(blFrts[i]);
            //
            // 			if(bluFt->IsSymm())
            // 			{
            // 				//need to remove segment
            // 			}
        }
    }

    if (blNod->getPerNode()) {
        if (blNod->getPerNode()->GetUpperNode()) {
            blNod->getPerNode()->GetUpperNode()->GetNeigFronts(blFrts, &nblFrts);
            if (nblFrts) {
                RmvUperNeigFronts(blNod->getPerNode());
            }
        }
    }
}

void BLMesh::UpdateSymnode()
{
    int i;

    for (i = 0; i < m_nNodes; i++) {
        if (m_pNodes[i].bsysm) {
            if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
                m_pNodes[i].coord[0] = m_sysValue;
            } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
                m_pNodes[i].coord[1] = m_sysValue;
            } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
                m_pNodes[i].coord[2] = m_sysValue;
            }

            // update boundary points coordinates
        }
    }
}

bool BLMesh::IsSymLine(int p1, int p2)
{
    return m_mapSymline[std::min(p1, p2)].find(std::max(p1, p2)) != m_mapSymline[std::min(p1, p2)].end();
}

void BLMesh::CalSymplnBdry(int *nBdry, int **pBdry, int ifc, bool add_symm)
{
    int i, j, pid1, pid2, pid3, pid4, conn[4], min, max, cnt;
    bool isedg = false;
    Elem *pSymBdry = nullptr;
    BLNode *pNod1, *pNod2, *pNodUp1, *pNodUp2;
    BLFront *blFront = nullptr, *blFrontUp = nullptr;

    std::set<std::array<int, 2>> sym_line;
    auto insert_sym = [&sym_line](int id1, int id2) {
        auto pair = std::array<int, 2>{std::min(id1, id2), std::max(id1, id2)};
        sym_line.insert(pair);
    };
    auto erase_sym = [&sym_line](int id1, int id2) {
        auto pair = std::array<int, 2>{std::min(id1, id2), std::max(id1, id2)};
        sym_line.erase(pair);
    };
    for (i = 0; i < m_nSymBdrys; i++) {
        pSymBdry = &m_pSymBdrys[i];
        int idx1 = pSymBdry->conn[0];
        int idx2 = pSymBdry->conn[1];

        if (std::find(m_pNodes[idx1].isymfc.begin(), m_pNodes[idx1].isymfc.end(), ifc) != m_pNodes[idx1].isymfc.end() ||
            std::find(m_pNodes[idx2].isymfc.begin(), m_pNodes[idx2].isymfc.end(), ifc) != m_pNodes[idx2].isymfc.end()) {
            continue;
        }

        if (pSymBdry->nconn >= 0) // wall boundaries
        {
            insert_sym(m_pSymBdrys[i].conn[0], m_pSymBdrys[i].conn[1]);
        }
    }

    int nCurSymBdry = m_nSymBdrys, nKeepSym = 0;
    int idx1, idx2;

    for (auto bdy_front : m_vBdyFront) {
        int nnode;
        BLNode *nodes[3];
        bdy_front->GetNodes(&nnode, nodes);
        for (int i = 0; i < 3; i++) {
            // if (m_pNodes[nodes[i]->GetNodIdx()].isymfc == ifc && m_pNodes[nodes[(i + 1) % 3]->GetNodIdx()].isymfc == ifc)
            //{
            int idx1 = nodes[i]->GetNodIdx();
            int idx2 = nodes[(i + 1) % 3]->GetNodIdx();

            if (std::find(m_pNodes[idx1].isymfc.begin(), m_pNodes[idx1].isymfc.end(), ifc) != m_pNodes[idx1].isymfc.end() ||
                std::find(m_pNodes[idx2].isymfc.begin(), m_pNodes[idx2].isymfc.end(), ifc) != m_pNodes[idx2].isymfc.end()) {

                pNod1 = nodes[i];
                pNod2 = nodes[(i + 1) % 3];
                erase_sym(pNod1->GetDecentID(), pNod2->GetDecentID());
                pNodUp1 = pNod1->GetUpperNode();
                pNodUp2 = pNod2->GetUpperNode();
                if (pNodUp1 != nullptr && pNodUp2 != nullptr) // so blFrontUp==nullptr
                {
                    pid1 = pNod1->GetNodIdx();
                    pid2 = pNod2->GetNodIdx();
                    pid3 = pNodUp1->GetNodIdx();
                    pid4 = pNodUp2->GetNodIdx();

                    // may be an symmetric line created when creating pyramids
                    bool isedg1 = IsSymLine(pid1, pid4);
                    bool isedg2 = IsSymLine(pid2, pid3);

                    if (isedg1 && !isedg2) {
                        conn[0] = pid1;
                        conn[1] = pid4;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);

                        conn[2] = pid2;
                        AddElem(3, conn, BLEntityTopology::TRIANGLE, ifc);

                        conn[0] = pid1;
                        conn[1] = pid3;
                        insert_sym(conn[0], conn[1]);
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                    } else if (!isedg1 && isedg2) {
                        conn[0] = pid2;
                        conn[1] = pid3;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);

                        conn[2] = pid1;
                        AddElem(3, conn, BLEntityTopology::TRIANGLE, ifc);

                        conn[0] = pid2;
                        conn[1] = pid4;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);
                    } else if (!isedg1 && !isedg2) {
                        conn[0] = pid1;
                        conn[1] = pid3;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);
                        conn[0] = pid1;
                        conn[1] = pid2;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);
                        conn[0] = pid2;
                        conn[1] = pid4;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);
                    } else {
                        // impossible
                        spdlog::info("Error: impossible case!\n");
                        throw(std::string("impossible case"));
                    }
                } else if (pNodUp1 == nullptr && pNodUp2 == nullptr) {
                    pid1 = pNod1->GetNodIdx();
                    pid2 = pNod2->GetNodIdx();
                    conn[0] = pid1;
                    conn[1] = pid2;
                    AddSymBdry(2, conn, BLEntityTopology::LINE);
                    insert_sym(conn[0], conn[1]);
                } else {
                    if (pNodUp1 == nullptr && pNodUp2 != nullptr) {
                        pid1 = pNod1->GetNodIdx();
                        pid2 = pNod2->GetNodIdx();
                        pid3 = pNodUp2->GetNodIdx();
                    } else if (pNodUp1 != nullptr && pNodUp2 == nullptr) {
                        pid1 = pNod2->GetNodIdx();
                        pid2 = pNod1->GetNodIdx();
                        pid3 = pNodUp1->GetNodIdx();
                    }

                    // may be an symmetric line created when creating pyramids
                    isedg = IsSymLine(pid1, pid3);

                    if (isedg) {
                        conn[0] = pid1;
                        conn[1] = pid3;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);

                        conn[2] = pid2;
                        AddElem(3, conn, BLEntityTopology::TRIANGLE, ifc);
                    } else { // I think it should not come to here too.Actually, it comes to here
                        conn[0] = pid1;
                        conn[1] = pid2;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);

                        conn[0] = pid2;
                        conn[1] = pid3;
                        AddSymBdry(2, conn, BLEntityTopology::LINE);
                        insert_sym(conn[0], conn[1]);
                    }
                }
            }
        }
    }

    for (auto i : m_pyramid_symline) {
        sym_line.insert(i);
    }

    cout << "Sym size=" << sym_line.size() << endl;

    /* remove some sym line if tow node in symm face inner the symm face */
    /*  1-------2----7---------
     *    \     /
     *	  \   /
     *       4
     *        \
     *          6
     *           \
     *
     *    remove (2,4)
     */
    if (add_symm) {
        removeNonManifoldPoint(sym_line);
    }

    *nBdry = sym_line.size();
    *pBdry = new int[2 * (*nBdry)];
    cnt = 0;
    for (auto i : sym_line) {

        (*pBdry)[cnt * 2 + 0] = i[0];
        (*pBdry)[cnt * 2 + 1] = i[1];
        cnt++;
    }
}

void BLMesh::RemvDeletedSymBdry()
{
    // remove boundary with deleted flag
    int cnt = 0, i;
    for (i = 0; i < m_nSymBdrys; i++) {
        if (m_pSymBdrys[i].nconn >= 0) {
            m_pSymBdrys[cnt].nconn = m_pSymBdrys[i].nconn;
            m_pSymBdrys[cnt].conn[0] = m_pSymBdrys[i].conn[0];
            m_pSymBdrys[cnt].conn[1] = m_pSymBdrys[i].conn[1];

            if (m_pSymBdrys[i].nconn < 0) {
                m_pSymBdrys[cnt].nconn = -abs(m_pSymBdrys[cnt].nconn);
            }

            cnt++;
        }
    }
    m_nSymBdrys = cnt;
}
void BLMesh::OutputSymplnBdry(string file, int nBdry, int *pBdry, int ifc)
{
    const char *filename = file.data();
    int i, idx;
    std::set<int> setpnt;
    std::set<int>::iterator sit;
    FILE *fout = nullptr;

    for (i = 0; i < nBdry; i++) {
        setpnt.insert(pBdry[2 * i + 0]);
        setpnt.insert(pBdry[2 * i + 1]);
    }
    map<int, int> ptmap;

    if (setpnt.size() != nBdry) {
        printf("Error:(%d, %d) %s %s %d\n", setpnt.size(), nBdry, __FILE__, __FUNCTION__, __LINE__);
        // exit(0);
    }

    fout = fopen(filename, "w");
    if (fout == nullptr) {
        throw(std::string("no authority to write file"));
    }
    // fprintf(fout, "%d %d 0 0 0 0 0\n", setpnt.size(), nBdry);
    fprintf(fout, "# vtk DataFile Version 2.0\nboundary layer mesh\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d double\n", setpnt.size());
    sit = setpnt.begin();
    i = 0;
    while (sit != setpnt.end()) {
        idx = *sit;
        if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
            fprintf(fout, "%f %f 0.0\n", m_pNodes[idx].coord[1], m_pNodes[idx].coord[2]);
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
            fprintf(fout, "%f %f 0.0\n", m_pNodes[idx].coord[0], m_pNodes[idx].coord[2]);
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
            fprintf(fout, "%f %f 0.0\n", m_pNodes[idx].coord[0], m_pNodes[idx].coord[1]);
        }
        ptmap[idx] = i++;
        ++sit;
    }
    fprintf(fout, "CELLS %d %d\n", nBdry, nBdry * 3);

    for (i = 0; i < nBdry; i++) {
        fprintf(fout, " 2 %d %d\n", ptmap[pBdry[i * 2 + 0]], ptmap[pBdry[i * 2 + 1]]);
    }
    fprintf(fout, "CELL_TYPES %d\n", nBdry);
    for (i = 0; i < nBdry; i++) {
        fprintf(fout, " 3\n");
    }

    fclose(fout);
    fout = nullptr;
}

/*
 * @
 */

void BLMesh::CalSymplnMsh(int nBdry, int *pBdry, int *nBdryElm, int **pBdryElm, double symval, int ifc)
{

    int i;
    int nbpt = 0, nbelm = 0, *belm = nullptr, *belmN = nullptr, *orientElm = nullptr;
    int npt = 0, nelm = 0, *elm = nullptr, *l_to_g = nullptr;
    double *bpt = nullptr, *pt = nullptr, *ptn = nullptr;

    std::set<int> spnt;
    std::set<int>::iterator sit;

    for (i = 0; i < nBdry; i++) {
        spnt.insert(pBdry[i * 2 + 0]);
        spnt.insert(pBdry[i * 2 + 1]);
    }
    nbpt = spnt.size();
    nbelm = nBdry;

    l_to_g = new int[nbpt];
    map<int, int> ptmap;
    bpt = new double[nbpt * 2];
    belm = new int[2 * nbelm];

    i = 0;
    sit = spnt.begin();
    while (sit != spnt.end()) {
        ptmap[*sit] = i;
        l_to_g[i] = *sit;
        if (*sit == 209636) {
            exit(-1);
        }
        ++i;
        ++sit;
    }

    for (i = 0; i < nbpt; i++) {
        if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
            bpt[i * 2 + 0] = m_pNodes[l_to_g[i]].coord[1];
            bpt[i * 2 + 1] = m_pNodes[l_to_g[i]].coord[2];
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
            bpt[i * 2 + 0] = m_pNodes[l_to_g[i]].coord[0];
            bpt[i * 2 + 1] = m_pNodes[l_to_g[i]].coord[2];
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
            bpt[i * 2 + 0] = m_pNodes[l_to_g[i]].coord[0];
            bpt[i * 2 + 1] = m_pNodes[l_to_g[i]].coord[1];
        }
    }
    for (i = 0; i < nbelm; i++) {
        belm[i * 2 + 0] = ptmap[pBdry[i * 2 + 0]];
        belm[i * 2 + 1] = ptmap[pBdry[i * 2 + 1]];
    }

    vector<std::array<int, 2>> loops;

    OrientLoop(nbpt, nbelm, bpt, belm, &belmN, &orientElm, loops);
#ifdef _DEBUG
    OutputFr2("testbdry.vtk", nbpt, nbelm, bpt, belmN);
#endif
//	OutputFr2("testbdry.vtk", nbpt, nbelm, bpt, belmN);
#ifdef USE_TRIANGLE

    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    Eigen::MatrixXd H;
    Eigen::MatrixXd H1;

    V.resize(nbpt, 2);
    E.resize(nbelm, 2);
    H.resize(loops.size(), 2);

    std::vector<double> point_size(nbpt, 0);
    for (int i = 0; i < nbpt; i++) {
        V(i, 0) = bpt[2 * i + 0];
        V(i, 1) = bpt[2 * i + 1];
    }

    for (int i = 0; i < nbelm; i++) {
        E(i, 0) = belmN[2 * i + 0];
        E(i, 1) = belmN[2 * i + 1];
        point_size[E(i, 0)] += (V.row(E(i, 0)) - V.row(E(i, 1))).norm() / 2;
        point_size[E(i, 1)] += (V.row(E(i, 0)) - V.row(E(i, 1))).norm() / 2;
    }

    bool eps = 1e-10;
    for (int i = 0; i < loops.size(); i++) {
        double l1 = (V(E(loops[i][0], 0), 0) - V(E(loops[i][0], 1), 0)) / 2;
        double l2 = (V(E(loops[i][0], 0), 1) - V(E(loops[i][0], 1), 1)) / 2;
        H(i, 0) = (V(E(loops[i][0], 0), 0) + V(E(loops[i][0], 1), 0)) / 2 - l2 * eps;
        H(i, 1) = (V(E(loops[i][0], 0), 1) + V(E(loops[i][0], 1), 1)) / 2 + l1 * eps;
    }

    //	cout << H;
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;

    igl::triangle::triangulate(V, E, H, "Yq32.5", V2, F2);
    Eigen::MatrixXd Vbe = V2;

    BLMESH::zju::Mesh mesh;
    // for (int i = 0; i < 10; i++)
    //	cout << V2(i, 0) << " ";
    // cout << endl;
    mesh.vertex = V2;
    mesh.topo = F2;
    mesh.vertex.conservativeResize(V2.rows(), 3);

    BLMESH::zju::SurfaceRemesh remesher(mesh);
    std::function<double(double, double, double)> sizingfunction = [point_size, V](double x, double y, double z) -> double {
        double ans = 1e9;
        for (int i = 0; i < V.rows(); i++) {
            double d = sqrt(pow(x - V(i, 0), 2) + pow(y - V(i, 1), 2));
            double size = 0.3 * d + point_size[i];
            ans = min(ans, size);
        }
        return ans;
    };
    std::function<std::array<double, 3>(double, double, double)> profunction = [](double x, double y, double z) {
        std::array<double, 3> pos;
        pos[2] = 0;
        pos[0] = x;
        pos[1] = y;
        return pos;
    };
    std::set<int> fixid;
    for (int i = 0; i < nbpt; i++) {
        fixid.insert(i);
    }
    remesher.set_pro_function(profunction);
    remesher.set_sizing_function(sizingfunction);
    remesher.set_fix_id(fixid);
    remesher.iso_remesh(25);
    V2 = mesh.vertex.leftCols(2);
    F2 = mesh.topo;
    Eigen::VectorXd C;
    Eigen::MatrixXi FF = F2;
    igl::bfs_orient(FF, F2, C);

    for (int i = 0; i < nbpt; i++) {
        if (V(i, 0) != V2(i, 0) || V(i, 1) != V2(i, 1)) {
            spdlog::info("the distance=" << (V.row(i) - V2.row(i)).norm());
            spdlog::info(V.row(i));
            spdlog::info(Vbe.row(i));
            spdlog::info(V2.row(i));
            V2.row(i) = V.row(i);
            // exit(-1);
        }
    }
    npt = V2.rows();
    pt = new double[2 * npt];
    for (int i = 0; i < npt; i++) {
        pt[2 * i + 0] = V2(i, 0);
        pt[2 * i + 1] = V2(i, 1);
    }
    nelm = F2.rows();
    elm = new int[3 * nelm];
    for (int i = 0; i < nelm; i++) {
        elm[3 * i + 0] = F2(i, 0);
        elm[3 * i + 1] = F2(i, 1);
        elm[3 * i + 2] = F2(i, 2);
    }
#else

        GenMesh2D(nbpt, nbelm, bpt, belmN, &npt, &pt, &nelm, &elm);
#endif

    cout << nelm << " " << npt << endl;

    ptn = new double[3 * npt];

    for (i = 0; i < npt; i++) {
        // #ifndef _NEW_SYMM
        if (m_sysPlane == SymmetryPlane::SYSMMETRY_X) {
            ptn[i * 3 + 0] = symval;
            ptn[i * 3 + 1] = pt[i * 2 + 0];
            ptn[i * 3 + 2] = pt[i * 2 + 1];
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Y) {
            ptn[i * 3 + 0] = pt[i * 2 + 0];
            ptn[i * 3 + 1] = symval;
            ptn[i * 3 + 2] = pt[i * 2 + 1];
        } else if (m_sysPlane == SymmetryPlane::SYSMMETRY_Z) {
            ptn[i * 3 + 0] = pt[i * 2 + 0];
            ptn[i * 3 + 1] = pt[i * 2 + 1];
            ptn[i * 3 + 2] = symval;
        }
        // #endif
    }

    // check orient
    BLVector centor(0, 0, 0);
    for (int i = 0; i < m_nNodes; i++) {
        centor += BLVector(m_pNodes[i].coord);
    }
    centor = centor / m_nNodes;
    double p[3];
    for (int k = 0; k < 3; k++) {
        p[k] = centor[k];
    }
    if (GEOM_FUNC::orient3d(ptn + 3 * elm[0], ptn + 3 * elm[1], ptn + 3 * elm[2], p) < 0) {

        for (int i = 0; i < nelm; i++) {
            swap(elm[3 * i + 0], elm[3 * i + 1]);
        }
        spdlog::info("symmetry plane swaped!");
    }

    //
    UpdateSymplnGrid(npt, nbpt, nelm, ptn, elm, &l_to_g, ifc);

#ifdef _DEBUG
    OutputPls("symmetry.pls", npt, nelm, ptn, elm);
#endif

    *nBdryElm = nelm;
    *pBdryElm = new int[nelm * 3];
    for (i = 0; i < nelm; i++) {
        (*pBdryElm)[i * 3 + 0] = l_to_g[elm[i * 3 + 0]];
        (*pBdryElm)[i * 3 + 1] = l_to_g[elm[i * 3 + 1]];
        (*pBdryElm)[i * 3 + 2] = l_to_g[elm[i * 3 + 2]];
    }

    // free
    if (belm) {
        delete[] belm;
    }
    if (belmN) {
        delete[] belmN;
    }
    if (orientElm) {
        delete[] orientElm;
    }
    if (elm) {
        delete[] elm;
    }
    if (l_to_g) {
        delete[] l_to_g;
    }
    if (bpt) {
        delete[] bpt;
    }
    if (pt) {
        delete[] pt;
    }
    if (ptn) {
        delete[] ptn;
    }
}

void DFS(std::map<int, std::set<int>> graph, vector<int> path, vector<vector<int>> &ans, vector<vector<int>> &failed)
{
    int current = path.back();
    if (graph.find(current) == graph.end()) {
        failed.push_back(path);
        return;
    }
    for (auto next : graph[current]) {
        auto graph_bak = graph;
        auto path_back = path;
        /* update graph */
        graph_bak[next].erase(current);
        graph_bak[current].erase(next);
        if (graph_bak[next].empty()) {
            graph_bak.erase(next);
        }
        if (graph_bak[current].empty()) {
            graph_bak.erase(current);
        }
        /* check is completed \& update path */
        if (next == path.front()) {
            ans.push_back(path);
            continue;
        } else {
            path_back.push_back(next);
        }
        /* find all path */
        DFS(graph_bak, path_back, ans, failed);
    }
}
/*
 * @date: 2021/11
 * @author： yhf
 * @brife: 该函数目的是修正dtiso2d的输入，不知为何输入里可能会存在多余的边。分为两步：
 * 1. DFS计算所有可能的环
 * 2. 对于有重叠的环，计算面积，取面积最大的环
 */
void BLMesh::removeNonManifoldPoint(std::set<std::array<int, 2>> &symlines)
{
    if (symlines.empty()) {
        return;
    }
    std::set<std::array<int, 2>> ans;
    std::map<int, std::set<int>> graph;
    for (auto i : symlines) {
        graph[i[0]].insert(i[1]);
        graph[i[1]].insert(i[0]);
    }
    while (!graph.empty()) {
        int start_point = (*graph.begin()).first;
        int oritation = m_pNodes[start_point].bsysm;
        int x = (oritation);
        int y = (oritation + 1) % 3;

        vector<int> path{start_point};
        vector<vector<int>> failed_path;
        vector<vector<int>> all_path;
        DFS(graph, path, all_path, failed_path);
        if (all_path.size() > 1) {
            cout << "Find " << all_path.size() << " path in loop" << endl;
        }
        double max_area = -1e20;
        vector<int> target;
        for (auto p : all_path) {
            double area = 0;
            for (int i = 0; i < p.size(); i++) {
                int node1 = p[(i - 1 + p.size()) % p.size()];
                int node2 = p[i];

                double *coord1 = m_pNodes[node1].coord;
                double *coord2 = m_pNodes[node2].coord;
                area += coord1[x] * coord2[y] - coord2[x] * coord1[y];
            }
            if (abs(area) > max_area) {
                max_area = abs(area);
                target = p;
            }
        }
        for (int i = 0; i < target.size(); i++) {
            int j = (target.size() + i - 1) % target.size();
            ans.insert(std::array<int, 2>{std::min(target[j], target[i]), std::max(target[j], target[i])});
        }
        for (auto i : all_path) {
            for (auto j : i) {
                graph.erase(j);
            }
        }
        for (auto i : failed_path) {
            for (auto j : i) {
                graph.erase(j);
            }
        }
    }

    symlines = ans;
}
std::function<double(std::array<double, 3>)> BLMesh::sizefuntion = nullptr;
void BLMesh::SetSizingFunction(std::function<double(std::array<double, 3>)> func) { sizefuntion = func; }

std::function<double(std::array<double, 3>)> &BLMesh::GetSizingFunction() { return sizefuntion; }

void BLMesh::GlobalOptmize()
{
    //	auto points=GetNarrowConstrainedPointidx();
    std::vector<BLVector> coordinate;
    std::vector<std::vector<int>> connector;

    for (int i = 0; i < m_nNodes; i++) {
        coordinate.push_back(BLVector(m_pNodes[i].coord));
    }
    for (int i = 0; i < m_nElems; i++) {
        vector<int> conn(m_pElems[i].nconn);
        if (m_pElems[i].topo == BLEntityTopology::PRISM) {
            conn[0] = m_pElems[i].conn[0];
            conn[1] = m_pElems[i].conn[2];
            conn[2] = m_pElems[i].conn[1];
            conn[3] = m_pElems[i].conn[3];
            conn[4] = m_pElems[i].conn[5];
            conn[5] = m_pElems[i].conn[4];
        } else if (m_pElems[i].topo == BLEntityTopology::PYRAMID) {
            conn[0] = m_pElems[i].conn[0];
            conn[1] = m_pElems[i].conn[2];
            conn[2] = m_pElems[i].conn[1];
            conn[3] = m_pElems[i].conn[3];
            conn[4] = m_pElems[i].conn[4];
        } else if (m_pElems[i].topo == BLEntityTopology::TETRAHEDRON) {
            conn[0] = m_pElems[i].conn[0];
            conn[1] = m_pElems[i].conn[1];
            conn[2] = m_pElems[i].conn[2];
            conn[3] = m_pElems[i].conn[3];
        }
        //	bool add = false;
        // if (points.find(conn[j]) != points.end())
        //	add = true;
        //	if (add)
        connector.push_back(conn);
    }

    int maxinteration = 10;
    //	while (MeshOptimize(coordinate, connector)) { maxinteration--; if (!maxinteration)break; }
    for (int i = 0; i < m_nNodes; i++) {
        m_pNodes[i].coord[0] = coordinate[i].x;
        m_pNodes[i].coord[1] = coordinate[i].y;
        m_pNodes[i].coord[2] = coordinate[i].z;
    }
}

std::vector<double> BLMesh::recommand_length_calculation()
{
    auto ans = std::vector<double>(m_nSurfNodes, std::numeric_limits<double>::max());
    std::set<int> setbdry;
    for (int i = 0; i < noutbdry; i++) {
        setbdry.insert(outbdry[3 * i + 0]);
        setbdry.insert(outbdry[3 * i + 1]);
        setbdry.insert(outbdry[3 * i + 2]);
    }
    for (auto i : setbdry) {
        double *coord = m_pNodes[i].coord;
        auto node = (BLNode *)m_pNodes[i].pointer;
        if (!node) {
            continue;
        }
        double *decent_coord = m_pNodes[node->GetDecentID()].coord;

        if (node->GetDecentID() == node->GetNodIdx()) {
            ans[node->GetNodIdx()] = -1;
            continue;
        }

        double length = (BLVector(decent_coord) - BLVector(coord)).magnitude();
        ans[node->GetDecentID()] = length;
    }
    return ans;
}
void BLMesh::resetZeroHeight(int NumNodes, BLNode **nodes)
{

    // Caculate initial boundary layer mesh height. The value is subject to MIN_LAYER,
    // Make sure that MIN_LAYER layer generate for every cell succeed if no intersection happen.

    if (cf.step_len * pow(cf.ratio1, MIN_LAYER) > ave_front_size) {
        cf.step_len = ave_front_size / pow(cf.ratio1, MIN_LAYER);
        for (int i = 0; i < NumNodes; i++) {
            BLNode *node_ptr = nodes[i];
            node_ptr->m_h0 = cf.step_len;
        }
    }

    for (int i = 0; i < NumNodes; i++) {
        BLNode *node_ptr = nodes[i];
        double node_size = 0;
        if (node_ptr->GetNeigFronts().empty()) {
            continue;
        }
        for (auto j : node_ptr->GetNeigFronts()) {
            node_size += j->GetFrontSize();
        }
        node_size /= node_ptr->GetNeigFronts().size();

        if (node_size / pow(cf.ratio1, MIN_LAYER) < node_ptr->m_h0) {
            node_ptr->m_h0 = node_size / pow(cf.ratio1, MIN_LAYER);
        }
    }

    // smooth m_h0 to avoid abrut change in neighbour front
    queue<BLNode *> Q;
    for (int i = 0; i < NumNodes; i++) {
        Q.push(nodes[i]);
    }

    while (!Q.empty()) {
        auto node = Q.front();
        Q.pop();
        ;
        if (node->GetNeigNods().empty()) {
            continue;
        }
        double min_height = std::numeric_limits<double>::max();
        double ave_height = 0;
        for (auto j : node->GetNeigNods()) {
            min_height = std::min(min_height, j->m_h0);
            ave_height += j->m_h0;
        }
        ave_height /= node->GetNeigNods().size();
        if (node->m_h0 > min_height * cf.ratio1 * cf.ratio2) {
            node->m_h0 = min_height * cf.ratio1 * cf.ratio2;
            for (auto &j : node->GetNeigNods()) {
                Q.push(j);
            }
        }
    }
}
