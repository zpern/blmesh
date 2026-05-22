#include "NormalSmoothStrategy.h"
#ifndef __APPLE__
#include <omp.h>
#include <spdlog/spdlog.h>
#endif
#include <queue>
NormalSmoothStrategy::NormalSmoothStrategy()
    : smooth_times(10)
    , node_array(nullptr)
{}

NormalSmoothStrategy::NormalSmoothStrategy(BLNode **node, MBLNode *pnodes, int num_front)
    : node_array(node)
    , pNodes(pnodes)
    , nFrtNods(num_front)
{}

void NormalSmoothStrategy::SmoothNormal()
{
    int i;
    tmp.clear();
    tmp.reserve(nFrtNods);

#ifdef USE_OPENMP
#pragma omp parallel for if (nFrtNods > 1000) schedule(guided, 10)
#endif

    for (int i = 0; i < nFrtNods; i++) {
        int a = node_array[i]->GetBeitaVisu(node_array[i]->GetNormal());
        if (node_array[i]->GetBeitaVisu(node_array[i]->GetNormal()) <= 1) {
            tmp.push_back(node_array[i]);
        }
    }

    vector<int> psize;
    for (int k = 0; k < 100; k++) {
        myvec.resize(nFrtNods, BLVector(-2, -1, -1));
        vector<BLNode *> next_tmp;

        int size = tmp.size();
        if (!size) {
            break;
        }
        psize.push_back(size);
        if (psize.size() > 12 && psize[psize.size() - 10] == size) {
            break;
        }

        for (i = 0; i < size; i++) {

            SmoothNormalOnce(tmp[i], i);
        }

        for (i = 0; i < size; i++) {
            if (isnan(myvec[i].x) || (myvec[i].x == -2)) {
                continue;
            }

            if (tmp[i]->GetBSys()) {
                int iNod = tmp[i]->GetNodIdx();
                auto ans = myvec[i];
                int decent_id = tmp[i]->GetDecentID();
                // if (decent_id = 6555) {
                //	spdlog::info("debug here");
                // }
                int nodeid = iNod;
                Eigen::RowVector3d start_point(pNodes[nodeid].coord[0], pNodes[nodeid].coord[1], pNodes[nodeid].coord[2]);
                Eigen::RowVector3d normal(ans[0], ans[1], ans[2]);
                std::vector<int> faceid = pNodes[nodeid].isymfc;
                if (faceid.size() == 1) {
                    faceid_to_sp[faceid[0]].adjustNormal(start_point, normal);
                } else {

                    for (int j = 0; j < 10; j++) {
                        faceid_to_sp[faceid[0]].adjustNormal(start_point, normal);
                        faceid_to_sp[faceid[1]].adjustNormal(start_point, normal);
                    }
                }
                BLVector n(normal(0), normal(1), normal(2));
                tmp[i]->SetNormal(n.normalized());
            } else {
                tmp[i]->SetNormal(myvec[i]);
            }
        }
        if (k < 3) {
            for (auto i : tmp) {
                next_tmp.push_back(i);
            }
        } else {
            for (auto i : tmp) {
                if (i->if_need_smooth) {
                    next_tmp.push_back(i);
                }
            }
        }

        tmp.swap(next_tmp);
    }

#ifdef SMOOTH_HOR_NORMAL

    for (int i = 0; i < nFrtNods; i++) {
        horsmooth(node_array[i]);
    }
#endif
}

void NormalSmoothStrategy::horsmooth(BLNode *blNod)
{
    if (!blNod->GetLowerNode()) {
        return;
    }
    double diff = blNod->GetNormal() * blNod->GetLowerNode()->GetNormal() - 1;
    blNod->SetNormal((blNod->GetNormal() * pow(2.0, diff) + blNod->GetLowerNode()->GetNormal() * (1 - pow(2.0, diff))));
}

void NormalSmoothStrategy::GetEnegy()
{
    std::vector<BLNode *> blNodes;
    int nNeigFrts, nNeigNods;
    energy = 0;
    double maxc = -1e10;
    double minc = 1e10;
    int i;
#ifndef _DEBUG

#endif
    for (i = 0; i < nFrtNods; i++) {
        BLNode *blNod = node_array[i];
        blNodes = blNod->GetNeigNods();
        nNeigNods = blNodes.size();
        for (int j = 0; j < nNeigNods; j++) {
            BLVector b = blNodes[j]->GetNormal();
            double c = b * blNod->GetNormal();
            int id = blNodes[j]->GetNodIdx();
            double distance = BLVector(pNodes[id].coord[0] - pNodes[blNod->GetNodIdx()].coord[0],
                                       pNodes[id].coord[1] - pNodes[blNod->GetNodIdx()].coord[1],
                                       pNodes[id].coord[2] - pNodes[blNod->GetNodIdx()].coord[2])
                                  .magnitude();
            double cost = acos(c); // / distance;
            maxc = max(cost, maxc);
            minc = min(cost, minc);
            energy += (unsigned long long)cost * cost;
        }
    }
    cout << " energy= " << energy << " max energy=" << maxc << " min cos=" << minc << endl;
}

NormalSmoothStrategy::~NormalSmoothStrategy() {}

void NormalSmoothStrategy::SetBlNode(BLNode *node)
{
    // blNod = node;
}

void NormalSmoothStrategy::SetSmoothTimes(int s_times) { smooth_times = s_times; }

void NormalSmoothStrategy::SetNumberofFrontNode(int num) { nFrtNods = num; }

void NormalSmoothStrategy::SetSymm(BLNode *blNod, int id)
{

    if (blNod->GetBSys()) {
        auto ans = blNod->GetNormal();
        int nodeid = id;
        Eigen::RowVector3d start_point(pNodes[nodeid].coord[0], pNodes[nodeid].coord[1], pNodes[nodeid].coord[2]);
        Eigen::RowVector3d normal(ans[0], ans[1], ans[2]);
        std::vector<int> faceid = pNodes[nodeid].isymfc;

        if (faceid.size() == 1) {
            faceid_to_sp[faceid[0]].adjustNormal(start_point, normal);
        } else {

            for (int j = 0; j < 1; j++) {
                faceid_to_sp[faceid[0]].adjustNormal(start_point, normal, true);
                faceid_to_sp[faceid[1]].adjustNormal(start_point, normal, true);
            }
        }

        BLVector n(normal(0), normal(1), normal(2));
        blNod->SetNormal(n.normalized());
    }
}

void NormalSmoothStrategy::SetFaceidSP(const std::map<int, TiGER::SymmetryPlane> &sp) { faceid_to_sp = sp; }
