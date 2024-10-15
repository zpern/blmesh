#include <spdlog/spdlog.h> 
 #include "IterationNormalSmoothStrategy.h"



IterationNormalSmoothStrategy::IterationNormalSmoothStrategy(BLNode **node, MBLNode* pnodes,int num_front) :NormalSmoothStrategy(node, pnodes, num_front)
{
}

void IterationNormalSmoothStrategy::SmoothNormal()
{
	/*
	BLNode* blNods[MAX_NCONN*2]; int nBLNod; bool flag = false;
	blNod->GetNeigNods(blNods,&nBLNod,flag);
	vector<double> weight(nBLNod, 1.0 / nBLNod);
	BLVector N_p(0, 0, 0);
	for (int i = 0; i < m_vecNeighFronts.size(); i++) {
		N_p += weight[i] * m_vecNeighFronts[i]->GetNormal();
	}
	N_p.normalize();
	int iter_count = 0;
	while (true) {
		if (iter_count > 200)
			break;
		iter_count++;
		vector<double> alpha(m_vecNeighFronts.size());
		double alpha_sum = 0.0, weight_sum = 0.0;
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			alpha[i] = acos(m_vecNeighFronts[i]->GetNormal()*N_p);
			alpha_sum += alpha[i];
		}
		if (alpha_sum < 1e-4 || isnan(alpha_sum))
			break;
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			weight[i] = weight[i] * alpha[i] / alpha_sum;
			weight_sum += weight[i];
		}
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			weight[i] /= weight_sum;
		}
		BLVector N_pnew(0, 0, 0);
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			N_pnew += weight[i] * m_vecNeighFronts[i]->GetNormal();
		}
		N_pnew.normalize();
		if ((N_p - N_pnew).magnitude() < 1e-4)
			break;
		double beita = 0.5;
		N_p = beita * N_p + (1.0 - beita)*N_pnew;
		N_p.normalize();
	}
	vec = N_p;
	if (bsys)
	{
		//printf("idx: %d, sym: %d\n", m_iNod, pNodes[m_iNod].symaxis);
		if (m_iSymAxis == 0)
			vec.x = 0;
		else if (m_iSymAxis == 1)
			vec.y = 0;
		else if (m_iSymAxis == 2)
			vec.z = 0;
		//vec.y = 0;
	}


	m_Normal = vec;
	beitaVisu = GetBeitaVisu(vec);
	return vec;*/
}


IterationNormalSmoothStrategy::~IterationNormalSmoothStrategy()
{
}
