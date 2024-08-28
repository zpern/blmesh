#include "BLNode.h"
#include "BLMesh.h"
#include "BLMesh_define.h"
#include "LightMemoryPool.h"
#include <spdlog/spdlog.h> 
 #include <set>

//#define X_PLANE_SYMMETRY
#define Y_PLANE_SYMMETRY
//#define Z_PLANE_SYMMETRY
//a callback function of threads
#ifdef  MEMORY_DEBUG
std::set<BLNode*> BLNode::my_vec = std::set<BLNode*>();
#endif
#ifdef USE_MEMORY_POOL
INIT(BLNode)
#endif

unsigned  fieldcal(void *para)
{
	FieldArg *fieldArg = (FieldArg *)para;
	return 0;
}

const BLVector& BLNode::GetNormal() const
{
	return m_Normal;
}

void BLNode::SetNormal(BLVector vec)
{
	m_Normal.x = vec.x;
	m_Normal.y = vec.y;
	m_Normal.z = vec.z;

	// 	if (m_blType == BLMType::blm2d)
	// 	{
	// 		m_uValue = sqrtf(vec.x*vec.x + vec.y*vec.y);
	// 	}
	// 	else
	// 	{
	// 		m_uValue = sqrtf(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
	// 	}
}
void BLNode::GetRealNeigNods(BLNode** blNods, int *nBLNod, bool flag) {
	int i, j, nSize, nNods;
	BLNode* blNodes[MAX_FRONT_NODES];
	std::set<BLNode*> setNods;
	std::set<BLNode*>::iterator sit;

	nSize = m_vecNeighFronts.size();
	for (i = 0; i < nSize; i++)
	{
		BLVector vecN;
		BLFront *blFront = m_vecNeighFronts[i];
		blFront->GetNodes(&nNods, blNodes);

		for (j = 0; j < nNods; j++)
		{
			//if(flag && !blNodes[j]->GetBSys())
			if (flag && bsys && !blNodes[j]->GetBSys())
				continue;
			else if (flag && !blNodes[j]->GetBdryPt())
				continue;

			if (blNodes[j] != this)
				setNods.insert(blNodes[j]);
		}
	}

	*nBLNod = setNods.size();
	sit = setNods.begin();
	i = 0;
	while (sit != setNods.end())
	{
		blNods[i++] = *sit;
		++sit;
	}
}
void BLNode::GetVirtualNeigNods(BLNode ** blNods, int * nBLNod, Node* pnode, bool flag)
{
	if (!GetVirtualFlag()) {
		*nBLNod = 0;
		return;
	}
	int i, j, nSize, nNods;
	BLNode* blNodes[MAX_FRONT_NODES];
	std::set<BLNode*> setNods;
	std::set<BLNode*>::iterator sit;

	nSize = m_vecNeighFronts.size();
	for (i = 0; i < nSize; i++)
	{
		BLVector vecN;
		BLFront *blFront = m_vecNeighFronts[i];
		blFront->GetNodes(&nNods, blNodes);

		for (j = 0; j < nNods; j++)
		{
			//if(flag && !blNodes[j]->GetBSys())
			if (flag && bsys && !blNodes[j]->GetBSys())
				continue;
			else if (flag && !blNodes[j]->GetBdryPt())
				continue;

			if (blNodes[j] != this) {
				BLVector v1(pnode[blNodes[j]->GetNodIdx()].coord[0], pnode[blNodes[j]->GetNodIdx()].coord[1], pnode[blNodes[j]->GetNodIdx()].coord[2]);
				BLVector v2(pnode[GetNodIdx()].coord[0], pnode[GetNodIdx()].coord[1], pnode[GetNodIdx()].coord[2]);
				if ((v1 - v2).magnitude() < 1e-7) {
					setNods.insert(blNodes[j]);
				}
			}
		}
	}

	*nBLNod = setNods.size();
	sit = setNods.begin();
	i = 0;
	while (sit != setNods.end())
	{

		blNods[i++] = *sit;
		++sit;
	}
}

BLNode::~BLNode(void)
{
	for (auto &i : m_vecNeighFronts)
		i = nullptr;
	for (auto &i : m_vecNeighNodes)
		i = nullptr;
#ifdef  MEMORY_DEBUG
	my_vec.erase(this);
#endif
	//m_vecNeighFronts.swap(std::vector <BLFront*>());
	//m_vecNeighNodes.swap(std::vector <BLNode*>());
}

void BLNode::AddDeleteCount()
{
	delete_count++;
	if (delete_count == m_vecNeighFronts.size()) {
		//if(!m_bBdry)
		delete this;
	}
}

void BLNode::CalConcaveConvex(Node * pNode)
{

}

BLVector BLNode::GetCoord(Node * pNode)
{
	return BLVector(pNode[GetNodIdx()].coord[0], pNode[GetNodIdx()].coord[1], pNode[GetNodIdx()].coord[2]);
}

vector<BLNode*>& BLNode::GetNeigNods()
{
	if (m_vecNeighNodes.empty()) {
		int i, j, nSize, count = 0;
		nSize = m_vecNeighFronts.size();
		m_vecNeighNodes.reserve(nSize);
		for (i = 0; i < nSize; i++)
		{
			for (j = 0; j < 3; j++)
			{
				if (m_vecNeighFronts[i]->m_pBLNods[j] == this) {
					m_vecNeighNodes.push_back(m_vecNeighFronts[i]->m_pBLNods[(j + 1) % 3]);
					break;
				}

			}
		}
	}

	return m_vecNeighNodes;
}


BLVector BLNode::GetNormal(Node* pNodes)
{
	BLVector vec(0.0, 0.0, 0.0);
	int i, nSize, nNods, iNod1, iNod2, iNod3;
	BLNode* blNodes[MAX_FRONT_NODES];
	double ratio;

	nSize = m_vecNeighFronts.size();

	for (i = 0; i < nSize; i++)
	{
		BLVector vecN;
		BLFront* blFront = m_vecNeighFronts[i];

		blFront->GetNodes(&nNods, blNodes);


		if (nNods == 3)
		{
			BLVector vectmp1, vectmp2;

			iNod1 = blNodes[0]->GetNodIdx();
			iNod2 = blNodes[1]->GetNodIdx();
			iNod3 = blNodes[2]->GetNodIdx();

			vectmp1.x = pNodes[iNod2].coord[0] - pNodes[iNod1].coord[0];
			vectmp1.y = pNodes[iNod2].coord[1] - pNodes[iNod1].coord[1];
			vectmp1.z = pNodes[iNod2].coord[2] - pNodes[iNod1].coord[2];

			vectmp2.x = pNodes[iNod3].coord[0] - pNodes[iNod1].coord[0];
			vectmp2.y = pNodes[iNod3].coord[1] - pNodes[iNod1].coord[1];
			vectmp2.z = pNodes[iNod3].coord[2] - pNodes[iNod1].coord[2];

			//inward normal
			//vecN = vectmp1^vectmp2;

			//outward normal
			vecN = vectmp2 ^ vectmp1;
		}
		else
		{
			;
		}
		vecN.normalize();
		blFront->SetNormal(vecN);

		vec.x += vecN.x;
		vec.y += vecN.y;
		vec.z += vecN.z;
	}

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

	vec.normalize();

	m_Normal.x = vec.x;
	m_Normal.y = vec.y;
	m_Normal.z = vec.z;

	return vec;
}
double BLNode::GetBeitaVisu(BLVector normal) {

	double minBeita = 1.0;
	for (int i = 0; i < m_vecNeighFronts.size(); i++)
	{
		BLFront* blFront = m_vecNeighFronts.at(i);
		minBeita = std::min(minBeita, blFront->GetNormal() * normal);
	}
	return minBeita;
}
double BLNode::GetBeitaVisu()
{
	return beitaVisu;
}
void BLNode::SetBeitaVisu(double val)
{
	beitaVisu = val;
}
///*直线与角（对角）求交，用于WF算法，使用克莱姆法则和行列式的叉乘求法*///


void lineAngleInter(BLVector A, BLVector C, BLVector I, BLVector B, BLVector P, BLVector &pos, bool &isInter, bool &rightDir, double &t) {
	I = -1 * I;
	double dm = (A * (C ^ I));
	if (abs(dm) < 1e-7) {
		isInter = false;//平行
		return;
	}
	double m = (B * (C ^I)) / dm;
	double n = (A * (B ^ I)) / dm;
	//cout << m << " " << n << endl;
	t = A * (C ^ B) / dm;
	if (m*n >= 0 && abs(t) < 50) {
		rightDir = (m > 0 || n > 0);
		isInter = true;
		pos = P + m * A + n * C;
		//cout << (m * A + n * C + t * I-B).magnitude()<<" ";
		return;
	}
	isInter = false;
	return;


}
BLVector BLNode::GetTheMostnNormalNormal(Node* pNodes) {
	assert(m_vecNeighFronts.size());
	vector<double> weight(m_vecNeighFronts.size(), 1.0 / m_vecNeighFronts.size());
	BLVector N_p(0, 0, 0);
	for (int i = 0; i < m_vecNeighFronts.size(); i++) {
		N_p += weight[i] * m_vecNeighFronts[i]->GetNormal();
	}
	N_p.normalize();
	int iter_count = 0;
	while (true) {
		if (iter_count > 100)
			break;
		iter_count++;
		vector<double> alpha(m_vecNeighFronts.size());
		double alpha_sum = 0.0, weight_sum = 0.0;
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			alpha[i] = acos(m_vecNeighFronts[i]->GetNormal()*N_p);
			alpha_sum += alpha[i];
		}
		if (alpha_sum < 1e-10 || std::isnan(alpha_sum))
			break;
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			weight[i] = weight[i] * alpha[i] / alpha_sum;
			weight_sum += weight[i];
		}
		/*normalize*/
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			weight[i] /= weight_sum;
		}
		BLVector N_pnew(0, 0, 0);
		for (int i = 0; i < m_vecNeighFronts.size(); i++) {
			N_pnew += weight[i] * m_vecNeighFronts[i]->GetNormal();
		}
		N_pnew.normalize();
		if ((N_p - N_pnew).magnitude() < 1e-10)
			break;
		double beita = ((iter_count + 50) / (iter_count + 100));
		N_p = beita * N_p + (1.0 - beita)*N_pnew;
		N_p.normalize();
	}
	return N_p;
}
BLVector BLNode::GetCenterNormal(Node* pNodes) {
	//1. find two face whose normal dot product is closest to -1
	int n0 = -1, n1 = -1;
	double min = 2;
	BLVector naive_normal = GetSimpleNormal(pNodes), vec;
	if (m_vecNeighFronts.size() < 2) {
		return naive_normal;
	}
	for (auto i = 0; i < m_vecNeighFronts.size(); i++) {
		for (auto j = i + 1; j < m_vecNeighFronts.size(); j++) {
			if (min > m_vecNeighFronts[i]->GetNormal() * m_vecNeighFronts[j]->GetNormal()) {
				n0 = i;
				n1 = j;
				min = m_vecNeighFronts.at(i)->GetNormal() * m_vecNeighFronts.at(j)->GetNormal();
			}
		}
	}
	BLVector I = (m_vecNeighFronts[n0]->GetNormal() ^ m_vecNeighFronts[n1]->GetNormal());
	I.normalize();
	BLVector B = (m_vecNeighFronts[n0]->GetNormal() + m_vecNeighFronts[n1]->GetNormal());
	B.normalize();

	std::vector<BLVector> interArray;
	std::vector<BLVector> interArrayTemp;
	double tmax = -1e10;
	double tmin = 1e10;
	BLVector pmax, pmin;
	for (int i = 0; i < m_vecNeighFronts.size(); i++)
	{
		BLVector blNorm;
		BLFront* blFront = m_vecNeighFronts[i];
		int numNod;
		BLNode* node[6];
		std::vector<BLVector> side;
		blFront->GetNodes(&numNod, node);
		int intTyp, intCod;
		double ans[3];
		for (auto k = 0; k < numNod; k++) {
			if (node[k]->GetNodIdx() != m_iNod) {
				BLVector bv(pNodes[node[k]->GetNodIdx()].coord[0] - pNodes[m_iNod].coord[0], pNodes[node[k]->GetNodIdx()].coord[1] - pNodes[m_iNod].coord[1], pNodes[node[k]->GetNodIdx()].coord[2] - pNodes[m_iNod].coord[2]);
				bv.normalize();
				side.push_back(bv);
			}
		}
		if (side.size() != 2) {
			cout << "front id=" << blFront->GetTriIdx() << endl;
			cout << "node id=" << m_iNod << endl;
			cout << "error neighbor vector" << endl;
		}

		BLVector pos;
		double t;
		bool isInter, rightDir;
		BLVector oldPos(pNodes[m_iNod].coord[0], pNodes[m_iNod].coord[1], pNodes[m_iNod].coord[2]);
		lineAngleInter(side.at(0), side.at(1), I, B, oldPos, pos, isInter, rightDir, t);
		if (isInter) {
			if (t > tmax && t > -0.10) {
				tmax = t;
				pmax = pos;
			}
			if (t < tmin && t < 0.10) {
				tmin = t;
				pmin = pos;
			}
		}
		if (isInter) {
			interArrayTemp.push_back(pos);
		}
	}
	if (tmax > -1e9) {
		interArray.push_back(pmax);
	}
	if (tmin < 1e9) {
		interArray.push_back(pmin);
	}
	//cout << interArray.size();
	if (interArray.size() == 0) {
		vec = B;
	}
	if (interArray.size() == 1) {
		BLVector C = interArray[0] - BLVector(pNodes[m_iNod].coord[0], pNodes[m_iNod].coord[1], pNodes[m_iNod].coord[2]);

		C.normalize();
		BLVector v1 = C + I;
		BLVector v2 = C - I;
		v1.normalize();
		v2.normalize();
		if ((v1)* naive_normal > (v2)* naive_normal)
			vec = v1;
		else
			vec = v2;
		//vec = I;
		vec.normalize();
	}
	if (interArray.size() == 2) {
		BLVector C = interArray[0] - BLVector(pNodes[m_iNod].coord[0], pNodes[m_iNod].coord[1], pNodes[m_iNod].coord[2]);
		BLVector D = interArray[1] - BLVector(pNodes[m_iNod].coord[0], pNodes[m_iNod].coord[1], pNodes[m_iNod].coord[2]);
		C.normalize();
		D.normalize();
		vec = C + D;
		vec.normalize();
	}
	return vec;
}
BLVector BLNode::GetSimpleNormal(Node* pNodes) {
	BLVector simple_normal;

	for (auto i : m_vecNeighFronts) {
		simple_normal += i->GetNormal();
	}
	simple_normal.normalize();
	return simple_normal;
}
BLVector BLNode::GetCirculeDotNormal(Node* pNodes) {

	int num_neigh_front = m_vecNeighFronts.size();
	if (num_neigh_front < 3 && num_neigh_front>11)
		return GetNaiveNormal(pNodes);
	for (int i = 0; i < num_neigh_front; i++) {
		for (int j = 0; j < num_neigh_front; j++) {
			for (int k = 0; k < num_neigh_front; k++) {

			}
		}
	}
}
BLVector BLNode::GetGeometryNormal(Node* pNodes) {
	BLVector ans(0, 1, 0);
	if (m_vecNeighFronts.size() == 1) {
		return m_vecNeighFronts[0]->GetNormal();
	}
	int size = m_vecNeighFronts.size();
	double min_bei_vise = -10;
	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			for (int k = j + 1; k < size; k++) {
				BLVector centor = solveCenterPointOfCircle(m_vecNeighFronts[i]->GetNormal(), m_vecNeighFronts[j]->GetNormal(), m_vecNeighFronts[k]->GetNormal());
				BLVector normal = centor.normalized();
				if (std::isnan(normal.x))
					continue;
				double d = GetBeitaVisu(normal);
				if (min_bei_vise < d) {
					min_bei_vise = d;
					ans = normal;
				}
			}
		}
	}
	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			BLVector centor = (m_vecNeighFronts[i]->GetNormal() + m_vecNeighFronts[j]->GetNormal());

			BLVector normal = centor.normalized();
			if (std::isnan(normal.x))
				continue;
			double d = GetBeitaVisu(normal);
			if (min_bei_vise < d) {
				min_bei_vise = d;
				ans = normal;
			}
		}
	}

	return ans;

}
BLVector BLNode::GetNaiveNormal(Node* pNodes) {
	BLVector naiveNormal;
	std::vector< std::vector<BLVector> > normgroup;
	int neigh_front_size = m_vecNeighFronts.size();
	std::vector<BLVector> norms(neigh_front_size);
	std::vector<bool> flag(neigh_front_size);
	double angle_eps = rad(25);
	int cnt = 0;
	for (int i = 0; i < m_vecNeighFronts.size(); i++)
	{
		BLVector blNorm;
		BLFront *blFront = m_vecNeighFronts[i];
		blNorm = blFront->GetNormal();

		norms[cnt] = blNorm;
		flag[cnt] = 0;
		cnt++;
	}

	//classify norms
	for (int i = 0; i < m_vecNeighFronts.size(); i++)
	{
		for (int k = 0; k < normgroup.size(); k++)
		{
			std::vector<BLVector> vtmp = normgroup[k];
			BLVector ntmp = vtmp[0];

			if (ntmp*norms[i] > cos(angle_eps))
			{
				normgroup[k].push_back(norms[i]);
				flag[i] = 1;
			}
		}

		if (!flag[i])
		{
			std::vector<BLVector> vtmp;
			vtmp.push_back(norms[i]);
			normgroup.push_back(vtmp);
			flag[i] = 1;
		}
	}

	m_nNormContri = normgroup.size();

	//update norms
	BLVector fnorm(0, 0, 0);
	for (int i = 0; i < m_nNormContri; i++)
	{
		BLVector ntmp(0, 0, 0);
		std::vector<BLVector> vtmp = normgroup[i];
		for (int k = 0; k < vtmp.size(); k++)
		{
			ntmp += vtmp[k];
		}
		ntmp.x /= vtmp.size();
		ntmp.y /= vtmp.size();
		ntmp.z /= vtmp.size();

		ntmp.normalize();

		naiveNormal += ntmp;
	}

	naiveNormal.normalize();
	if (naiveNormal.magnitude() < 0.9) {
		if (m_nNormContri)
			naiveNormal = normgroup[0][0];
		else {
			naiveNormal = BLVector(1.0, 0, 0);
		}
		naiveNormal.normalize();
	}
	return naiveNormal;
}
BLVector BLNode::GetNormal(Node* pNodes, int type)
{
	BLVector ans;
	double max_cost = -1;

	double cost;
	constexpr double throud = 30 * PI / 180.0;
	constexpr double throud2 = 10 * PI / 180.0;
	constexpr double throud3 = 1 * PI / 180.0;
	while (true)
	{
		BLVector tmp;



		if (HasLowerNode()) {
			tmp = GetLowerNode()->GetNormal();
			cost = GetBeitaVisu(tmp);
			if (cost > max_cost) {
				max_cost = cost;
				ans = tmp;
			}
			if (cost > throud) {
				beitaVisu = cost;
				break;
			}
		}


		tmp = GetSimpleNormal(pNodes);
		cost = GetBeitaVisu(tmp);
		if (cost > max_cost) {
			max_cost = cost;
			ans = tmp;
		}
		if (cost > throud) {
			beitaVisu = cost;
			break;
		}


		tmp = GetNaiveNormal(pNodes);
		cost = GetBeitaVisu(tmp);
		if (cost > max_cost) {
			max_cost = cost;
			ans = tmp;
		}
		if (cost > throud) {
			beitaVisu = cost;
			break;
		}

		tmp = GetCenterNormal(pNodes);
		cost = GetBeitaVisu(tmp);
		if (cost > max_cost) {
			max_cost = cost;
			ans = tmp;
		}
		if (cost > throud2) {
			beitaVisu = cost;
			break;
		}

		tmp = GetGeometryNormal(pNodes);
		cost = GetBeitaVisu(tmp);
		if (cost > max_cost) {
			max_cost = cost;
			ans = tmp;
		}
		if (cost > throud3) {
			beitaVisu = cost;
			break;
		}

		//tmp = GetTheMostnNormalNormal(pNodes);
		//cost = GetBeitaVisu(tmp);
		//if (cost > max_cost) {
		//	max_cost = cost;
		//	ans = tmp;
		//}
		//if (cost > throud3) {
		//	beitaVisu = cost;
		//	break;
		//}



		beitaVisu = max_cost;
		break;
	}

	if (GetBSys())
	{
		int nodeid = m_iNod;
		Eigen::RowVector3d start_point(pNodes[nodeid].coord[0], 
			pNodes[nodeid].coord[1], pNodes[nodeid].coord[2]);
		Eigen::RowVector3d normal(ans[0], ans[1], ans[2]);
		int faceid = pNodes[nodeid].isymfc;
		
		//ans[GetSymAxis()] = 0;
		//vec.y = 0;
	}
	m_Normal = ans;
	static double throuhold = sin(60 * PI / 180.0);
	if (beitaVisu < throuhold) {
		auto nei_fronts = GetNeigFronts();
		for (auto &i : nei_fronts) {
			i->SetiSSimple(false);
		}
	}
	return ans;
}




BLVector BLNode::GetHeight()
{
	double d= GetHeightLength();
	if (std::isnan(d)) {
		spdlog::info("warning~ not a number");
	}
	return m_Normal * d;
}


BLVector BLNode::GetHeight(double aver)
{
	BLVector vec(0, 0, 0), normal(0, 0, 0);
	int idx;

	double uVal, lVal, uLen, gradient, dis, disbase;
	double alpha = 1.0, betal = 1.0, gama, ratio = 1.0, uRatio;

	BLFront *blFront = m_vecNeighFronts[0];


	normal.x = m_Normal.x;
	normal.y = m_Normal.y;
	normal.z = m_Normal.z;

	idx = m_iNod;

	gradient = normal.magnitude();
	alpha = normal.x / gradient;
	betal = normal.y / gradient;
	gama = normal.z / gradient;

	ratio = pow(1.3, m_ilayNum);
	//dis = (uVal - 0.95*uVal)*ratio/gradient;
	//dis = /*0.002*/0.036/*0.72*/*(m_ilayNum+1);

	dis = m_h0 * ratio;
	//dis *= exp(-uVal);
	dis = (1 - (uVal - aver) / aver)*dis;


	vec.x = dis * alpha;
	vec.y = dis * betal;
	vec.z = dis * gama;

	return vec;
}

double BLNode::GetHeightLength()
{
	if (expect_height_ < 0) {
		double uVal, lVal, uLen, gradient, dis, disbase;
		double ratio = 1.0;
		if (m_ilayNum > cf.layer_ratio)
		{
			ratio = pow(cf.ratio1, cf.layer_ratio);
			ratio *= pow(cf.ratio2, m_ilayNum - cf.layer_ratio);
		}
		else
			ratio = pow(cf.ratio1, m_ilayNum);
		if (adjacent) {
			auto nodes = this->GetNeigFronts();
			double length = 0;
			for (auto n : nodes) {
				double l = n->GetFrontSize();
				length += l;
			}
			length /= nodes.size();

			m_h0 = length;
		}
		
		expect_height_ = (m_h0 * ratio);
	}
#ifdef CHANGE_STEP_BY_DISTANCE
	return expect_height_ * GetDistanceRatio()*(1 + GetHightRatio());
#else
	return expect_height_ * (1 + GetHightRatio())*(1 + GetFixedHightRatio());
#endif
}

void BLNode::SetStopFlag(bool status)
{

	if (getPerNode())
		getPerNode()->m_bStopPropagate = status;
	m_bStopPropagate = status;
}
