#include "SimpleNormalSmoothStrategy.h"
#include <spdlog/spdlog.h> 
 #include <set>
using namespace std;


SimpleNormalSmoothStrategy::SimpleNormalSmoothStrategy(BLNode **node, Node* pnodes,int num_front):NormalSmoothStrategy(node,pnodes, num_front)
{

	
}

void SimpleNormalSmoothStrategy::SmoothNormalOnce(BLNode *blNod, int iq)
{
	//if (blNod->GetVirtualFlag())
	//	return;
	blNod->if_need_smooth =false;
	constexpr double thredhold =25;
	double min_front_size=std::numeric_limits<double>::max();


	int nNeigFrts, nNeigNods;
	int id = blNod->GetNodIdx();
	BLVector centor(pNodes[id].coord);;
#ifdef SMOOTH_FRONT_SIZE
	BLVector current_node = blNod->GetHeight();// ;
	 double size = current_node.magnitude();

	BLVector x(1.0, 0.0, 0.0);
	if (abs(current_node[1]) < 1e-8) {
		x = BLVector(0.0, 1.0, 0.0);
	}

	BLVector z = current_node.normalized();
	BLVector y = x ^ z;
	y.normalize();
	x = z ^ y;
	x.normalize();
#endif
	const std::vector<BLNode*>& blNodes = blNod->GetNeigNods();
	nNeigNods = blNodes.size();
	if (nNeigNods==0)
		return;


	BLVector sum(0, 0, 0);
	double height = 0;
	int count = 0;
	const std::vector <BLFront*>& blfronts=blNod->GetNeigFronts();
	// 周围点数应该等于周围前沿数  
	for (int i = 0; i < nNeigNods; i++) {
		min_front_size =min(min_front_size, blfronts[i]->GetMinFrontSize());
	}
	//frontsize /= nNeigFrts;
	BLVector desentNode (pNodes[blNod->GetDecentID()].coord);
	double height1 = (desentNode-centor).magnitude();

	/*高度宽度比**/
	double ratio = height1 / min_front_size;
	static double lengths[200];
	//std::vector<double> lengths;
	double length_sum = 0;
	for (int i = 0; i < nNeigNods; i++) {
		int idx = blNodes[i]->GetNodIdx();
		BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
		double length = (coord1 - centor).magnitude2();
		lengths[i]=length;
		length_sum += length;
	}
	length_sum /= nNeigNods;

	for (int i = 0; i < nNeigNods; i++) {
		count++;
		height += blNodes[i]->GetHightRatio();
		int id = blNodes[i]->GetNodIdx();
		
		double sp = blNodes[i]->GetBeitaVisu()*2 / (PI);

		int idx=blNodes[i]->GetNodIdx();
		BLVector coord1(pNodes[idx].coord[0],pNodes[idx].coord[1],pNodes[idx].coord[2]);
		double length = lengths[i] / length_sum;
		double angle=(coord1 - centor).normalized()*blNodes[i]->GetNormal();

		sum += blNodes[i]->GetNormal() * (atan(0.01*pow(length, 5*angle)/(sp*sp)));

#ifdef SMOOTH_FRONT_SIZE
		coord = (coord - centor + blNodes[i]->GetHeight()) / max(0.005, size);
		neigh_coord.x += coord * x;
		neigh_coord.y += coord * y;
		neigh_coord.z += coord * z;
#endif
	}
	sum.normalize();
	if (blNod->GetLowerNode()) {
		sum += blNod->GetLowerNode()->GetNormal();
		sum.normalize();
	}
#ifdef SMOOTH_FRONT_SIZE
	if(count)
		height /= count;
	else {
		height = blNod->GetHightRatio();
	}
	
	double fai, theta;
	if (abs(neigh_coord.x) < 1e-8) {
		if (neigh_coord.y>0)
			fai = PI / 2;
		else
			fai = -1 * PI / 2;
	}
	else {
		fai = atan(neigh_coord.y / neigh_coord.x);
	}
	if (fai > PI / 2) {
		fai -= PI;
	}

	double t = neigh_coord.x*cos(fai) + neigh_coord.y*sin(fai);
	//cout << fai << endl;
	if (abs(neigh_coord.z) < 1e-8) {
		if (t > 0) {
			theta = PI / 2;

		}
		else {
			theta = -1 * PI / 2;
		}
	}
	else {
		theta = atan(t / neigh_coord.z);
	}
	if (theta > PI / 2) {
		theta -= PI;
	}
	//cout << neigh_coord.x << " " << neigh_coord.y << " " << neigh_coord.z << endl;
	//cout << fai << "fai " << theta << " theta" << endl;
	theta = min(angle - PI / 6, theta);
	BLVector norm0(0, 0, 0);
	norm0 = norm0 + sin(theta)*cos(fai)*x;
	norm0 = norm0 + sin(theta)*sin(fai)*y;
	norm0 = norm0 + cos(theta)* z;
	//cout << endl << "id=" << blNod->GetNodIdx() << endl;

	norm0.normalize();
	if(blNod->GetVirtualFlag())
		norm0 = (norm0*0.05 +sum + blNod->GetNormal()).normalized();
	else
		norm0 = (norm0*0.1 + sum + blNod->GetNormal()).normalized();
#else
	double r_borm;
	r_borm = (pow(1.7,  ratio) - 1)*15 ;
	BLVector norm0=(r_borm *sum+ blNod->GetNormal());
	norm0.normalize();
	int t = 0;
#endif
	
	constexpr double output_thredhold = thredhold * 0.8;
	
	while (blNod->GetBeitaVisu(norm0) * 180 / PI < output_thredhold) {
		if (t > 10) {
			return;
		}
		norm0 = (norm0 + 0.7*blNod->GetNormal());
		norm0.normalize();
		t++;
	}
	SetSymm(blNod, blNod->GetNodIdx());

	if (norm0*  blNod->GetNormal() < 0.9985) {
		for (int i = 0; i < nNeigNods; i++) {
			blNodes[i]->if_need_smooth = true;

		}
		blNod->if_need_smooth = true;
			
	}
	myvec[iq] = norm0;
//	blNod->SetNormal(norm0);
}


SimpleNormalSmoothStrategy::~SimpleNormalSmoothStrategy()
{
}
