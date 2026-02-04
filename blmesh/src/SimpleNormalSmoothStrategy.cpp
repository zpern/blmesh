#include "SimpleNormalSmoothStrategy.h"
#include <spdlog/spdlog.h> 
 #include <set>
using namespace std;
#define OLD

SimpleNormalSmoothStrategy::SimpleNormalSmoothStrategy(BLNode **node, MBLNode* pnodes,int num_front):NormalSmoothStrategy(node,pnodes, num_front)
{

	
}
#ifdef OLD
void SimpleNormalSmoothStrategy::SmoothNormalOnce(BLNode *blNod, int iq)
{
	//if (blNod->GetVirtualFlag())
	//	return;
	blNod->if_need_smooth =false;
	constexpr double thredhold =35;
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
    BLVector neigh_coord;
#endif
	const std::vector<BLNode*>& blNodes = blNod->GetNeigNods();
	nNeigNods = blNodes.size();
	if (nNeigNods==0)
		return;


	BLVector sum(0, 0, 0);
	double height = 0;
	int count = 0;
	const std::vector <BLFront*>& blfronts=blNod->GetNeigFronts();
    nNeigFrts = blfronts.size();
	// 周围点数应该等于周围前沿数  
	for (int i = 0; i < nNeigFrts; i++) {
		min_front_size =min(min_front_size, blfronts[i]->GetMinFrontSize());
	}
	//frontsize /= nNeigFrts;
	BLVector descentNode (pNodes[blNod->GetDecentID()].coord);
	double height1 = (descentNode-centor).magnitude();

	/*高度宽度比**/
	double ratio = height1 / min_front_size;
    std::vector<double> lengths(nNeigNods);
	double length_sum = 0;
	for (int i = 0; i < nNeigNods; i++) {
		int idx = blNodes[i]->GetNodIdx();
		BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
		double length = (coord1 - centor).magnitude2();
		lengths[i]=length;
		length_sum += length;
	}
    double length_ave = length_sum / nNeigNods;

	for (int i = 0; i < nNeigNods; i++) {
		count++;
		height += blNodes[i]->GetHightRatio();
		
		double sp = blNodes[i]->GetBeitaVisu()*2 / (PI);

		int idx=blNodes[i]->GetNodIdx();
		BLVector coord1(pNodes[idx].coord[0],pNodes[idx].coord[1],pNodes[idx].coord[2]);
		double length = lengths[i]/length_ave;
		double angle=(coord1 - centor).normalized()*blNodes[i]->GetNormal();

		double x = pow(length, 2 + 2 * angle) / (sp * sp);
		double weight = 0.30 + 0.5 * x / (1 + x); // 保证权重在 [0.5, 1) 区间
		sum += blNodes[i]->GetNormal() * weight;
		//sum += blNodes[i]->GetNormal() * (atan(0.01*pow(length, 5*angle)/(sp*sp)));

#ifdef SMOOTH_FRONT_SIZE
		coord1 = (coord1 - centor + blNodes[i]->GetHeight()) / max(0.005, size);
		neigh_coord.x += coord1 * x;
		neigh_coord.y += coord1 * y;
		neigh_coord.z += coord1 * z;
#endif
	}
	sum.normalize();
	if (blNod->GetLowerNode()) {
		sum += blNod->GetLowerNode()->GetNormal();
		sum.normalize();
	}
#ifdef SMOOTH_FRONT_SIZE
	if(count)
		double height_ave =height / count;
	else {
		double height_ave = blNod->GetHightRatio();
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
	//cout << neigh_coord.x << " " << neigh_coord.y << " " << neigh_coord.z << endl;
	//cout << fai << "fai " << theta << " theta" << endl;
	//theta = std::min(angle - PI / 6, theta);
	BLVector norm0(0, 0, 0);
	norm0 = norm0 + sin(theta)*cos(fai)*x;
	norm0 = norm0 + sin(theta)*sin(fai)*y;
	norm0 = norm0 + cos(theta)* z;
	//cout << endl << "id=" << blNod->GetNodIdx() << endl;

	norm0.normalize();
	if(blNod->GetVirtualFlag())
		norm0 = (norm0*0.1 +sum + blNod->GetNormal()).normalized();
	else
		norm0 = (norm0*0.2 + sum + blNod->GetNormal()).normalized();
#else
	double r_borm;
	r_borm = (pow(1.7,  ratio) - 1)*15 ;
	BLVector norm0=(r_borm *sum+ blNod->GetNormal());
	norm0.normalize();
#endif
	
	constexpr double output_thredhold = thredhold * 0.8;
    int n = 0;
	while (blNod->GetBeitaVisu(norm0) * 180 / PI < output_thredhold) {
		if (n > 10) {
			return;
		}
		norm0 = (norm0 + 0.7*blNod->GetNormal());
		norm0.normalize();
		n++;
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
#else
void SimpleNormalSmoothStrategy::SmoothNormalOnce(BLNode * blNod, int iq)
    {
        // if (blNod->GetVirtualFlag())
        //	return;
        blNod->if_need_smooth = false;
        constexpr double thredhold = 35;
        double min_front_size = std::numeric_limits<double>::max();

        int nNeigFrts, nNeigNods;
        int id = blNod->GetNodIdx();
        int deid = blNod->GetDecentID();
        BLVector centor(pNodes[id].coord);
        BLVector centor1(pNodes[deid].coord);
        auto eq = [](double a, double b, double eps = 1e-2) { return std::abs(a - b) <= eps; };

#ifdef SMOOTH_FRONT_SIZE
        BLVector current_node = blNod->GetHeight(); // ;
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
        BLVector neigh_coord;
#endif
        const std::vector<BLNode *> &blNodes = blNod->GetNeigNods();
        nNeigNods = blNodes.size();
        if (nNeigNods == 0) {
            return;
        }

        BLVector sum(0, 0, 0);
        double height = 0;
        int count = 0;
        const std::vector<BLFront *> &blfronts = blNod->GetNeigFronts();
        // 周围点数应该等于周围前沿数,对称面存在时则少一
        nNeigFrts = (int)blfronts.size();
        for (int i = 0; i < nNeigFrts; i++) {
            min_front_size = std::min(min_front_size, blfronts[i]->GetMinFrontSize());
        }
        // frontsize /= nNeigFrts;
        BLVector descentNode(pNodes[blNod->GetDecentID()].coord);
        double height1 = (descentNode - centor).magnitude();

        /*高度宽度比**/
        double ratio = height1 / min_front_size;
        std::vector<double> lengths(nNeigNods);
        double length_sum = 0;
        for (int i = 0; i < nNeigNods; i++) {
            int idx = blNodes[i]->GetNodIdx();
            BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
            double length = (coord1 - centor).magnitude2();
            lengths[i] = length;
            length_sum += length;
        }
        double length_ave = length_sum / nNeigNods;

        for (int i = 0; i < nNeigNods; i++) {
            count++;
            height += blNodes[i]->GetHightRatio();

            double sp = blNodes[i]->GetBeitaVisu() * 2 / (PI);
            int idx = blNodes[i]->GetNodIdx();
            BLVector coord1(pNodes[idx].coord[0], pNodes[idx].coord[1], pNodes[idx].coord[2]);
            double length = length_ave / lengths[i];
            double angle = abs((coord1 - centor).normalized() * blNodes[i]->GetNormal());

            double x = pow(length, 3 + 2 * angle) / (sp * sp);
            double weight = 0.5 + x / (1 + x); // 保证权重在 [0.5, 1) 区间
            sum += blNodes[i]->GetNormal() * weight;
            // sum += blNodes[i]->GetNormal() * (atan(0.01*pow(length, 5*angle)/(sp*sp)));

#ifdef SMOOTH_FRONT_SIZE
            coord1 = (coord1 - centor + blNodes[i]->GetHeight()) / max(0.005, size);
            neigh_coord.x += coord1 * x;
            neigh_coord.y += coord1 * y;
            neigh_coord.z += coord1 * z;
#endif
        }
        sum.normalize();
        if (blNod->GetLowerNode()) {
            sum += blNod->GetLowerNode()->GetNormal();
            sum.normalize();
        }
#ifdef SMOOTH_FRONT_SIZE
        if (count) {
            double height_ave = height / count;
        } else {
            double height_ave = blNod->GetHightRatio();
        }

        double fai, theta;
        if (abs(neigh_coord.x) < 1e-8) {
            if (neigh_coord.y > 0) {
                fai = PI / 2;
            } else {
                fai = -1 * PI / 2;
            }
        } else {
            fai = atan(neigh_coord.y / neigh_coord.x);
        }

        double t = neigh_coord.x * cos(fai) + neigh_coord.y * sin(fai);
        // cout << fai << endl;
        if (abs(neigh_coord.z) < 1e-8) {
            if (t > 0) {
                theta = PI / 2;

            } else {
                theta = -1 * PI / 2;
            }
        } else {
            theta = atan(t / neigh_coord.z);
        }
        // cout << neigh_coord.x << " " << neigh_coord.y << " " << neigh_coord.z << endl;
        // cout << fai << "fai " << theta << " theta" << endl;
        // theta = std::min(angle - PI / 6, theta);
        BLVector norm0(0, 0, 0);
        norm0 = norm0 + sin(theta) * cos(fai) * x;
        norm0 = norm0 + sin(theta) * sin(fai) * y;
        norm0 = norm0 + cos(theta) * z;
        // cout << endl << "id=" << blNod->GetNodIdx() << endl;

        norm0.normalize();
        if (blNod->GetVirtualFlag()) {
            norm0 = (norm0 * 0.1 + sum + blNod->GetNormal()).normalized();
        } else {
            norm0 = (norm0 * 0.2 + sum + blNod->GetNormal()).normalized();
        }
#else
        double r_borm;
        r_borm = r_borm = 5 / (1.0 + 3 * ratio);
        BLVector norm0 = (r_borm * sum + blNod->GetNormal());
        norm0.normalize();
#endif

        constexpr double output_thredhold = thredhold * 0.8;
        int n = 0;
        while (blNod->GetBeitaVisu(norm0) * 180 / PI < output_thredhold) {
            if (n > 10) {
                return;
            }
            norm0 = (norm0 + 0.7 * blNod->GetNormal());
            norm0.normalize();
            n++;
        }
        SetSymm(blNod, blNod->GetNodIdx());

        if (norm0 * blNod->GetNormal() < 0.9985) {
            for (int i = 0; i < nNeigNods; i++) {
                blNodes[i]->if_need_smooth = true;
            }
            blNod->if_need_smooth = true;
        }
        myvec[iq] = norm0;
        //	blNod->SetNormal(norm0);
    }
#endif

SimpleNormalSmoothStrategy::~SimpleNormalSmoothStrategy()
{
}
