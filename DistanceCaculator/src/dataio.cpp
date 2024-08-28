#include"dataio.h"
#include <spdlog/spdlog.h> 
 #include "kdtree.h"

void rayTriIntersect(Ray r, Triangle t, bool& b, BLVector& interPos, bool& inBou, float &tt, float &u, float &v,double length) {
	if (!rayAABBIntersect(r, t.box,length)) {
		b = false;
		return;
	}
	BLVector e1 = t.pos[1] - t.pos[0];
	BLVector e2 = t.pos[2] - t.pos[0];
	BLVector q = BLVector::crossProduct(r.oritation, e2);
	float a = BLVector::dotProduct(q, e1);
	if (abs(a)<float(1e-9)) {
		b = false;
		return;
	}
	float f = float(1.0) / a;
	BLVector s = r.startPos - t.pos[0];
	u = f * BLVector::dotProduct(q, s);
	if (u<float(0.0)) {
		b = false;
		return;
	}
	BLVector rp = BLVector::crossProduct(s, e1);
	v = f * BLVector::dotProduct(rp, r.oritation);
	if (v<float(0.0f) || u + v>float(1.0f)) {
		b = false;
		return;
	}
	tt = f * BLVector::dotProduct(rp, e2);
	if (tt<0) {
		b = false;
		return;
	}
	b = true;
	interPos = r.startPos + tt * r.oritation;
	inBou = false;
	/*for test*/
	if (u<1e-2 || v<1e-2 || u>1 - 1e-2 || v>1 - 1e-2 || u + v>1 - 1e-2)
		inBou = true;
	return;
}
void rayTriFirstIntersect(Ray r, bool &b, BLVector &interPos, Triangle* &t, bool &inBou, float &u, float &v,TNode* head, kdTree& kdt,double length)
{
	bool success;
	bool inBouIn = false;
	b = false;
	float tt, um, vm;
	BLVector pos;
	Triangle* bestTri;
	BLVector bestPos;
	float min = float(1e10);
	std::vector<Triangle*> vt = kdt.getTriangleMayInterWithRay(head, r,length);
	std::map<int, int> iToi;
	//cout << "size=" << vt.size() << endl;
    vector<pair<double, int> > idx_array(vt.size());
	for (int i = 0; i < vt.size(); i++) {
		idx_array[i].first = ((vt[i]->box.min + vt[i]->box.max) / 2 - r.startPos).magnitude();
		idx_array[i].second = i;
	}
	sort(idx_array.begin(), idx_array.end());


	for (unsigned int i = 0; i< vt.size(); i++) {
		 
		if (iToi.find(vt[idx_array[i].second]->index) != iToi.end())
			continue;
		else {
			iToi[vt[idx_array[i].second]->index] = 1;
		}
//		rayTriIntersect(r, *vt[idx_array[i].second], success, pos, inBouIn, tt, um, vm);
		inBou = inBou || inBouIn;
		if (success) {
			if ((pos - r.startPos).length()<min) {
				min = (pos - r.startPos).length();
				bestTri = vt[idx_array[i].second];
				bestPos = pos;
				u = um;
				v = vm;
				b = true;
			}
			break;
		}
	}
	t = bestTri;
	interPos = bestPos;
	return;
}
map<int, double> GetDistanceRatio(map<int, pair<BLVector, BLVector> >& point, vector<vector<int> >& graph)
{
	cout << "begin to get the distance of every node" << endl;;
	map<int, double> dis;
	TNode* root = new TNode();
	Triangle* tri_array = new Triangle[graph.size()];
	AABB space;
	space.max = BLVector(-1e20,-1e20,-1e20);
	space.min = BLVector(1e20, 1e20, 1e20);
	for (int i = 0; i < graph.size(); i++) {
		tri_array[i].box.min = BLVector(1e20, 1e20, 1e20);
		tri_array[i].box.max = BLVector(-1e20, -1e20, -1e20);
		for (int j = 0; j < 3; j++) {
			tri_array[i].pos[j] = point[graph[i][j]].first;
			for (int k = 0; k < 3; k++) {
				space.max[k]= max(space.max[k], tri_array[i].pos[j][k]);
				space.min[k]=min(space.min[k], tri_array[i].pos[j][k]);
				tri_array[i].box.max[k] = max(tri_array[i].box.max[k], tri_array[i].pos[j][k]);
				tri_array[i].box.min[k] = min(tri_array[i].box.min[k], tri_array[i].pos[j][k]);
			}
		}
		tri_array[i].index = i;
	}
	double length = (space.max - space.min).magnitude() / 10;
	root->dim = 3;
	root->max = space.max+BLVector(1,1,1);
	root->min = space.min - BLVector(1, 1, 1);
	kdTree kdt;
	vector<Triangle*> tri_vec;
	for (int i = 0; i < graph.size(); i++) {
		tri_vec.push_back(tri_array+i);
	}

	for (auto i : point) {
		dis[i.first] = -1;
	}
	cout << "begin to creat kdtree" << endl;
	if (!kdt.buildTree(root, tri_vec, 0)) {
		cout << "build kdtree failed" << endl;
		return dis;
	}
	cout << "creat kdtree done" << endl;
	for (auto i : point) {
		//cout << i.first << endl;
		Ray r;
		r.startPos = i.second.first+ i.second.second*1e-7;
		r.oritation = i.second.second;
		bool isinter,inBou; BLVector interpos;
		Triangle* inter_tri;
		float u, v;
		rayTriFirstIntersect(r, isinter, interpos, inter_tri, inBou, u, v, root, kdt, length);
		if (isinter) {
			dis[i.first] = (interpos - i.second.first).magnitude();
		}
	}




	/*delete some memory*/
	for (auto &i : tri_vec) {
		i = nullptr;
	}
	delete[]tri_array;
	return dis;
}

bool rayAABBIntersect(Ray r,AABB b, double length){
    // sp 线起点
    // amin amax  表示 AABB包围盒
        static const double EPS = 1e-9f;

        // 光线方向
        double d[3];
        double sp[3];
        double amin[3];
        double amax[3];
        amin[0]=b.min.x;
        amin[1]=b.min.y;
        amin[2]=b.min.z;
        amax[0]=b.max.x;
        amax[1]=b.max.y;
        amax[2]=b.max.z;
        double tmax,tmin;
        sp[0]=r.startPos.x;
        sp[1]=r.startPos.y;
        sp[2]=r.startPos.z;
        d[0] = r.oritation.x;
        d[1] = r.oritation.y;
        d[2] = r.oritation.z;
        // 因为是线段 所以参数t取值在0和1之间
        tmin = 0.0;
        tmax = length;

        for (int i = 0; i < 3; i++)
        {
            // 如果光线某一个轴分量为 0，且在包围盒这个轴分量之外，那么直接判定不相交
            if (fabsf(d[i]) < EPS)
            {
                if (sp[i] < amin[i] || sp[i] > amax[i])
                    return false;
            }
            else
            {
                const double ood = 1.0f / d[i];
                // 计算参数t 并令 t1为较小值 t2为较大值
                double t1 = (amin[i] - sp[i]) * ood;
                double t2 = (amax[i] - sp[i]) * ood;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > tmin) tmin = t1;
                if (t2 < tmax) tmax = t2;
                // 判定不相交
                if (tmin > tmax) return false;
            }
        }
        return true;
}



