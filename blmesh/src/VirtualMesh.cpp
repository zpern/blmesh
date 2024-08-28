#ifndef _VIRTUAL_MESH_H_
#define _VIRTUAL_MESH_H_
#include <spdlog/spdlog.h> 
 #include "VirtualMesh.h"
#include "BLVector.h"

bool VM::checkOri() {

	std::vector < std::map<int, std::vector<int>> > edge_map(nSN0);
	for (int i = 0; i < pnSNO; i++) {
		std::array<int, 3> triangle;


		for (int j = 0; j < 3; j++) {
			triangle[j] = ppnSFFmO[3 * i + j];
		}
		for (int j = 0; j < 3; j++) {
			int e1 = triangle[j];
			int e2 = triangle[(j + 1) % 3];
			if (e1 > e2)  swap(e1, e2);
			edge_map[e1][e2].push_back(i);
			edge_map[e2][e1].push_back(i);
		}
	}
	int bad_ori = 0;
	for (int i = 0; i < edge_map.size();i++) {
		for (auto j : edge_map[i]) {
			int e1 = i;
			int e2 = j.first;
			int k1 = j.second[0];
			int k2 = j.second[1];
			std::vector<int> triangle1(2);
			std::vector<int> triangle2(2);

			for (int k = 0; k < 3; k++) {
				if (ppnSFFmO[3 * k1 + k] != e1 && ppnSFFmO[3 * k1 + k] != e2) {
					triangle1[0] = ppnSFFmO[3 * k1 + (k + 1) % 3];
					triangle1[1] = ppnSFFmO[3 * k1 + (k + 2) % 3];
				}
			}		

			for (int k = 0; k < 3; k++) {
				if (ppnSFFmO[3 * k2 + k] != e1 && ppnSFFmO[3 * k2 + k] != e2) {
					triangle2[0] = ppnSFFmO[3 * k2 + (k + 1) % 3];
					triangle2[1] = ppnSFFmO[3 * k2 + (k + 2) % 3];
				}
			}

			if (triangle1[0] != triangle2[1] || triangle1[1] != triangle2[0]) {
				cout << triangle1[0] << " " << triangle2[1] << " " << triangle1[1] << " " << triangle2[0] << " " << endl;
				bad_ori++;
			}

		}
	}
	if(bad_ori)
	cout <<"NUm.Bad_ori"<< bad_ori << endl;
	return !bad_ori;
}
bool VM::checkManifold()
{
	std::vector < std::map<int, std::vector<int>> > edge_map(nSN0);
	for (int i = 0; i < pnSNO; i++) {
		std::array<int, 3> triangle;


		for (int j = 0; j < 3; j++) {
			triangle[j] = ppnSFFmO[3 * i + j];
		}
		for (int j = 0; j < 3; j++) {
			int e1 = triangle[j];
			int e2 = triangle[(j + 1) % 3];
			if (e1 > e2)  swap(e1, e2);
			edge_map[e1][e2].push_back(i);
			edge_map[e2][e1].push_back(i);
		}
	}
	bool manifold = true;
	for (int i = 0; i < edge_map.size();i++) {
		for (auto j : edge_map[i])
			if (j.second.size() != 2) {
				if (i < j.first) {
					//cout << "non-manifold find,  point id=(" << i << "," << j.first << ")" << " face_connnect= " << j.second.size() << endl;
					//cout << "Connected face:";
					//for (auto k : j.second) {
					//	cout << k << " ";
					//}
					//cout << endl;
				}
				manifold = false;
			}
	}

	return manifold;
}
//#pragma optimize("",off)
void VM::convert2Manifold()
{
	//saveVTK("test_non_manifold.vtk");

	if (checkManifold()) {
		return;
	}
	else {
		spdlog::info("find non_manifold edge,try to fix it by generate virtual point" );
	}
	std::vector<std::map<int, std::vector<int> > > edge_map(nSN0);

	std::vector<std::set<int> > point_map(nSN0);

	for (int i = 0; i < pnSNO; i++) {
		std::array<int, 3> triangle;


		for (int j = 0; j < 3; j++) {
			triangle[j] = ppnSFFmO[3 * i + j];
		}
		for (int j = 0; j < 3; j++) {
			point_map[triangle[j]].insert(i);
			int e1 = triangle[j];
			int e2 = triangle[(j + 1) % 3];
			if (e1 > e2)  std::swap(e1, e2);
			edge_map[e1][e2].push_back(i);
			edge_map[e2][e1].push_back(i);
		}
	}
	std::set<std::array<int, 2>> non_manifold_edges;
	std::set<int> non_manifold_points;
	for (int i = 0; i < edge_map.size(); i++) {
		for (auto j : edge_map[i]) {
			if (j.second.size() != 2) {
				if (j.second.size() % 2 != 0) {
					spdlog::info(j.second.size());
					throw std::logic_error("can not handle the complex case in dealing with non-manifold edges");
					
				}
				non_manifold_edges.insert(std::array<int, 2>{i, j.first});
				non_manifold_points.insert(i);
				non_manifold_points.insert(j.first);

			}
		}
	}

	auto normal = [](int triangle_id, int* ppnSFFmO, double* ppdSNC0) {
		BLVector n1(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id]);
		BLVector n2(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 1]);
		BLVector n3(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 2]);

		return ((n3 - n2) ^ (n1 - n2)).normalized();

	};
	//�ҵ���Աߵ��������߼�����

	auto centor = [normal](int triangle_id, int* ppnSFFmO, double* ppdSNC0,int pid1,int pid2) {
		BLVector n1(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id]);
		BLVector n2(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 1]);
		BLVector n3(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 2]);
		BLVector n4(ppdSNC0 + 3*pid2);
		BLVector n5(ppdSNC0 + 3*pid1);
		BLVector e = n4 - n5;
		BLVector center = (n1 + n2 + n3) / 3;
		
		BLVector nor = -1*normal(triangle_id, ppnSFFmO,ppdSNC0);
		BLVector v=nor^e;
		
		v.normalize();
	//	cout << e << endl;
		double oritation = v * (center - n4);
		oritation = oritation / abs(oritation);


		return oritation* v;

	};
	// ע�������תԭ��ֻ�����ڱ߽������
	// �����ǲ��ü���׼��ȥ�жϵģ��������������������³��
	auto filp = [edge_map, centor, normal](int triangle_id, int point_id, int centor_point_id, int* ppnSFFmO, double* ppdSNC0) {
		for (auto j : edge_map[centor_point_id]) {
			if (j.first == point_id) {
				if (j.second.size() == 2) {
					return j.second[0] + j.second[1] - triangle_id;
				}
				else {
					double max_value = -1e10;
					int max_id = -1;
					for (auto k : j.second) {
						if (k != triangle_id) {
							BLVector c1 = centor(k, ppnSFFmO, ppdSNC0, point_id, centor_point_id);
							BLVector c2 = centor(triangle_id, ppnSFFmO, ppdSNC0, point_id, centor_point_id);
							BLVector d1 = c1-c2;
							BLVector d2 = -1*normal(triangle_id, ppnSFFmO, ppdSNC0);
							BLVector d3 = -1*normal(k, ppnSFFmO, ppdSNC0);
							// һ�������lengthӦ��ԽСԽ��
							double length = d1.length();
							d1.normalize();
							d2.normalize();
							d3.normalize();
							double obj_value = (d1 * d2- d3*d1)/sqrt(length);
							// ������ȷ
							if((d1 * d2)*(d3*d1)<0)
							if (obj_value > max_value) {
								max_value = obj_value;
								max_id = k;
							}
						}
					}
					return max_id;
				}
			}
			
		}
		return -2;
	};

	std::map<int, int> new_triangle_point_id;
	vector<int> new_point;//�µ��id���±�+ pnSNOΪid,�����ֵΪ�ο����id
	for (auto p : non_manifold_points) {
		std::set<int> triangle_aside_byP = point_map[p];

		
		while (triangle_aside_byP.size()) {
			
			//��״��ת����
			std::set<int> isread;

			//flip
			int first_triangle = *triangle_aside_byP.begin();
			int first_edge_point;
			for (int k = 0; k < 3; k++) {
				if (ppnSFFmO[3 * first_triangle + k] == p) {
					first_edge_point = ppnSFFmO[3 * first_triangle + (k + 1) % 3];
				}
			}


			int next_triangle = -1;
			int next_edge_point = first_edge_point;
			int pre_triangle = first_triangle;
			int loop = 0;
			while (next_triangle != first_triangle) {

				loop++;
				if (loop > 50) {
					cout  << next_triangle << " " ;
				}
				if (loop > 100) {

					
					cout << "Too much iterration in flip function!" << endl;;
					cout << "The out-most mesh is saved as too_much_iteration.vtk" << endl;
					cout << "Point index = " << p << endl;
					saveVTK("too_much_iteration.vtk");

					throw std::logic_error("can not convert model to non manifold!");
					
				}

				isread.insert(pre_triangle);
				next_triangle = filp(pre_triangle, next_edge_point, p, ppnSFFmO, ppdSNC0);
				next_edge_point = ppnSFFmO[3 * next_triangle + 2] + ppnSFFmO[3 * next_triangle] + ppnSFFmO[3 * next_triangle + 1] - next_edge_point - p;
				pre_triangle = next_triangle;
			}
			int new_id = new_point.size() + nSN0;
			new_point.push_back(p);

			for (auto i : isread) {
				for (int j = 0; j < 3; j++) {
					if (ppnSFFmO[3 * i + j] == p) {
						new_triangle_point_id[3 * i + j] = new_id;
					}
				}
			}


			for (auto i : isread)
				triangle_aside_byP.erase(i);
		}
	}
	double* new_coord_array = new double[nSN0 *3 + 3 * new_point.size()];
	for (int i = 0; i < 3 * nSN0; i++) {
		new_coord_array[i] = ppdSNC0[i];
	}
	for (int i = 0; i < new_point.size(); i++) {
		for (int k = 0; k < 3; k++)
			new_coord_array[3 * nSN0 + 3 * i + k] = ppdSNC0[3 * new_point[i] + k];
	}
	for (auto i : new_triangle_point_id) {
		ppnSFFmO[i.first] = i.second;
	}
	delete ppdSNC0;
	ppdSNC0 = new_coord_array;
	nSN0 += new_point.size();
	cout << "Add " << new_point.size() << " points" << endl;
	if (checkManifold()) {
		cout << "succeed" << endl;
		rmoveFreeNodeInTop();
	//	saveVTK("toplayer.vtk");
		return;

	}
	else {
		saveVTK("After.vtk");
		throw std::logic_error("impossible case!");
	}

	if (checkOri()) {
		cout << "succeed" << endl;

	}
	else {
		throw std::logic_error("error oritation");
	}
	return;
}
//#pragma optimize("",on)
//�ϲ�������
void VM::rmoveFreeNodeInTop() {
	vector<int> use_count(nSN0,0);
	vector<int> lgmap(nSN0,0);
	for (int i = 0; i < 3*pnSNO; i++) {
		use_count[ppnSFFmO[i]]++;
	}
	int valid_point_ptr=0;
	for (int i = 0; i < nSN0; i++) {
		if (use_count[i]) {
			for (int k = 0; k < 3; k++) {
				ppdSNC0[3 * valid_point_ptr + k] = ppdSNC0[3 * i + k];
			}
			lgmap[i] = valid_point_ptr;
			valid_point_ptr++;
		}
	}
	for (int i = 0; i < 3 * pnSNO; i++) {
		ppnSFFmO[i]=lgmap[ppnSFFmO[i]];
	}
	if (nSN0 != valid_point_ptr) {
		cout <<"Remove "<< nSN0 - valid_point_ptr <<"bare point!"<< endl;
	}
	nSN0 = valid_point_ptr;
}
//#pragma optimize("",off)
void VM::convert2NonManifold()
{

	std::map<std::array<double,3>,int> first_id_map;
	std::map<int, int> duplicated_node;
	for (int i = 0; i < nSN0; i++) {
		std::array<double, 3> coordinate;
		coordinate[0] = *(ppdSNC0 + 3 * i);
		coordinate[1] = *(ppdSNC0 + 3 * i+1);
		coordinate[2] = *(ppdSNC0 + 3 * i+2);
		if (first_id_map.find(coordinate) == first_id_map.end()) {
			first_id_map[coordinate] = i;
		}
		else {
			duplicated_node[i] = first_id_map[coordinate];
		}
	}
	//û����֤idһ�������
	std::vector<int> offset(nSN0,0);
	for (auto i : duplicated_node) {
		offset[i.first] = i.first - i.second;
	}
	int diff=0;
	for (int i = 0; i < nSN0; i++) {
		if (offset[i]) {
			offset[i] += diff;
			diff += 1;

		}
		else
			offset[i] = diff;
	}
	
	for (int i = 0; i < nSN0; i++) {
		*(ppdSNC0 + 3 * (i - offset[i]))=*(ppdSNC0 + 3 * i);
		*(ppdSNC0 + 3 * (i - offset[i])+1) = *(ppdSNC0 + 3 * i + 1);
		*(ppdSNC0 + 3 * (i - offset[i])+2) = *(ppdSNC0 + 3 * i + 2);
	}
	nSN0 -= diff;
	spdlog::info("remove {}  duplicated nodes diff={}",duplicated_node.size(),diff);

	for (int i = 0; i < pnSNO; i++) {
		for (int j = 0; j < 3; j++) {
			int temp = ppnSFFmO[3 * i + j];
			ppnSFFmO[3*i+j] = ppnSFFmO[3 * i + j]-offset[ppnSFFmO[3 * i + j]];
			if (ppnSFFmO[3 * i + j] >= nSN0|| ppnSFFmO[3 * i + j]<0) {
				cout << "Doing removal of duplicated point,point id=" << temp << endl;
				cout << "offset=" << offset[ppnSFFmO[3 * i + j]] << endl;
				cout << "new id=" << temp - offset[ppnSFFmO[3 * i + j]] << endl;
				cout << "Max_point=" << nSN0 << endl;
				cout << "convert to manifold failed!" << endl;
				throw std::logic_error("impossible case!");
			}
		}
	}

	checkManifold();
}

//#pragma optimize("",on)

void VM::saveVTK(std::string filename)
{
	//return 0;
	int i, j, npt, nelm, sidx;
	npt = nSN0;
	//nelm = m_nElems;
	nelm = pnSNO;

	std::ofstream fout(filename);

	fout<< "# vtk DataFile Version 3.0" << endl;
	fout<<"Really cool data"<<endl;
	fout<<"ASCII" << endl;
	fout<<"DATASET UNSTRUCTURED_GRID" << endl;
	fout<<"POINTS " << npt<<"double"<<endl;

	for (i = 0; i < npt; i++)
	{
		fout <<  ppdSNC0[3 * i + 0] << " " << ppdSNC0[3 * i + 1] << " " << ppdSNC0[3 * i + 2] << endl;;
	}

	sidx = nelm;
	//sidx = 0;
	//sidx = 2914;
		//nelm = sidx + 100;
	fout << "CELLS " << sidx << " " << sidx * 4 << endl;;
	for (i = 0; i < sidx; i++)
		//for (i=sidx; i<nelm; i++)
	{
		fout << "3 " << ppnSFFmO[3 * i + 0] << " " << ppnSFFmO[3 * i + 1] << " " << ppnSFFmO[3 * i + 2] << endl;;
	}

	fout<<"CELL_TYPES "<< sidx << endl;

	for (i = 0; i < sidx; i++)
		fout<< 5 << endl;;
	//for testing
	fout<<"CELL_DATA "<< nelm << endl;
	fout<<"NORMALS vectors double"<< endl;

	auto normal = [](int triangle_id, int* ppnSFFmO, double* ppdSNC0) {
		BLVector n1(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id]);
		BLVector n2(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 1]);
		BLVector n3(ppdSNC0 + 3 * ppnSFFmO[3 * triangle_id + 2]);

		return ((n3 - n2) ^ (n1 - n2)).normalized();



	};
	for (i = 0; i < nelm; i++)
	{
		BLVector N = normal(i, ppnSFFmO, ppdSNC0);
		fout << N[0] << " " << N[1] << " " << N[2] << endl;
	}

	return ;
}
#endif
