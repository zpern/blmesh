// OctreeAgent.cpp: implementation of the OctreeAgent class.


#include <spdlog/spdlog.h> 
 #include "OctreeAgent.h"
#include <vector>

//OctreeAgent::OctreeAgent(int* pEle, double* pCoord)
OctreeAgent::OctreeAgent(DynamicArray<SearchTriangle>&  tri_elm, Node* pNod):pEle(tri_elm)
{
	//this->pCoord = pCoord;
	this->pNod = pNod;

}

OctreeAgent::~OctreeAgent()
{

}

BLVector OctreeAgent::getnormal(int data)
{
	BLVector pos[3];
	for (int i = 0; i < 3; i++) {
		for(int j=0;j<3;j++)
			pos[i][j]=pNod[pEle[data][i]].coord[j];
	}
	return BLVector::crossProduct(pos[0] - pos[1],pos[0]-pos[2]);
}

std::array<BLVector, 3> OctreeAgent::getcoord(int data)
{
	std::array<BLVector, 3> ans;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			ans[i][j] = pNod[pEle[data][i]].coord[j];
	}
	return ans;
}
std::array<double, 3> OctreeAgent::getcoord(int data,int dim)
{
	std::array<double, 3> ans;
	for (int i = 0; i < 3; i++) {
			ans[i] = pNod[pEle[data][i]].coord[dim];
	}
	return ans;
}
void OctreeAgent::divideData(const OCCUBE& cube, const std::vector<int>& vecData, std::vector<int>(& outData)[8])
{
	int idx, cubeIdx;

	//changed start
	int i,j;
	OCCUBE subcube[OC_SUBNODES];
	
	//divide cube to eight sub-cubes
	CUBECOORD coorcenter, len;
	coorcenter = (cube.lower + cube.upper)/2;
	len = (cube.upper - cube.lower)/2; 
	
	//subcube[UNE]
	subcube[UNE].lower = coorcenter;
	subcube[UNE].upper = coorcenter + len;
	
	//subcube[UNW]
	subcube[UNW].lower = coorcenter;
	subcube[UNW].lower.x -= len.x;
	subcube[UNW].upper = coorcenter;
	subcube[UNW].upper.y += len.y;
	subcube[UNW].upper.z += len.z;
	
	//subcube[USW]
	subcube[USW].lower = coorcenter;
	subcube[USW].lower.x -= len.x;
	subcube[USW].lower.y -= len.y;
	subcube[USW].upper = coorcenter;
	subcube[USW].upper.z += len.z;
	
	//subcube[USE]
	subcube[USE].lower = coorcenter;
	subcube[USE].lower.y -= len.y;
	subcube[USE].upper = coorcenter;
	subcube[USE].upper.x += len.x;
	subcube[USE].upper.z += len.z;
	
	//subcube[LNE]
	subcube[LNE].lower = coorcenter;
	subcube[LNE].lower.z -= len.z;
	subcube[LNE].upper = coorcenter;
	subcube[LNE].upper.x += len.x;
	subcube[LNE].upper.y += len.y;
	
	//subcube[LNW]
	subcube[LNW].lower = coorcenter;
	subcube[LNW].lower.x -= len.x;
	subcube[LNW].lower.z -= len.z;
	subcube[LNW].upper = coorcenter;
	subcube[LNW].upper.y += len.y;
	
	//subcube[LSW]
	subcube[LSW].lower = coorcenter - len;
	subcube[LSW].upper = coorcenter;
	
	//subcube[LSE]
	subcube[LSE].lower = coorcenter;
	subcube[LSE].lower.y -= len.y;
	subcube[LSE].lower.z -= len.z;
	subcube[LSE].upper = coorcenter;
	subcube[LSE].upper.x += len.x;


	for (i=0; i<vecData.size(); i++)
	{
		/*
		idx = vecData[i];

		cubeIdx = patchCube(idx, cube);

		if( cubeIdx != -1)
			outData[cubeIdx].push_back(idx);
		else
			outData[8].push_back(idx);
		*/

		idx = vecData[i];
		for (j=0; j<OC_SUBNODES; j++)
		{
			if( patchInCube(idx, subcube[j]) )
				outData[j].push_back(idx);
		}
	}
}

int OctreeAgent::patchCube(int patchIdx,const OCCUBE& cube)
{
	int i;
	OCCUBE subcube[OC_SUBNODES];
	
	//divide cube to eight sub-cubes
	CUBECOORD coorcenter, len;
	coorcenter = (cube.lower + cube.upper)/2;
	len = (cube.upper - cube.lower)/2; 
	
	//subcube[UNE]
	subcube[UNE].lower = coorcenter;
	subcube[UNE].upper = coorcenter + len;
	
	//subcube[UNW]
	subcube[UNW].lower = coorcenter;
	subcube[UNW].lower.x -= len.x;
	subcube[UNW].upper = coorcenter;
	subcube[UNW].upper.y += len.y;
	subcube[UNW].upper.z += len.z;
	
	//subcube[USW]
	subcube[USW].lower = coorcenter;
	subcube[USW].lower.x -= len.x;
	subcube[USW].lower.y -= len.y;
	subcube[USW].upper = coorcenter;
	subcube[USW].upper.z += len.z;
	
	//subcube[USE]
	subcube[USE].lower = coorcenter;
	subcube[USE].lower.y -= len.y;
	subcube[USE].upper = coorcenter;
	subcube[USE].upper.x += len.x;
	subcube[USE].upper.z += len.z;
	
	//subcube[LNE]
	subcube[LNE].lower = coorcenter;
	subcube[LNE].lower.z -= len.z;
	subcube[LNE].upper = coorcenter;
	subcube[LNE].upper.x += len.x;
	subcube[LNE].upper.y += len.y;
	
	//subcube[LNW]
	subcube[LNW].lower = coorcenter;
	subcube[LNW].lower.x -= len.x;
	subcube[LNW].lower.z -= len.z;
	subcube[LNW].upper = coorcenter;
	subcube[LNW].upper.y += len.y;
	
	//subcube[LSW]
	subcube[LSW].lower = coorcenter - len;
	subcube[LSW].upper = coorcenter;
	
	//subcube[LSE]
	subcube[LSE].lower = coorcenter;
	subcube[LSE].lower.y -= len.y;
	subcube[LSE].lower.z -= len.z;
	subcube[LSE].upper = coorcenter;
	subcube[LSE].upper.x += len.x;
	

	for (i=0; i<OC_SUBNODES; i++)
	{
		if( patchInCube(patchIdx, subcube[i]) )
			return i;
	}

	return -1;  //the patch is in cube, but not in any sub-cube
}

/*
bool OctreeAgent::patchInCube(int patchIdx, const OCCUBE& cube)
{
	int i, ptIdx;
	bool bIn = true;

	for (i=0; i<3; i++)
	{
		ptIdx = pEle[patchIdx*4 + i];
		Coord coord(pCoord[ptIdx*3 + 0], pCoord[ptIdx*3 + 1], pCoord[ptIdx*3 + 2]);

		bIn = bIn && (coord <= cube.upper) && (coord >= cube.lower);
	}

	return bIn;
}*/
bool OctreeAgent::LineInCube(BLVector start, BLVector end, const OCCUBE& cube) {

	static const double EPS = 1e-7;
	double d[3];
	double sp[3];
	double amin[3];
	double amax[3];
	amin[0] = cube.lower.x;
	amin[1] = cube.lower.y;
	amin[2] = cube.lower.z;
	amax[0] = cube.upper.x;
	amax[1] = cube.upper.y;
	amax[2] = cube.upper.z;
	double tmax, tmin;
	sp[0] = start.x;
	sp[1] = start.y;
	sp[2] = start.z;
	d[0] = end.x-start.x;
	d[1] = end.y-start.y;
	d[2] = end.z-start.z;
	// 因为是线段 所以参数t取值在0和1之间
	tmin = 1e20;
	tmax = -1e20;

	for (int i = 0; i < 3; i++)
	{
		// 如果光线某一个轴分量为 0，且在包围盒这个轴分量之外，那么直接判定不相交
		if (fabs(d[i]) < EPS)
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

			// 判定不相交
			if (t1 < 0&&t2<0) return false;
			if (t1 > 1&&t2>1) return false;
		}
	}

	return true;
}
/*三角形完全在cube里*/
inline bool OctreeAgent::patchTotalInCube(int patchIdx, const OCCUBE& cube)
{
	
	int  ptIdx;
	if ( patchIdx<0) {
		//cout << "warning, patchIdx error patch idx=" <<patchIdx<<" elmsize="<<numElem<< endl;
		return true;
	}
	std::array<double, 6> box;
	getBox(patchIdx,box);
	return box[0] < cube.upper.x&&box[1] < cube.upper.y&&box[2] < cube.upper.z&&\
		box[3] > cube.lower.x&&box[4] > cube.lower.y&&box[5] > cube.lower.z;


}
bool OctreeAgent::LineTotalInCube(BLVector start, BLVector end, const OCCUBE & cube)
{
	if (start.x < cube.lower.x || start.y < cube.lower.y || start.z < cube.lower.z)
		return false;
	if (start.x > cube.upper.x || start.y > cube.upper.y || start.z > cube.upper.z)
		return false;
	if (end.x < cube.lower.x || end.y < cube.lower.y || end.z < cube.lower.z)
		return false;
	if (end.x > cube.upper.x || end.y > cube.upper.y || end.z > cube.upper.z)
		return false;
	return true;
}
bool OctreeAgent::patchInCube(int patchIdx, const OCCUBE& cube)
{
	int i = 0, ptIdx;// = pEle[patchIdx * 4 + 0];
//	bool bIn = true;
	int s = 0;

	//

	//CUBECOORD coord_max(pNod[ptIdx].coord[0], pNod[ptIdx].coord[1], pNod[ptIdx].coord[2]);
	//CUBECOORD coord_min(pNod[ptIdx].coord[0], pNod[ptIdx].coord[1], pNod[ptIdx].coord[2]);
	//ptIdx = pEle[patchIdx * 4 + 1];

	//coord_max.x = max(pNod[ptIdx].coord[0], coord_max.x);
	//coord_max.y = max(pNod[ptIdx].coord[1], coord_max.y);
	//coord_max.z = max(pNod[ptIdx].coord[2], coord_max.z);
	//coord_min.x = min(pNod[ptIdx].coord[0], coord_min.x);
	//coord_min.y = min(pNod[ptIdx].coord[1], coord_min.y);
	//coord_min.z = min(pNod[ptIdx].coord[2], coord_min.z);

	//ptIdx = pEle[patchIdx * 4 + 2];

	//coord_max.x = max(pNod[ptIdx].coord[0], coord_max.x);
	//coord_max.y = max(pNod[ptIdx].coord[1], coord_max.y);
	//coord_max.z = max(pNod[ptIdx].coord[2], coord_max.z);
	//coord_min.x = min(pNod[ptIdx].coord[0], coord_min.x);
	//coord_min.y = min(pNod[ptIdx].coord[1], coord_min.y);
	//coord_min.z = min(pNod[ptIdx].coord[2], coord_min.z);


	//Coord dis = (cube.lower + cube.upper)  - (coord_max + coord_min) ;
	//Coord maxd = (cube.upper - cube.lower) + (coord_max - coord_min) ;
	//return (abs(dis.x) <= maxd.x&&abs(dis.y) <= maxd.y&&abs(dis.z) <= maxd.z);

	//这里修改了,因为该函数的性能对整体性能影响较大，因此做了深度优化，采用优化写法

	for (int i = 0; i < 3; i++) {
		ptIdx = pEle[patchIdx][i];
		if (pNod[ptIdx].coord[0] > cube.lower.x)
			s |= 1;
		if (pNod[ptIdx].coord[0] <= cube.upper.x)
			s |= 2;

		if (pNod[ptIdx].coord[1] > cube.lower.y)
			s |=4;
		if (pNod[ptIdx].coord[1] <= cube.upper.y)
			s |=8;

		if (pNod[ptIdx].coord[2] > cube.lower.z)
			s |= 16;
		if (pNod[ptIdx].coord[2] <= cube.upper.z)
			s |= 32;
	}

	return s == 63;
	
}

void OctreeAgent::getBox(const int & index, std::array<double, 6>& box)
{
	for (int i = 0; i < 3; i++) {
		box[i] =  pNod[pEle[index ][0]].coord[i];
		box[3 + i] = pNod[pEle[index ][0]].coord[i];
	}
	for (int i = 0; i < 3; i++) {
		box[i] = max(box[i], pNod[pEle[index ][1]].coord[i]);
		box[3 + i] = min(box[3 + i], pNod[pEle[index ][1]].coord[i]);
		box[i] = max(box[i], pNod[pEle[index ][2]].coord[i]);
		box[3 + i] = min(box[3 + i], pNod[pEle[index ][2]].coord[i]);
	}
	// TODO: 在此处插入 return 语句
}
