#include <spdlog/spdlog.h> 
 #include "optdihangle.h"
#include "iso3d_define.h"

float getAngleSin(MeshElment me)
{
	ElemQual *eleQual = nullptr;
	eleQual = me.elmQual;
	int val = eleQual->iReserved;
	int idx = getPriority(val);
	
	if(eleQual)
		return (eleQual->sinValue * (1+0.25*idx));
	else
		return -1.0;
}

/* ---------------------------------------------------------------------------------------
 * 调用getAngleSin，比较ang1和ang2
 * -------------------------------------------------------------------------------------*/
int compareDiheAng(MeshElment me1, MeshElment me2)
{
	return 1;
}

//setter function
void setAcute(int &val, int acute)
{
	assert(acute == 0 || acute == 1);
	val &= (~0x1);
	val |= acute;
}

void setPriority(int &val, int priority)
{
	assert(priority >= 0 && priority <= 7);
	val &= (~(0x7 << 1));
	val |= (priority << 1);
}

void setAngleCode(int &val, int angCode)
{
	assert(angCode >= 0 && angCode <= 5);
	val &= (~(0x7 << 4));
	val |= (angCode << 4);
}

void setHeapIndex(int &val, int idx)
{
	assert(idx >= 0);
	val &= (0x7F); //127
	val |= (idx << 7);
}

//getter function
int getAcute(int val)
{
	return (val & 0x1);
}

int getPriority(int val)
{
	return ( (val >> 1) & 0x7 );
}

int getAngleCode(int val)
{
	return ( (val >> 4) & 0x7 );
}

int getHeapIndex(int val)
{
	return (val >> 7);
}

int getEdge(MeshElment mshElm, int ep[2])
{
	int code;
	code = getAngleCode(mshElm.elmQual->iReserved);
	ep[0] = tetra_pnt[code][1];
	ep[1] = tetra_pnt[code][2];

	return 1;
}

int getEdgeFace(MeshElment mshElm, int fc[2])
{
	int code;
	code = getAngleCode(mshElm.elmQual->iReserved);
	fc[0] = tetra_pnt[code][0];
	fc[1] = tetra_pnt[code][3];

	return 1;
}

//ActElemHeap
/************************************************************************/
/* 含有n个节点的堆的高度为lgn，故插入一个节点的时间复杂度为O（lgn）     */
/************************************************************************/
int ActElemHeap::insertElem(MeshElment elm)
{
	if(nHeapSize == nMaxHeapSize)
		allocHeap();

	m_actElmHeap[nHeapSize] = elm;
	setHeapIndex(m_actElmHeap[nHeapSize].elmQual->iReserved, nHeapSize+1);
	minHeapifyUp(nHeapSize);
	nHeapSize++;

#if 0
   if(nHeapSize % 10 == 0 )
   {
	   printfinfo();
	   removeElem(3);
	   printfinfo();
   }
#endif

	return 0;
}


//idx为从1开始
void ActElemHeap::removeElem(int idx)
{
	assert(idx >= 0);
	if (nHeapSize <= 0)
		return;
	//index设为无效值
	setHeapIndex(m_actElmHeap[idx].elmQual->iReserved, 0);
	m_actElmHeap[idx] = m_actElmHeap[nHeapSize-1];
	//更新index
	setHeapIndex(m_actElmHeap[idx].elmQual->iReserved, idx+1);
	nHeapSize--;

	//
	minHeapifyDown(idx);
}

MeshElment ActElemHeap::getMinElem()
{
	MeshElment elm;
	if(isEmpty())
		return elm;

	elm = m_actElmHeap[0];
	
	return elm;
}

void ActElemHeap::rmMinElem()
{
	if(isEmpty())
		return ;
	//设置index
//	if (m_actElmHeap[0].oldElem == 8386)
//	{
//		spdlog::info("8386.ireserved = {}.\n", 0);
//	}
	setHeapIndex(m_actElmHeap[0].elmQual->iReserved, 0);
	//if (nHeapSize > 1)
	//{
		m_actElmHeap[0] = m_actElmHeap[nHeapSize-1];
		//更新index
//		if (m_actElmHeap[0].oldElem == 8386)
//		{
//			spdlog::info("8386.ireserved = {}; nHeapSize = {}\n", 1, nHeapSize - 1);
//		}
		setHeapIndex(m_actElmHeap[0].elmQual->iReserved, 1);
	//}
	nHeapSize--;

	minHeapifyDown(0);
}

//从上往下调整
void ActElemHeap::minHeapifyDown(int beg)
{
	int i, j;
	MeshElment temp;

	i = beg;
	j = 2*i+1;	//左子女编号
	temp = m_actElmHeap[i];

	while(j < nHeapSize)
	{
		if ( j < nHeapSize-1 && getAngleSin(m_actElmHeap[j]) > getAngleSin(m_actElmHeap[j+1]) )
			j++;
		if(getAngleSin(temp) < getAngleSin(m_actElmHeap[j]))
			break;
		else
		{
			m_actElmHeap[i] = m_actElmHeap[j];
			//更新index
			setHeapIndex(m_actElmHeap[i].elmQual->iReserved, i+1);
			i = j;
			j = 2*j+1;
		}
	}

	m_actElmHeap[i] = temp;
	//更新index
	setHeapIndex(m_actElmHeap[i].elmQual->iReserved, i+1);
}

//从下往上调整
void ActElemHeap::minHeapifyUp(int beg)
{
	int i, j;
	MeshElment temp;

	j = beg;
	i = (j-1)/2;	//父亲节点
	temp = m_actElmHeap[j];
	while(j>0)
	{
		if(getAngleSin(m_actElmHeap[i]) < getAngleSin(temp))
			break;
		else
		{
			m_actElmHeap[j] = m_actElmHeap[i];
			//更新index
			setHeapIndex(m_actElmHeap[j].elmQual->iReserved, j+1);
			j = i;
			i = (i-1)/2;
		}
	}
	m_actElmHeap[j] = temp;
	//更新index
	setHeapIndex(m_actElmHeap[j].elmQual->iReserved, j+1);
}

void ActElemHeap::addPriority(int idx)
{
	int prio, index;

	index = getHeapIndex(m_actElmHeap[idx].elmQual->iReserved);
	prio = getPriority(m_actElmHeap[idx].elmQual->iReserved);
	assert(index == idx+1);
	if(prio == 6)
	{
		//spdlog::info("%lf\n", m_actElmHeap[idx].elmQual->sinValue);
		removeElem(idx);
	}
	else
	{
		setPriority(m_actElmHeap[idx].elmQual->iReserved, prio+1);
		minHeapifyDown(idx);
	}
}

void ActElemHeap::printfinfo()
{
	MeshElment temp;
	//temp = getMinAngle();
	//do 
	//{
		
	//	spdlog::info("value: %lf\n", getAngleSin(temp));

	//	temp = getMinAngle();
	//}// while (temp);
	for (int i=0; i<nHeapSize; i++)
	{
		temp = m_actElmHeap[i];
		printf("index: %d, value: %lf\n", getHeapIndex(temp.elmQual->iReserved), temp.elmQual->sinValue);
	}
	spdlog::info("\n");
}


int ActElemHeap::heapSize()
{ 
	return nHeapSize; 
}