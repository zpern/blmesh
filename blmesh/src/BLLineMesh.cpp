#include "BLLineMesh.h"
#include <spdlog/spdlog.h> 
 #include <cmath>
#include <vector>
#include <string>

BLLineMesh::BLLineMesh(void)
{
}


BLLineMesh::~BLLineMesh(void)
{
}

#if 0
void BLLineMesh::GenMeshLine(int nbpt, double *bpt, double *spa, int *npt,
	                         double **pt, int *nElm, int **pElm)
{
	int nsg, tnsg, i, j, k, idx, cnt, istart;
	double len, aspa;

	if(nbpt%2 != 0)
	{
		spdlog::info("Error: wrong number of points!\n");
		exit(0);
	}
	cnt = nbpt/2;
	tnsg = 0;
	for (idx=0; idx<cnt; idx++)
	{
		aspa = (spa[idx*2+0]+spa[idx*2+1])/2;
		len = fabs(bpt[idx*2+0] - bpt[idx*2+1]);
		nsg = (len/aspa + 0.5);
		tnsg += nsg;
	}

	*pt = new double[2*tnsg + cnt];
	*pElm = new int[tnsg*2];

	for (i=0; i<nbpt; i++)
		(*pt)[i] = bpt[i];

	j = nbpt;
	k = 0;
	for (idx=0; idx<cnt; idx++)
	{
		aspa = (spa[idx*2+0]+spa[idx*2+1])/2;
		len = fabs(bpt[idx*2+0] - bpt[idx*2+1]);
		nsg = (len/aspa + 0.5);
		istart = j;
		//right now only uniform spacing is considered

		if (bpt[2*idx+0] < bpt[2*idx+1])
		{
			for (i=1; i<nsg; i++)
				(*pt)[j++] = bpt[2*idx+0] + i*aspa;

			for (i=1; i<nsg-1; i++)
			{
				(*pElm)[2*k+0] = istart+i-1;
				(*pElm)[2*k+1] = istart+i;
				k++;
			}

			if(nsg = 1)
			{
				(*pElm)[2*k+0] = 2*idx+0;
				(*pElm)[2*k+1] = 2*idx+1;
				k++;
			}
			else
			{
				(*pElm)[2*k+0] = 2*idx+0;
				(*pElm)[2*k+1] = istart;
				k++;
				(*pElm)[2*k+0] = j-1;
				(*pElm)[2*k+1] = 2*idx+1;
				k++;
			}
		} 
		else
		{
			for (i=1; i<nsg; i++)
				(*pt)[j++] = bpt[2*idx+1] + i*aspa;

			for (i=1; i<nsg-1; i++)
			{
				(*pElm)[2*k+0] = istart+i-1;
				(*pElm)[2*k+1] = istart+i;
				k++;
			}

			if(nsg = 1)
			{
				(*pElm)[2*k+0] = 2*idx+0;
				(*pElm)[2*k+1] = 2*idx+1;
				k++;
			}
			else
			{
				(*pElm)[2*k+0] = 2*idx+1;
				(*pElm)[2*k+1] = istart;
				k++;

				(*pElm)[2*k+0] = j-1;
				(*pElm)[2*k+1] = 2*idx+0;
				k++;
			}
		}
	}
}
#else
void BLLineMesh::GenMeshLine(int nbpt, double *bpt, double *spa, int *npt,
	double **pt, int *nElm, int **pElm)
{
	int nsg, tnsg, i, j, k, idx, cnt, istart;
	double len, aspa;

	if(nbpt%2 != 0){

		throw(std::string("Error: wrong number of points!"));
	}
	cnt = nbpt/2;
	tnsg = 0;
	for (idx=0; idx<cnt; idx++)
	{
		aspa = (spa[idx+0]+spa[idx+cnt])/2;
		len = fabs(bpt[idx+0] - bpt[idx+cnt]);
		nsg = (len/aspa + 0.5);
		tnsg += nsg;
	}

	*pt = new double[tnsg + cnt];
	*pElm = new int[tnsg*2];
	*npt = tnsg+cnt;
	*nElm = tnsg;

	for (i=0; i<nbpt; i++)
		(*pt)[i] = bpt[i];

	j = nbpt;
	k = 0;
	for (idx=0; idx<cnt; idx++)
	{
		aspa = (spa[idx+0]+spa[idx+cnt])/2;
		len = fabs(bpt[idx+0] - bpt[idx+cnt]);
		nsg = (len/aspa + 0.5);
		istart = j;
		//right now only uniform spacing is considered

		if (bpt[idx+0] < bpt[idx+cnt])
		{
			for (i=1; i<nsg; i++)
				(*pt)[j++] = bpt[idx+0] + i*aspa;

			for (i=1; i<nsg-1; i++)
			{
				(*pElm)[2*k+0] = istart+i-1;
				(*pElm)[2*k+1] = istart+i;
				k++;
			}

			if(nsg == 1)
			{
				(*pElm)[2*k+0] = idx+0;
				(*pElm)[2*k+1] = idx+cnt;
				k++;
			}
			else
			{
				(*pElm)[2*k+0] = idx+0;
				(*pElm)[2*k+1] = istart;
				k++;
				(*pElm)[2*k+0] = j-1;
				(*pElm)[2*k+1] = idx+cnt;
				k++;
			}
		} 
		else
		{
			for (i=1; i<nsg; i++)
				(*pt)[j++] = bpt[idx+cnt] + i*aspa;

			for (i=1; i<nsg-1; i++)
			{
				(*pElm)[2*k+0] = istart+i-1;
				(*pElm)[2*k+1] = istart+i;
				k++;
			}

			if(nsg == 1)
			{
				(*pElm)[2*k+0] = idx+0;
				(*pElm)[2*k+1] = idx+cnt;
				k++;
			}
			else
			{
				(*pElm)[2*k+0] = idx+cnt;
				(*pElm)[2*k+1] = istart;
				k++;

				(*pElm)[2*k+0] = j-1;
				(*pElm)[2*k+1] = idx+0;
				k++;
			}
		}
	}
}
#endif