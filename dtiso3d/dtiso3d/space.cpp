/* ----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------
 *
 * 三维各向同性Delaunay网格生成器 (版本号：0.3)
 * 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 *
 * 陈建军 中国 浙江大学工程与科学计算研究中心
 * 版权所有	  2005年9月15日
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 10, 26
 * 
 * 联系方式
 *   电话：+86-571-87953165
 *   传真：+86-571-87953167
 *   邮箱：zdchenjj@yahoo.com.cn
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 *
 * ------------------------------------------------------------------------------
 * ------------------------------------------------------------------------------*/

#include <spdlog/spdlog.h> 
 #include "space.h"
#include <stdio.h>
#include <assert.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

/* 清空背景网格对象 free the object of background mesh */
int freeBGMesh(BGMesh *pBGMesh)
{
	if (pBGMesh)
	{
		if (pBGMesh->pNodes)
		{
			free(pBGMesh->pNodes);
			pBGMesh->pNodes = nullptr;
		}
		if (pBGMesh->pElems)
		{
			free(pBGMesh->pElems);
			pBGMesh->pElems = nullptr;
		}
		pBGMesh->nNodes = 0;
		pBGMesh->nElems = 0;

		return 1;
	}

	return 0;
}

/* 清空源对象 free the object of source */
int freeSource(Source *pSource)
{
	if (pSource)
	{
		if (pSource->pPntS)
		{
			free(pSource->pPntS);
			pSource->pPntS = nullptr;
		}
		if (pSource->pLinS)
		{
			free(pSource->pLinS);
			pSource->pLinS = nullptr;
		}
		if (pSource->pTriS)
		{
			free(pSource->pTriS);
			pSource->pTriS = nullptr;
		}
		pSource->nPntS = 0;
		pSource->nLinS = 0;
		pSource->nTriS = 0;

		return 1;
	}

	return 0;
}

/*
/*
 * 读取BA3文件
 * read a .ba3 file
 */
int readBA3(const char *fname, BGMesh *pBGMesh, Source *pSource)
{
	FILE* fp = nullptr;
	INTEGER iTok = 0;
	INTEGER i,j,k;
	char cTok[512];

	assert(fname);
	assert(pBGMesh);
	assert(pSource);

	fp = fopen(fname, "r");
	if (!fp)
	{
		spdlog::info("cannot read file %s\n", fname);
		return 0;
	}

	freeBGMesh(pBGMesh);
	freeSource(pSource);

	fgets(cTok, 512, fp);
	fgets(cTok, 512, fp);
	sscanf(cTok, "%d%d%d%d%d\n", &pBGMesh->nNodes, &pBGMesh->nElems,
		&pSource->nPntS, &pSource->nLinS, &pSource->nTriS);
	
	if (pBGMesh->nNodes > 0)
		pBGMesh->pNodes = (oldMBLNode *)malloc(sizeof(oldMBLNode)*pBGMesh->nNodes);
	else
		pBGMesh->pNodes = nullptr;

	if (pBGMesh->nElems > 0)
		pBGMesh->pElems = (oldElem *)malloc(sizeof(oldElem)*pBGMesh->nElems);
	else
		pBGMesh->pElems = nullptr;

	if (pSource->nPntS > 0)
		pSource->pPntS = (PntSource*)malloc(sizeof(PntSource)*pSource->nPntS);
	else
		pSource->pPntS = nullptr;

	if (pSource->nLinS > 0)
		pSource->pLinS = (LinSource*)malloc(sizeof(LinSource)*pSource->nLinS);
	else
		pSource->pLinS = nullptr;

	if (pSource->nTriS > 0)
		pSource->pTriS = (TriSource*)malloc(sizeof(TriSource)*pSource->nTriS);
	else
		pSource->pTriS = nullptr;

	if ((!pBGMesh->pNodes && pBGMesh->nNodes > 0) || 
		(!pBGMesh->pElems && pBGMesh->nElems > 0) ||
		(!pSource->pPntS && pSource->nPntS > 0) ||
		(!pSource->pLinS && pSource->nLinS > 0) ||
		(!pSource->pTriS && pSource->nTriS > 0))
	{
		fclose(fp);
		freeBGMesh(pBGMesh);
		freeSource(pSource);
		spdlog::info("not enough memory \n");

		return -1;
	}
	
	for (i = 0; i < pBGMesh->nNodes; i++)
	{
		fscanf(fp, "%d", &iTok);
		for (j = 0; j < DIM; j++)
		{
#ifdef _DOUBLE
			fscanf(fp, "%lf", &((pBGMesh->pNodes)[i].pt)[j]);
#else
			fscanf(fp, "%f", &((pBGMesh->pNodes)[i].pt)[j]);
#endif
		}
		fgets(cTok, 512, fp);
		fgets(cTok, 512, fp);
		fgets(cTok, 512, fp);
		fgets(cTok, 512, fp);
	}
	
	
	for (i = 0; i < pBGMesh->nElems; i++)
	{
		fscanf(fp, "%d", &iTok);
		for (j = 0; j <= DIM; j++)
		{
			fscanf(fp, "%d", &((pBGMesh->pElems)[i].form)[j]);
			((pBGMesh->pElems)[i].form)[j]--;
		}
	}
	fgets(cTok,512,fp);

	fgets(cTok,512,fp);
	for (i = 0; i < pSource->nPntS; i++)
	{
		fgets(cTok,512,fp); 
		for (j = 0; j < DIM; j++)
		{
#ifdef _DOUBLE
			fscanf(fp, "%lf", &((pSource->pPntS)[i].cen)[j]);
#else
			fscanf(fp, "%f", &((pSource->pPntS)[i].cen)[j]);
#endif
		}

#ifdef _DOUBLE
		fscanf(fp, "%lf%lf%lf", 
			&((pSource->pPntS)[i].space), 
			&((pSource->pPntS)[i].in_rad),
			&((pSource->pPntS)[i].ou_rad));
#else
		fscanf(fp, "%f%f%f", 
			&((pSource->pPntS)[i].space), 
			&((pSource->pPntS)[i].in_rad),
			&((pSource->pPntS)[i].ou_rad));
#endif
		if ((pSource->pPntS)[i].in_rad > (pSource->pPntS)[i].ou_rad)
		{
			fclose(fp);
			free(pBGMesh);
			free(pSource);
			spdlog::info("Inner radius of a point source is greater than its outer radius\n");

			return -1;
		}

		fgets(cTok,512,fp); 
	}

	fgets(cTok,512,fp);
	for (i = 0; i< pSource->nLinS; i++)
	{
		fgets(cTok,512,fp);
		for (j = 0; j < 2; j++)
		{
			for (k = 0; k < DIM; k++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &(((pSource->pLinS)[i].src)[j].cen)[k]);
#else
				fscanf(fp, "%f", &(((pSource->pLinS)[i].src)[j].cen)[k]);
#endif
			}

#ifdef _DOUBLE
			fscanf(fp, "%lf%lf%lf", 
				&(((pSource->pLinS)[i].src)[j].space), 
				&(((pSource->pLinS)[i].src)[j].in_rad),
				&(((pSource->pLinS)[i].src)[j].ou_rad));
#else
			fscanf(fp, "%f%f%f", 
				&(((pSource->pLinS)[i].src)[j].space), 
				&(((pSource->pLinS)[i].src)[j].in_rad),
				&(((pSource->pLinS)[i].src)[j].ou_rad));
#endif
			if (((pSource->pLinS)[i].src)[j].in_rad > ((pSource->pLinS)[i].src)[j].ou_rad)
			{
				fclose(fp);
				free(pBGMesh);
				free(pSource);
				spdlog::info("Inner radius of a point source of a line source is greater than its outer radius\n");

				return -1;
			}

			fgets(cTok,512,fp);
		}
	}
	
	fgets(cTok,512,fp);
	for (i = 0; i< pSource->nTriS; i++)
	{
		fgets(cTok,512,fp);
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < DIM; k++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &(((pSource->pTriS)[i].src)[j].cen)[k]);
#else
				fscanf(fp, "%f", &(((pSource->pTriS)[i].src)[j].cen)[k]);
#endif
			}

#ifdef _DOUBLE
			fscanf(fp, "%lf%lf%lf", 
				&(((pSource->pTriS)[i].src)[j].space), 
				&(((pSource->pTriS)[i].src)[j].in_rad),
				&(((pSource->pTriS)[i].src)[j].ou_rad));
#else
			fscanf(fp, "%f%f%f", 
				&(((pSource->pTriS)[i].src)[j].space), 
				&(((pSource->pTriS)[i].src)[j].in_rad),
				&(((pSource->pTriS)[i].src)[j].ou_rad));
#endif
			if (((pSource->pTriS)[i].src)[j].in_rad > ((pSource->pTriS)[i].src)[j].ou_rad)
			{
				fclose(fp);
				free(pBGMesh);
				free(pSource);
				spdlog::info("Inner radius of a point source of a Triangle source is greater than its outer radius\n");

				return -1;
			}

			fgets(cTok,512,fp);
		}
	}

	fclose(fp);
	return 1;
}

/*
 * 写BA3文件
 * write a .ba3 file
 */
int writeBA3(const char *fname, BGMesh *pBGMesh, Source *pSource)
{
	FILE* fp = nullptr;
	INTEGER iTok = 0;
	INTEGER i,j,k;

	assert(fname);
	assert(pBGMesh);
	assert(pSource);

	fp = fopen(fname, "w");
	if (!fp)
	{
		spdlog::info("cannot read file %s\n", fname);
		return 0;
	}

	fprintf(fp, "* background mesh\n");

	fprintf(fp, "%d %d %d %d %d\n", pBGMesh->nNodes, pBGMesh->nElems,
		pSource->nPntS, pSource->nLinS, pSource->nTriS);
	
	for (i = 0; i < pBGMesh->nNodes; i++)
	{
		fprintf(fp, "%d\t", i + 1);
		for (j = 0; j < DIM; j++)
		{
			fprintf(fp, "%f\t", ((pBGMesh->pNodes)[i].pt)[j]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "1.0 0.0 0.0 2\n");
		fprintf(fp, "0.0 1.0 0.0 2\n");
		fprintf(fp, "0.0 0.0 1.0 2\n");
	}
	
	
	for (i = 0; i < pBGMesh->nElems; i++)
	{
		fprintf(fp, "%d\t", i + 1);
		for (j = 0; j <= DIM; j++)
		{
			fprintf(fp, "%d\t", ((pBGMesh->pElems)[i].form)[j] + 1);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "* points\n");
	for (i = 0; i < pSource->nPntS; i++)
	{
		fprintf(fp, "%d\n", i + 1);
		for (j = 0; j < DIM; j++)
		{
			fprintf(fp, "%f\t", ((pSource->pPntS)[i].cen)[j]);
		}

		fprintf(fp, "%f\t%f\t%f\n", 
			((pSource->pPntS)[i].space), 
			((pSource->pPntS)[i].in_rad),
			((pSource->pPntS)[i].ou_rad));
	}

	fprintf(fp, "* lines\n");
	for (i = 0; i< pSource->nLinS; i++)
	{
		fprintf(fp, "%d\n", i + 1);
		for (j = 0; j < 2; j++)
		{
			for (k = 0; k < DIM; k++)
			{
				fprintf(fp, "%f\t", (((pSource->pLinS)[i].src)[j].cen)[k]);
			}
			fprintf(fp, "%f\t%f\t%f\n", 
				(((pSource->pLinS)[i].src)[j].space), 
				(((pSource->pLinS)[i].src)[j].in_rad),
				(((pSource->pLinS)[i].src)[j].ou_rad));
		}
	}

	fprintf(fp, "* Triangles\n");
	for (i = 0; i< pSource->nTriS; i++)
	{
		fprintf(fp, "%d\n", i + 1);
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < DIM; k++)
			{
				fprintf(fp, "%f\t", (((pSource->pTriS)[i].src)[j].cen)[k]);
			}
			fprintf(fp, "%f\t%f\t%f\n", 
				(((pSource->pTriS)[i].src)[j].space), 
				(((pSource->pTriS)[i].src)[j].in_rad),
				(((pSource->pTriS)[i].src)[j].ou_rad));
		}
	}

	fclose(fp);

	return 1;
}

