#include <stdio.h>
#include <spdlog/spdlog.h> 
 #include "IntIntMap.h"
#include "iso3d_utility.h"
#include "locsmoothing.h"
#include "optdihangle.h"
#include "skirtpolylist.h"
#include "myvector.h"

LocSmoother::LocSmoother(SPRImpl* aSprimpl) 
	:sprimpl(aSprimpl)
{
	g_pLocSmoothingElems = nullptr;
	g_pLocSmoothingNodes = nullptr;
	g_pLocSmoothingQuals = nullptr;
	g_nLocSmoothingElems = 0;
	g_nLocSmoothingNodes = 0;
	g_nLocSmoothingCentNode = 0;
}

/* 初始化函数 */
int LocSmoother::setMesh_LocSmoothing(divide *divi, REAL *quals)
{
	int i, j, k, l, m;
	INTEGER iia, iib, iic, iid;
	INTEGER jja, jjb, jjc, jjd; 
	INTEGER iLoc;
	INTEGER eGlob1, eGlob2;
	INTEGER iIdx;
	int nRet = 0;

	g_pLocSmoothingNodes = (oldMBLNode*)malloc(sizeof(oldMBLNode)*sprimpl->num_vertices);
	g_pLocSmoothingElems = (oldElem*)malloc(sizeof(oldElem)*divi->num);
	g_pLocSmoothingQuals = (TetraElemQual*)malloc(sizeof(TetraElemQual)*divi->num);
	if (!g_pLocSmoothingNodes || !g_pLocSmoothingElems || !g_pLocSmoothingQuals)
	{
		nRet = -1;
		goto FAIL;
	}
	g_nLocSmoothingElems = divi->num;
	g_nLocSmoothingNodes = sprimpl->num_vertices;
	memset(g_pLocSmoothingNodes, 0, sizeof(oldMBLNode)*sprimpl->num_vertices);
	memset(g_pLocSmoothingElems, 0, sizeof(oldElem)*divi->num);
	memset(g_pLocSmoothingQuals, 0, sizeof(TetraElemQual)*divi->num);

	/* 将Node信息拷贝到网格中 */
	for (i = 0; i < sprimpl->num_vertices; i++)
	{
		memcpy(g_pLocSmoothingNodes[i].pt, sprimpl->vertices[i], sizeof(APOINT));
	}

	/* 将结果加入到单元数组中 */
	for (i = 0; i < divi->num; i++)
	{
		/* 不能利用中间空位 */
		iLoc = i;//m_nElems++;

		memset(&g_pLocSmoothingElems[iLoc], 0, sizeof(oldElem));

		for (j = 0; j <= DIM; j++)
		{
			iIdx = divi->te[i][j];
			g_pLocSmoothingElems[iLoc].form[j] = iIdx;
			g_pLocSmoothingElems[iLoc].neig[j] = NULL_NEIG;
			g_pLocSmoothingNodes[iIdx].iReserved = iLoc;
		}
	}

	/* 更新内部的相邻关系 */
	for (i = 0; i < divi->num; i++)
		for (j = i + 1; j < divi->num; j++)
		{
			eGlob1 = i;
			eGlob2 = j;

			for (k = 0; k <= 3; k++)
			{
				if (g_pLocSmoothingElems[eGlob1].neig[k] >= 0)
					continue;

				DNC(k, iia, iib, iic, iid)
					iib = divi->te[i][iib];
				iic = divi->te[i][iic];
				iid = divi->te[i][iid];

				for (m = 0; m <= 3; m++)
				{
					if (g_pLocSmoothingElems[eGlob2].neig[m] >= 0)
						continue;

					DNC(m, jja, jjb, jjc, jjd)
					jjb = divi->te[j][jjb];
					jjc = divi->te[j][jjc];
					jjd = divi->te[j][jjd];

					if ((iib == jjb || iib == jjc || iib == jjd) &&
						(iic == jjb || iic == jjc || iic == jjd) &&
						(iid == jjb || iid == jjc || iid == jjd))
					{
						g_pLocSmoothingElems[eGlob1].neig[k] = eGlob2;
						g_pLocSmoothingElems[eGlob2].neig[m] = eGlob1;

						break;
					}
				}
			}
		}

	/* 记录单元质量信息 */
	for (i = 0; i < divi->num; i++)
	{
		g_pLocSmoothingQuals[i].sinValue = quals[i];
	}
//	tetEleQual.iReserved = 0;
//	::setAngleCode(tetEleQual.iReserved , iEdge);
//	::setAcute(tetEleQual.iReserved , bAcute ? 1 : 0);


	goto SUCCESS;

FAIL:
	free_LocSmoothing();

SUCCESS:

	return nRet;
}

int LocSmoother::setMesh_LocSmoothing(INTEGER elems[][4], int numOfElems, REAL *quals)
{
	int i, j, k, l, m;
	INTEGER iia, iib, iic, iid;
	INTEGER jja, jjb, jjc, jjd; 
	INTEGER iLoc;
	INTEGER eGlob1, eGlob2;
	INTEGER iIdx;
	int nRet = 0;
#if 1
	int cellIdx, faceNdIdx1, faceNdIdx2, faceNdIdx3, minFaceNdIdx;
	int faceIt, faceTail, faceData, faceElem, faceCode, iNode;
	int nodeFacesData[MAX_SPHERE_SIZE], nodeFacesHash[MAX_SPHERE_SIZE], numOfFaces = 0;
#endif

	g_pLocSmoothingNodes = (oldMBLNode*)malloc(sizeof(oldMBLNode)*sprimpl->num_vertices);
	g_pLocSmoothingElems = (oldElem*)malloc(sizeof(oldElem)*numOfElems);
	g_pLocSmoothingQuals = (TetraElemQual*)malloc(sizeof(TetraElemQual)*numOfElems);
	if (!g_pLocSmoothingNodes || !g_pLocSmoothingElems || !g_pLocSmoothingQuals)
	{
		nRet = -1;
		goto FAIL;
	}
	g_nLocSmoothingElems = numOfElems;
	g_nLocSmoothingNodes = sprimpl->num_vertices;
	memset(g_pLocSmoothingNodes, 0, sizeof(oldMBLNode)*sprimpl->num_vertices);
	memset(g_pLocSmoothingElems, 0, sizeof(oldElem)*numOfElems);
	memset(g_pLocSmoothingQuals, 0, sizeof(TetraElemQual)*numOfElems);

	/* 将Node信息拷贝到网格中 */
	for (i = 0; i < sprimpl->num_vertices; i++)
	{
		memcpy(g_pLocSmoothingNodes[i].pt, sprimpl->vertices[i], sizeof(APOINT));
	}

#if 0
	/* 将结果加入到单元数组中 */
	for (i = 0; i < numOfElems; i++)
	{
		/* 不能利用中间空位 */
		iLoc = i;//m_nElems++;
		for (j = 0; j <= DIM; j++)
		{
			iIdx = elems[i][j];
			g_pLocSmoothingElems[iLoc].form[j] = iIdx;
			g_pLocSmoothingElems[iLoc].neig[j] = NULL_NEIG;
			g_pLocSmoothingNodes[iIdx].iReserved = iLoc;
		}
	}

	/* 更新内部的相邻关系 */
	for (i = 0; i < numOfElems; i++)
		for (j = i + 1; j < numOfElems; j++)
		{
			eGlob1 = i;
			eGlob2 = j;

			for (k = 0; k <= 3; k++)
			{
				if (g_pLocSmoothingElems[eGlob1].neig[k] >= 0)
					continue;

				DNC(k, iia, iib, iic, iid)
					iib = elems[i][iib];
				iic = elems[i][iic];
				iid = elems[i][iid];

				for (m = 0; m <= 3; m++)
				{
					if (g_pLocSmoothingElems[eGlob2].neig[m] >= 0)
						continue;

					DNC(m, jja, jjb, jjc, jjd)
					jjb = elems[j][jjb];
					jjc = elems[j][jjc];
					jjd = elems[j][jjd];

					if ((iib == jjb || iib == jjc || iib == jjd) &&
						(iic == jjb || iic == jjc || iic == jjd) &&
						(iid == jjb || iid == jjc || iid == jjd))
					{
						g_pLocSmoothingElems[eGlob1].neig[k] = eGlob2;
						g_pLocSmoothingElems[eGlob2].neig[m] = eGlob1;

						break;
					}
				}
			}
		}
#else
	/* 通过哈希表，更快地更新相邻关系 */
	for (cellIdx = 0; cellIdx < numOfElems; cellIdx++)
		for (i = 0; i <= DIM; i++) 
		{
			faceNdIdx1 = elems[cellIdx][(i+1)%(DIM+1)];
			faceNdIdx2 = elems[cellIdx][(i+2)%(DIM+1)];
			faceNdIdx3 = elems[cellIdx][(i+3)%(DIM+1)];
			minFaceNdIdx = faceNdIdx1 < faceNdIdx2 ? faceNdIdx1 : faceNdIdx2;
			minFaceNdIdx = minFaceNdIdx < faceNdIdx3 ? minFaceNdIdx : faceNdIdx3;

			/* 求得包含节点的面片 */
			faceTail = faceIt = g_pLocSmoothingNodes[minFaceNdIdx].iReserved;
			while (faceIt > 0) 
			{
				faceData = nodeFacesData[faceIt];
				faceElem = faceData >> 2;
				if (faceElem != cellIdx)
				{
					faceCode = faceData & 0x3;
					/* 确定faceElem是否包含这个面 */
					for (m = 1; m <= DIM; m++)
					{
						iNode = elems[faceElem][(faceCode + m)%(DIM+1)];
						if (iNode != faceNdIdx1 && iNode != faceNdIdx2 && iNode != faceNdIdx3)
							break;
					}
					if (m > DIM)
					{/* get the neighbour */
						g_pLocSmoothingElems[cellIdx].neig[i] = faceElem;
						assert(g_pLocSmoothingElems[faceElem].neig[faceCode] == NULL_NEIG);
						g_pLocSmoothingElems[faceElem].neig[faceCode] = cellIdx;

						/* detach the face */
						faceIt == g_pLocSmoothingNodes[minFaceNdIdx].iReserved ?
							g_pLocSmoothingNodes[minFaceNdIdx].iReserved = nodeFacesHash[faceIt] :
							nodeFacesHash[faceTail] = nodeFacesHash[faceIt];
						break;
					}
				}
				faceTail = faceIt;
				faceIt = nodeFacesHash[faceIt];
			}
			if (faceIt <= 0)
			{/* 没有找到 */
				numOfFaces++;
				if (numOfFaces >= MAX_SPHERE_SIZE)
				{
					spdlog::info("Error in setMesh_LocSmoothing(...).\n");
					spdlog::info("Please increase MAX_SPHERE_SIZE:{}.\n", MAX_SPHERE_SIZE);
					exit(1);
				}
				g_pLocSmoothingElems[cellIdx].neig[i] = NULL_NEIG;
				nodeFacesData[numOfFaces] = (cellIdx << 2) | i;
				nodeFacesHash[numOfFaces] = 0;
				if (faceTail <= 0)
				{/* 第一个 */
					g_pLocSmoothingNodes[minFaceNdIdx].iReserved = numOfFaces;
				}
				else
				{	
					nodeFacesHash[faceTail] = numOfFaces;
				}
			}
		}

	/* 将结果加入到单元数组中 */
	for (i = 0; i < numOfElems; i++)
	{
		/* 不能利用中间空位 */
		iLoc = i;//m_nElems++;
		for (j = 0; j <= DIM; j++)
		{
			iIdx = elems[i][j];
			g_pLocSmoothingElems[iLoc].form[j] = iIdx;
			g_pLocSmoothingNodes[iIdx].iReserved = iLoc;
		}
	}
#endif
	/* 记录单元质量信息 */
	for (i = 0; i < numOfElems; i++)
	{
		g_pLocSmoothingQuals[i].sinValue = quals[i];
	}
//	tetEleQual.iReserved = 0;
//	::setAngleCode(tetEleQual.iReserved , iEdge);
//	::setAcute(tetEleQual.iReserved , bAcute ? 1 : 0);


	goto SUCCESS;

FAIL:
	free_LocSmoothing();

SUCCESS:

	return nRet;
}
int LocSmoother::setCentNode_LocSmoothing(int iCentNod)
{
	assert(iCentNod >= 0 && iCentNod < g_nLocSmoothingNodes);
	g_nLocSmoothingCentNode = iCentNod;
	return 0;
}

int LocSmoother::getMeshQuality_LocSmoothing(double *elemQuality, int *size)
{
	int i, j, m;
	REAL *verts[4];
	REAL thisqual;

	assert(size);
	*size = g_nLocSmoothingElems;

	if (elemQuality)
	{
		for (i = 0; i < *size; i++)
		{
			elemQuality[i] = g_pLocSmoothingQuals[i].sinValue;
#if 0
			for (m = 0; m <= DIM; m++)
				verts[m] = g_pLocSmoothingNodes[g_pLocSmoothingElems[i].form[m]].pt;

//			thisqual = tetquality(verts[0], verts[1], verts[3], verts[2], QUALWARPEDMINSINE);
			thisqual = tetquality(verts[0], verts[1], verts[3], verts[2], QUALVLRMS3RATIO);
			
			if ((elemQuality[i] == 0.0 && fabs(thisqual) >= 1.0e-3) || 
			(elemQuality[i] != 0.0 && fabs((thisqual - elemQuality[i])/elemQuality[i]) >= 1.0e-3))
			{
				spdlog::info("In  ::getMeshQuality_LocSmoothing(...)\n");
				printf("Some errors may happen when evaluating the element quality: tetEleQual.sinValue = %f; worstqual = %f.\n", 
				elemQuality[i], thisqual);
				spdlog::info("Loc. information.\n");
				spdlog::info("Elem. Idx: {}.\n", i);
				printf("Node. Idx: %d, %d, %d, %d.\n", g_pLocSmoothingElems[i].form[0], g_pLocSmoothingElems[i].form[1],
					g_pLocSmoothingElems[i].form[2], g_pLocSmoothingElems[i].form[3]);
				spdlog::info("Node. Coords.:\n");
				for (m = 0; m <= DIM; m++)
				{
					spdlog::info("#Node {}:\n%24.20g\n%24.20g\n%24.20g\n", m+1, verts[m][0], verts[m][1], verts[m][2]);
				}
				exit(1);
			}
#endif
		}
	}

	return *size;
}

 void LocSmoother::free_LocSmoothing()
 {
	if (g_pLocSmoothingElems)
	{
		free(g_pLocSmoothingElems);
		g_pLocSmoothingElems = nullptr;
	}
	if (g_pLocSmoothingNodes)
	{
		free(g_pLocSmoothingNodes);
		g_pLocSmoothingNodes = nullptr;
	}
	if (g_pLocSmoothingQuals)
	{
		free(g_pLocSmoothingQuals);
		g_pLocSmoothingQuals = nullptr;
	}

	g_nLocSmoothingElems = 0;
	g_nLocSmoothingNodes = 0;
	g_nLocSmoothingCentNode = 0;
 }

/* 网格质量优化代码，以下函数将网格优化定义在一个局部网格上 */
/* given a set of tets incident to a vertex, and the quality
   of the worst quality function that varies with that vertex,
   compute the active set A of quality functions very near
   the worst */
void LocSmoother::getactiveset(Sphere sph,
                  int nSph,
				  REAL quals[],
				  REAL qualgrads[][3],
                  REAL activegrads[][3],
                  int *numactive,
                  REAL worstqual,
				  int qualmeasure)
{
    int i, j;
    const REAL ACTIVESETFACTOR = 1.03;
	int QUAL_VALUES_PER_ELE = 1;

	if (qualmeasure == QUALMINSINE || qualmeasure == QUALWARPEDMINSINE)
		QUAL_VALUES_PER_ELE = 6;

    /* reset number of active gradients to zero */
    *numactive = 0;
    
    /* if we are including surface quadrics, give them a chance to
       enter the active set */
    assert(nSph > 0);

	worstqual = *std::min_element(quals, quals+ nSph * QUAL_VALUES_PER_ELE);

	for (i = 0; i < nSph*QUAL_VALUES_PER_ELE; i++)
	{
		if (quals[i] < worstqual * ACTIVESETFACTOR)
		{
			memcpy(activegrads[*numactive], qualgrads[i], sizeof(REAL) * 3);
			(*numactive)++;
		}
	}


#if 0
    for (i = 0; i < nSph; i++)
    {
		for (j=0; j<QUAL_VALUES_PER_ELE; j++)
        {
			if (quals[6*QUAL_VALUES_PER_ELE+j] < worstqual * ACTIVESETFACTOR)
			{
				memcpy(activegrads[*numactive], qualgrads[6*i+j], sizeof(REAL)*3);
				(*numactive)++;
			}
         }
    }
#else

	for (i = 0; i < nSph*QUAL_VALUES_PER_ELE; i++)
    {
		if (quals[i] < worstqual * ACTIVESETFACTOR)
		{
			memcpy(activegrads[*numactive], qualgrads[i], sizeof(REAL)*3);
			(*numactive)++;
		}
	}
#endif
    
	if (*numactive == 0)
	{
		spdlog::info("didn't find any active quality functions.\n");
		spdlog::info("nSph = {}.\n", nSph);
		
	}


		//if (nSph > 0)
		//{
		//	int iCentNode = 0;
		//	for (i = 0; i <= DIM; i++)
		//	{
		//		iCentNode = g_pLocSmoothingElems[sph[0]].form[i];
		//		for (j = 1; j < nSph; j++)
		//			if (g_pLocSmoothingElems[sph[0]].form[0] != iCentNode &&
		//				g_pLocSmoothingElems[sph[0]].form[1] != iCentNode &&
		//				g_pLocSmoothingElems[sph[0]].form[2] != iCentNode &&
		//				g_pLocSmoothingElems[sph[0]].form[3] != iCentNode)
		//			{
		//				break;
		//			}
		//		if (j >= nSph)
		//		{
		//			spdlog::info("Smoothing node {}.\n", iCentNode);
		//			break;
		//		}
		//	}
		//}
  //  }
    /* we must have at least found the worst */
	/* assert(*numactive > 0); */
}

/* returns the basis B of S union M. S is a set of points known
   to be in the basis */
void LocSmoother::findbasis(REAL S[][3],
               REAL M[][3],
               REAL B[][3],
               int sizeS,
               int sizeM,
               int *sizeB)
{
    REAL p[3],q[3],r[3],s[3];        /* vertices */
    REAL s1[3], t1[3], d1[3], d2[3]; /* temporary vertices */
    REAL origin[3] = {0.0, 0.0, 0.0};
    REAL localS[4][3];               /* for passing to recursive calls */
    REAL localM[MAX_SPHERE_SIZE][3]; /* for passing to recursive colls */
    int m;

    assert(sizeM > 0);
    
   
    /* we assume that M was passed to us shuffled, so that taking
       the last element is equivalent to removing a random one. */
	memcpy(p, M[sizeM-1], sizeof(REAL)*3);

    /* if M has only one element */
    if (sizeM == 1)
    {
        /* and S has no elements */
        if (sizeS == 0)
        {           
            /* then the single element in M must be the
               entire basis, just send back M */
            /*pointarraycopy(M, B, 1); */
			memcpy(B, M, sizeof(REAL)*3);
            *sizeB = 1;
            return;
        }
        
        
        /* otherwise, because we assume the last element
           we just removed is not part of the basis, assign
           the basis to be the elements of S */
        /*pointarraycopy(S, B, sizeS); */
		memcpy(B, S, sizeof(REAL)*3*sizeS);
        *sizeB = sizeS;
    }
    /* M has more than one element. Throw one out (p), and look for the
       basis assuming p is not part of it. */
    else
    {   
        /* make a new copy of M minus the last element */
       /* pointarraycopy(M, localM, sizeM-1);
        pointarraycopy(S, localS, sizeS); */
		memcpy(localM, M, sizeof(REAL)*3*(sizeM-1));
		memcpy(localS, S, sizeof(REAL)*3*sizeS);
        findbasis(localS, localM, B, sizeS, sizeM-1, sizeB);
    }
    
    /* now the we have determined the basis without p, we need to 
       go back and check whether p actually is part of the basis. */
    
    switch (*sizeB)
    {
        /* if the returned basis has just one point q, we just need to check
           whether p is closer to the origin than q */
        case 1:
            /* fetch the actual coordinates from the mesh */
            memcpy(q, B[0], sizeof(REAL)*3);
			
            /* compute the vector from q to p */
            for (m = 0; m < DIM; m++)
				d1[m] = p[m] - q[m];
            
            /* check the sign of the dot product. >=0 means p doesn't
               improve the basis */
            if (q[0]*d1[0] + q[1]*d1[1] + q[2]*d1[2] >= 0)
                /* this is a good B, send it back!*/
                return;
            break;
            
        /* check whether p improves the basis using math I don't understand */        
        case 2:
            /* fetch coordinates from the mesh */
			memcpy(q, B[0], sizeof(REAL)*3);
			memcpy(r, B[1], sizeof(REAL)*3);

            /* compute vector s from r to p & vector t from r to q */
            for (m = 0; m < DIM; m++)
			{
				s1[m] = p[m] - r[m];
				t1[m] = q[m] - r[m];
			}
            /* now a couple of cross products */
            /* cross(s1, t1, d1); */
			d1[0] = s1[1]*t1[2] - s1[2]*t1[1];
			d1[1] = s1[2]*t1[0] - s1[0]*t1[2];
			d1[2] = s1[0]*t1[1] - s1[1]*t1[0];
			/* cross(r, t1, d2); */
			d2[0] = r[1]*t1[2] - r[2]*t1[1];
			d2[1] = r[2]*t1[0] - r[0]*t1[2];
			d2[2] = r[0]*t1[1] - r[1]*t1[0];

            /* can p improve the basis? */
            if (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2] >= 0)
                /* nope! send back B as is. */
                return;
                        
            break;
        case 3:
            /* fetch coordinates from the mesh */
			memcpy(q, B[0], sizeof(REAL)*3);
			memcpy(r, B[1], sizeof(REAL)*3);
			memcpy(s, B[2], sizeof(REAL)*3);
            
            /* does p improve the basis? */
         //   if (orient3d(&behave, p, q, r, s) * orient3d(&behave, origin, q, r, s) <= 0)
			if (GEOM_FUNC::orient3d(p, q, r, s) * GEOM_FUNC::orient3d(origin, q, r, s) <= 0)
            {
                /* nope! send back B as is. */
                return;
            }
            break;
        default:
            /* B has size of 4, and any basis of this size is optimal */
            return;
            break;
    }
    
    /* if we have made it this far, we know that p actually is a part of
       any basis of S union M */
       
    /* if p was the last element of M, or if S already has three other basis
       points, we're done and just need to send back S union p. */
    if ((sizeM == 1) || (sizeS == 3))
    {
        /* p in basis and it's the last possible point, return S U p\n") */
        /* the final size of B is the size of S + 1 */
        *sizeB = sizeS + 1;
        /* copy S into B */
        //pointarraycopy(S, B, sizeS);
		memcpy(B, S, sizeof(REAL)*3*sizeS);
        /* and add p at the end */
		//vcopy(p, B[*sizeB - 1]);
		memcpy(B[*sizeB - 1], p, sizeof(REAL)*3);
        return;
    }
    else
    {/* there may be more basis points to find! move p from M to the known
       basis point set S, and go again */
        /* p in basis, more points to check. Moving p from M to S\n");
            spdlog::info("here's how everything looked before the move:\n");
            printbasisarrays(S,M,B,sizeS,sizeM,sizeB,p); */

        /* create the new S */
        //pointarraycopy(S, localS, sizeS);
		memcpy(localS, S, sizeof(REAL)*3*sizeS);
        /* add p to the end of it */
        //vcopy(p, localS[sizeS]);
        memcpy(localS[sizeS], p, sizeof(REAL)*3);

        /* create the new M, leaving off the last element */
        //pointarraycopy(M, localM, sizeM-1);
        memcpy(localM, M, sizeof(REAL)*3*(sizeM-1));

        /* find any basis points remaining in M */
        findbasis(localS, localM, B, sizeS+1, sizeM-1, sizeB);
        
        return;
    }        
}

/* finds the point on the convex hull of P nearest
   the origin */
void LocSmoother::minconvexhullpoint(REAL P[][3],
                        int sizeP,
                        REAL nearest[])
{
    REAL B[4][3];                /* the basis for the convex hull point */
    int sizeB=0;                   /* size of the basis */
    REAL empty[4][3];            /* empty set for known basis */
    REAL *p, *q, *r;             /* basis points */
    REAL pmq[3];                 /* p minus q */
    REAL c, d, l;                /* scalar factors */
    REAL s[3], t[3], s2[3], t2[3];/* temporary points */
    REAL sxt[3], sxr[3], rxt[3], temp[3], sxt_len_sq, sxr_len_sq, rxt_len_sq; /* temporary cross products */
    int m;
	const REAL NEARESTMIN = 1.0e-13; /* if closest point computation has a factor smaller than this, make it the origin */
    
    assert(sizeP > 0);
    
    /* find a basis for the minimum point on the convex hull */
    findbasis(empty, P, B, 0, sizeP, &sizeB);
        
    switch(sizeB)
    {
        /* if the basis is just a single point, return that point */
        case 1:
            memcpy(nearest, B[0], sizeof(REAL)*3);
            return;
        /* for two points, find the closest point to the origin on
           the line between the two points */
        case 2:
            p = B[0];
            q = B[1];
			for (m = 0; m < DIM; m++)
				pmq[m] = p[m] - q[m];
    
            /*
              nearest = q - dot(q,p-q)/(length(p-q)^2) * (p-q)
            */
            l = sqrt(pmq[0]*pmq[0] + pmq[1]*pmq[1] + pmq[2]*pmq[2]);
            
            /* if these points are the same, just return one of them */
            if (l == 0.0)
            {
				memcpy(nearest, B[0], sizeof(REAL)*3);
                return;
            }
            
            c = (q[0]*pmq[0] + q[1]*pmq[1] + q[2]*pmq[2])/(l*l);//dot(q,pmq) / (l * l);
			for (m = 0; m < DIM; m++)
			{
				nearest[m] = c * pmq[m];
				nearest[m] = q[m] - nearest[m];
			}
            /*vscale(c,pmq,nearest);
              vsub(q,nearest,nearest); */
            
            return; 
        /* for three points, find the point closest to the origin
           on the triangle that the they form */
        case 3:
            p = B[0];
            q = B[1];
            r = B[2];
            
            /*vsub(p,r,s);
            vsub(q,r,t);*/
			for (m = 0; m < DIM; m++)
			{
				s[m] = p[m] - r[m];
				t[m] = q[m] - r[m];
			}
			/* cross(s,t,sxt); */
            sxt[0] = s[1]*t[2] - s[2]*t[1];
			sxt[1] = s[2]*t[0] - s[0]*t[2];
			sxt[2] = s[0]*t[1] - s[1]*t[0];
			/* cross(r,t,rxt); */
			rxt[0] = r[1]*t[2] - r[2]*t[1];
			rxt[1] = r[2]*t[0] - r[0]*t[2];
			rxt[2] = r[0]*t[1] - r[1]*t[0];
            /*cross(s,r,sxr);*/
            sxr[0] = s[1]*r[2] - s[2]*r[1];
			sxr[1] = s[2]*r[0] - s[0]*r[2];
			sxr[2] = s[0]*r[1] - s[1]*r[0];

            /* if any of these cross products is really tiny, give up
               and return the origin */
            if (sqrt(sxt_len_sq = sxt[0]*sxt[0] + sxt[1]*sxt[1] + sxt[2]*sxt[2]) < NEARESTMIN ||
				sqrt(rxt_len_sq = rxt[0]*rxt[0] + rxt[1]*rxt[1] + rxt[2]*rxt[2]) < NEARESTMIN ||
				sqrt(sxr_len_sq = sxr[0]*sxr[0] + sxr[1]*sxr[1] + sxr[2]*sxr[2]) < NEARESTMIN)
            {
                nearest[0] = nearest[1] = nearest[2] = 0.0;
                return;
            }
            
           /* c = dot(sxt,rxt) / dot(sxt,sxt);
              d = dot(sxt,sxr) / dot(sxt,sxt); */
			c = (sxt[0]*rxt[0] + sxt[1]*rxt[1] +sxt[2]*rxt[2]) / sxt_len_sq;
            d = (sxt[0]*sxr[0] + sxt[1]*sxr[1] +sxt[2]*sxr[2]) / sxt_len_sq;

			/* 
			 vscale(c,s,s2);
             vscale(d,t,t2);
			 vsub(r,s2,temp);
             vsub(temp,t2,nearest);
			 */
			for (m = 0; m < DIM; m++)
			{
				s2[m] = c * s[m];
				t2[m] = d * t[m];
				temp[m] = r[m] - s2[m];
				nearest[m] = temp[m] - t2[m];
			}
           
            return;
            /* if the basis has four points, they must enclose the origin
           so just return the origin. */
        case 4:
            nearest[0] = 0.0;
            nearest[1] = 0.0;
            nearest[2] = 0.0;
            return;
        default:
            printf("Error, basis size %d is bogus, dying\n",sizeB);
            exit(1);
            break;
    }
}

/* for our initial step size, we use the distance
   to the next intersection with another quality funtion.
   this is the point at which the other quality function
   becomes the worst. we use a single-term taylor expansion
   to approximate all of the quality functions as lines,
   so we'll have to do a line search to find our ultimate
   step size. */
REAL LocSmoother::getinitialalpha(INTEGER iNod,
                     Sphere sph,
					 int nSph,
					 REAL quals[],
				     REAL qualgrads[][3],
                     REAL d[3],
                     REAL r,
                     REAL worstqual,
					 int qualmeasure)
{
    int i,j;
	const REAL HUGEFLOAT = 1.0e100;
	const REAL RATEEPSILON = 1.0e-6;
    REAL alpha = HUGEFLOAT;
	/* rate must be worse by this much to reset alpha */
	
    REAL newalpha;
    REAL rate;
    
    /* if we are including surface quadrics add check for it */
    assert(nSph > 0);
    /* if (improvebehave.usequadrics && hasquadric(vtx) && ((struct vertextype *) arraypoolforcelookup(&vertexinfo, vtx))->kind != FREEVERTEX) */
    
    for (i = 0; i < nSph; i++)
    {
		switch (qualmeasure)
        {
		case QUALWARPEDMINSINE:
        case QUALMINSINE:
			for (j=0; j<6; j++)
			{
				/* if this function improves more slowly
					than any in the active set, then it might
					end up as the objective. */
				//rate = dot(d,incidenttets[i].sinegrad[j]);
				rate = d[0]*qualgrads[6*i + j][0] + d[1]*qualgrads[6*i + j][1] + d[2]*qualgrads[6*i + j][2]; 
				if (rate + RATEEPSILON < r)
				{
					/* compute the approximation of when this
						function will become the objective */
					newalpha = (quals[6*i + j] - worstqual) / (r - rate);
                
					/* if this is smaller than our current step size,
						use it for the step size */
					if (newalpha < alpha)
					{
						alpha = newalpha;
					}
				}
			}
			break;
		case QUALVLRMS3RATIO:
			 /* if this function improves more slowly
                   than any in the active set, then it might
                   end up as the objective. */
//                rate = dot(d,incidenttets[i].vlrms3rgrad);
			rate = d[0]*qualgrads[i][0] + d[1]*qualgrads[i][1] + d[2]*qualgrads[i][2]; 
            if (rate + RATEEPSILON < r)
            {
                /* compute the approximation of when this
                    function will become the objective */
                newalpha = (quals[i] - worstqual) / (r - rate);

                /* if this is smaller than our current step size,
                    use it for the step size */
                if (newalpha < alpha)
                {
                    alpha = newalpha;
                }
            }
            break;
		default:
			spdlog::info("Unsupported quality measure: {}.\n", qualmeasure);
			exit(1);
			break;
        }
    }
   
    if (alpha < 0.0)
		alpha = 0.0;

    return alpha;
}

/* given two values a and b and their gradients, compute the 
   gradient of their product grad(a*b) */
void LocSmoother::gradproduct(REAL a,
                 REAL b, 
                 REAL grada[3],
                 REAL gradb[3],
                 REAL prod[3])
{
    prod[0] = grada[0] * b + gradb[0] * a;
    prod[1] = grada[1] * b + gradb[1] * a;
    prod[2] = grada[2] * b + gradb[2] * a;
}

/* given two values top and bottom and their gradients, compute the 
   gradient of their quotient grad(top / bottom) */
void LocSmoother::gradquotient(REAL top,
                  REAL bot, 
                  REAL gradtop[3],
                  REAL gradbot[3],
                  REAL quot[3])
{
    REAL denom = bot * bot;
    quot[0] = (bot * gradtop[0] - top * gradbot[0]) / denom;
    quot[1] = (bot * gradtop[1] - top * gradbot[1]) / denom;
    quot[2] = (bot * gradtop[2] - top * gradbot[2]) / denom;
}

/* compute Z, a quantity associated with circumradius computation
   TODO this code is lifted from Jonathan's tetcircumcenter computation
   in primitives.c */
REAL LocSmoother::getZ(REAL *tetorg,
          REAL *tetdest,
          REAL *tetfapex,
          REAL *tettapex)
{
    REAL xot, yot, zot, xdt, ydt, zdt, xft, yft, zft;
    REAL otlength, dtlength, ftlength;
    REAL xcrossdf, ycrossdf, zcrossdf;
    REAL xcrossfo, ycrossfo, zcrossfo;
    REAL xcrossod, ycrossod, zcrossod;
    REAL xct, yct, zct;

    /* Use coordinates relative to the apex of the tetrahedron. */
    xot = tetorg[0] - tettapex[0];
    yot = tetorg[1] - tettapex[1];
    zot = tetorg[2] - tettapex[2];
    xdt = tetdest[0] - tettapex[0];
    ydt = tetdest[1] - tettapex[1];
    zdt = tetdest[2] - tettapex[2];
    xft = tetfapex[0] - tettapex[0];
    yft = tetfapex[1] - tettapex[1];
    zft = tetfapex[2] - tettapex[2];
    /* Squares of lengths of the origin, destination, and face apex edges. */
    otlength = xot * xot + yot * yot + zot * zot;
    dtlength = xdt * xdt + ydt * ydt + zdt * zdt;
    ftlength = xft * xft + yft * yft + zft * zft;
    /* Cross products of the origin, destination, and face apex vectors. */
    xcrossdf = ydt * zft - yft * zdt;
    ycrossdf = zdt * xft - zft * xdt;
    zcrossdf = xdt * yft - xft * ydt;
    xcrossfo = yft * zot - yot * zft;
    ycrossfo = zft * xot - zot * xft;
    zcrossfo = xft * yot - xot * yft;
    xcrossod = yot * zdt - ydt * zot;
    ycrossod = zot * xdt - zdt * xot;
    zcrossod = xot * ydt - xdt * yot;
    
    /* Calculate offset (from apex) of circumcenter. */
    xct = (otlength * xcrossdf + dtlength * xcrossfo + ftlength * xcrossod);
    yct = (otlength * ycrossdf + dtlength * ycrossfo + ftlength * ycrossod);
    zct = (otlength * zcrossdf + dtlength * zcrossfo + ftlength * zcrossod);
        
    /* Calculate the length of this vector, which is Z */
    return sqrt(xct * xct + yct * yct + zct * zct);
}

/* get the information about this tet needed for non-smooth
   optimization of the current quality measure */
void LocSmoother::getoptinfo(INTEGER iNod,
                 INTEGER iElem,
                 REAL *qual,
                 REAL qualgrad[][3],
                 REAL *volume,
                 REAL volumegrad[3],
                 int qualmeasure)
{
	REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};       /* the vertices of the tet */
    REAL edgelength[3][4];  /* the lengths of each of the edges of the tet */
    REAL edgegrad[3][4][3]; /* the gradient of each edge length wrt vtx1 */
    REAL facenormal[4][3];  /* the normals of each face of the tet */
    REAL facearea[4];       /* areas of the faces of the tet */
    REAL facegrad[4][3];    /* the gradient of each of the face areas wrt vtx1 */
    int i, j, k, l, m;           /* loop indices */
    int edgecount=0;          /* keep track of current edge */
    REAL ejk[3];            /* vector representing edge from j to k */
    REAL ejl[3];            /* vector representing edge from j to l */
    REAL t[3];
    REAL u[3];
    REAL v[3];
    REAL e1[3] = {0.0, 0.0, 0.0};
    REAL e2[3] = {0.0, 0.0, 0.0};
    
    /* temporary variables */
    REAL diff[3];
    REAL term1[3];
    REAL term2[3];
    REAL factor;
    REAL c;
    REAL top, bot;
    REAL gradtop[3];
    REAL gradbot[3];
    REAL gradquot[3];
    
    /* radius ratio vars */
    REAL zValue;
    REAL twooverZ;
    REAL gradZtfac, gradZufac, gradZvfac;
    REAL gradZt[3];
    REAL gradZu[3];
    REAL gradZv[3];
    REAL gradZ[3];
    REAL faceareasum;
    REAL rootZareasum;
    REAL facegradsum[3];
    REAL vdott;
    REAL tdotu;
    REAL udotv;
    REAL tlen2;
    REAL ulen2;
    REAL vlen2;
    REAL uminusv[3];
    REAL vminust[3];
    REAL tminusu[3];
    REAL umvlen2;
    REAL vmtlen2;
    REAL tmulen2;
    REAL normfac = sqrt(3.0) * 6.0;
    REAL rnrrgradterm1[3], rnrrgradterm2[3];
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    
    /* V / lrms^3 ratio vars */
    REAL edgelengthsum = 0.0;
    REAL lrms;
    REAL gradlrms[3];
    REAL vlrmsterm1[3];
    REAL vlrmsterm2[3];
	INTEGER iElemNd;
	int ia, ib, ic, id;
	REAL *erray[][2] = {{nullptr, nullptr}, {u, v}, {t, u}, {v, t}};

    /* get tet vertices, 注意绕向：将iNod放置到第0个顶点位置 */
    for (m = 0; m <= DIM; m++)
	{
		iElemNd = g_pLocSmoothingElems[iElem].form[m];
		if (iElemNd == iNod)
		{

//			memcpy(point[0], g_pLocSmoothingNodes[iElemNd].pt, sizeof(APOINT));
			point[0] = g_pLocSmoothingNodes[iElemNd].pt;
			DNC(m, ia, ib, ic, id);
//			memcpy(point[1], g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[ib]].pt, sizeof(APOINT));
//			memcpy(point[2], g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[id]].pt, sizeof(APOINT));
//			memcpy(point[3], g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[ic]].pt, sizeof(APOINT));
			point[1] = g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[ib]].pt;
			point[2] = g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[id]].pt;
			point[3] = g_pLocSmoothingNodes[g_pLocSmoothingElems[iElem].form[ic]].pt;
			break;
		}
	}
	assert(m <= DIM);
    /* set some vectors */
    for (m = 0; m < DIM; m++)
	{
		t[m] = point[1][m] - point[0][m];
		v[m] = point[2][m] - point[0][m];
		u[m] = point[3][m] - point[0][m];
	}

    /* calculate the volume*6 of the tetrahedron using orientation */
    *volume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]) / 6.0;
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        /* for each vertex i in the loop, the ith face is the face
           opposite i, so that face's normal is found by taking the
           cross product of two edges of the opposite face */
        /* TODO implement this cross product with Orient2D calls? */
        
        /* one edge on opposite face */
        //vsub(point[k], point[j], ejk);
        
        /* another edge originating from the same vertex */
        //vsub(point[l], point[j], ejl);
        for (m = 0; m < DIM; m++)
		{
			ejk[m] = point[k][m] - point[j][m];
			ejl[m] = point[l][m] - point[j][m];
		}

        /* compute a normal vector to this face */
        //cross(ejk, ejl, facenormal[i]);
        facenormal[i][0] = ejk[1]*ejl[2] - ejk[2]*ejl[1];
		facenormal[i][1] = ejk[2]*ejl[0] - ejk[0]*ejl[2];
		facenormal[i][2] = ejk[0]*ejl[1] - ejk[1]*ejl[0];

        /* if i=0, this is also the gradient of the volume * 6
           with respect to vertex 0 */
        if (i == 0)
        {
            volumegrad[0] = -facenormal[i][0] / 6.0;
            volumegrad[1] = -facenormal[i][1] / 6.0;
            volumegrad[2] = -facenormal[i][2] / 6.0;
        }
        
        /* compute (2 *area)^2 for this face */
        facearea[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        /* now get the real area */
        facearea[i] = sqrt(facearea[i]) / 2.0;
            
        /* compute the gradient of the area for this face */
        if (i == 0)
        {
            /* this face doesn't include vtx1, gradient is zero */
            facegrad[i][0] = 0.0;
            facegrad[i][1] = 0.0;
            facegrad[i][2] = 0.0;
        }
        else
        {
            assert(facearea[i] > 0);
            /* gradient scaled by the face's area */
            factor = 1.0 / (4.0 * facearea[i]);
            
            /* handle each face separately */
			memcpy(e1, erray[i][0], sizeof(REAL)*3);
			memcpy(e2, erray[i][1], sizeof(REAL)*3);
            
            /* find the vector from elk to elj */
            for (m = 0; m < DIM; m++)
				diff[m] = e1[m] - e2[m];
            
            /* compute first term of gradient */
            c = e2[0]*diff[0] + e2[1]*diff[1] + e2[2]*diff[2];
            term1[0] = c * e1[0];
            term1[1] = c * e1[1];
            term1[2] = c * e1[2];
            
            /* compute the second term */
            c = e1[0]*diff[0] + e1[1]*diff[1] + e1[2]*diff[2];
            term2[0] = c * e2[0];
            term2[1] = c * e2[1];
            term2[2] = c * e2[2];
            
            /* now, combine the terms, scaled with the 1/4A */
            facegrad[i][0] = factor * (term1[0] - term2[0]);
            facegrad[i][1] = factor * (term1[1] - term2[1]);
            facegrad[i][2] = factor * (term1[2] - term2[2]);
        }
        
            
        /* compute edge lengths for quality measures that need them */
        if (qualmeasure == QUALMINSINE ||
            qualmeasure == QUALWARPEDMINSINE ||
            qualmeasure == QUALVLRMS3RATIO)
        {
            for (j = i + 1; j < 4; j++) {
            
                /* e1 is edge from point i to point j */
                for (m = 0; m < DIM; m++)
					e1[m] = point[j][m] - point[i][m];
				edgelength[i][j] = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
            
                /* also compute the gradient of the length of this edge */
            
                /* if vtx1 isn't one of the edge's endpoints, the gradent is zero */
                if (i != 0)
                {
                    edgegrad[i][j][0] = 0.0;
                    edgegrad[i][j][1] = 0.0;
                    edgegrad[i][j][2] = 0.0;
                }
                /* otherwise, it's in the negative direction of this edge,
                   and scaled by edge length */
                else
                { 
                    factor = -1.0 / edgelength[i][j];
					for (m = 0; m < DIM; m++)
						edgegrad[i][j][m] = factor * e1[m];
                 }
            }
        }
    }
    
    /* if the quality measure is minimum sine */
    if (qualmeasure == QUALMINSINE || qualmeasure == QUALWARPEDMINSINE)
    {
		edgecount = 0;
        /* for each edge in the tetrahedron */
        for (i = 0; i < 3; i++) {
            for (j = i + 1; j < 4; j++) {
                k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
                l = 6 - i - j - k;
            
                /* compute the sine of this dihedral angle */
                qual[edgecount] = (3 * (*volume) * edgelength[i][j]) / (2 * facearea[k] * facearea[l]);
                
                /* if we are warping the minimum sine */
                if (qualmeasure == QUALWARPEDMINSINE)
                {
                    /* and this is an obtuse angle */
                    if (facenormal[k][0]*facenormal[l][0] + facenormal[k][1]*facenormal[l][1] + facenormal[k][2]*facenormal[l][2] > 0)
                    {                     
                        /* scale the sin down by WARPFACTOR */
                        qual[edgecount] *= WRAPPED_SINE_RATIO;
                    }
                }
            
                /* compute the gradient of the sine
               
                   we need the gradient of this expression:
               
                   3 * V * lij
                   ------------
                   2 * Ak * Al
               
                   so, find the gradient of the top product, the bottom product, then the quotient
                */
                top = (*volume) * edgelength[i][j];
                bot = facearea[k] * facearea[l];
            
                /* find gradient of top */
                gradproduct(*volume, edgelength[i][j], volumegrad, edgegrad[i][j], gradtop);
            
                /* find gradient of bottom */
                gradproduct(facearea[k], facearea[l], facegrad[k], facegrad[l], gradbot);
            
                /* now, find the gradient of the quotient */
                gradquotient(top, bot, gradtop, gradbot, gradquot);
            
                /* now scale with constant factor */
                c = 1.5;//3.0 / 2.0;
                qualgrad[edgecount][0] = c * gradquot[0];
                qualgrad[edgecount][1] = c * gradquot[1];
                qualgrad[edgecount][2] = c * gradquot[2];
            
                edgecount++;
            }
        }
    }
    
    /* compute stuff for radius ratio */
    if (qualmeasure == QUALRADIUSRATIO)
    {
        /* compute intermediate quantity Z */
        zValue = getZ(point[0], point[1], point[2], point[3]);
        
        twooverZ = 2.0 / zValue;
        
        /* some dot products */
        vdott = v[0]*t[0] + v[1]*t[1] + v[2]*t[2];//dot(v, t);
        tdotu = t[0]*u[0] + t[1]*u[1] + t[2]*u[2];//dot(t, u);
        udotv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];//dot(u, v);
        
        /* some vector lengths */
		for (m = 0; m < DIM; m++)
		{
			uminusv[m] = u[m] - v[m];//vsub(u, v, uminusv);
			vminust[m] = v[m] - t[m];//vsub(v, t, vminust);
			tminusu[m] = t[m] - u[m];//vsub(t, u, tminusu);
		}
        tlen2 = (t[0] * t[0]) + (t[1] * t[1]) + (t[2] * t[2]);
        ulen2 = (u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]);
        vlen2 = (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]);
        umvlen2 = (uminusv[0] * uminusv[0]) + (uminusv[1] * uminusv[1]) + (uminusv[2] * uminusv[2]);
        vmtlen2 = (vminust[0] * vminust[0]) + (vminust[1] * vminust[1]) + (vminust[2] * vminust[2]);
        tmulen2 = (tminusu[0] * tminusu[0]) + (tminusu[1] * tminusu[1]) + (tminusu[2] * tminusu[2]);
        
        /* compute Z's gradient */
        gradZtfac = twooverZ * 
                 (
                    (ulen2 * vdott - vlen2 * tdotu) * (ulen2 - vlen2) - 
                    (ulen2 * vlen2 + tlen2 * udotv) * (umvlen2)
                 );
        gradZufac = twooverZ * 
                 (
                    (vlen2 * tdotu - tlen2 * udotv) * (vlen2 - tlen2) - 
                    (vlen2 * tlen2 + ulen2 * vdott) * (vmtlen2)
                 );
        gradZvfac = twooverZ * 
                 (
                    (tlen2 * udotv - ulen2 * vdott) * (tlen2 - ulen2) - 
                    (tlen2 * ulen2 + vlen2 * tdotu) * (tmulen2)
                 );
        
        for (m = 0; m < DIM; m++)
		{
			/* compute t, u, v components of gradient */
			gradZt[m] = gradZtfac*t[m];//vscale(gradZtfac, t, gradZt);
			gradZu[m] = gradZufac*u[m];//vscale(gradZufac, u, gradZu);
			gradZv[m] = gradZvfac*v[m];//vscale(gradZvfac, v, gradZv);
            
			/* add the components together to form grad(Z) */
			gradZ[m] = gradZt[m] + gradZu[m];//vadd(gradZt, gradZu, gradZ);
			gradZ[m] += gradZv[m];			 //vadd(gradZv, gradZ, gradZ);
		}
        
        /* compute sqrt (Z * (sum of face areas)) */
        faceareasum = facearea[0] +
                      facearea[1] +
                      facearea[2] +
                      facearea[3];
        rootZareasum = sqrt(zValue * faceareasum);
        
        /* set the actual root normalized radius ratio */
        *qual = (normfac * (*volume)) / rootZareasum;
        
        assert(*qual >= 0.0);
        
        
		for (m = 0; m < DIM; m++)
		{
			/* compute the first term */
			rnrrgradterm1[m] = 1.0/rootZareasum * volumegrad[m];//vscale((1.0 / rootZareasum), volumegrad, rnrrgradterm1);
			/* sum of face gradients */
			facegradsum[m] = zValue *(facegrad[0][m] + facegrad[1][m] + facegrad[2][m] + facegrad[3][m]);
			gradZ[m] *= faceareasum;
			/* compute the second term */
			rnrrgradterm2[m] = ((*volume) / (2 * (rootZareasum * rootZareasum * rootZareasum))) * (facegradsum[m] + gradZ[m]);
			qualgrad[0][m] = normfac*(rnrrgradterm1[m] - rnrrgradterm2[m]);
		}
    }
    
    /* if the quality measure is volume to edge length ratio */
    if (qualmeasure == QUALVLRMS3RATIO)
    {
        /* for each edge in the tetrahedron */
        for (i = 0; i < 3; i++) 
        {
            for (j = i + 1; j < 4; j++) 
            {
                k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
                l = 6 - i - j - k;
                
                /* accumulate edge length sum */
                edgelengthsum += edgelength[i][j] * edgelength[i][j];
            }
        }
        
        /* compute the root mean square */
        lrms = sqrt((1.0 / 6.0) * edgelengthsum);
        
        normfac = 6.0 * sqrt(2.0);
        
        /* compute the raw ratio */
        *qual = (normfac * (*volume)) / (lrms * lrms * lrms);
               
        /* compute gradient of lrms */
		for (m = 0; m < DIM; m++)
		{
			/* compute the terms of the gradient of the ratio */
			vlrmsterm1[m] = (1.0 / (lrms * lrms * lrms)) * volumegrad[m];
			gradlrms[m] = (-1.0/(6.0*lrms)) * (t[m] + u[m] + v[m]);
			vlrmsterm2[m] = (3.0 * (*volume)) / (lrms * lrms * lrms * lrms) * gradlrms[m];
			/* add terms and normalize */
			qualgrad[0][m] = normfac * (vlrmsterm1[m] - vlrmsterm2[m]);
		}
    }
}

/* compute the (square) of the minimum sine
   of all the dihedral angles in the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL LocSmoother::minsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
	REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
    REAL facenormal[4][3]; /* the normals of each face of the tet */
    REAL dx, dy, dz;       /* intermediate values of edge lengths */
    REAL facearea2[4];     /* areas of the faces of the tet */
    REAL pyrvolume;        /* volume of tetrahedron */
    REAL sine2, minsine2;  /* the sine (squared) of the dihedral angle */
    int i, j, k, l;          /* loop indices */
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const REAL HUGEFLOAT = 1.0e100;

    /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
    /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
				/* Encountered zero-area face */
                sine2 = 0.0;
            }
            
            /* update minimum sine */
            if (sine2 < minsine2)
            {
                minsine2 = sine2;
            }
        }
    }
    
    return sqrt(minsine2) * pyrvolume;
}

REAL LocSmoother::radtodeg(REAL inangle)
{
    return (inangle * 180) / PI;
}

/* compute the minimum or maximum angle of the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL LocSmoother::minmaxangle(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4, bool max)
{
	REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
    REAL facenormal[4][3]; /* the normals of each face of the tet */
    REAL dx, dy, dz;       /* intermediate values of edge lengths */
    REAL pyrvolume;        /* volume of tetrahedron */
    int i, j, k, l;          /* loop indices */
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	const REAL HUGEFLOAT = 1.0e100;
    REAL minangle = HUGEFLOAT;
    REAL maxangle = 0.0;
    REAL angle, tantheta;
    REAL dotproduct;
    
   /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
    /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
	/* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }

    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the tangent of the angle using the tangent formula:

               tan(theta_ij) = - 6 * V * l_ij
                               --------------
                                dot(n_k, n_l)

               because this formula is accurate in the entire range.
            */
            dotproduct = facenormal[k][0]*facenormal[l][0] + facenormal[k][1]*facenormal[l][1] + facenormal[k][2]*facenormal[l][2];
            
            if (dotproduct != 0.0)
            {
                tantheta = (-pyrvolume * sqrt(edgelength[i][j])) / dotproduct;
            
                /* now compute the actual angle */
                angle = atan(tantheta);
            }
            else
            {
                angle = PI / 2.0;
            }
            
            /* adjust angle for sign of dot product */
            if (dotproduct > 0)
            {
                angle += PI;
            }

            /* make negative angles positive */
            if (angle < 0)
            {
                angle += 2.0 * PI;
            }
            
            if (dotproduct == 0.0) angle = PI / 2.0;
            
            if (angle < minangle) minangle = angle;
            if (angle > maxangle) maxangle = angle;
        }
    }
    
    /*
    assert(radtodeg(maxangle) <= 180.0);
    assert(minangle >= 0.0);
    */
    
    if (max) return radtodeg(maxangle);
    return radtodeg(minangle);
}

/* warp the sine of the dihedral angle to penalize obtuse angles more than acute */
REAL LocSmoother::warpsine(REAL sine)
{
    return sine * WRAPPED_SINE_RATIO; 
}

/* compute the (square) of the minimum sine
   of all the dihedral angles in the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL LocSmoother::warpedminsine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
    REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
    REAL facenormal[4][3]; /* the normals of each face of the tet */
    REAL dx, dy, dz;       /* intermediate values of edge lengths */
    REAL facearea2[4];     /* areas of the faces of the tet */
    REAL pyrvolume;        /* volume of tetrahedron */
    REAL sine2, minsine2;  /* the sine (squared) of the dihedral angle */
    int i, j, k, l;          /* loop indices */
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	const REAL HUGEFLOAT = 1.0e100; /* 这个值要尽量大，否则不保险 */

    /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
   /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    minsine2 = HUGEFLOAT;     /* start with absurdly big value for sine */
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
               /* Encountered zero-area face */
                sine2 = 0.0;
            }
            
            /* check whether this angle is obtuse */
            if (facenormal[k][0] * facenormal[l][0] + facenormal[k][1] * facenormal[l][1] + facenormal[k][2] * facenormal[l][2] > 0)
            {
                /* if so, warp it down */
                sine2 = warpsine(sqrt(sine2));
                sine2 *= sine2;
            }
            
            /* update minimum sine */
            if (sine2 < minsine2)
            {
                minsine2 = sine2;
            }
        }
    }
     
    return sqrt(minsine2) * pyrvolume;
}

/* compute the (square) of the minimum sine
   of all the dihedral angles in the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL LocSmoother::minsineandedgeratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
    REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL edgelength[3][4];  /* the lengths of each of the edges of the tet */
    REAL facenormal[4][3];  /* the normals of each face of the tet */
    REAL dx, dy, dz;        /* intermediate values of edge lengths */
    REAL facearea2[4];      /* areas of the faces of the tet */
    REAL pyrvolume;         /* volume of tetrahedron */
    REAL sine2, minsine2;   /* the sine (squared) of the dihedral angle */
    REAL minsine;
	const REAL HUGEFLOAT = 1.0e100;
    REAL shortest = HUGEFLOAT;
    REAL longest = 0.0;     /* shortest and longest edges in teh tet */
    REAL edgeratio;         /* ratio of shortest to longest edge */
    int i, j, k, l;           /* loop indices */
	const REAL SINEEQUILATERAL = 0.94280903946;
    
   /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
   /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);

	/* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
            
            /* keep track of longest and shortest edge */
            if (edgelength[i][j] > longest) longest = edgelength[i][j];
            if (edgelength[i][j] < shortest) shortest = edgelength[i][j];
        }
    }
    
    minsine2 = 10.0e10;     /* start with absurdly big value for sine */
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
                /* Encountered zero-area face */
                sine2 = 0.0;
            }
            /* update minimum sine */
            if (sine2 < minsine2)
            {
                minsine2 = sine2;
            }
        }
    }
    
    /* edge ratio, scaled down for parity with sin of equilateral tet's dihedrals */
    edgeratio = sqrt(shortest / longest) * SINEEQUILATERAL;
    minsine = sqrt(minsine2) * pyrvolume;
      
    if (edgeratio < minsine) 
    {
        return edgeratio;
    }
    return minsine;
}

/* compute the mean of the sines
   of all the dihedral angles in the tet defined
   by the four vertices (vtx1, vtx2, vtx3, vtx4)
*/
REAL LocSmoother::meansine(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
    REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
    REAL facenormal[4][3]; /* the normals of each face of the tet */
    REAL dx, dy, dz;       /* intermediate values of edge lengths */
    REAL facearea2[4];     /* areas of the faces of the tet */
    REAL pyrvolume;        /* volume of tetrahedron */
    REAL sine2;            /* the sine (squared) of the dihedral angle */
    REAL sinesum=0.0;      /* the accumulating sum of the sines */
    int i, j, k, l;          /* loop indices */
    
    /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
   /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
	/* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }

    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        
        /* compute edge lengths (squared) */
        for (j = i + 1; j < 4; j++) {
            dx = point[i][0] - point[j][0];
            dy = point[i][1] - point[j][1];
            dz = point[i][2] - point[j][2];
            edgelength[i][j] = dx * dx + dy * dy + dz * dz;
        }
    }
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
            k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
            l = 6 - i - j - k;
            
            /* compute the expression for minimum sine, squared, over 4 
               The reason it's over 4 is because the area values we have
               are actually twice the area squared */
            /* if either face area is zero, the sine is zero */
            if (facearea2[k] > 0 && facearea2[l] > 0)
            {
                sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
            }
            else
            {
                /* Encountered zero-area face */
                sine2 = 0.0;
            }
            
            /* accumulate sine */
            sinesum += sqrt(sine2);
        }
    }
    
    /* average sine */
    return (sinesum / 6.0) * pyrvolume;
}

/* the inradius to circumradius ratio */
REAL LocSmoother::radiusratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
    REAL *point[4] = {nullptr, nullptr, nullptr, nullptr};      /* tet vertices */
    REAL facenormal[4][3]; /* the normals of each face of the tet */
    REAL facearea2[4];     /* areas of the faces of the tet */
    REAL pyrvolume;        /* volume of tetrahedron */
    REAL zValue;                /* quantity needed for circumradius */
    REAL facesum=0.0;       /* sum of the areas of the faces */
    int i, j, k, l;          /* loop indices */
    REAL sign;
    REAL qual;
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    
    /* get tet vertices */
	point[0] = vtx1;
	point[1] = vtx2;
	point[2] = vtx3;
	point[3] = vtx4;
    
   /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }
    
    /* for each vertex/face of the tetrahedron */
    for (i = 0; i < 4; i++) {
        j = (i + 1) & 3;
        if ((i & 1) == 0) {
            k = (i + 3) & 3;
            l = (i + 2) & 3;
        } else {
            k = (i + 2) & 3;
            l = (i + 3) & 3;
        }
        
        /* compute the normal for each face */
        facenormal[i][0] =
            (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
            (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
        facenormal[i][1] =
            (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
            (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
        facenormal[i][2] =
            (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
            (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);
            
        /* compute (2 *area)^2 for this face */
        facearea2[i] = facenormal[i][0] * facenormal[i][0] +
            facenormal[i][1] * facenormal[i][1] +
            facenormal[i][2] * facenormal[i][2];
        facesum += sqrt(facearea2[i]) * 0.5;
    }
    
    /* compute Z */
    zValue = getZ(point[0], point[1], point[2], point[3]);
    
    /* now we are ready to compute the radius ratio, which is
    
       (108 * V^2) / Z (A1 + A2 + A3 + A4)
    
       (use 3 instead of 108 because pyrvolume = 6V)
    */
    /* use sqrt for now... */
    sign = (pyrvolume < 0.0) ? -1.0 : 1.0;
    
    qual = sign * sqrt((3.0 * pyrvolume * pyrvolume) / (zValue * facesum));
    
    return qual;
}

/* compute the ratio of the tet volume to the cube of
   the rms edge length */
REAL LocSmoother::vlrms3ratio(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4)
{
    REAL *point[4] = { vtx1, vtx2, vtx3, vtx4};      /* tet vertices */
   // REAL edgelength[3][4]; /* the lengths of each of the edges of the tet */
    REAL dx, dy, dz;       /* intermediate values of edge lengths */
    REAL pyrvolume;        /* volume of tetrahedron */
    int i, j, k, l;          /* loop indices */
    REAL edgelengthsum = 0.0;
    REAL lrms;             /* root mean squared of edge length */
    REAL qual;
    REAL E[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    

    
   /* calculate the volume*6 of the tetrahedron */
    pyrvolume = (REAL) GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);
    
    /* if the volume is zero, the quality is zero, no reason to continue */
    if (pyrvolume <= 0.0)
    {
        return 0.0;
    }
    
    
    /* for each edge in the tetrahedron */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 4; j++) {
			dx = point[i][0] - point[j][0];
			dy = point[i][1] - point[j][1];
			dz = point[i][2] - point[j][2];
            
            edgelengthsum += dx * dx + dy * dy + dz * dz;
        }
    }
    
    /* compute the root mean square */
	const double SQ_2 = sqrt(2.0);
	const double LRM =1.0 / (6*sqrt(6.0));
    /* compute the normalized ratio of volume to lrms^3 */
    
    return (SQ_2* pyrvolume) / (LRM*edgelengthsum*sqrt(edgelengthsum));
}

REAL LocSmoother::tetquality(APOINT vtx1, APOINT vtx2, APOINT vtx3, APOINT vtx4, int qualmeasure)
{
    REAL quality = 0.0; /* the quality of this tetrahedron */
    
    switch (qualmeasure)
    {
        case QUALMINSINE:
            quality = minsine(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALMEANSINE:
            quality = meansine(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALMINSINEANDEDGERATIO:
            quality = minsineandedgeratio(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALRADIUSRATIO:
            quality = radiusratio(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALVLRMS3RATIO:
            quality = vlrms3ratio(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALWARPEDMINSINE:
            quality = warpedminsine(vtx1, vtx2, vtx3, vtx4);
            break;
        case QUALMINANGLE:
            quality = minmaxangle(vtx1, vtx2, vtx3, vtx4, false);
            break;
        case QUALMAXANGLE:
            quality = minmaxangle(vtx1, vtx2, vtx3, vtx4, true);
            break;
        default:
            spdlog::info("I don't know what quality measure {} is. Dying...\n", qualmeasure);
            exit(1);
    }
    
    return quality;
}

/* find the best step to take to improve all of the quality
   functions affected by a vertex vtx in the search direction
   d with an expected rate of improvement r */
void LocSmoother::nonsmoothlinesearch(INTEGER iNod,
	                     Sphere sph,
						 int nSph,
                         REAL d[],
                         REAL inworstqual,
						 REAL *ouworstqual,
                         REAL *alpha,
                         REAL r,
						 int qualmeasure)
{
	const REAL MINSTEPSIZE = 1.0e-5;/* minimum step size */
	const int MAXLINEITER = 50; /* maximum iterations in non-smooth line search */
	const REAL HUGEFLOAT = 1.0e100;
    int numiter = 0;      /* number of iterations */
    REAL v[3], *verts[4];            /* the vertex to be modified */
    REAL origvertex[3]; /* to save the original vertex position */
    REAL offset[3];     /* the offset to move the vertex, alpha * d */
    REAL worstqual;     /* the current worst quality */
    REAL origworstqual; /* the original worst quality */
    REAL thisqual;    /* the quality of the current tet */
    REAL oldworstqual;  /* the worst quality at the last step */
    int i, m, k;                /* loop index */
    INTEGER iElemNd;

    /* save the original worst quality */
    origworstqual = oldworstqual = inworstqual;
    *ouworstqual = 0.0;

    /* fetch the original vertex coordinates from the mesh */
	memcpy(v, g_pLocSmoothingNodes[iNod].pt, sizeof(APOINT));
	memcpy(origvertex, g_pLocSmoothingNodes[iNod].pt, sizeof(APOINT));
     
    /* keep trying until alpha gets too small or we exceed maximum
       number of iterations */
    while ((*alpha > MINSTEPSIZE) && (numiter < MAXLINEITER))
	{   
        /* compute the offset from original vertex positon,
           alpha * d &  move the vertex */
        //vscale(*alpha, d, offset);
        //vadd(v, offset, v);
        for (m = 0; m < DIM; m++)
		{
			offset[m] = (*alpha) * d[m];
			v[m] += offset[m];
		}

        /* recompute all of the quality functions affected
           by v's position, taking note of the smallest one */
        worstqual = HUGEFLOAT; 
        for (i = 0; i < nSph; i++)
        {
			/* 注意，调整iNod的坐标至v (注意：4个顶点的顺序很重要）*/
			for (m = 0, k = 0; m <= DIM; m++)
			{
				iElemNd = g_pLocSmoothingElems[sph[i]].form[m];
				verts[k++] = (iElemNd != iNod) ? g_pLocSmoothingNodes[iElemNd].pt : v;
			}
			/* 节点顺序要反一反 */
			thisqual = tetquality(verts[0], verts[1], verts[3], verts[2], qualmeasure);
            
            /* is this the worst quality we've seen? */
            if (thisqual < worstqual) worstqual = thisqual;
        }
          
         /* if this is not the first iteration, and
           we did better on the last iteration, use
           the step size from the previous iteration */
        if ((oldworstqual > origworstqual) && (oldworstqual > worstqual))
        {
            /* last step did better than current alpha, return it */;
            *alpha = (*alpha) * 2;
            assert(*alpha > 0.0);
            *ouworstqual = oldworstqual;

	//		spdlog::info("Half the alpha values.\n");
            /* put vertex back where it started */
			//vcopy(origvertex, v);
			//memcpy(v, origvertex, sizeof(REAL)*3);
			//memcpy(g_pLocSmoothingNodes[iNod].pt, origvertex, sizeof(REAL)*3);
            return;
        }
        
        /* if we have succeeded in gaining 90% of the expected
           improvement, accept this initial step size */
        if ((worstqual - origworstqual) > (0.9 * (*alpha) * r))
        {
            /* put vertex back where it started */
            //vcopy(origvertex, v);
            //memcpy(v, origvertex, sizeof(REAL)*3);
			//memcpy(g_pLocSmoothingNodes[iNod].pt, origvertex, sizeof(REAL)*3);
			*ouworstqual = worstqual;
		//	spdlog::info("Gaining 90%% of the expected improvement.\n");
            return;
        }
        
        /* this alpha isn't working, going to half it...*/
        /* put vertex back where it started */
        //vcopy(origvertex, v);
        memcpy(v, origvertex, sizeof(REAL)*3);
		//memcpy(g_pLocSmoothingNodes[iNod].pt, origvertex, sizeof(REAL)*3);

        /* cut alpha down by half and try again */
        *alpha = (*alpha) / 2.0;
        
        /* save the worst quality from this step */
        oldworstqual = worstqual;

		numiter++;
    }
    
    /* no positive alpha could be found that improved things... give up and return zero */
    *alpha = 0.0;
}

/* 
 * Smart smoothing combining smart Laplacian method & nonsmooth optimisation method 
 * return true if vertex is moved successfully 
 */
bool LocSmoother::combinedSmoothing(INTEGER iNod,
	           Sphere sph,
			   int nSph,
			   bool bSmartLap,
			   REAL inworstqual,
			   REAL qualthreshold,
			   REAL *outworstqual,
			   REAL *outquals,
			   int qualmeasure)
{
    REAL *v = nullptr;             /* the vertex to be altered */
	REAL *verts[4] = {nullptr, nullptr, nullptr, nullptr};
    REAL origpoint[3];   /* to save the original vertex location */
    int numactive;       /* number of quality functions in the active set */
    int i, j, m, k;      /* loop index */
    int numiter = 0;     /* number of optimization iterations */
    REAL worstqual;      /* the numerical value of the worst quality function */
    REAL thisqual;       /* quality of the current tet */
    REAL oldworstqual, ouworstqual;   /* the numerical value of the worst quality function at the last step */
    REAL improvement;    /* how much we've improved this step */
    REAL d[3];           /* direction to move vertex */
    REAL dlength;        /* length of d */
    REAL r;              /* expected rate of improvement */
    REAL newr;           /* temp r var */
    REAL alpha;          /* step size */
    REAL newalpha;       /* candidate new step size */
    REAL rate;           /* the rate a particular function will improve in the direction d */
    REAL change[3];      /* change we'll make to the point */
    struct vertextype *vinfo;
    REAL qe;
    REAL allgrads[3];
    /* the gradients in the active set */
    REAL quals[MAX_SPHERE_SIZE*6], qualgrads[MAX_SPHERE_SIZE*6][3], volumes[MAX_SPHERE_SIZE], volumegrads[MAX_SPHERE_SIZE][3], activegrads[MAX_SPHERE_SIZE][3];
	const REAL HUGEFLOAT = std::numeric_limits<REAL>::max();
	const REAL DEPSILON = 1.0e-5;
	const int MAXSMOOTHITER = 50;				/* maximum iterations of non smooth optimization */
	const REAL MINSMOOTHITERIMPROVE = 1.0e-5;	/* minimum quality improvement in a smoothing step */
	TetraElemQual tetEleQual;
	REAL q;
    int QUAL_VALUES_PER_ELE = 1;
	REAL avg[DIM], movement[DIM];
	int num_incident_vtx = 0;
	INTEGER iElem, iSecond, iElemNd;
#ifdef _USING_STD_LIB
	std::map<INTEGER, INTEGER> mapG2L;
	std::map<INTEGER, INTEGER>::iterator mapIt;
#else
	IntIntMap mapG2L;
	int value;
#endif
	bool bRet = false;

	if (bSmartLap)
	{
		avg[0] = avg[1] = avg[2] = 0.0;
		for (i = 0, num_incident_vtx = 0; i < nSph; i++)
		{
			iElem = sph[i];
			for (m = 0; m <= DIM; m++)
			{

				iSecond = g_pLocSmoothingElems[iElem].form[m];
				if (iSecond != iNod)
				{
#ifdef _USING_STD_LIB
					mapIt = mapG2L.find(iSecond);
					if (mapIt == mapG2L.end())
#else
					if (!mapG2L.get(iSecond, &value))
#endif
					{
						avg[0] += g_pLocSmoothingNodes[iSecond].pt[0];
						avg[1] += g_pLocSmoothingNodes[iSecond].pt[1];
						avg[2] += g_pLocSmoothingNodes[iSecond].pt[2];
#ifdef _USING_STD_LIB
						mapG2L[iSecond] = num_incident_vtx;
#else
						mapG2L.put(iSecond, num_incident_vtx);
#endif
						
						num_incident_vtx++;
					}
				}
			}
		}
		
		for (m = 0; m < DIM; m++)
		{
			avg[m] /= num_incident_vtx;
		    avg[m] = 2.0 * g_pLocSmoothingNodes[iNod].pt[m] - avg[m];
		}
		worstqual = HUGEFLOAT; 
        for (i = 0; i < nSph; i++)
        {
			/* 注意，调整iNod的坐标至v (注意：4个顶点的顺序很重要）*/
			for (m = 0, k = 0; m <= DIM; m++)
			{
				iElemNd = g_pLocSmoothingElems[sph[i]].form[m];
				verts[k++] = (iElemNd != iNod) ? g_pLocSmoothingNodes[iElemNd].pt : avg;
			}
			/* 节点顺序要反一反 */
			thisqual = tetquality(verts[0], verts[1], verts[3], verts[2], qualmeasure);
            
			outquals[i] = thisqual;

			if (thisqual <= 0.0 || thisqual <= inworstqual)
				break;
            /* is this the worst quality we've seen? */
            if (thisqual < worstqual) worstqual = thisqual;
        }
		
		//if (g_nLocSmoothingCentNode == -100)
		//{
		//	spdlog::info("i == {}; nSph = {}.\n", i, nSph);
		//	for (m = 0; m < nSph; m++)
		//	{
		//		spdlog::info("outquals[{}] = {}.\n", m, outquals[m]);
		//	}
		//}

		if (i >= nSph)
		{
			assert(worstqual > inworstqual && worstqual > 0.0);
			memcpy(g_pLocSmoothingNodes[iNod].pt, avg, sizeof(REAL)*3);
			*outworstqual = worstqual;
			/* 不管后面优化是否成功，都需要返回TRUE */
			bRet = true;

			/* 别忘了这行 (节点位置更新了，需要更新单元质量数据 */
			for (j = 0; j < nSph; j++)
			{
				g_pLocSmoothingQuals[sph[j]].iReserved = 0;
				g_pLocSmoothingQuals[sph[j]].sinValue = outquals[j];
				//updateMeshQual(nSph, sph, outquals);
			}
			if (worstqual >= qualthreshold)
				return bRet;
		} 
	}

	if (qualmeasure == QUALMINSINE || qualmeasure == QUALWARPEDMINSINE)
		QUAL_VALUES_PER_ELE = 6;

    /* record the node to be altered */
    v = g_pLocSmoothingNodes[iNod].pt;
    
    /* record the original position of the vertex */
    memcpy(origpoint, v, sizeof(APOINT));
    
    /* find the worst quality of all incident tets */
    worstqual = HUGEFLOAT; 
    for (i = 0; i < nSph; i++)
    {
		tetEleQual = g_pLocSmoothingQuals[sph[i]];
		thisqual = tetEleQual.sinValue; 
        
		if (thisqual < worstqual) worstqual = thisqual;
    }
    if (worstqual <= 0.0)
    {
        printf("Tangled element is found, with quality = %g\n", worstqual);
    }
    assert(worstqual > 0.0);
    
    *outworstqual = worstqual;
    
    allgrads[0] = allgrads[1] = allgrads[2] = 0.0;

    /* Search the active set of quality functions that are
       nearly as bad as f (the quality difference with the worst value ranges within 3% */
    for (i = 0; i < nSph; i++)
    {
        /* compute gradient info for this tet */
        getoptinfo(iNod, sph[i], &quals[QUAL_VALUES_PER_ELE*i], &qualgrads[QUAL_VALUES_PER_ELE*i], &volumes[i], volumegrads[i], qualmeasure);
#ifdef _DEBUG
		tetEleQual = g_pLocSmoothingQuals[sph[i]];
		thisqual = HUGEFLOAT;
		for (j = 0; j < QUAL_VALUES_PER_ELE; j++)
		{
			if (quals[QUAL_VALUES_PER_ELE*i+j] < thisqual)
				thisqual = quals[QUAL_VALUES_PER_ELE*i+j];
		}
		if ((tetEleQual.sinValue == 0.0 && fabs(thisqual) >= 1.0e-3) || 
			(tetEleQual.sinValue != 0.0 && fabs((thisqual - tetEleQual.sinValue)/tetEleQual.sinValue) >= 1.0e-3))
		{
			spdlog::info("In  ::combinedSmoothing(...)\n");
			printf("Some errors may happen when evaluating the element quality: tetEleQual.sinValue = %f; worstqual = %f.\n", 
			tetEleQual.sinValue, thisqual);
			spdlog::info("Loc Elem Idx:{}.\n", i);
			for (m = 0; m <= DIM; m++)
			{
				printf("#%d\t:(%f,%f,%f).\n", g_pLocSmoothingElems[i].form[m],
					g_pLocSmoothingNodes[g_pLocSmoothingElems[i].form[m]].pt[0],
					g_pLocSmoothingNodes[g_pLocSmoothingElems[i].form[m]].pt[1],
					g_pLocSmoothingNodes[g_pLocSmoothingElems[i].form[m]].pt[2]);
			}
//			exit(1);
		}
#endif
    }

    getactiveset(sph, nSph, quals, qualgrads, activegrads, &numactive, worstqual, qualmeasure);
	if (numactive <= 0)
	{
		spdlog::info("Unexpected error: empty active set (1).\n");
		printf("A possible error case happens when the worst quality of an element is very small, e.g., smaller than 1.0e-15,\n");
		spdlog::info("and round-off errors of float points dominate the value of the worst quality.\n");
		goto FINISH;
	}
	assert(numactive > 0);

#if 0
	if (numactive <= 0)
	{
		spdlog::info("Unexpected error: numactive == 0.\n");
		goto FINISH;
	}
#endif
    
    /* d <- point on the convex hull of all gradients nearest origin */
    minconvexhullpoint(activegrads, numactive, d);
    
    /* if d is the origin, we can't improve this vertex with smoothing. */
    dlength = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (dlength < DEPSILON)
	{
	//	spdlog::info("Unexpected error: dlength < %e\n", DEPSILON);
        goto FINISH;
	}
       
    /* otherwise, it's time to do some smoothing! */
    do
    {
        /* 比较STELLA和OPT-MS的代码，我们发现STELLA中搜索方向未做归一化。参照OPT-MS的作法，我们在这对
		 * 搜索方向做归一化
		 * 我倾向认为OPT-MS的做法是对的，具体的几何意义需要再发觉
		 */
		for (m = 0; m < DIM; m++)
			d[m] /= dlength;
		 
		/* find r, the expected rate of improvement. r is computed as the minimum
           dot product between the improvment direction d and all of the active 
           gradients, so it's like "how fast do we move in the direction of the
           gradient least favored by d." */
        
        /* start with an absurdly big r */
        r = HUGEFLOAT;
        /* find the smallest dot product */
        for (i=0; i<numactive; i++)
        {
            newr = d[0]*activegrads[i][0] + d[1]*activegrads[i][1] + d[2]*activegrads[i][2];
            
            if (newr <= 0.0)
            {   
                /* if we have already moved the vertex some, this is still a success */
                if ((origpoint[0] != v[0]) || (origpoint[1] == v[1]) || (origpoint[2] == v[2]))
                {
					/* record this change */   
					*outworstqual = worstqual;    

					//if (g_nLocSmoothingCentNode == -100)
					//{
					//	spdlog::info("newr <= 0.0; bRet = {}\n", bRet);
					//}
					//for (m = 0; m < nSph; m++)
					//{
					//	spdlog::info("quals[{}] = {}.\n", m, quals[m]);
					//}

					/* 别忘了更新网格质量 */
					goto FINISH;
                }
				else
				{
		//			spdlog::info("Unexpected error: newr <= 0.0\n");
					goto FINISH;
				}
            }
            assert(newr > 0.0);
            
            if (newr < r)
            {
                r = newr;
            }
        }
        
        /* save the worst quality from the previous step */
        oldworstqual = worstqual;
        
        /* initialize alpha to the nearest intersection with another
           quality function */
        alpha = getinitialalpha(iNod, sph, nSph, quals, qualgrads, d, r, worstqual, qualmeasure);
		assert(alpha >= 0.0);

        /* if we didn't find a limit for alpha above, at least limit it
           so that we don't invert any elements. */
        if (alpha == HUGEFLOAT)
        {
            /* using volume gradients to pick alpha...*/
            for (i = 0; i < nSph; i++)
            {
                /* if moving in the direction d will decrease Sph
                   this element's volume */
                rate = d[0]*volumegrads[i][0] + d[1]*volumegrads[i][1] + d[2]*volumegrads[i][2];
                if (rate < 0.0) 
                {
                    newalpha = -volumes[i] / (2.0 * rate);
                    
                    /* if this is smaller than the current step size,
                       use it */
                    if (newalpha < alpha)
                    {
                        alpha = newalpha;
                    }
                }
            }
        }
        
        /* do normal line search */
        nonsmoothlinesearch(iNod, sph, nSph, d, worstqual, &ouworstqual, &alpha, r, qualmeasure);
        assert(alpha >= 0.0);
   ///     spdlog::info("iNod={}: alpha={}; worstqual = {}; ouworstqual = {}\n", iNod, alpha, worstqual, ouworstqual);

		if (alpha > 0.0)
		{
			/* move vertex in direction d step size alpha */
			bRet = true; /* 点的位置已改变，需要返回TRUE */

			for (m = 0; m < DIM; m++)
			{
				change[m] = alpha*d[m];
				v[m] += change[m]; //?
			}
         
			/* recompute quality information */
			oldworstqual = worstqual;
			worstqual = ouworstqual;
            /* how much did we improve this step? */
			improvement = worstqual - oldworstqual;
			assert(improvement >= 0.0);
			
			/* recompute the active set */
			for (i = 0; i < nSph; i++)
			{
				/* compute gradient info for this tet */
				getoptinfo(iNod, sph[i], &quals[QUAL_VALUES_PER_ELE*i], &qualgrads[QUAL_VALUES_PER_ELE*i], &volumes[i], volumegrads[i], qualmeasure);
			}

			getactiveset(sph, nSph, quals, qualgrads, activegrads, &numactive, worstqual, qualmeasure);
#if 0
			if (numactive <= 0)
			{
				bool bAcute;
				float qualsTemp[6];
				spdlog::info("Unexpected error: empty active set (2).\n");
				spdlog::info("worstqual = {}.\n", worstqual);
				/* recompute the active set */
				for (i = 0; i < nSph; i++)
				{
					thisqual = HUGEFLOAT;
					for (j = 0; j < 6; j++)
						if (thisqual > quals[6*i+j])
							thisqual = quals[6*i+j];
					for (m = 0; m <= DIM; m++)
						verts[m] = g_pLocSmoothingNodes[g_pLocSmoothingElems[sph[i]].form[m]].pt;
					if (i+1 == 12)
					{
						spdlog::info("0: %24.20g %24.20g %24.20g\n", verts[0][0], verts[0][1], verts[0][2]);
						spdlog::info("1: %24.20g %24.20g %24.20g\n", verts[1][0], verts[1][1], verts[1][2]);
						spdlog::info("2: %24.20g %24.20g %24.20g\n", verts[2][0], verts[2][1], verts[2][2]);
						spdlog::info("3: %24.20g %24.20g %24.20g\n", verts[3][0], verts[3][1], verts[3][2]);
					}
					printf("elem %d: %f %f %f", i+1, thisqual, 
						tetquality(verts[0], verts[1], verts[3], verts[2], QUALWARPEDMINSINE), 
						tetqual(
						g_pLocSmoothingElems[sph[i]].form[0], g_pLocSmoothingElems[sph[i]].form[1],
						g_pLocSmoothingElems[sph[i]].form[2], g_pLocSmoothingElems[sph[i]].form[3],
						qualsTemp, bAcute, QUALTYPE::DIHEDRAL_ANGLE_SINE));

//					for (j = 0; j < 6; j++)
//						spdlog::info("{} ", quals[6*i+j]);
					spdlog::info("\n");
				}	
				exit(1);
			}
#endif
			assert(numactive > 0);
        
			/* d <- point on the convex hull of all gradients nearest origin */
			minconvexhullpoint(activegrads, numactive, d);
        
			numiter++;
        
			dlength = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		}
	} 
	while (alpha > 0.0 && dlength > DEPSILON && 
           numiter < MAXSMOOTHITER && improvement > MINSMOOTHITERIMPROVE);
    
FINISH:
	/* record this change */
//    memcpy(g_pLocSmoothingNodes[iNod].pt, v, sizeof(APOINT));     
	if (bRet)
	{
		*outworstqual = worstqual;  
    
		if (qualmeasure == QUALMINSINE || qualmeasure == QUALWARPEDMINSINE)
		{
			for (i = 0; i < nSph; i++)
			{
				outquals[i] = HUGEFLOAT;
				for (j = 0; j < 6; j++)
				{
					if (quals[6*i+j] < outquals[i])
						outquals[i] = quals[6*i+j];
				}
			}
		}
		else
			memcpy(outquals, quals, sizeof(REAL)*nSph);


		/* 别忘了这行 (节点位置更新了，需要更新单元质量数据 */
		for (j = 0; j < nSph; j++)
		{
			g_pLocSmoothingQuals[sph[j]].iReserved = 0;
			g_pLocSmoothingQuals[sph[j]].sinValue = outquals[j];
			//updateMeshQual(nSph, sph, outquals);
		}
	}

    return bRet;
}