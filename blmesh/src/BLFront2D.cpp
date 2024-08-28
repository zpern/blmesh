#include <spdlog/spdlog.h> 
 #include "BLFront2D.h"
#include "BLFront.h"


BLFront2D::BLFront2D(void)
{
	m_blType = BLMType::blm2d;
	m_nNods = 2;

}


BLFront2D::~BLFront2D(void)
{
}


BLVector BLFront2D::GetNormal()
{
	return m_normal;
}

void BLFront2D::GetNodes(int* nNods, BLNode** blNods)
{
	int i = 0;
	*nNods = 2;

	for (i=0; i<*nNods; i++)
	{
		blNods[i] = m_pBLNods[i];
	}
}

void BLFront2D::ReplaceFrontNods(BLNode* blNod, BLNode* blNodNew)
{
	int i;

	for (i=0; i<m_nNods; i++)
	{
		if(blNod == m_pBLNods[i])
			m_pBLNods[i] = blNodNew;
	}
}

bool BLFront2D::IncludeNode(BLNode* blNod)
{
	int i;

	for (i=0; i<m_nNods; i++)
	{
		if(blNod == m_pBLNods[i])
			return true;
	}

	return false;
}
