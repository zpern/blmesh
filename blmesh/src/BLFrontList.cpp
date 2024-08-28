#include <spdlog/spdlog.h> 
 #include "BLFrontList.h"


BLFrontList::BLFrontList(void)
{
	m_fHead = NULL;
	m_fCurr = NULL;
}


BLFrontList::~BLFrontList(void)
{
}


void BLFrontList::AddFront(BLFront* blFront)
{
	if (blFront == NULL)
	{
		throw(string("Error: null front!"));
	}
	if (m_fHead == NULL)
	{
		m_fHead = blFront;
		m_fCurr = blFront;
	}
	else
	{
		m_fCurr->m_pNxtFront = blFront;
		blFront->m_pPrevFront = m_fCurr;

		m_fCurr = blFront;
		if(m_fCurr == NULL)
		{
			throw(string("Error: null front!"));
		}
	}
}

/*
*@author yhf
*/
void BLFrontList::DeleteFront(BLFront* blFront)
{
	if (blFront == m_fHead) {
		m_fHead = m_fHead->m_pNxtFront;
		if (m_fHead != nullptr)
		{
			m_fHead->m_pPrevFront = nullptr;
			//m_fCurr = m_fHead;
		}
	}
	if (m_fCurr == blFront)
	{
		if (blFront->m_pNxtFront)
			m_fCurr = blFront->m_pNxtFront;
		else
			m_fCurr = blFront->m_pPrevFront;
	}
	if(blFront->m_pNxtFront)
		blFront->m_pNxtFront->m_pPrevFront = blFront->m_pPrevFront;
	if (blFront->m_pPrevFront)
		blFront->m_pPrevFront->m_pNxtFront= blFront->m_pNxtFront;
	blFront->m_pNxtFront = NULL;
	blFront->m_pPrevFront = NULL;

}