#ifndef __BLMesh_BLFrontList_h__
#define	__BLMesh_BLFrontList_h__

#include "BLFront.h"

class BLFrontList
{
public:
	BLFrontList(void);
	~BLFrontList(void);

public:
	BLFront* GetNextFront()
	{
		BLFront* pFront = m_fCurr;
		m_fCurr = m_fCurr->m_pNxtFront;

		return pFront;
	}

	bool HasNextFront()
	{
		//return m_fCurr->m_pNxtFront != NULL;
		return m_fCurr != nullptr;
	}

	void RestoreFront()
	{ 
		m_fCurr = m_fHead; 
	}

	void AddFront(BLFront* blFront);

	void DeleteFront(BLFront* blFront);

private:
	BLFront* m_fHead;
	BLFront* m_fCurr;
};

#endif