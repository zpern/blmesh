#pragma once
#include "BLFront.h"
class BLFront2D : public BLFront
{
public:
	BLFront2D(void);
	~BLFront2D(void);

public:
	virtual BLVector GetNormal();
	virtual void GetNodes(int* nNods, BLNode** blNods);
	virtual void ReplaceFrontNods(BLNode* blNod, BLNode* blNodNew);
	virtual bool IncludeNode(BLNode* blNod);

	virtual void GetNeigbourFronts(int *nblFront, BLFront** blFronts){}
	virtual void AddNeigbourFronts(BLFront* blFront){}
	virtual void RmvNeigbourFronts(BLFront* blFront){}
	virtual void SetSTriIdx(int nidx, int idx, int itri){}
	virtual int GetSTriIdx(int nidx, int idx){return 0;}
	virtual int GetThirdNodIdx(BLNode* blNod1, BLNode* blNod2) { return 0;}
};

