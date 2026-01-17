#ifndef _DISTANCE_CALCULATION_H_ 
#define  _DISTANCE_CALCULATION_H_
#include "../include/BLMesh.h"

class DistanceCalculator {
public:
	DistanceCalculator();
	~DistanceCalculator();
	void ReadInput(INPUTFORMAT input,ConfigArgc cf);
	std::vector<double> getHeightProjection();
protected:
	void PreGenerate();
	void Calculate();
protected:
	static bool isexist;
	BLMesh* blm_;
	std::vector<double> recommand_array_;// output array,  vary from 0-1, just recommand value
};


#endif