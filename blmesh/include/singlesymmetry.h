#ifndef _SYMM_SIMPLE_H
#define _SYMM_SIMPLE_H
#include <vector>
#include <array>
constexpr double eps = 1e-10;
enum SymmType{
X_PLUS=1,
Y_PLUS=2,
Z_PLUS=3,
X_MINUS=4,
Y_MINUS=5,
Z_MINUS=6
};
struct SymmetryPlaneV{
	double value;
	SymmType st;
};
class SymmetryJudger {

public:
	SymmetryJudger() {}
	//true = in, false = out
	bool judge(double* coordinate) {
		bool in = true;
		for (auto p : planes_) {
			double value = p.value;
			int index = ((int)p.st - 1) % 3;
			double my_value = coordinate[index];
			if (p.st <= 3) {
				if (my_value  < value + eps)
					in = false;
			}
			else {
				if (my_value  > value -eps)
					in = false;
			}
			
		}
		return in;
	}
	void addPlane(SymmetryPlaneV sp) {
		planes_.push_back(sp);
	}
private:
	std::vector<SymmetryPlaneV> planes_;
};




#endif