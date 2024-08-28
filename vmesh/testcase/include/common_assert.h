#ifndef _COMMON_
#define _COMMON_
#include "Eigen/Dense"
//#include "igl/copyleft/tetgen/tetrahedralize.h"
#define DELETE_ARRAY(X) if(X){delete[] X;X=nullptr;}
static int CheckClosed(const std::vector<std::vector<double> > & V, const std::vector<std::vector<int> > & F) {
	std::vector<std::vector<double > >  TV;
	std::vector<std::vector<int > >  TT;
	std::vector<std::vector<int> >  TF;
	std::string op = "p";
	return 0;
//	return igl::copyleft::tetgen::tetrahedralize(V, F, op, TV, TT, TF);
}
static int CheckClosed(int pnSN0,int pnSEO,int* ppnSFFmO,double* ppdSNC0) {
	std::vector<std::vector<double> >  V(pnSN0, std::vector<double>(3));
	std::vector<std::vector<int> > F(pnSEO, std::vector<int>(3));
	for (int i = 0; i < pnSN0; i++) {
		for (int j = 0; j < 3; j++) {
			V[i][j] = ppdSNC0[3 * i + j];
		}
	}
	for (int i = 0; i < pnSEO; i++) {
		for (int j = 0; j < 3; j++) {
			F[i][j] = ppnSFFmO[3 * i + j];
		}
	}
	return CheckClosed(V, F);
}

#endif