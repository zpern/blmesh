#pragma once
#ifndef SPERICAL_VORONOI_GEOMETRY_H_
#define SPERICAL_VORONOI_GEOMETRY_H_

#include <functional>
#include <memory>

#include "NormalSmoothStrategy.h"
#include "BLVector.h"
#include "BLMesh_define.h"


namespace VORONOISMOOTHING {


	using DofDefinition = std::pair<double, double>;


	/*!
	*
	*
	* \param p1 must be normal
	* \param p2 must be normal
	* \return arc length between p1 and p2
	*/
	inline double PointGeodesicLength(const BLVector& p1, const BLVector& p2) {
		return acos(p1*p2)/PI;
	}



	// a linear or non linear function mappring the any domain to [0,1)
	class WeightFunc {
	public:
		/**
		 * .
		 *
		 * \param val the input parameter, must within DofDefinition
		 * \return return a double value within [0,1)
		 */
		WeightFunc() {};
		virtual double operator()(const double& val) const = 0;
		virtual void setExpoent(const double& val) = 0;
		virtual	double getExpoent() const =0;
		virtual DofDefinition getDomainOfDefinition()=0;
	};

	///// a linear weight function, the upperbound of domian of definition equal to PI/2*upperbound
	//template<int upperbound>
	//class LinearWeightFunc : public WeightFunc {
	//public:
	//	double operator()(const double& val) const {
	//		return val / upperbound * 2 / PI;
	//	}
	//
	//	double getExpoent() const { return 1.0; }

	//};


	/// a non-linear weight function, the upperbound of domian of definition equal to PI/2*upperbound
	template<int upperbound>
	class ExpWightFunc :public WeightFunc {
	public:
		ExpWightFunc() {};
		double operator()(const double& val) const {
			return pow(val / upperbound * 2 / PI, exponent_index_);
		}
		DofDefinition getDomainOfDefinition() {
			DofDefinition def;
			def.first = 0;
			def.second = upperbound * PI / 2;
			return def;
		}

		void setExpoent(const double& val) { exponent_index_ = val; }
		double getExpoent() const { return exponent_index_; }
	public:
		double exponent_index_;
	};
	ExpWightFunc<1> ;
	ExpWightFunc<2>;


	struct VoronoiPoint {
		BLVector pos;
		std::shared_ptr<WeightFunc> func;
	};


	/**
	 ******************************** height ratio to exp ***********************************
	 * most of time, we have a higher tolerance for low height ratio layers,
	 * which means higher input
	 */
	static	double ratio2Exponent(const double& val) {
	//	return 0.1;
		double p = val;
		if (p > 0.3)
			p = 0.3;
		return 0.05* 1.0 / p;
	}


	/**
	 ******************************** geometry function ************************************* 
	 *.
	 */
	
	static BLVector Quaternions(const double& angle, const BLVector& input, const BLVector& axis) {
		BLVector vrot=cos(angle)*input+(1-cos(angle))*(input*axis)*axis+sin(angle)*(input^axis);
		return vrot;
	}


	static BLVector CentorOfThreePoint(VoronoiPoint p1, VoronoiPoint p2, VoronoiPoint p3,double& cost) {
		if (p1.func->getExpoent() == p2.func->getExpoent() &&p2.func->getExpoent() == p3.func->getExpoent()) {
			return solveCenterPointOfCircle(p1.pos, p2.pos, p3.pos);
		}
		else {
			BLVector bp[3];
			bp[0] = p1.pos;
			bp[1] = p2.pos;
			bp[2] = p3.pos;

			short it = 10;
			while (it--) {
				BLVector np=(bp[0]+bp[1]+bp[2]).normalized();
				
				double cost1 = p1.func->operator()(PointGeodesicLength(np,p1.pos));
				double cost2 = p2.func->operator()(PointGeodesicLength(np, p2.pos));
				double cost3 = p3.func->operator()(PointGeodesicLength(np, p3.pos));

				if (cost1 < cost2&&cost1 < cost3) {
					bp[0] = np;
				}
				else if (cost2 < cost3) {
					bp[1] = np;
				}
				else {
					bp[2] = np;
				}
				cost = (cost1 + cost2 + cost3) / 3;

			}
			return (bp[0] + bp[1] + bp[2]).normalized();
		}
		
		return BLVector(1, 0, 0);
	}

	static BLVector CentorOfTwoPoint(VoronoiPoint start, VoronoiPoint end,double& cost) {
		double s = 0;
		double e = 1;



		if (start.func == end.func) {
			return (start.pos + end.pos).normalized();
		}
		//	PointGeodesicLength(start.pos, end.pos);

		double l = PointGeodesicLength(start.pos, end.pos);

		double e1=	start.func->getExpoent();
		double e2 = end.func->getExpoent();

		double u1 = start.func->getDomainOfDefinition().second;
		double u2 = end.func->getDomainOfDefinition().second;
		double t = e1 / e2; //e1*log(x/u1)=e2*log((l-x)/u2)

		/** *************************Newton's Method******************* */

		//double x = l/2;

		vector<double> target_start{1e20*l,l-1e20*l,l/2};
		vector<double> target_start_save{ 1e20*l,l - 1e20*l,l / 2 };
		double x = target_start.back();
		target_start.pop_back();


		int max_iteration = 50;
		while (max_iteration--) {

			double f = pow(x / u1, e1) - pow((l - x) / u2, e2);
		
			double df = e1 * pow(x / u1, e1 - 1) / u1 + e2 * pow((l - x) / u2, e2 - 1) / u2;

			if (abs(f / df) < 1e-20||abs(f)<1e-30)
				break;
		//	cout<< x<<" "<< f << " " <<  df << " " << endl;
			x = x - f / df;
			if (x<0||x>l||std::isnan(x)) {
				if (target_start.empty()) {
					if (e1 < e2)
						x = 0;
					else
						x = l;
					break;
				}
				x = target_start.back();
				target_start.pop_back();
				max_iteration = 50;
				
			}
		}
		if (abs(pow(x / u1, e1) - pow((l - x) / u2, e2)) > 1e-10) {
			return BLVector(10.0, 0, 0);
		}
		//if (!max_iteration) {
		//	cout << "error=" << << endl;
		//}
		//else if(max_iteration!=49)
		//	cout << max_iteration << endl;
		target_start_save.push_back(x);

		double min_cost = 1000;
		for (auto i : target_start_save) {
		
			if ((pow(i / u1, e1) + pow((l - i) / u2, e2)) / 2 < min_cost) {
				min_cost = (pow(i / u1, e1) + pow((l - i) / u2, e2)) / 2;
				x = i;
			}
		}
	//	std::max_element(target_start_save.begin(), target_start_save.end(), []() {})
		cost = min_cost;

	//	cout << x << " " << l<<" "<<u1<<" "<<u2 <<" "<<e1<<" "<<e2<< endl;
	//	cout << Quaternions(x, start.pos, -(start.pos^end.pos).normalized())  << (start.pos + end.pos).normalized() << endl;
		return Quaternions(x,start.pos,-(start.pos^end.pos).normalized());


	}






}
#endif
