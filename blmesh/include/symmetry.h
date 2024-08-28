#pragma once
#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_
#include <Eigen/Dense>
#include "igl/AABB.h"
#include "binary_tree.hpp"
namespace TiGER {
	class SymmetryPlane {
	public:
		SymmetryPlane() {};
		void init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
			tree.init(V, F);
		}
		Eigen::MatrixXd V_;
		Eigen::MatrixXi F_;
		SymmetryPlane(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
			tree.init(V, F);
		}
		Eigen::RowVector3d project(const Eigen::RowVector3d& vec) {
			Eigen::RowVector3d ans;
			tree.projection(vec, ans);
			if (std::isnan(ans(0)))
				throw std::runtime_error("invalid backgroundmesh");
			return ans;
		}
		/**
		* @brife adjust a normal that started from point
		* try to project the normal
		*/
		void adjustNormal(const Eigen::RowVector3d& point, Eigen::RowVector3d& normal) {
			Eigen::RowVector3d endpoint = point + normal;
			Eigen::RowVector3d nep=project(endpoint);
			normal = nep - point;
		}
		void connect(const SymmetryPlane& sp) {

		}
		TiGER::common::BinaryTree<Eigen::Matrix3d> tree;
	};
}
#endif //!_SYMMETRY_H_