#pragma once
#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_
#include <Eigen/Dense>
#include "igl/AABB.h"
#include "igl/boundary_loop.h"
#include "igl/project_to_line_segment.h"
#include "igl/avg_edge_length.h"
#include "binary_tree.hpp"
namespace TiGER {
	class SymmetryPlane {
	public:
		SymmetryPlane() {};
        enum SType {
            curved_face = 0,
            x = 1,
            y = 2,
            z = 3
        } stype;
		void init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
			tree.init(V, F);
            V_ = V;
            F_ = F;
            reference_length = 0.02 * igl::avg_edge_length(V_, F_);
			igl::boundary_loop(F,L_);
            judgyType(V_);
		}
		
        Eigen::MatrixXd V_;
		Eigen::MatrixXi F_;
		Eigen::VectorXi L_;
        double reference_length = 1;
		
        SymmetryPlane(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
			tree.init(V, F);
            V_ = V;
            F_ = F;
            reference_length=0.02*igl::avg_edge_length(V_,F_);
            igl::boundary_loop(F, L_);
            judgyType(V_);
		}
		
        Eigen::RowVector3d project(const Eigen::RowVector3d& vec) {
			Eigen::RowVector3d ans;
			tree.projection(vec, ans);
			if (std::isnan(ans(0)))
				throw std::runtime_error("invalid backgroundmesh");
			return ans;
		}
		/**
		* @brife a normal that started from point
		* try to project the normal
		*/

        void adjustNormal(const Eigen::RowVector3d& point, Eigen::RowVector3d& normal, bool boundary = false) {
            Eigen::RowVector3d endpoint = point + normal;
            Eigen::RowVector3d nep;
            if (boundary) {
                nep= projectToLoop(V_, L_, endpoint);
                normal = nep - point;
            }
            else {
                nep = project(endpoint);
                normal = nep - point;
                if (stype == SType::x) {
                    normal(0) = 0; 
                }
                if (stype == SType::y) {
                    normal(1) = 0;
                }
                if (stype == SType::z) {
                    normal(2) = 0;
                }
            }
            
            normal.normalize();
		}
        void adjustNormalSecond(const Eigen::RowVector3d& point, Eigen::RowVector3d& normal, bool boundary = false)
        {
            Eigen::RowVector3d endpoint = point + normal;
            Eigen::RowVector3d nep;
            if( boundary )
            {
                nep = projectToLoop(V_, L_, endpoint);
                normal = nep - point;
            }
            else
            {
                nep = project(endpoint);
                normal = nep - point;
                if( stype == SType::x )
                {
                    normal(0) = 0;
                }
                if( stype == SType::y )
                {
                    normal(1) = 0;
                }
                if( stype == SType::z )
                {
                    normal(2) = 0;
                }
            }

            normal.normalize();
        }

		void connect(const SymmetryPlane& sp) {

		}
		TiGER::common::BinaryTree<Eigen::Matrix3d> tree;
		



        Eigen::RowVector3d closestPointOnSegment(const Eigen::RowVector3d& A,
            const Eigen::RowVector3d& B,
            const Eigen::RowVector3d& P) {
            Eigen::RowVector3d AP = P - A;  // 向量AP
            Eigen::RowVector3d AB = B - A;  // 向量AB
            double AB_squared = AB.squaredNorm();  // 线段AB的长度的平方

            if (AB_squared == 0.0) return A;  // A和B是同一个点

            // 投影点的参数化表示中的系数
            double t = AP.dot(AB) / AB_squared;
            t = std::max(0.0, std::min(1.0, t));  // 将t限制在[0, 1]范围内

            return A + t * AB;  // 计算投影点
        }

        // 将点P投影到由点集V_和索引集L_定义的线段上
        Eigen::RowVector3d projectToLoop(const Eigen::MatrixXd& V_,
            const Eigen::VectorXi& L_,
            const Eigen::RowVector3d& P) {
            double min_distance = std::numeric_limits<double>::max();
            Eigen::RowVector3d closest_point;

            // 遍历每条线段
            for (int i = 0; i < L_.size() - 1; ++i) {
                int index_A = L_(i);
                int index_B = L_(i + 1);

                Eigen::RowVector3d A = V_.row(index_A);
                Eigen::RowVector3d B = V_.row(index_B);

                Eigen::RowVector3d current_closest = closestPointOnSegment(A, B, P);
                double current_distance = (P - current_closest).norm();

                // 更新最近点
                if (current_distance < min_distance) {
                    min_distance = current_distance;
                    closest_point = current_closest;
                }
            }

            // 检查loop是否闭合，如果闭合，还需计算最后一个点和第一个点组成的线段
            if (L_(0) != L_(L_.size() - 1)) {
                Eigen::RowVector3d A = V_.row(L_(L_.size() - 1));
                Eigen::RowVector3d B = V_.row(L_(0));
                Eigen::RowVector3d current_closest = closestPointOnSegment(A, B, P);
                double current_distance = (P - current_closest).norm();

                if (current_distance < min_distance) {
                    min_distance = current_distance;
                    closest_point = current_closest;
                }
            }

            return closest_point;
        }
        
        void judgyType(const Eigen::MatrixXd& V_) {
            if (V_.rows() == 0)
                return;
            double eps = reference_length * 1e-1;
            Eigen::RowVector3d maxV(V_.row(0));
            Eigen::RowVector3d minV(V_.row(0));
            for (int i = 0; i < V_.rows(); i++) {
                for (int j = 0; j < 3; j++) {
                    maxV(j) = std::max(maxV(j), V_(i,j));
                    minV(j) = std::min(minV(j), V_(i, j));
                }
            }

            bool found_SType=false;
            for (int j = 0; j < 3; j++) {
                if (maxV(j) - minV(j) < eps) {
                    stype = SType(j + 1);
                    found_SType = true;
                    break;
                }
            }
            if( !found_SType )
            {
                stype = SType(0);
            }

        }
	};
}
#endif //!_SYMMETRY_H_