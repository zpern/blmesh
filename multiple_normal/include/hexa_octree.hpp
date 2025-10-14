#ifndef SRC_HEXA_COMMON_BINARY_TREE_HPP_
#define SRC_HEXA_COMMON_BINARY_TREE_HPP_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <array>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#ifdef max
#undef max
#undef min
#endif

#define DIM DerivedV::ColsAtCompileTime

#ifdef max
#undef max
#undef min
#endif

namespace HEXA {

    /**
     * @brief Axis-aligned spatial tree for all simplices of dimensionality DerivedV::ColsAtCompileTime.
     *
     * For example:
     * - simplex: triangle, tree: octree
     * - simplex: edge,     tree: quadtree
     *
     * @tparam DerivedV The simplex, described by an m*n matrix, where m is the number of points and
     *                  n = DerivedV::ColsAtCompileTime.
     */
    template <typename DerivedV>
    class BinaryTree {
    public:
        using BoxVector = Eigen::Matrix<double, 1, DerivedV::ColsAtCompileTime>;
        using Vertices = Eigen::Matrix<double, -1, DerivedV::ColsAtCompileTime>;
        using Topos = Eigen::Matrix<int, -1, DerivedV::RowsAtCompileTime>;
        using BoundingBox = std::array<double, 2 * DerivedV::ColsAtCompileTime>;
        using HashArray = std::array<size_t, 1 << DerivedV::ColsAtCompileTime>;

        BinaryTree() {}

        // Overload the assignment operator
        BinaryTree& operator=(const BinaryTree& other) {
            // Guard self-assignment
            if (this == &other) {
                return *this;
            }
            // Copy simple fields
            hash_array = other.hash_array;
            tree_node_array_ = other.tree_node_array_;
            tree_data_array_ = other.tree_data_array_;
            diag_len_of_box_ = other.diag_len_of_box_;

            // Return the existing object so we can chain this operator
            return *this;
        }

        explicit BinaryTree(const std::array<double, DerivedV::ColsAtCompileTime * 2>& initial_box) {
            tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(), initial_box);
            generateChildHash();
            diag_len_of_box_ = 0;
            for (int i = 0; i < DIM; i++) {
                diag_len_of_box_ += (initial_box[i] - initial_box[DIM + i]) *
                    (initial_box[i] - initial_box[DIM + i]);
            }
            diag_len_of_box_ = std::sqrt(diag_len_of_box_);
        }

        explicit BinaryTree(const Vertices& V, const Topos& F) {
            std::array<double, 2 * DerivedV::ColsAtCompileTime> box_array;
            for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                box_array[i] = V(0, i);
                box_array[i + DerivedV::ColsAtCompileTime] = V(0, i);
            }
            for (size_t i = 1; i < V.rows(); i++) {
                for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
                    box_array[k] = std::min(box_array[k], V(i, k));
                    box_array[k + DerivedV::ColsAtCompileTime] =
                        std::max(box_array[k + DerivedV::ColsAtCompileTime], V(i, k));
                }
            }

            init(box_array);
            size_t index = std::numeric_limits<size_t>::max();
            for (int i = 0; i < F.rows(); i++) {
                DerivedV M;
                for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
                    M.row(k) = V.row(F(i, k));
                }
                index = insert(M, i, index);
            }
        }

        explicit BinaryTree(const DerivedV& box) {
            std::array<double, DerivedV::ColsAtCompileTime * 2> initial_box;
            for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                initial_box[i] = box(0, i);
                initial_box[i + DerivedV::ColsAtCompileTime] = box(1, i);
            }
            tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(), initial_box);
            generateChildHash();
            for (int i = 0; i < DIM; i++) {
                diag_len_of_box_ += (box[i] - box[DIM + i]) * (box[i] - box[DIM + i]);
            }
            diag_len_of_box_ = std::sqrt(diag_len_of_box_);
        }

        explicit BinaryTree(BoxVector& v1, BoxVector& v2) {
            std::array<double, DerivedV::ColsAtCompileTime * 2> initial_box;
            for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                initial_box[i] = v1(i);
                initial_box[i + DerivedV::ColsAtCompileTime] = v2(i);
            }
            tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(), initial_box);
            generateChildHash();
            for (int i = 0; i < DIM; i++) {
                diag_len_of_box_ += (initial_box[i] - initial_box[DIM + i]) *
                    (initial_box[i] - initial_box[DIM + i]);
            }
            diag_len_of_box_ = std::sqrt(diag_len_of_box_);
        }

        struct Data {
            Data() : simplex(DerivedV(DerivedV())), local_hint(0) {}

            explicit Data(const DerivedV& s)
                : simplex(s),
                local_hint(0) {  // Set local hint to 0 since the most common ancestor is the root.
                for (int i = 0; i < DIM; i++) {
                    box[i] = s(0, i);
                    box[i + DIM] = s(0, i);
                }
                for (int j = 1; j < s.rows(); j++) {
                    for (int i = 0; i < DIM; i++) {
                        box[i] = std::min(box[i], s(j, i));
                        box[i + DIM] = std::max(box[i + DIM], s(j, i));
                    }
                }
            }
            DerivedV simplex;
            BoundingBox box;
            size_t local_hint;
            size_t id;
        };

        struct TreeNode {
            TreeNode() = default;
            TreeNode(const size_t& current_id, const size_t& father_id, const BoundingBox& bbox)
                : id(current_id), fatherid(father_id), depth(0) {
                setBox(bbox);
            }
            size_t id;                  ///< Index of the current node
            std::size_t fatherid;       ///< Index of the parent node
            std::vector<size_t> children;  ///< Indices of the children nodes in the tree
            std::vector<size_t> datas;     ///< Indices of the data stored in this node
            BoundingBox box;
            BoundingBox coordcentor;
            size_t depth;

            inline void setBox(const BoundingBox& bbox) {
                box = bbox;
                for (int i = 0; i < DIM; i++) {
                    coordcentor[i] = (box[i] + box[DIM + i]) * 0.5;
                }
            }
            inline bool isLeafNode() const { return children.empty(); }
        };

        // Initialize tree using vertices and topologies
        void init(const Vertices& V, const Topos& F) {
            std::array<double, 2 * DerivedV::ColsAtCompileTime> box_array;
            for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                box_array[i] = V(0, i);
                box_array[i + DerivedV::ColsAtCompileTime] = V(0, i);
            }
            for (size_t i = 1; i < V.rows(); i++) {
                for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
                    box_array[k] = std::min(box_array[k], V(i, k));
                    box_array[k + DerivedV::ColsAtCompileTime] =
                        std::max(box_array[k + DerivedV::ColsAtCompileTime], V(i, k));
                }
            }

            init(box_array);
            size_t index = std::numeric_limits<size_t>::max();
            for (int i = 0; i < F.rows(); i++) {
                DerivedV M;
                for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
                    M.row(k) = V.row(F(i, k));
                }
                index = insert(M, i, index);
            }
        }

        // Initialize tree using vertices only (each vertex is inserted as a simplex)
        void init(const Vertices& V) {
            std::array<double, 2 * DerivedV::ColsAtCompileTime> box_array;
            for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                box_array[i] = V(0, i);
                box_array[i + DerivedV::ColsAtCompileTime] = V(0, i);
            }
            for (size_t i = 1; i < V.rows(); i++) {
                for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
                    box_array[k] = std::min(box_array[k], V(i, k));
                    box_array[k + DerivedV::ColsAtCompileTime] =
                        std::max(box_array[k + DerivedV::ColsAtCompileTime], V(i, k));
                }
            }

            init(box_array);
            size_t index = std::numeric_limits<size_t>::max();
            for (int i = 0; i < V.rows(); i++) {
                DerivedV M;
                M.row(0) = V.row(i);
                index = insert(M, i, index);
            }
        }

        // Retrieve hints from all data entries
        void get_hint(std::vector<int>& hints) {
            hints.reserve(tree_data_array_.size());
            for (int i = 0; i < tree_data_array_.size(); ++i) {
                hints.push_back(tree_data_array_[i].local_hint);
            }
        }

        // Initialize tree using a bounding box
        void init(const BoundingBox& box) {
            if (tree_node_array_.size()) {
                throw std::runtime_error("The binary tree has already been initialized.");
                return;
            }
            tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(), box);
            generateChildHash();
            diag_len_of_box_ = 0;
            for (int i = 0; i < DIM; i++) {
                diag_len_of_box_ += (box[i] - box[DIM + i]) * (box[i] - box[DIM + i]);
            }
            diag_len_of_box_ = std::sqrt(diag_len_of_box_);
        }

        /**
         * @brief Insert a simplex into the tree.
         *
         * @param M The simplex (n x dim matrix, where n is the number of points).
         * @param simplex_id The ID for the simplex.
         * @param hint_id A hint index for insertion (default is max size_t).
         * @return The index of the inserted element.
         * @note Assumes that M intersects with the box of node[id].
         */
        size_t insert(const DerivedV& M, const size_t& simplex_id,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) {
            Data data(M);
            data.id = simplex_id;
            size_t index = tree_data_array_.size();
            tree_data_array_.push_back(data);
            if (tree_data_array_.size() <= hint_id) {
                insert(simplex_id, 0);
            }
            else {
                size_t id = tree_data_array_[hint_id].local_hint;
                while (id != 0 && !beFullyIn(tree_node_array_[id].box, data.box)) {
                    id = tree_node_array_[id].fatherid;
                }
                insert(simplex_id, id);
            }
            return index;
        }

        /**
         * @brief Reinsert an updated simplex into the tree.
         *
         * @param M The simplex.
         * @param simplex_id The ID of the simplex.
         * @param hint_id A hint index for reinsertion.
         * @return The simplex_id.
         * @note Assumes that M intersects with the box of node[id].
         */
        size_t reinsert(const DerivedV& M, const size_t& simplex_id,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) {
            Data data(M);
            data.id = simplex_id;
            size_t index = simplex_id;
            tree_data_array_[simplex_id] = data;
            if (tree_data_array_.size() <= hint_id) {
                insert(simplex_id, 0);
            }
            else {
                size_t id = tree_data_array_[hint_id].local_hint;
                while (id != 0 && !beFullyIn(tree_node_array_[id].box, data.box)) {
                    id = tree_node_array_[id].fatherid;
                }
                insert(simplex_id, id);
            }
            return simplex_id;
        }

        // Remove a simplex from the tree using recursion.
        void remove(const size_t& simplex_id,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) {
            if (tree_data_array_.size() <= hint_id) {
                removeRecursion(simplex_id, 0);
            }
            else {
                size_t id = tree_data_array_[hint_id].local_hint;
                while (!beFullyIn(tree_node_array_[id].box, tree_data_array_[simplex_id].box)) {
                    id = tree_node_array_[id].fatherid;
                }
                removeRecursion(simplex_id, id);
            }
            // Optionally mark data as deleted.
            // tree_data_array_.erase(simplex_id);
        }

        // Query the tree for simplices intersecting with the given bounding box.
        void query(const BoundingBox& box, std::vector<std::pair<DerivedV, BoundingBox>>& result,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) const {
            std::vector<size_t> simplex_id_result;
            query(box, simplex_id_result, hint_id);
            if (DerivedV::RowsAtCompileTime > 1) {
                std::sort(simplex_id_result.begin(), simplex_id_result.end());
                auto index_it = std::unique(simplex_id_result.begin(), simplex_id_result.end());
                simplex_id_result.erase(index_it, simplex_id_result.end());
            }
            result.reserve(simplex_id_result.size());
            for (auto& id : simplex_id_result) {
                result.push_back(std::make_pair(tree_data_array_[id].simplex,
                    tree_data_array_[id].box));
            }
        }

        // Query the tree for simplices within a radius of a given point.
        void query(const BoxVector& query_xyz, const double& radius,
            std::vector<size_t>& result,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) const {
            BoundingBox box;
            for (int i = 0; i < DerivedV::ColsAtCompileTime; ++i) {
                box[i] = query_xyz[i] - radius;
                box[i + DerivedV::ColsAtCompileTime] = query_xyz[i] + radius;
            }
            result.clear();
            query(box, result, hint_id);
            if (DerivedV::RowsAtCompileTime > 1) {
                std::sort(result.begin(), result.end());
                auto index_it = std::unique(result.begin(), result.end());
                result.erase(index_it, result.end());
            }
        }

        // Query the tree for simplices intersecting with the given bounding box.
        void query(const BoundingBox& box, std::vector<size_t>& result,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) const {
            if (tree_data_array_.empty()) return;
            if (tree_data_array_.size() <= hint_id) {
                queryRecursion(box, 0, result);
            }
            else {
                size_t id = tree_data_array_[hint_id].local_hint;
                while (id != 0 && !beFullyIn(tree_node_array_[id].box, box)) {
                    id = tree_node_array_[id].fatherid;
                }
                queryRecursion(box, id, result);
            }
            if (DerivedV::RowsAtCompileTime > 1) {
                std::sort(result.begin(), result.end());
                auto index_it = std::unique(result.begin(), result.end());
                result.erase(index_it, result.end());
            }
        }

        // Perform a projection query on the tree.
        size_t projection(const BoxVector& xyz, BoxVector& ans_xyz,
            const double reference_box_len = -1,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) {
            size_t simplex_id;
            double distance;
            Eigen::Vector3d barycentric_coordinates;  // TODO: Adjust for dimensions (e.g., Vector2d)
            return queryNearestTriangle(xyz, simplex_id, distance, ans_xyz,
                barycentric_coordinates, reference_box_len,
                hint_id);
        }

        /**
         * @brief Query the nearest triangle to a given point.
         *
         * @param xyz Query point.
         * @param simplex_id (Output) The id of the nearest simplex.
         * @param distance (Output) The distance to the nearest simplex.
         * @param closet_point (Output) The closest point on the simplex.
         * @param barycentric_coordinates (Output) The barycentric coordinates of the closest point.
         * @param reference_box_len Optional reference box length to start the query.
         * @param hint_id A hint index for the query.
         * @return The id of the nearest simplex.
         */
        size_t queryNearestTriangle(
            const BoxVector& xyz, size_t& simplex_id, double& distance,
            BoxVector& closet_point, Eigen::Vector3d& barycentric_coordinates,
            const double reference_box_len = -1,
            const size_t& hint_id = (std::numeric_limits<size_t>::max)()) {
            if (tree_data_array_.empty() || tree_node_array_.empty()) {
                throw std::runtime_error("Tree is empty!");
                return 0;
            }

            distance = std::numeric_limits<double>::max();
            double query_length = reference_box_len;
            if (query_length <= 0) {
                query_length = diag_len_of_box_ * 1e-4;
            }
            assert(query_length > 0);

            std::vector<size_t> result;
            while (distance > 2 * query_length) {
                query_length *= 2;
                std::array<double, 6> box;
                for (int i = 0; i < 3; i++) {
                    box[i] = xyz(0, i) - query_length;
                    box[3 + i] = xyz(0, i) + query_length;
                }
                result.clear();
                this->query(box, result, hint_id);
                for (auto tsimplex_id : result) {
                    double dis;
                    BoxVector tcloset_point;
                    Eigen::Vector3d tbarycentric_coordinates;
                    Eigen::Matrix<double, 1, 3> P;
                    P(0, 0) = xyz.x();
                    P(0, 1) = xyz.y();
                    P(0, 2) = xyz.z();

                    this->point_simplex_distance(P, tree_data_array_[tsimplex_id].simplex,
                        dis, tcloset_point,
                        tbarycentric_coordinates);
                    if (distance > dis) {
                        distance = dis;
                        barycentric_coordinates = tbarycentric_coordinates;
                        simplex_id = tsimplex_id;
                        closet_point = tcloset_point;
                    }
                    if (distance == 0) break;
                }
                if (result.size() > 0 && distance > 0) {
                    query_length = distance;
                    std::array<double, 6> box;
                    for (int i = 0; i < 3; i++) {
                        box[i] = xyz(0, i) - query_length;
                        box[3 + i] = xyz(0, i) + query_length;
                    }
                    result.clear();
                    this->query(box, result, hint_id);
                    for (auto tsimplex_id : result) {
                        double dis;
                        BoxVector tcloset_point;
                        Eigen::Vector3d tbarycentric_coordinates;
                        Eigen::Matrix<double, 1, 3> P;
                        P(0, 0) = xyz.x();
                        P(0, 1) = xyz.y();
                        P(0, 2) = xyz.z();

                        this->point_simplex_distance(P, tree_data_array_[tsimplex_id].simplex,
                            dis, tcloset_point,
                            tbarycentric_coordinates);
                        if (distance > dis) {
                            distance = dis;
                            barycentric_coordinates = tbarycentric_coordinates;
                            simplex_id = tsimplex_id;
                            closet_point = tcloset_point;
                        }
                    }
                    break;
                }
            }
            return simplex_id;
        }

        // Test function to print leaf node information.
        void test() {
            for (auto i : tree_node_array_) {
                if (i.isLeafNode()) {
                    std::cout << "id=" << i.id << " father id=" << i.fatherid
                        << " size=" << i.datas.size() << std::endl;
                }
            }
            std::cout << "hello world!" << std::endl;
        }

        // Get triangle instances stored in the tree.
        std::vector<std::vector<DerivedV> > getTriInstance() {
            std::vector<std::vector<DerivedV> > ans;
            for (int i = 0; i < tree_node_array_.size(); i++) {
                if (tree_node_array_[i].isLeafNode() && tree_node_array_[i].datas.size()) {
                    ans.push_back(std::vector<DerivedV>());
                    for (int j = 0; j < tree_node_array_[i].datas.size(); j++) {
                        ans.back().push_back(tree_data_array_[tree_node_array_[i].datas[j]].simplex);
                    }
                }
            }
            return ans;
        }

        // Get instances with their IDs and bounding boxes.
        std::vector<std::vector<std::pair<size_t, BoundingBox> > > getInstance() {
            std::vector<std::vector<std::pair<size_t, BoundingBox> > > ans;
            for (int i = 0; i < tree_node_array_.size(); i++) {
                if (tree_node_array_[i].isLeafNode() && tree_node_array_[i].datas.size()) {
                    ans.push_back(std::vector<std::pair<size_t, BoundingBox> >());
                    for (int j = 0; j < tree_node_array_[i].datas.size(); j++) {
                        ans.back().push_back(std::make_pair(
                            tree_node_array_[i].datas[j],
                            tree_data_array_[tree_node_array_[i].datas[j]].box));
                    }
                }
            }
            return ans;
        }

        // Get instances with their simplices and bounding boxes.
        std::vector<std::vector<std::pair<DerivedV, BoundingBox>>> getInstance2() {
            std::vector<std::vector<std::pair<DerivedV, BoundingBox>>> ans;
            for (int i = 0; i < tree_node_array_.size(); i++) {
                if (tree_node_array_[i].isLeafNode() && tree_node_array_[i].datas.size()) {
                    ans.push_back(std::vector<std::pair<DerivedV, BoundingBox>>());
                    for (int j = 0; j < tree_node_array_[i].datas.size(); j++) {
                        ans.back().push_back(std::make_pair(
                            tree_data_array_[tree_node_array_[i].datas[j]].simplex,
                            tree_data_array_[tree_node_array_[i].datas[j]].box));
                    }
                }
            }
            return ans;
        }

        // Compute the distance between a point and a simplex defined by vertices and elements (for triangles).
        template <typename Derivedp, typename DerivedPV, typename DerivedEle,
            typename Derivedsqr_d, typename Derivedc, typename Derivedb>
        void point_simplex_distance(const Eigen::MatrixBase<Derivedp>& p,
            const std::vector<DerivedPV>& V,
            const std::vector<DerivedEle>& Ele,
            const typename DerivedEle::Index primitive,
            Derivedsqr_d& sqr_d,
            Eigen::MatrixBase<Derivedc>& c,
            Eigen::PlainObjectBase<Derivedb>& bary) {
            typedef typename Derivedp::Scalar Scalar;
            typedef Eigen::Matrix<Scalar, 1, Derivedp::ColsAtCompileTime> Vector;
            typedef Vector Point;
            typedef Eigen::Matrix<typename Derivedb::Scalar, 1, 3> BaryPoint;

            const auto& Dot = [](const Point& a, const Point& b) -> Scalar {
                return a.dot(b);
                };
            // Real-time collision detection, Ericson, Chapter 5.
            const auto& ClosestBaryPtPointTriangle =
                [&Dot](Point p, Point a, Point b, Point c,
                    BaryPoint& bary_out) -> Point {
                        // Check if p is in the vertex region outside a.
                        Vector ab = b - a;
                        Vector ac = c - a;
                        Vector ap = p - a;
                        Scalar d1 = Dot(ab, ap);
                        Scalar d2 = Dot(ac, ap);
                        if (d1 <= 0.0 && d2 <= 0.0) {
                            bary_out << 1, 0, 0;
                            return a;
                        }
                        // Check if p is in the vertex region outside b.
                        Vector bp = p - b;
                        Scalar d3 = Dot(ab, bp);
                        Scalar d4 = Dot(ac, bp);
                        if (d3 >= 0.0 && d4 <= d3) {
                            bary_out << 0, 1, 0;
                            return b;
                        }
                        // Check if p is in the edge region of ab; if so, return its projection onto ab.
                        Scalar vc = d1 * d4 - d3 * d2;
                        if (a != b) {
                            if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
                                Scalar v = d1 / (d1 - d3);
                                bary_out << 1 - v, v, 0;
                                return a + v * ab;
                            }
                        }
                        // Check if p is in the vertex region outside c.
                        Vector cp = p - c;
                        Scalar d5 = Dot(ab, cp);
                        Scalar d6 = Dot(ac, cp);
                        if (d6 >= 0.0 && d5 <= d6) {
                            bary_out << 0, 0, 1;
                            return c;
                        }
                        // Check if p is in the edge region of ac.
                        Scalar vb = d5 * d2 - d1 * d6;
                        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
                            Scalar w = d2 / (d2 - d6);
                            bary_out << 1 - w, 0, w;
                            return a + w * ac;
                        }
                        // Check if p is in the edge region of bc.
                        Scalar va = d3 * d6 - d5 * d4;
                        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
                            Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                            bary_out << 0, 1 - w, w;
                            return b + w * (c - b);
                        }
                        // p is inside the face region; compute q through barycentric coordinates (u,v,w).
                        Scalar denom = 1.0 / (va + vb + vc);
                        Scalar v = vb * denom;
                        Scalar w = vc * denom;
                        bary_out << 1.0 - v - w, v, w;
                        return a + ab * v + ac * w;
                };

            assert(Ele[0].cols() <= 3 && "Only simplices up to triangles are considered");
            assert((Derivedb::RowsAtCompileTime == 1 || Derivedb::ColsAtCompileTime == 1) &&
                "bary must be an Eigen Vector or RowVector");
            assert(((Derivedb::RowsAtCompileTime == -1 || Derivedb::ColsAtCompileTime == -1) ||
                (Derivedb::RowsAtCompileTime == Ele[0].rows() || Derivedb::ColsAtCompileTime == Ele[0].rows())) &&
                "bary must be dynamic or match the size of Ele.cols()");

            BaryPoint tmp_bary;
            c = (const Derivedc&)ClosestBaryPtPointTriangle(
                p, V[Ele[primitive](0)],
                V[Ele[primitive](1 % Ele[0].cols())],
                V[Ele[primitive](2 % Ele[0].cols())],
                tmp_bary);
            bary.resize(Derivedb::RowsAtCompileTime == 1 ? 1 : Ele[0].cols(),
                Derivedb::ColsAtCompileTime == 1 ? 1 : Ele[0].cols());
            bary.head(Ele[0].cols()) = tmp_bary.head(Ele[0].cols());
            sqr_d = (p - c).norm();
        }

        // Compute the distance between a point and a simplex defined by a matrix of vertices.
        template <typename Derivedp, typename DerivedPV, typename Derivedsqr_d,
            typename Derivedc, typename Derivedb>
        void point_simplex_distance(const Eigen::MatrixBase<Derivedp>& p,
            const DerivedPV& V, Derivedsqr_d& sqr_d,
            Eigen::MatrixBase<Derivedc>& c,
            Eigen::PlainObjectBase<Derivedb>& bary) {
            typedef typename Derivedp::Scalar Scalar;
            typedef Eigen::Matrix<Scalar, 1, Derivedp::ColsAtCompileTime> Vector;
            typedef Vector Point;
            typedef Eigen::Matrix<typename Derivedb::Scalar, 1, 3> BaryPoint;

            const auto& Dot = [](const Point& a, const Point& b) -> Scalar {
                return a.dot(b);
                };
            // Real-time collision detection, Ericson, Chapter 5.
            const auto& ClosestBaryPtPointTriangle =
                [&Dot](const Point p, const Point a, const Point b, const Point c,
                    BaryPoint& bary_out) -> Point {
                        Vector ab = b - a;
                        Vector ac = c - a;
                        Vector ap = p - a;
                        Scalar d1 = Dot(ab, ap);
                        Scalar d2 = Dot(ac, ap);
                        if (d1 <= 0.0 && d2 <= 0.0) {
                            bary_out << 1, 0, 0;
                            return a;
                        }
                        Vector bp = p - b;
                        Scalar d3 = Dot(ab, bp);
                        Scalar d4 = Dot(ac, bp);
                        if (d3 >= 0.0 && d4 <= d3) {
                            bary_out << 0, 1, 0;
                            return b;
                        }
                        Scalar vc = d1 * d4 - d3 * d2;
                        if (a != b) {
                            if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
                                Scalar v = d1 / (d1 - d3);
                                bary_out << 1 - v, v, 0;
                                return a + v * ab;
                            }
                        }
                        Vector cp = p - c;
                        Scalar d5 = Dot(ab, cp);
                        Scalar d6 = Dot(ac, cp);
                        if (d6 >= 0.0 && d5 <= d6) {
                            bary_out << 0, 0, 1;
                            return c;
                        }
                        Scalar vb = d5 * d2 - d1 * d6;
                        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
                            Scalar w = d2 / (d2 - d6);
                            bary_out << 1 - w, 0, w;
                            return a + w * ac;
                        }
                        Scalar va = d3 * d6 - d5 * d4;
                        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
                            Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                            bary_out << 0, 1 - w, w;
                            return b + w * (c - b);
                        }
                        Scalar denom = 1.0 / (va + vb + vc);
                        Scalar v = vb * denom;
                        Scalar w = vc * denom;
                        bary_out << 1.0 - v - w, v, w;
                        return a + ab * v + ac * w;
                };

            assert(DerivedV::RowsAtCompileTime <= 3 && "Only simplices up to triangles are considered");
            assert((Derivedb::RowsAtCompileTime == 1 || Derivedb::ColsAtCompileTime == 1) &&
                "bary must be an Eigen Vector or RowVector");
            assert(((Derivedb::RowsAtCompileTime == -1 || Derivedb::ColsAtCompileTime == -1) ||
                (Derivedb::RowsAtCompileTime == DerivedV::RowsAtCompileTime ||
                    Derivedb::ColsAtCompileTime == DerivedV::RowsAtCompileTime)) &&
                "bary must be dynamic or match the size of the simplex");

            BaryPoint tmp_bary;
            c = (const Derivedc&)ClosestBaryPtPointTriangle(
                p, V.row(0),
                V.row(1 % DerivedV::RowsAtCompileTime),
                V.row(2 % DerivedV::RowsAtCompileTime),
                tmp_bary);
            bary.resize(Derivedb::RowsAtCompileTime == 1 ? 1 : DerivedV::RowsAtCompileTime,
                Derivedb::ColsAtCompileTime == 1 ? 1 : DerivedV::RowsAtCompileTime);
            bary.head(DerivedV::RowsAtCompileTime) = tmp_bary.head(DerivedV::RowsAtCompileTime);
            sqr_d = (p - c).norm();
        }

    protected:
        inline size_t getHashNum(const Data& data, const BoundingBox& box_centor) const {
            return getHashNum(data.box, box_centor);
        }

        inline size_t getHashNum(const BoundingBox& box, const BoundingBox& box_centor) const {
            size_t hash_num = 0;
            size_t hash_id = 1;
            for (int i = 0; i < DIM; i++) {
                //std::cout << box[i] << "  " << box_centor[i] << std::endl;
                if (box[i] < box_centor[i]) {
                    hash_num |= hash_id;
                }
                hash_id <<= 1;
                if (box[i + DIM] >= box_centor[i]) {
                    hash_num |= hash_id;
                }
                hash_id <<= 1;
            }
            return hash_num;
        }

        // Recursively query the tree for all data indices whose bounding boxes intersect with 'box'
        void queryRecursion(const BoundingBox& box, const size_t& node_id,
            std::vector<size_t>& result) const {
            const auto& coord_center = tree_node_array_[node_id].coordcentor;
            if (tree_node_array_[node_id].isLeafNode()) {
                for (auto i : tree_node_array_[node_id].datas) {
                    const auto& current_box = tree_data_array_[i].box;
                    size_t j = 0;
                    for (; j < DerivedV::ColsAtCompileTime; j++) {
                        if (current_box[j] > box[j + DIM] || current_box[j + DIM] < box[j])
                            break;
                    }
                    if (j == DerivedV::ColsAtCompileTime)
                        result.push_back(i);
                }
            }
            else {
                for (int i = 0; i < hash_array.size(); i++) {
                    auto data_hash = getHashNum(box, coord_center);
                    if (!((~data_hash) & hash_array[i])) {  // if intersected
                        queryRecursion(box, tree_node_array_[node_id].children[i], result);
                    }
                }
            }
        }

        // Recursively remove a simplex from the tree.
        void removeRecursion(const size_t& simplex_id, const size_t& node_id) {
            Data& data = tree_data_array_[simplex_id];
            const auto& coord_center = tree_node_array_[node_id].coordcentor;
            if (tree_node_array_[node_id].isLeafNode()) {
                for (auto it = tree_node_array_[node_id].datas.begin();
                    it != tree_node_array_[node_id].datas.end(); it++) {
                    if (*it == simplex_id) {
                        tree_node_array_[node_id].datas.erase(it);
                        break;
                    }
                }
            }
            else {
                for (int i = 0; i < hash_array.size(); i++) {
                    auto data_hash = getHashNum(data, coord_center);
                    if (!((~data_hash) & hash_array[i])) {  // if intersected
                        removeRecursion(simplex_id, tree_node_array_[node_id].children[i]);
                    }
                }
            }
        }

    protected:
        // Insert a simplex (by id) into the tree node with given node_id.
        void insert(const size_t& simplex_id, const size_t& node_id) {
            if (tree_node_array_[node_id].isLeafNode()) {
                tree_node_array_[node_id].datas.push_back(simplex_id);

                // Decompose the node if too many data entries and the maximum tree height is not reached.
                if (tree_node_array_[node_id].datas.size() > kMaxAllowedNodeDataSize &&
                    tree_node_array_[node_id].depth < kMaxAllowedTreeHeight &&
                    tree_node_array_.size() < tree_data_array_.size() / 2) {
                    size_t child_id = tree_node_array_.size();
                    std::vector<size_t> data_hashs;
                    for (auto j = 0; j < tree_node_array_[node_id].datas.size(); j++) {
                        data_hashs.push_back(
                            getHashNum(tree_data_array_[tree_node_array_[node_id].datas[j]],
                                tree_node_array_[node_id].coordcentor));
                    }
                    for (int i = 0; i < hash_array.size(); i++) {
                        tree_node_array_[node_id].children.push_back(child_id + i);
                        tree_node_array_.push_back(
                            TreeNode(child_id + i, node_id, BoundingBox()));
                        tree_node_array_.back().depth = tree_node_array_[node_id].depth + 1;
                    }
                    for (int i = 0; i < hash_array.size(); i++) {
                        BoundingBox box;
                        size_t hash_num = 1;
                        for (int j = 0; j < DIM; j++) {
                            if (hash_array[i] & hash_num) {  // construct the box for the child
                                box[j] = tree_node_array_[node_id].box[j];
                                box[j + DIM] = tree_node_array_[node_id].coordcentor[j];
                            }
                            else {
                                box[j] = tree_node_array_[node_id].coordcentor[j];
                                box[j + DIM] = tree_node_array_[node_id].box[j + DIM];
                            }
                            hash_num <<= 2;
                        }
                        // Set the bounding box of the child node.
                        tree_node_array_[child_id + i].setBox(box);
                        // Insert data into the child node if their bounding boxes intersect.
                        for (auto j = 0; j < tree_node_array_[node_id].datas.size(); j++) {
                            if (!((~data_hashs[j]) & hash_array[i])) {
                                insert(tree_node_array_[node_id].datas[j], child_id + i);
                            }
                        }
                    }
                    tree_node_array_[node_id].datas.clear();
                }
            }
            else {  // If the node has children
                auto data_hash = getHashNum(tree_data_array_[simplex_id].box,
                    tree_node_array_[node_id].coordcentor);
                for (int i = 0; i < hash_array.size(); i++) {
                    if (!((~data_hash) & hash_array[i])) {  // if intersected
                        if (beFullyIn(tree_node_array_[tree_node_array_[node_id].children[i]].box,
                            tree_data_array_[simplex_id].box)) {
                            tree_data_array_[simplex_id].local_hint =
                                tree_node_array_[node_id].children[i];
                        }
                        insert(simplex_id, tree_node_array_[node_id].children[i]);
                    }
                }
            }
        }

    public:
        /**
         * @brief Generate the hash array for all octants (or quadrants).
         *
         * For example, in a quadtree:
         * - 0101 represents the southwestern quadrant,
         * - 1010 represents the northeastern quadrant.
         *
         * In an octree, similar logic applies.
         * This function is called once during tree construction.
         */
        void generateChildHash() {
            hash_array.fill(0);
            for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                size_t current_size = 1 << i;
                for (int j = 0; j < current_size; j++) {
                    hash_array[j + current_size] = hash_array[j] | (1 << (2 * i + 1));
                    hash_array[j] = hash_array[j] | (1 << (2 * i));
                }
            }
        }

        /**
         * @brief Check whether box1 completely contains box2.
         */
        inline bool beFullyIn(const BoundingBox& box1, const BoundingBox& box2) const {
            for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
                if (box1[i] > box2[i] || box1[i + DIM] < box2[i + DIM])
                    return false;
            }
            return true;
        }

        // Write the tree to a VTK file.
        void WriteTreeToVTK(const std::string& filename) {
            std::cout << "Writing VTK file: " << filename << "\n";
            FILE* vtkFile = fopen(filename.c_str(), "w");
            if (!vtkFile) {
                perror("Unable to open file for writing");
                return;
            }

            // Write VTK header.
            fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
            fprintf(vtkFile, "Octree Simplex Data\n");
            fprintf(vtkFile, "ASCII\n");
            fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");

            // Use a set to avoid duplicate data entries.
            std::set<size_t> uniqueDataIndices;
            for (const auto& node : tree_node_array_) {
                uniqueDataIndices.insert(node.datas.begin(), node.datas.end());
            }

            // Compute total number of points.
            size_t totalPoints = 0;
            for (auto dataIndex : uniqueDataIndices) {
                const auto& data = tree_data_array_[dataIndex];
                totalPoints += data.simplex.rows();
            }

            // Write point data.
            fprintf(vtkFile, "POINTS %zu float\n", totalPoints);
            for (auto dataIndex : uniqueDataIndices) {
                const auto& data = tree_data_array_[dataIndex];
                for (int i = 0; i < data.simplex.rows(); ++i) {
                    fprintf(vtkFile, "%lf %lf %lf\n", data.simplex(i, 0),
                        data.simplex(i, 1),
                        DerivedV::ColsAtCompileTime == 3 ? data.simplex(i, 2) : 0.0);
                }
            }

            // Write cell information.
            fprintf(vtkFile, "CELLS %zu %zu\n", uniqueDataIndices.size(),
                uniqueDataIndices.size() * 4);  // Each cell: number of vertices + vertex indices.
            size_t currentPointIdx = 0;
            for (auto dataIndex : uniqueDataIndices) {
                fprintf(vtkFile, "3 %zu %zu %zu\n", currentPointIdx, currentPointIdx + 1,
                    currentPointIdx + 2);
                currentPointIdx += 3;
            }

            // Write cell types.
            fprintf(vtkFile, "CELL_TYPES %zu\n", uniqueDataIndices.size());
            for (size_t i = 0; i < uniqueDataIndices.size(); ++i) {
                fprintf(vtkFile, "5\n");  // VTK_TRIANGLE
            }

            fclose(vtkFile);
        }

    public:
        size_t kMaxAllowedTreeHeight = 15;
        size_t kMaxAllowedNodeDataSize = 80;

        HashArray hash_array;
        std::vector<TreeNode> tree_node_array_;
        std::vector<Data> tree_data_array_;
        double diag_len_of_box_ = 0.0;
    };

}  // namespace HEXA

#endif  // SRC_TiGER_COMMON_BINARY_TREE_HPP_
