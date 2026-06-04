#include "vmesh.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

	// 本地 PGrid_ALM0 调试程序使用的输入数据。
	// 这里不依赖 TiGER 大程序中的 SurfaceMesh / HybridParameters，
	// 只保留 API_Gen_Boundary_ALM_Mesh 真正需要的输入参数。
	struct BoundaryLayerInput
	{
		// -------- 边界层控制参数 --------
		int nLN = 5;									// 边界层层数，对应 BoundaryMeshing.txt 中的 nLN:
		double dLen = 0.01;								// 第一层高度，对应 dLen:
		double dRto = 1.2;								// 层高增长率，对应 dRto:
		std::vector<int> layer_vec;						// 分面层数覆盖参数；PGrid_ALM0 调试文件通常为空。
		std::vector<double> length_vec;					// 分面第一层高度覆盖参数；PGrid_ALM0 调试文件通常为空。
		std::vector<double> ratio_vec;					// 分面增长率覆盖参数；PGrid_ALM0 调试文件通常为空。
		double bisostop = 1.0;							// 各向同性停止高度，对应 bisostop:
		int max_layer_diff = 6;							// 相邻面允许的最大层数差。
		double max_ratio_diff = 0.25;					// 相邻层高比例限制，对应 max_ratio_diff:
		std::vector<double> max_skewness = {0.95, 1.0}; // 最大偏斜度，[棱柱, 金字塔]。
		std::vector<double> max_orth = {0.95, 1.0};		// 最大非正交限制，[棱柱, 金字塔]。
		double clearance = -1;							// 边界层留空距离。
		bool use_multiple_normals = false;				// 是否启用多法向处理，对应 b_use_multiple_normals:

		// -------- 表面网格输入参数 --------
		std::map<int, int> face_types; // 几何面 id -> 边界类型。0=远场，1=物面，2=对称面，3=不推层面，4=周期面，5=相邻网格面。
		std::vector<double> coords;	   // 表面节点坐标，按 x,y,z 连续展开。
		std::vector<int> tris;		   // 表面三角形连接关系，按 v0,v1,v2 连续展开。
		std::vector<int> face_ids;	   // 每个三角形所在的几何面 id。
	};

	// vmesh API 返回的数组由 API 内部使用 new[] 分配。
	// 本地调试程序结束前统一用 delete[] 释放，并把指针置空。
	template <typename T>
	void releaseArray(T *&ptr)
	{
		delete[] ptr;
		ptr = nullptr;
	}

	// 严格检查 BoundaryMeshing.txt 的关键字顺序。
	// 本地测试入口必须和 TiGER/test/test_api/test_ALM.hpp::readBoundary 的 PGrid_ALM0 格式保持一致。
	void requireKey(const std::string &actual, const char *expected)
	{
		if (actual != expected)
		{
			throw std::runtime_error("Expected " + std::string(expected) + ", got " + actual);
		}
	}

	// 读取严格 PGrid_ALM0 格式的 BoundaryMeshing.txt。
	//
	// 文件头部必须按以下顺序出现：
	//   nLN:
	//   dLen:
	//   dRto:
	//   max_prism_skewness:
	//   max_pyramid_skewness:
	//   max_ratio_diff:
	//   bisostop:
	//   b_use_multiple_normals:
	//   boundary_info:
	//   PointsAndCells:
	//
	// 文件中的点 id 和三角形 id 都从 1 开始；面内点id从0开始，这里会按 id 放入连续数组。
	BoundaryLayerInput readBoundaryMeshing(const std::string &path)
	{
		std::ifstream file(path);
		if (!file)
		{
			throw std::runtime_error("Cannot open input file: " + path);
		}

		BoundaryLayerInput input;
		std::string key;

		file >> key >> input.nLN;
		requireKey(key, "nLN:");

		file >> key >> input.dLen;
		requireKey(key, "dLen:");

		file >> key >> input.dRto;
		requireKey(key, "dRto:");

		file >> key >> input.max_skewness[0];
		requireKey(key, "max_prism_skewness:");

		file >> key >> input.max_skewness[1];
		requireKey(key, "max_pyramid_skewness:");

		file >> key >> input.max_ratio_diff;
		requireKey(key, "max_ratio_diff:");

		file >> key >> input.bisostop;
		requireKey(key, "bisostop:");

		int multiple_normal = 0;
		file >> key >> multiple_normal;
		requireKey(key, "b_use_multiple_normals:");
		input.use_multiple_normals = multiple_normal != 0;

		file >> key;
		requireKey(key, "boundary_info:");

		// boundary_info 段格式：face_id face_type，直到 PointsAndCells。
		while (file >> key)
		{
			if (key == "PointsAndCells:")
			{
				break;
			}

			int face_id = 0;
			try
			{
				face_id = std::stoi(key);
			}
			catch (...)
			{
				throw std::runtime_error("Bad face id in boundary_info: " + key);
			}

			int face_type = 0;
			file >> face_type;
			if (!file)
			{
				throw std::runtime_error("Bad boundary_info entry for face " + std::to_string(face_id));
			}
			input.face_types[face_id] = face_type;
		}

		if (key != "PointsAndCells:")
		{
			throw std::runtime_error("Expected PointsAndCells:");
		}

		int tri_count = 0;
		int point_count = 0;
		file >> tri_count >> point_count;
		if (tri_count <= 0 || point_count <= 0)
		{
			throw std::runtime_error("Invalid PointsAndCells size in: " + path);
		}

		input.coords.resize(static_cast<std::size_t>(point_count) * 3);
		for (int i = 0; i < point_count; ++i)
		{
			int point_id = 0;
			double x = 0.0;
			double y = 0.0;
			double z = 0.0;
			file >> point_id >> x >> y >> z;
			if (point_id < 1 || point_id > point_count)
			{
				throw std::runtime_error("Point id out of range: " + std::to_string(point_id));
			}

			const int offset = 3 * (point_id - 1);
			input.coords[offset + 0] = x;
			input.coords[offset + 1] = y;
			input.coords[offset + 2] = z;
		}

		input.tris.resize(static_cast<std::size_t>(tri_count) * 3);
		input.face_ids.resize(tri_count);
		for (int i = 0; i < tri_count; ++i)
		{
			int tri_id = 0;
			int v0 = 0;
			int v1 = 0;
			int v2 = 0;
			int face_id = 0;
			file >> tri_id >> v0 >> v1 >> v2 >> face_id;
			if (tri_id < 1 || tri_id > tri_count)
			{
				throw std::runtime_error("Triangle id out of range: " + std::to_string(tri_id));
			}

			const int offset = 3 * (tri_id - 1);
			input.tris[offset + 0] = v0;
			input.tris[offset + 1] = v1;
			input.tris[offset + 2] = v2;
			input.face_ids[tri_id - 1] = face_id;
		}

		return input;
	}

	// 根据 vmesh 单元类型返回该单元的节点数。
	// 体网格连接关系是混合单元的扁平数组，写 VTK 时需要靠这个函数推进 offset。
	int nodesPerCell(int type)
	{
		switch (type)
		{
		case TiGER::TETRAHEDRON:
			return 4;
		case TiGER::PRISM:
			return 6;
		case TiGER::PYRAMID:
			return 5;
		default:
			throw std::runtime_error("Unsupported element type: " + std::to_string(type));
		}
	}

	// 将 vmesh 的单元类型转换成 VTK legacy unstructured grid 的单元类型编号。
	int vtkCellType(int type)
	{
		switch (type)
		{
		case TiGER::TETRAHEDRON:
			return 10;
		case TiGER::PRISM:
			return 13;
		case TiGER::PYRAMID:
			return 14;
		default:
			throw std::runtime_error("Unsupported element type: " + std::to_string(type));
		}
	}

	// 将 API 返回的边界层体网格写成 VTK 文件，便于用 ParaView 检查棱柱/金字塔/四面体结果。
	void writeVolumeVtk(const std::string &path, const double *coords, int point_count, const int *cells, const int *cell_types, int cell_count)
	{
		std::ofstream out(path);
		if (!out)
		{
			throw std::runtime_error("Cannot write output file: " + path);
		}

		int connectivity_size = 0;
		for (int i = 0; i < cell_count; ++i)
		{
			connectivity_size += nodesPerCell(cell_types[i]) + 1;
		}

		out << "# vtk DataFile Version 3.0\n";
		out << "vmesh PGrid_ALM0 volume mesh\n";
		out << "ASCII\n";
		out << "DATASET UNSTRUCTURED_GRID\n";
		out << "POINTS " << point_count << " double\n";
		out << std::setprecision(17);
		for (int i = 0; i < point_count; ++i)
		{
			out << coords[3 * i + 0] << " " << coords[3 * i + 1] << " " << coords[3 * i + 2] << "\n";
		}

		out << "CELLS " << cell_count << " " << connectivity_size << "\n";
		int offset = 0;
		for (int i = 0; i < cell_count; ++i)
		{
			const int n = nodesPerCell(cell_types[i]);
			out << n;
			for (int j = 0; j < n; ++j)
			{
				out << " " << cells[offset + j];
			}
			out << "\n";
			offset += n;
		}

		out << "CELL_TYPES " << cell_count << "\n";
		for (int i = 0; i < cell_count; ++i)
		{
			out << vtkCellType(cell_types[i]) << "\n";
		}
	}

	// 将 API 返回的边界层顶面三角网格写成 VTK 文件。
	// 在 TiGER 大程序流程里，这个顶面网格会继续传给四面体网格生成器。
	void writeTopSurfaceVtk(const std::string &path, const double *coords, int point_count, const int *tris, int tri_count, const double *point_sizing)
	{
		std::ofstream out(path);
		if (!out)
		{
			throw std::runtime_error("Cannot write output file: " + path);
		}

		out << "# vtk DataFile Version 3.0\n";
		out << "vmesh PGrid_ALM0 top surface for tetra meshing\n";
		out << "ASCII\n";
		out << "DATASET UNSTRUCTURED_GRID\n";
		out << "POINTS " << point_count << " double\n";
		out << std::setprecision(17);
		for (int i = 0; i < point_count; ++i)
		{
			out << coords[3 * i + 0] << " " << coords[3 * i + 1] << " " << coords[3 * i + 2] << "\n";
		}

		out << "CELLS " << tri_count << " " << tri_count * 4 << "\n";
		for (int i = 0; i < tri_count; ++i)
		{
			out << "3 " << tris[3 * i + 0] << " " << tris[3 * i + 1] << " " << tris[3 * i + 2] << "\n";
		}

		out << "CELL_TYPES " << tri_count << "\n";
		for (int i = 0; i < tri_count; ++i)
		{
			out << "5\n";
		}

		if (point_sizing != nullptr)
		{
			out << "POINT_DATA " << point_count << "\n";
			out << "SCALARS sizing double 1\n";
			out << "LOOKUP_TABLE default\n";
			for (int i = 0; i < point_count; ++i)
			{
				out << point_sizing[i] << "\n";
			}
		}
	}

} // namespace

int main(int argc, char *argv[])
{
	// 命令行参数：
	//   argv[1]：严格 PGrid_ALM0 格式的 BoundaryMeshing.txt。
	//   argv[2]：输出体网格 VTK 文件路径。
	//   argv[3]：输出边界层顶面 VTK 文件路径。
	const std::string input_path = argc > 1 ? argv[1] : "BoundaryMeshing.txt";
	const std::string volume_output = argc > 2 ? argv[2] : "boundary_layer.vtk";
	const std::string top_output = argc > 3 ? argv[3] : "top_layer.vtk";

	// -------- API 输出：体网格 --------
	double *volume_coords = nullptr;  // 体网格节点坐标，长度 = volume_point_count * 3。
	int volume_point_count = 0;		  // 体网格节点数量。
	int *volume_cells = nullptr;	  // 体网格单元连接关系，按单元类型顺序连续存储节点 id。
	int *volume_cell_types = nullptr; // 体网格单元类型：11=四面体，13=棱柱，14=金字塔。
	int volume_cell_count = 0;		  // 体网格单元数量。

	// -------- API 输出：边界层顶面网格，后续四面体剖分使用 --------
	double *top_coords = nullptr;	// 顶面网格节点坐标，长度 = top_point_count * 3。
	int top_point_count = 0;		// 顶面网格节点数量。
	int top_tri_count = 0;			// 顶面三角形数量。
	int *top_cell_types = nullptr;	// 顶面单元类型数组。
	int *top_tris = nullptr;		// 顶面三角形连接关系，长度 = top_tri_count * 3。
	int *local_to_global = nullptr; // 顶面局部点 id 到体网格全局点 id 的映射。
	double *point_sizing = nullptr; // 顶面节点尺寸场，后续四面体网格可用它控制尺寸。

	// -------- API 输出：生成体网格的边界面信息 --------
	int boundary_face_count = 0;  // 生成体网格的边界面数量。
	int *boundary_mesh = nullptr; // 边界面连接关系，每 4 个 id 为一组；三角形第 4 位为 -1。
	int *boundary_face = nullptr; // 每个边界面对应的原始几何面 id。

	try
	{
		BoundaryLayerInput input = readBoundaryMeshing(input_path);

		const int status = TiGER::API_Gen_Boundary_ALM_Mesh(
			// -------- 输入：表面网格 --------
			input.coords.data(),					   // pdSNC：表面节点坐标数组。
			static_cast<int>(input.coords.size() / 3), // nSN：表面节点数量。
			input.tris.data(),						   // pnSFFm：表面三角形连接关系。
			nullptr,								   // pnSFTp：表面单元类型；当前只处理三角形，可传空。
			input.face_ids.data(),					   // pnSFPt：每个三角形所在的几何面 id。
			static_cast<int>(input.face_ids.size()),   // nSF：表面三角形数量。
			input.face_types,						   // pnFT：几何面 id 到边界类型的映射。

			// -------- 输入：边界层参数 --------
			input.nLN,			  // 边界层层数。
			input.layer_vec,	  // 分面层数覆盖。
			input.max_layer_diff, // 相邻面最大层数差。
			input.max_ratio_diff, // 相邻层高比例限制。
			input.dLen,			  // 第一层高度。
			input.length_vec,	  // 分面第一层高度覆盖。
			input.dRto,			  // 层高增长率。
			input.ratio_vec,	  // 分面增长率覆盖。
			input.max_skewness,	  // 棱柱/金字塔最大偏斜度。
			input.max_orth,		  // 最大非正交限制。
			input.bisostop,		  // 各向同性停止高度。
			input.clearance,	  // 边界层留空距离。

			// -------- 输出：体网格 --------
			&volume_coords,
			&volume_point_count,
			&volume_cells,
			&volume_cell_types,
			&volume_cell_count,

			// -------- 输出：边界层顶面网格 --------
			&top_coords,
			&top_point_count,
			&top_tri_count,
			&top_cell_types,
			&top_tris,
			&local_to_global,
			&point_sizing,

			// -------- 输出：边界面信息 --------
			&boundary_face_count,
			&boundary_mesh,
			&boundary_face,

			// -------- 控制参数 --------
			true,						// b_have_pyramid：是否允许生成金字塔过渡单元。
			input.use_multiple_normals, // b_use_multiple_normals：是否启用多法向。
			false,						// b_output_io_file：是否让 API 额外写调试文件。
			"virtualmesh",				// filename：API 内部调试输出名/目录。
			std::array<double, 12>());	// per_matrix：周期面控制矩阵，本地测试使用默认空矩阵。

		if (status != 0)
		{
			std::cerr << "API_Gen_Boundary_ALM_Mesh returned " << status << "\n";
			return status;
		}

		writeVolumeVtk(volume_output, volume_coords, volume_point_count, volume_cells, volume_cell_types, volume_cell_count);
		writeTopSurfaceVtk(top_output, top_coords, top_point_count, top_tris, top_tri_count, point_sizing);

		std::cout << "Input: " << input_path << "\n";
		std::cout << "Wrote " << volume_output << "\n";
		std::cout << "Wrote " << top_output << "\n";
		std::cout << "Volume points=" << volume_point_count << ", cells=" << volume_cell_count << "\n";
		std::cout << "Top surface points=" << top_point_count << ", tris=" << top_tri_count << "\n";
		std::cout << "Boundary faces=" << boundary_face_count << "\n";
	}
	catch (const std::exception &e)
	{
		std::cerr << "vmesh main failed: " << e.what() << "\n";
		return 1;
	}

	releaseArray(volume_coords);
	releaseArray(volume_cells);
	releaseArray(volume_cell_types);
	releaseArray(top_coords);
	releaseArray(top_cell_types);
	releaseArray(top_tris);
	releaseArray(local_to_global);
	releaseArray(point_sizing);
	releaseArray(boundary_mesh);
	releaseArray(boundary_face);

	return 0;
}
