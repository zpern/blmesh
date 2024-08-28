

#include <spdlog/spdlog.h> 
 #include "SurfaceRemesh.h"
#include <iostream>

using namespace BLMESH::zju;



SurfaceRemesh::SurfaceRemesh(Mesh &xyz_mesh):xyz_mesh_(xyz_mesh)
{
  
  mesh = new PolyMesh();
  
  for(int i = 0; i < xyz_mesh_.vertex.rows(); ++i)
  {
    mesh->addVertex(
      xyz_mesh_.vertex(i, 0),
      xyz_mesh_.vertex(i, 1),
      xyz_mesh_.vertex(i, 2));
  }

  std::vector<size_t> faceVec;
  for(int i = 0; i < xyz_mesh_.topo.rows(); ++i)
  {
    faceVec.push_back(xyz_mesh_.topo(i, 0));
    faceVec.push_back(xyz_mesh_.topo(i, 1));
    faceVec.push_back(xyz_mesh_.topo(i, 2));
    mesh->addPolyFace(faceVec);
    faceVec.clear();
  }
}

void SurfaceRemesh::split_long_edges()
{
	int NE = mesh->numEdges();
	for (int i = 0; i < NE; i++)
	{		
		MEdge* e = mesh->edge(i);
		MHalfedge* he = e->halfEdge();
		MVert* p0 = he->fromVertex();
		MVert* p1 = he->toVertex();
		if (mesh->isBoundary(p0) && mesh->isBoundary(p1)) continue;

		double size_of_p0=	sizingfunction_(p0->position()(0), p0->position()(1), p0->position()(2));
		double size_of_p1=	sizingfunction_(p1->position()(0), p1->position()(1), p1->position()(2));

		double len = e->length();
		len = std::sqrt(pow(p0->x() - p1->x(), 2) + pow(p0->y() - p1->y(), 2));
		if (2*len > size_of_p0+size_of_p1)	mesh->splitEdgeTriangle(e);
	}
}

void SurfaceRemesh::collapse_short_edges()
{
	int NE = mesh->numEdges();
	for (int i = NE-1; i >=0; i--)
	{
		if (i > mesh->numEdges() - 1)	continue;
		MEdge* e = mesh->edge(i);
		MHalfedge* he = e->halfEdge();
		MVert* p0 = he->fromVertex();
		MVert* p1 = he->toVertex();
		if (!mesh->is_collapse_ok(he)) continue;
		if (mesh->isBoundary(p0) || mesh->isBoundary(p1)) continue;
		double size_of_p0 = sizingfunction_(p0->position()(0), p0->position()(1), p0->position()(2));
		double size_of_p1 = sizingfunction_(p1->position()(0), p1->position()(1), p1->position()(2));
		double len = e->length();
		if (len < size_of_p0/2||len<size_of_p1/2)
		{
			bool is_collapse = true;
		
			if (is_collapse)	mesh->collapseTriangle(he);
		}
	}
}

void SurfaceRemesh::equalize_valences()
{
	std::vector<int> target_valence;
	int deviation_pre, deviation_post;
	for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
	{
		if (mesh->isBoundary(*v_it))	target_valence.push_back(4);
		else target_valence.push_back(6);
	}
	for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it)) continue;
		MHalfedge* he1 = (*e_it)->halfEdge();
		MVert* v1 = he1->fromVertex();
		MVert* v2 = he1->toVertex();

		MHalfedge* he2 = (*e_it)->halfEdge()->next();
		MVert* v3 = he2->toVertex();
		MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next();
        MVert* v4 = he3->toVertex();
		
		deviation_pre = abs(int(mesh->valence(v1) - target_valence[v1->index()]))
			+ abs(int(mesh->valence(v2) - target_valence[v2->index()]))
			+ abs(int(mesh->valence(v3) - target_valence[v3->index()]))
			+ abs(int(mesh->valence(v4) - target_valence[v4->index()]));
		deviation_post = abs(int(mesh->valence(v1) - 1 - target_valence[v1->index()]))
			+ abs(int(mesh->valence(v2) - 1 - target_valence[v2->index()]))
			+ abs(int(mesh->valence(v3) + 1 - target_valence[v3->index()]))
			+ abs(int(mesh->valence(v4) + 1 - target_valence[v4->index()]));
		if (deviation_pre> deviation_post)
			mesh->flipEdgeTriangle(*e_it);
	}
}

void SurfaceRemesh::tangential_relaxation()
{
	mesh->updateMeshNormal();
	for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
	{
		if(fix_id_.find((*v_it)->index())!= fix_id_.end()) continue;
		if (mesh->isBoundary(*v_it)) continue;
		double count = 0.0;
		vec3d p = (*v_it)->position();
		vec3d q(0.0, 0.0, 0.0);
		std::vector<vec3d> vs;
		for (VertexVertexIter vv_it = mesh->vv_iter(*v_it); vv_it.isValid(); ++vv_it)
		{
			vs.push_back((*vv_it)->position());
			q += (*vv_it)->position();
			++count;
		}
		q /= count;
		vec3d n = (*v_it)->normal();
		n.normalize();
		vec3d npos = q + (n.dot(p - q)) * n;
		std::array<double, 3> pos = projectfunction_(npos.x(), npos.y(), npos.z());
		npos.x() = pos[0];
		npos.y() = pos[1];
		npos.z() = pos[2];
		//npos.z() = 0;
		
		double min1 = 1e10;
		double min2 = 1e10;
		for (int k = 0; k < vs.size(); k++) {
			vec3d v1 = npos;
			vec3d v2 = vs[k];
			vec3d v3 = vs[(k + 1) % vs.size()];
			vec3d v4 = p;

			double detleft1 = (v1(0) - v3(0)) * (v2(1) - v3(1));
			double detright1 = (v1(1) - v3(1)) * (v2(0) - v3(0));
			double det1 = detleft1 - detright1;
			min1 = std::min(min1, -det1);

			double detleft2 = (v4(0) - v3(0)) * (v2(1) - v3(1));
			double detright2 = (v4(1) - v3(1)) * (v2(0) - v3(0));
			double det2 = detleft2 - detright2;
			min2 = std::min(min2, -det2);
		}
		if (min1 > min2) {
	//		spdlog::info(min1);
			(*v_it)->setPosition(npos);
		}

	}
}

double SurfaceRemesh::calculate_target_lenth_shortest()
{
  double edge_min = std::numeric_limits<double>::max();
  double edge_max = std::numeric_limits<double>::min();

	for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		edge_min =  std::min(edge_min, (*e_it)->length());
		edge_max =  std::max(edge_max, (*e_it)->length());
	}
  double target_edge_length = edge_min * 1.5;
  target_edge_length = std::min(target_edge_length, edge_max);
	return target_edge_length;
}

double SurfaceRemesh::calculateTargetEdgeLength()
{

	double target_edge_length = 0.0;
	for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		target_edge_length += (*e_it)->length();
	}
	target_edge_length /= mesh->numEdges();
	return target_edge_length;
}

void SurfaceRemesh::iso_remesh(const int & iternum)
{
  assert(iternum > 0);
  double target_edge_length, high, low;
	target_edge_length = calculate_target_lenth_shortest();
	high = 4.0 / 3.0 * target_edge_length;
	low = 4.0 / 5.0 * target_edge_length;
 
	for (int i = 0; i < std::min(100, iternum); i++)
	{
		split_long_edges();
//		collapse_short_edges();
		equalize_valences();
		tangential_relaxation();
	}

	for (int i = 0; i < std::min(100, iternum); i++)
	{
		equalize_valences();
		tangential_relaxation();
	}


  Eigen::MatrixXd V_out;
  Eigen::MatrixXi T_out;

  std::vector<vec3d> V_out_vec;
  std::vector<vec3i> T_out_vec;

  for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
	{
		MVert* vert = *v_it;
		vec3d p(vert->x(), vert->y(), vert->z());
		V_out_vec.push_back(p);
	}

	for (FaceIter f_it = mesh->polyfaces_begin();f_it != mesh->polyfaces_end(); ++f_it)
	{
		std::vector<int> index;
		for (FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.isValid(); ++fv_it)
		{
			MVert* fv = *fv_it;
			index.push_back(fv->index());
		}

		vec3i temp(index[0], index[1], index[2]);
    T_out_vec.push_back(temp);
	}

  V_out.resize(V_out_vec.size(), 3);
  T_out.resize(T_out_vec.size(), 3);

  for(int i = 0; i < V_out_vec.size(); ++i)
  {
    V_out.row(i) = V_out_vec[i];
  }
  for(int i = 0; i < T_out_vec.size(); ++i)
  {
    T_out.row(i) = T_out_vec[i];
  }
  
  this->xyz_mesh_.vertex = V_out;
  this->xyz_mesh_.topo = T_out;
}



