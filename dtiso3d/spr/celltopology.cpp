/* ----------------------------------------------------------------------------
*
* 多学科应用模拟的赋能环境
* Enabling Environment for Multi-displinary Application Simulations
*
* 陈建军 中国 浙江大学工程与科学计算研究中心
* 版权所有	  2007年10月11日
* Chen Jianjun  Center for Engineering & Scientific Computation,
* Zhejiang University, P. R. China
* Copyright reserved, Dec. 25, 2008
* 
* 联系方式 (For further information, please conctact)
*   电话 (Tel)：+86-571-87953165
*   传真 (Fax)：+86-571-87953167
*   邮箱 (Mail)：chenjj@zju.edu.cn
*
* 文件名称 (File Name)：celltopology.cpp
* 初始版本 (Initial Version): V1.0
* 功能介绍 (Function Introduction：
		*     实现了查询单元拓扑信息的方法和数据结构
		*     Implement a set of data structures & algorithms for 
		*     cell topologies query
		* 
		*
		* -----------------------------修改记录 (Revision Record)------------------------
		* 修改者 (Revisor):
* 修改日期 (Revision Date):
* 当前版本 (Current Version):
* 修改介绍 (Revision Introduction):
* ------------------------------------------------------------------------------
* ------------------------------------------------------------------------------*/

#include <spdlog/spdlog.h> 
 #include "celltopology.h"

#include <string.h>
#include <assert.h>

namespace  EEMAS {

const TopologyInfo TopologyInfo::instance;


TopologyInfo::TopologyInfo()
{
  memset( dimMap, 0, sizeof(dimMap) );
  memset( adjMap, 0, sizeof(adjMap) );
  memset( edgeMap, 0, sizeof(edgeMap) );
  memset( faceMap, 0, sizeof(faceMap) );
  memset( vertAdjMap, 0, sizeof(vertAdjMap) );
  
  dimMap[NODE]			= 0;
  dimMap[LINE]			= 1;
  dimMap[POLYGON ]      = 2;
  dimMap[TRIANGLE]      = 2;
  dimMap[QUADRILATERAL] = 2;
  dimMap[POLYHEDRON]    = 3;
  dimMap[TETRAHEDRON]   = 3;
  dimMap[HEXAHEDRON]    = 3;
  dimMap[PRISM]         = 3;
  dimMap[PYRAMID]       = 3;
  dimMap[SEPTAHEDRON]   = 3;
  
  adjMap[LINE][0] = 2;

  adjMap[TRIANGLE][0] = 3;
  adjMap[TRIANGLE][1] = 3;
    
  adjMap[QUADRILATERAL][0] = 4;
  adjMap[QUADRILATERAL][1] = 4;

  adjMap[TETRAHEDRON][0] = 4;
  adjMap[TETRAHEDRON][1] = 6;
  adjMap[TETRAHEDRON][2] = 4;
  
  adjMap[HEXAHEDRON][0] = 8;
  adjMap[HEXAHEDRON][1] = 12;
  adjMap[HEXAHEDRON][2] = 6;
  
  adjMap[PRISM][0] = 6;
  adjMap[PRISM][1] = 9;
  adjMap[PRISM][2] = 5;
  
  adjMap[PYRAMID][0] = 5;
  adjMap[PYRAMID][1] = 8;
  adjMap[PYRAMID][2] = 5;
  
  adjMap[SEPTAHEDRON][0] = 7;
  adjMap[SEPTAHEDRON][1] = 11;
  adjMap[SEPTAHEDRON][2] = 6;  /* See description in TSTT mesh interface doc */

  int side;
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TRIANGLE-FIRST_FACE][side][0] = side;
    edgeMap[TRIANGLE-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[QUADRILATERAL-FIRST_FACE][side][0] = side;
    edgeMap[QUADRILATERAL-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side -3 ;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = 3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = side;
  }
  for (side = 8; side < 12; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = 4+(side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side;
    edgeMap[PRISM-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side - 3;
    edgeMap[PRISM-FIRST_FACE][side][1] = side;
  }
  for (side = 6; side < 9; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side-3;
    edgeMap[PRISM-FIRST_FACE][side][1] = 3+(side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side - 4;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = 4;
  }
  
  for (side = 0; side < 3; ++side)
  {
    faceMap[TETRAHEDRON-FIRST_VOL][side][0] = 3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[TETRAHEDRON-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][3] = 3;
  }
  faceMap[TETRAHEDRON-FIRST_VOL][3][0] = 3;
  faceMap[TETRAHEDRON-FIRST_VOL][3][1] = 2;
  faceMap[TETRAHEDRON-FIRST_VOL][3][2] = 1;
  faceMap[TETRAHEDRON-FIRST_VOL][3][3] = 0;

  for (side = 0; side < 4; ++side)
  {
    faceMap[HEXAHEDRON-FIRST_VOL][side][0] = 4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[HEXAHEDRON-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][3] = 4+(side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][4] = side + 4;
  }
  faceMap[HEXAHEDRON-FIRST_VOL][4][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][4][1] = 3;
  faceMap[HEXAHEDRON-FIRST_VOL][4][2] = 2;
  faceMap[HEXAHEDRON-FIRST_VOL][4][3] = 1;
  faceMap[HEXAHEDRON-FIRST_VOL][4][4] = 0;
  faceMap[HEXAHEDRON-FIRST_VOL][5][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][1] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][2] = 5;
  faceMap[HEXAHEDRON-FIRST_VOL][5][3] = 6;
  faceMap[HEXAHEDRON-FIRST_VOL][5][4] = 7;
  
  for (side = 0; side < 4; ++side)
  {
    faceMap[PYRAMID-FIRST_VOL][side][0] = 3;
    faceMap[PYRAMID-FIRST_VOL][side][1] = side;
    faceMap[PYRAMID-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[PYRAMID-FIRST_VOL][side][3] = 4;
  }
  faceMap[PYRAMID-FIRST_VOL][4][0] = 4;
  faceMap[PYRAMID-FIRST_VOL][4][1] = 3;
  faceMap[PYRAMID-FIRST_VOL][4][2] = 2;
  faceMap[PYRAMID-FIRST_VOL][4][3] = 1;
  faceMap[PYRAMID-FIRST_VOL][4][4] = 0;

  for (side = 0; side < 3; ++side)
  {
    faceMap[PRISM-FIRST_VOL][side][0] = 4;
    faceMap[PRISM-FIRST_VOL][side][1] = side;
    faceMap[PRISM-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][3] = 3+(side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][4] = side + 3;
  }
  faceMap[PRISM-FIRST_VOL][3][0] = 3;
  faceMap[PRISM-FIRST_VOL][3][1] = 2;
  faceMap[PRISM-FIRST_VOL][3][2] = 1;
  faceMap[PRISM-FIRST_VOL][3][3] = 0;
  faceMap[PRISM-FIRST_VOL][4][0] = 3;
  faceMap[PRISM-FIRST_VOL][4][1] = 3;
  faceMap[PRISM-FIRST_VOL][4][2] = 4;
  faceMap[PRISM-FIRST_VOL][4][3] = 5;
  
  int i;
  for (i = 0; i < 3; ++i)
  {
    vertAdjMap[TRIANGLE-FIRST_FACE][i][0] = 2;
    vertAdjMap[TRIANGLE-FIRST_FACE][i][1] = (i+1)%3;
    vertAdjMap[TRIANGLE-FIRST_FACE][i][2] = (i+2)%3;
  }
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][0] = 2;
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][2] = (i+3)%4;
  }
  
  
  unsigned tet_corner_data[] = { 1, 2, 3, 
                                 0, 3, 2,
                                 3, 0, 1,
                                 2, 1, 0 };
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[TETRAHEDRON-FIRST_FACE][i][0] = 3;
    for (unsigned j = 0; j < 3; ++j)
      vertAdjMap[TETRAHEDRON-FIRST_FACE][i][j+1] = tet_corner_data[3*i+j];
  }
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[PYRAMID-FIRST_FACE][i][0] = 3;
    vertAdjMap[PYRAMID-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[PYRAMID-FIRST_FACE][i][2] = (i+3)%4;
    vertAdjMap[PYRAMID-FIRST_FACE][i][3] = 4;
  }
  vertAdjMap[PYRAMID-FIRST_FACE][4][0] = 4;
  for (i = 0; i < 4; i++)
    vertAdjMap[PYRAMID-FIRST_FACE][4][i+1] = 3 - i;
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][0] = 3;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][2] = (i+3)%4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][3] = i+4;
  }
  for (i = 4; i < 8; ++i)
  {
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][0] = 3;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][1] = (i+3)%4+4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][2] = (i+1)%4+4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][3] = i-4;
  }
  
  for (i = 0; i < 3; ++i)
  {
    vertAdjMap[PRISM-FIRST_FACE][i][0] = 3;
    vertAdjMap[PRISM-FIRST_FACE][i][1] = (i+1)%3;
    vertAdjMap[PRISM-FIRST_FACE][i][2] = (i+2)%3;
    vertAdjMap[PRISM-FIRST_FACE][i][3] = i+3;
  }
  for (i = 3; i < 6; ++i)
  {
    vertAdjMap[PRISM-FIRST_FACE][i][0] = 3;
    vertAdjMap[PRISM-FIRST_FACE][i][1] = (i+2)%3+3;
    vertAdjMap[PRISM-FIRST_FACE][i][2] = (i+1)%3+3;
    vertAdjMap[PRISM-FIRST_FACE][i][3] = i-3;
  }
  
    // Build reverse vertex-vertex adjacency index map
  const EntityTopology types[] = { TRIANGLE, 
                                   QUADRILATERAL, 
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM, 
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (i = 0; i < num_types; ++i)
  {
    const unsigned num_vert = corners( types[i] );
    for (unsigned v = 0; v < num_vert; ++v)
    {
      unsigned num_v_adj;
      const unsigned* v_adj = adjacent_vertices( types[i], v, num_v_adj );
      unsigned* reverse = revVertAdjIdx[types[i]-FIRST_FACE][v];
      reverse[0] = num_v_adj;
      
      for (unsigned j = 0; j < num_v_adj; ++j)
      {
        unsigned num_j_adj, k;
        const unsigned* j_adj = adjacent_vertices( types[i], v_adj[j], num_j_adj );
        for (k = 0; k < num_j_adj && j_adj[k] != v; ++k);
        assert( k < num_j_adj ); // If this fails, vertAdjMap is corrupt!
        reverse[j+1] = k;
      }
    }
  }
}

int TopologyInfo::higher_order( EntityTopology topo,
                                     unsigned num_nodes,
                                     bool& midedge,
                                     bool& midface,
                                     bool& midvol)
{
  midedge = midface = midvol = false;
  if (topo >= MIXED || num_nodes < instance.adjMap[topo][0])
  {
      return 2;	/* INVALID_ARG */
  }
  
  unsigned dim = instance.dimMap[topo];
  assert( num_nodes >= instance.adjMap[topo][0] );
  unsigned nodes = num_nodes - instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  if (edges && nodes >= edges)
  {
    nodes -= edges;
    midedge = true;
  }
  if (faces && nodes >= faces)
  {
    nodes -= faces;
    midface = true;
  }
  if (1 == nodes)
  {
    if (2 == dim)
    {
      nodes -= 1;
      midface = true;
    }
    else if(3 == dim)
    {
      nodes -= 1;
      midvol = true;
    }
  }
  
  if (nodes)
  {
    return 4;	/* E_FAIL */
  }

  return 0;	/* S_OK */
}


int TopologyInfo::edge_vertices(EntityTopology topo,
                                unsigned edge, 
                                const unsigned** ev_array)
{
  if (!ev_array ||
	  topo < (EntityTopology)FIRST_FACE || 
      topo > (EntityTopology)LAST_VOL || 
      edge >= edges( topo ) )
  {
    topo = (EntityTopology)FIRST_FACE;
    edge = 0;

	return 2;	/* E_INVALID_ARG */
  }
  
  *ev_array = instance.edgeMap[topo-FIRST_FACE][edge];

  return 0;		/* S_OK */
}

int TopologyInfo::face_vertices( EntityTopology topo,
                                 unsigned face,
                                 unsigned& length,
								 const unsigned** fv_array)
{
  if (!fv_array ||
	  topo < (EntityTopology)FIRST_VOL || 
      topo > (EntityTopology)LAST_VOL || 
      face >= faces( topo ) )
  {
    topo = (EntityTopology)FIRST_VOL;
    face = 0;

	return 2;	/* E_INVALID_ARG */
  }
  
  length = instance.faceMap[topo-FIRST_VOL][face][0];
  *fv_array = instance.faceMap[topo-FIRST_VOL][face] + 1;

  return 0;
}


int TopologyInfo::side_vertices( EntityTopology topo,
                                 unsigned dim,
                                 unsigned side,
                                 unsigned& count_out,
                                 const unsigned** sv_array)
{
  static const unsigned all[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  int nErrCode = 0;

  count_out = 0;
  if (!sv_array)
  {
	  return 2;	/* E_INVALID_ARG */
  }
  *sv_array = NULL;

  if (dim != 0 && dim == dimension(topo))
  {
    count_out = corners( topo );
    *sv_array = all;
  }
  else if (dim == 1)
  {
    count_out = 2;
    nErrCode = edge_vertices( topo, side, sv_array);
  }
  else if( dim == 2)
  {
    nErrCode = face_vertices( topo, side, count_out, sv_array);
  } 
  else
  {
   	return 2;
  }    
  
  return nErrCode;
}
      


int TopologyInfo::side_number(EntityTopology topo,
                              unsigned num_nodes,
                              unsigned node_index,
                              unsigned& side_dim_out,
                              unsigned& side_num_out)
{
  if (topo >= (EntityTopology)MIXED || num_nodes < instance.adjMap[topo][0])
  {
    return 2; /* E_INVALID_ARG */
  }
  
  unsigned nodes = instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  side_num_out = node_index;

  if (side_num_out < nodes)
  {
    side_dim_out = 0;
    return 0;	/* S_OK */
  }
  num_nodes -= nodes;
  side_num_out -= nodes;
  
  if (edges && num_nodes >= edges)
  {
    if (side_num_out < edges)
    {
      side_dim_out = 1;
      return 0;	/* S_OK */
    }
    num_nodes -= edges;
    side_num_out -= edges;
  }
  if (faces && num_nodes >= faces)
  {
    if (side_num_out < faces)
    {
      side_dim_out = 2;
      return 0;	/* S_OK */
    }
    num_nodes -= faces;
    side_num_out -= faces;
  }
  if (side_num_out == 0)
  {
    side_dim_out = instance.dimMap[topo];
    side_num_out = 0;
    return 0;	/* S_OK */
  }
  
  return 2;	/* S_OK */
}
  
  

const unsigned* TopologyInfo::adjacent_vertices( EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_adj_out )
{
  const unsigned count = corners( topo );
  if (!count || index >= count)
  {
    num_adj_out = 0;
    return 0;
  }
  
  const unsigned* vect = instance.vertAdjMap[topo-FIRST_FACE][index];
  num_adj_out = vect[0];
  return vect + 1;
}

const unsigned* TopologyInfo::reverse_vertex_adjacency_offsets(
                                              EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_adj_out )
{
  const unsigned count = corners( topo );
  if (!count || index >= count)
  {
    num_adj_out = 0;
    return 0;
  }
  
  const unsigned* vect = instance.revVertAdjIdx[topo-FIRST_FACE][index];
  num_adj_out = vect[0];
  return vect + 1;
}



} //namepsace EEMAS
