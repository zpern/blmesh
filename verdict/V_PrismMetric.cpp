
/*
 *
 * WedgeMetric.cpp contains quality calculations for wedges
 *
 * Copyright (C) 2003 Sandia National Laboratories <cubit@sandia.gov>
 *
 * This file is part of VERDICT
 *
 * This copy of VERDICT is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * VERDICT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


#include <spdlog/spdlog.h> 
 #include "verdict.h"
#include "VerdictVector.hpp"
#include <memory.h> 

/*
  the wedge element

        4
        ^
       / \  
      / | \   
     / /1\ \       
   5/_______\6
    | /   \ |
    |/_____\|
   2         3
  
*/

/*!
  
  calculate the volume of a wedge

  this is done by dividing the wedge into 3 tets
  and summing the volume of each tet

*/

#if 0
double* tet_volume_4(int i0, int i1, int i2, int i3, VERDICT_REAL coordinates[][3], double v4[][3]) {

	for (int i = 0; i < 3; i++){
		v4[0][i] = coordinates[i0][i];
		v4[1][i] = coordinates[i1][i];
		v4[2][i] = coordinates[i2][i];
		v4[3][i] = coordinates[i3][i];
	}

	return *v4;
}
#endif

C_FUNC_DEF VERDICT_REAL v_prism_volume( int num_nodes, VERDICT_REAL coordinates[][3] )
{

  double volume = 0;
  VerdictVector side1, side2, side3;

  if ( num_nodes == 6 )
  {
    // divide the prism into 3 tets and calculate each volume

    side1.set( coordinates[1][0] - coordinates[0][0],
               coordinates[1][1] - coordinates[0][1],
               coordinates[1][2] - coordinates[0][2]);

    side2.set( coordinates[2][0] - coordinates[0][0],
               coordinates[2][1] - coordinates[0][1],
               coordinates[2][2] - coordinates[0][2]);


    side3.set( coordinates[3][0] - coordinates[0][0],
               coordinates[3][1] - coordinates[0][1],
               coordinates[3][2] - coordinates[0][2]);

    volume = side3 % (side1 * side2) / 6;

    side1.set( coordinates[4][0] - coordinates[1][0],
               coordinates[4][1] - coordinates[1][1],
               coordinates[4][2] - coordinates[1][2]);

    side2.set( coordinates[5][0] - coordinates[1][0],
               coordinates[5][1] - coordinates[1][1],
               coordinates[5][2] - coordinates[1][2]);


    side3.set( coordinates[3][0] - coordinates[1][0],
               coordinates[3][1] - coordinates[1][1],
               coordinates[3][2] - coordinates[1][2]);

    volume += side3 % (side1 * side2) / 6;

    side1.set( coordinates[5][0] - coordinates[1][0],
               coordinates[5][1] - coordinates[1][1],
               coordinates[5][2] - coordinates[1][2]);

    side2.set( coordinates[2][0] - coordinates[1][0],
               coordinates[2][1] - coordinates[1][1],
               coordinates[2][2] - coordinates[1][2]);


    side3.set( coordinates[3][0] - coordinates[1][0],
               coordinates[3][1] - coordinates[1][1],
               coordinates[3][2] - coordinates[1][2]);

    volume += side3 % (side1 * side2) / 6;
  }

  return (VERDICT_REAL)volume;
}

/*!
  the equiangle skewness of a wedge
*/
C_FUNC_DEF VERDICT_REAL v_prism_eqangle_skew(int num_nodes, VERDICT_REAL coordinates[][3])
{
	if (num_nodes != 6)
		return -1;

	double tri012[][3] = { {coordinates[0][0], coordinates[0][1], coordinates[0][2]},
					       {coordinates[1][0], coordinates[1][1], coordinates[1][2]},
					       {coordinates[2][0], coordinates[2][1], coordinates[2][2]} };
	double tri345[][3] = { {coordinates[3][0], coordinates[3][1], coordinates[3][2]},
						   {coordinates[4][0], coordinates[4][1], coordinates[4][2]},
					       {coordinates[5][0], coordinates[5][1], coordinates[5][2]} };
	double quad0143[4][3] = { {coordinates[0][0], coordinates[0][1], coordinates[0][2]},
						      {coordinates[1][0], coordinates[1][1], coordinates[1][2]},
						      {coordinates[4][0], coordinates[4][1], coordinates[4][2]},
						      {coordinates[3][0], coordinates[3][1], coordinates[3][2]} };
	double quad0253[4][3] = { {coordinates[0][0], coordinates[0][1], coordinates[0][2]},
							  {coordinates[2][0], coordinates[2][1], coordinates[2][2]},
							  {coordinates[5][0], coordinates[5][1], coordinates[5][2]},
						      {coordinates[3][0], coordinates[3][1], coordinates[3][2]} };
	double quad1254[4][3] = { {coordinates[1][0], coordinates[1][1], coordinates[1][2]},
							  {coordinates[2][0], coordinates[2][1], coordinates[2][2]},
							  {coordinates[5][0], coordinates[5][1], coordinates[5][2]},
							  {coordinates[4][0], coordinates[4][1], coordinates[4][2]} };

	double skew[5];
	skew[0] = v_tri_eqangle_skew(3, tri012);
	skew[1] = v_tri_eqangle_skew(3, tri345);
	skew[2] = v_quad_eqangle_skew(4, quad0143);
	skew[3] = v_quad_eqangle_skew(4, quad0253);
	skew[4] = v_quad_eqangle_skew(4, quad1254);

	double skew_max = skew[0];
	for (int i = 1; i < 5; i++)
	{
		if (skew[i] > skew_max)
			skew_max = skew[i];
	}

	return skew_max;
}

C_FUNC_DEF void v_prism_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
    unsigned int metrics_request_flag, PrismMetricVals *metric_vals )
{
  memset( metric_vals, 0, sizeof(PrismMetricVals) );

  // calculate volume
  if(metrics_request_flag & V_WEDGE_VOLUME)
    metric_vals->volume = v_prism_volume(num_nodes, coordinates);

  // calculate the equiangle skewness, added by BHWang
  if (metrics_request_flag & V_WEDGE_ALL)
  {
	  metric_vals->eqangle_skew = v_prism_eqangle_skew(num_nodes, coordinates);
  }
}

