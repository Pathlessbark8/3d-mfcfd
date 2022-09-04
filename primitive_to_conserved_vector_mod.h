/*  
	MFCFD is a 3D Computational Fluid Dynamics Solver based off q-LSKUM
    Copyright (C) 2022 Dhruv Saxena
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once
//
//
	#include "data_structure_mod.h"
//
//
//
		void primitive_to_conserved(int k,double *U) 
//
			{
//		
			
			double temp1, temp2, temp3;
//		
//			
			U[0] = point.prim[0][k];
			temp1 = point.prim[0][k]*point.prim[1][k];
			temp2 = point.prim[0][k]*point.prim[2][k];
			temp3 = point.prim[0][k]*point.prim[3][k];
//			U[4] = 2.50*point.prim[4][k] + 0.50*U[0]*(temp1*temp1 + temp2*temp2 + temp3*temp3)
//
//			The momentum components in the rotational frame ..
//
			U[1] = temp1*point.tan1[0][k] + temp2*point.tan1[1][k] + temp3*point.tan1[2][k];
			U[2] = temp1*point.tan2[0][k] + temp2*point.tan2[1][k] + temp3*point.tan2[2][k];
			U[3] = temp1*point.nor[0][k] + temp2*point.nor[1][k] + temp3*point.nor[2][k];
//
			temp1 = point.prim[1][k]*point.prim[1][k] + point.prim[2][k]*point.prim[2][k] + point.prim[3][k]*point.prim[3][k];
			U[4] = 2.50*point.prim[4][k] + 0.50*U[0]*temp1;
//	
		}
//		



__device__ void primitive_to_conserved_cuda(points &point,int k,double *U) 
//
			{
//		
			
			double temp1, temp2, temp3;
//		
//			
			U[0] = point.prim[0][k];
			temp1 = point.prim[0][k]*point.prim[1][k];
			temp2 = point.prim[0][k]*point.prim[2][k];
			temp3 = point.prim[0][k]*point.prim[3][k];
//
//			The momentum components in the rotational frame ..
//
			U[1] = temp1*point.tan1[0][k] + temp2*point.tan1[1][k] + temp3*point.tan1[2][k];
			U[2] = temp1*point.tan2[0][k] + temp2*point.tan2[1][k] + temp3*point.tan2[2][k];
			U[3] = temp1*point.nor[0][k] + temp2*point.nor[1][k] + temp3*point.nor[2][k];
//
			temp1 = point.prim[1][k]*point.prim[1][k] + point.prim[2][k]*point.prim[2][k] + point.prim[3][k]*point.prim[3][k];
			U[4] = 2.50*point.prim[4][k] + 0.50*U[0]*temp1;
//	
		}

__device__ void primitive_to_conserved_multi_nccl(splitPoints *splitPoint,int k,double *U) 
//
			{
//		
			
			double temp1, temp2, temp3;
//		
//			
			U[0] = splitPoint[k].prim[0];
			temp1 = splitPoint[k].prim[0]*splitPoint[k].prim[1];
			temp2 = splitPoint[k].prim[0]*splitPoint[k].prim[2];
			temp3 = splitPoint[k].prim[0]*splitPoint[k].prim[3];
//
//			The momentum components in the rotational frame ..
//
			U[1] = temp1*splitPoint[k].tan1[0] + temp2*splitPoint[k].tan1[1] + temp3*splitPoint[k].tan1[2];
			U[2] = temp1*splitPoint[k].tan2[0] + temp2*splitPoint[k].tan2[1] + temp3*splitPoint[k].tan2[2];
			U[3] = temp1*splitPoint[k].nor[0] + temp2*splitPoint[k].nor[1] + temp3*splitPoint[k].nor[2];
//
			temp1 = splitPoint[k].prim[1]*splitPoint[k].prim[1] + splitPoint[k].prim[2]*splitPoint[k].prim[2] + splitPoint[k].prim[3]*splitPoint[k].prim[3];
			U[4] = 2.50*splitPoint[k].prim[4] + 0.50*U[0]*temp1;
//	
		}
