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
