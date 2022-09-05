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
			U[0] = point.prim[k][0];
			temp1 = point.prim[k][0]*point.prim[k][1];
			temp2 = point.prim[k][0]*point.prim[k][2];
			temp3 = point.prim[k][0]*point.prim[k][3];
//			U[4] = 2.50*point.prim[k][4] + 0.50*U[0]*(temp1*temp1 + temp2*temp2 + temp3*temp3)
//
//			The momentum components in the rotational frame ..
//
			U[1] = temp1*point.tan1[k][0] + temp2*point.tan1[k][1] + temp3*point.tan1[k][2];
			U[2] = temp1*point.tan2[k][0] + temp2*point.tan2[k][1] + temp3*point.tan2[k][2];
			U[3] = temp1*point.nor[k][0] + temp2*point.nor[k][1] + temp3*point.nor[k][2];
//
			temp1 = point.prim[k][1]*point.prim[k][1] + point.prim[k][2]*point.prim[k][2] + point.prim[k][3]*point.prim[k][3];
			U[4] = 2.50*point.prim[k][4] + 0.50*U[0]*temp1;
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
			U[0] = point.prim[k][0];
			temp1 = point.prim[k][0]*point.prim[k][1];
			temp2 = point.prim[k][0]*point.prim[k][2];
			temp3 = point.prim[k][0]*point.prim[k][3];
//
//			The momentum components in the rotational frame ..
//
			U[1] = temp1*point.tan1[k][0] + temp2*point.tan1[k][1] + temp3*point.tan1[k][2];
			U[2] = temp1*point.tan2[k][0] + temp2*point.tan2[k][1] + temp3*point.tan2[k][2];
			U[3] = temp1*point.nor[k][0] + temp2*point.nor[k][1] + temp3*point.nor[k][2];
//
			temp1 = point.prim[k][1]*point.prim[k][1] + point.prim[k][2]*point.prim[k][2] + point.prim[k][3]*point.prim[k][3];
			U[4] = 2.50*point.prim[k][4] + 0.50*U[0]*temp1;
//	
		}
