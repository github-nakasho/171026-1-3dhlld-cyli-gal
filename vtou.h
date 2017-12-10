

/***********************************************************************
 *
 *	vtou.h 
 *
 *	convert from V (primitive variables) to U (conserved variables).
 *
 *
 *	i, j, k : loop index of r, phi, z
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	b2=bx*bx+by*by+bz*bz : b**2
 *
 *	
 *	2013 Apr. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int VtoU(void)
{
	
	int i, j, k;
	double v2, b2;
	
	
#pragma omp parallel for private(i, j, k, v2, b2)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				v2=V[1][i][j][k]*V[1][i][j][k]
						+V[2][i][j][k]*V[2][i][j][k]
						+V[3][i][j][k]*V[3][i][j][k];
				
				b2=V[4][i][j][k]*V[4][i][j][k]
						+V[5][i][j][k]*V[5][i][j][k]
						+V[6][i][j][k]*V[6][i][j][k];
				
				
				//-----V--->U-----
				
				U[0][i][j][k]=V[0][i][j][k];
				U[1][i][j][k]=V[0][i][j][k]*V[1][i][j][k];
				U[2][i][j][k]=x[i]*V[0][i][j][k]*V[2][i][j][k];
				U[3][i][j][k]=V[0][i][j][k]*V[3][i][j][k];
				U[4][i][j][k]=V[4][i][j][k];
				U[5][i][j][k]=V[5][i][j][k];
				U[6][i][j][k]=V[6][i][j][k];
				U[7][i][j][k]=0.5*(V[0][i][j][k]*v2+b2)+oneogm1*V[7][i][j][k];
				U[8][i][j][k]=V[8][i][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}

