

/***********************************************************************
 *
 *	utov.h 
 *
 *	convert from U or U1 (conserved variables) 
 *																	to V (primitive variables).
 *
 *	
 *	i, j, k : loop index of r, phi, z
 *	onerho=1.0/rho
 *	b2=br*br+bphi*bphi+bz*bz : B**2
 *	r2v2=rvr*rvr+amz*amz/r**2+rvz*rvz : (rho*v)**2
 *	dummyU : U or U1
 *	dummyV : V
 *
 *
 *	2012 Nov. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int UtoV(double dummyU[][ixmax][jymax][kzmax], 
				 double dummyV[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	double oneorho, b2, r2v2;
	
	
#pragma omp parallel for private(i, j, k, oneorho, b2, r2v2)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				oneorho=1.0/dummyU[0][i][j][k];
				
				b2=dummyU[4][i][j][k]*dummyU[4][i][j][k]
						+dummyU[5][i][j][k]*dummyU[5][i][j][k]
						+dummyU[6][i][j][k]*dummyU[6][i][j][k];
				
				r2v2=dummyU[1][i][j][k]*dummyU[1][i][j][k]
							+dummyU[2][i][j][k]*dummyU[2][i][j][k]*oneox[i]*oneox[i]
							+dummyU[3][i][j][k]*dummyU[3][i][j][k];
				
				
				dummyV[0][i][j][k]=dummyU[0][i][j][k];
				dummyV[1][i][j][k]=dummyU[1][i][j][k]*oneorho;
				dummyV[2][i][j][k]=dummyU[2][i][j][k]*oneorho*oneox[i];
				dummyV[3][i][j][k]=dummyU[3][i][j][k]*oneorho;
				dummyV[4][i][j][k]=dummyU[4][i][j][k];
				dummyV[5][i][j][k]=dummyU[5][i][j][k];
				dummyV[6][i][j][k]=dummyU[6][i][j][k];
				dummyV[7][i][j][k]=gm1*(dummyU[7][i][j][k]
																-0.5*(r2v2*oneorho+b2));
				dummyV[8][i][j][k]=dummyU[8][i][j][k];
				
				
			}
		}
	}
	
	
	return 0;
	
}




