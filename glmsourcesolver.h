

/***********************************************************************
 *
 *	glmsourcesolver.h
 *
 *	GLM-MHD source solver by using an operator-splitting approach.
 *
 *
 *	i, j, k : loop index of r, phi, z
 *
 *	
 *	2012 Nov. 06 : combine half step & full step.
 *	2012 Nov. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int GLMSourceSolver(double dummyU[][ixmax][jymax][kzmax], 
										double dummyV[][ixmax][jymax][kzmax])
{		

	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){		
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				//-----psi^(n+1/2)=e^(-dt/2*ch/cr)*psi^(n)-----
				
				dummyU[8][i][j][k]=psicoef*dummyU[8][i][j][k];
				dummyV[8][i][j][k]=psicoef*dummyV[8][i][j][k];
				/*
				dummyU[8][i][j][k]*=psicoef;
				dummyV[8][i][j][k]*=psicoef;
				*/
				
			}
		}
	}
	
	
	return 0;

}


