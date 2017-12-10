

/***********************************************************************
 *
 *	perturbation.h
 *
 *	add perturbation to (magneto-)hydrostatistical state
 *
 *	
 *	i, j, k : loop index of r, phi, z
 *	amp : perturbation amplitude
 *
 *	2013 Apr. 13 : coded by Sho Nakamura (Tohoku Univ.)
 *
 **********************************************************************/



int Perturbation(void)
{
	
	int i, j, k;
	double amp=0.01;
	
	srand((unsigned)time(NULL));
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){

				U[2][i][j][k]*=1.0+amp*(2.0*rand()/(RAND_MAX+1.0)-1.0);
				
				
			}
		}
	}
	
	printf("perturbation set complete.\n");
	
	return 0;
}


