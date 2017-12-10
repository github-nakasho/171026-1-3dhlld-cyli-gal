

/***********************************************************************
 *
 *	cartesiangrid.h
 *	
 *	making Cartesian coordinate grid.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	oneminlength=1.0/minlength
 *
 *	
 *	2012 Oct. 03 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int CartesianGrid(void)
{

	int i, j, k;
	double oneominlength;
	
	
	//-----x, no stretch renght-----
	
	dx[0]=dx0;
	x[0]=x0-dx0;
	
	for(i=0; i<ixmax; i++){
		
		x[i+1]=x[i]+dx[i];
		dx[i+1]=dx[i];
		
		
	}
	

	//-----y, no stretch renght-----
	
	dy=dy0;
	y[0]=y0-dy0;
	
	for(j=0; j<jymax; j++) y[j+1]=y[j]+dy;
	

	//-----z, no stretch renght-----
	
	dz[0]=dz0;
	z[0]=z0-dz0;
	
	for(k=0; k<kzmax; k++){
		
		z[k+1]=z[k]+dz[k];
		dz[k+1]=dz[k];
		
		
	}	
	
	
	oneominminlength=0.0;
	
	
#pragma omp parallel for private(i, j, k, oneominlength) shared(oneominminlength)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				//-----searching minimum length of [i][j][k] cell-----
				
				minlength[i][j][k]=min(dx[i], min(dy, dz[k]));
				
				
				//-----searching minimum length of all cell-----
				
				oneominlength=1.0/minlength[i][j][k];
				
				
#pragma omp flush(oneominminlength)
				if(oneominlength>oneominminlength)	oneominminlength=oneominlength;
				
				
			}
		}
	}

	
	minminlength=1.0/oneominminlength;
	
	
	return 0;

}



