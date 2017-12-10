

/***********************************************************************
 *
 *	yboundary.h
 *
 *	set j=0, j=jymax1 boundary condition
 *
 *	
 *	i, k : loop index of x, z
 *	m : U, U1 vector index
 *	dummyU : U or U1 (conserved variables)
 *
 *	
 *	2013 Apr. 12 : add initializeboundary.
 *	2012 Dec. 06 : separate left and right side condition 
 *									of simulation region .
 *	2012 Nov. 06 use pointer & reduce useless functions.
 *	2012 Nov. 03 add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 add periodic boundary.
 *	2012 Oct. 05 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 


//============================================================
//	j=0 initialize boundary condition
//============================================================


int YLeftInitializeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
			
				dummyU[m][i][0][k]=Uinitial[m][i][0][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}


//============================================================
//	j=jymax1 initialize boundary condition
//============================================================


int YRightInitializeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	

#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
			
				dummyU[m][i][jymax1][k]=Uinitial[m][i][jymax1][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}



//============================================================
//	j=0 periodic boundary condition
//============================================================


int YLeftPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
				
				dummyU[m][i][0][k]=dummyU[m][i][jymax1-1][k];
				
				
			}
    }
    
	}
	
	
	return 0;
	
}



//============================================================
//	j=jymax1 periodic boundary condition
//============================================================


int YRightPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
			
				dummyU[m][i][jymax1][k]=dummyU[m][i][1][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}



//============================================================
//	j=0 free boundary condition
//============================================================


int YLeftFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
				
				dummyU[m][i][0][k]=dummyU[m][i][1][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}



//============================================================
//	j=jymax1 free boundary condition
//============================================================


int YRightFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, k;
	int m;
	
	
#pragma omp parallel for private(i, k, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(k=0; k<kzmax; k++){
				
				dummyU[m][i][jymax1][k]=dummyU[m][i][jymax1-1][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}





