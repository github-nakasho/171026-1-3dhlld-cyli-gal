

/***********************************************************************
 *
 *	zboundary.h
 *
 *	set k=0, k=kz boundary condition
 *
 *
 *	i, j : loop index of x, y
 *	m : U, U1 vector index
 *	dummyU : U or U1 (conserved variables)
 *
 *
 *	2013 Apr. 13 : add ZLeftVZAntisymBSymBoundary, ZRightAbsorbBoundary.
 *	2012 Dec. 17 : add fixed boundary.
 *	2012 Dec. 06 : separate left and right side condition 
 *									of simulation region .
 *	2012 Nov. 06 : use pointer & reduce useless functions.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 


//============================================================
//	k=kzmax1 absorbing boundary condition
//============================================================


int ZRightAbsorbBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	int m;
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<9; m++){
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
        for(k=kzmax1-10; k<kzmax; k++){
				
					dummyU[m][i][j][k]-=(z_absorb_factor[k-(kzmax1-10)]
                               *(dummyU[m][i][j][k]-Uinitial[m][i][j][k]));
				
				
				}
			}
		}
    
	}
	
	
	return 0;
	
}



//============================================================
//	k=0 vz:antisymmetry, B:symmetry boundary condition
//============================================================


int ZLeftVzAntisymBSymBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	
	
#pragma omp parallel for private(i, j)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			
			dummyU[0][i][j][0]=dummyU[0][i][j][2];
			dummyU[1][i][j][0]=dummyU[1][i][j][2];
			dummyU[2][i][j][0]=dummyU[2][i][j][2];
			dummyU[3][i][j][0]=-dummyU[3][i][j][2];
			dummyU[4][i][j][0]=dummyU[4][i][j][2];
			dummyU[5][i][j][0]=dummyU[5][i][j][2];
			dummyU[6][i][j][0]=dummyU[6][i][j][2];
			dummyU[7][i][j][0]=dummyU[7][i][j][2];
			dummyU[8][i][j][0]=dummyU[8][i][j][2];
				
			dummyU[3][i][j][1]=0.0;
				
			
		}
	}
	
	
	return 0;
	
}



//============================================================
//	k=0 periodic boundary condition
//============================================================


int ZLeftPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
			
				dummyU[m][i][j][0]=dummyU[m][i][j][kzmax1-1];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}


//============================================================
//	k=kzmax1 periodic boundary condition
//============================================================


int ZRightPeriodicBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
				
				dummyU[m][i][j][kzmax1]=dummyU[m][i][j][1];
				
				
			}
    }
    
  }
	
	
	return 0;
	
}



//============================================================
//	k=0 free boundary condition 
//============================================================


int ZLeftFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
   
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
			
        dummyU[m][i][j][0]=dummyU[m][i][j][1];
				
      
			}
    }
    
  }
	
	
	return 0;
	
}



//============================================================
//	k=kzmax1 free boundary condition 
//============================================================


int ZRightFreeBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	int m;
	
	
#pragma omp parallel for private(i, j, m)
    for(m=0; m<9; m++){
        for(i=0; i<ixmax; i++){
            for(j=0; j<jymax; j++){
			
				
				dummyU[m][i][j][kzmax1]=dummyU[m][i][j][kzmax1-1];
				
			}
			
		}
	}
	
	
	return 0;
	
}



//============================================================
//	k=0 fixed boundary condition 
//============================================================


int ZLeftFixedBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
			
				dummyU[m][i][j][0]=Uinitial[m][i][j][0];
				
        
      }
		}
    
	}
	
	
	return 0;
	
}


//============================================================
//	k=kzmax1 fixed boundary condition 
//============================================================


int ZRightFixedBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, m;
	
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax; i++){
      for(j=0; j<jymax; j++){
			
				dummyU[m][i][j][kzmax1]=Uinitial[m][i][j][kzmax1];
			
        
			}
		}
    
	}
	
	
	return 0;
	
}


//============================================================
//	k=0 reflected boundary condition 
//============================================================


int ZLeftReflectedBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	
	
#pragma omp parallel for private(i, j)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			
			dummyU[0][i][j][0]=dummyU[0][i][j][1];
			dummyU[1][i][j][0]=dummyU[1][i][j][1];
			dummyU[2][i][j][0]=dummyU[2][i][j][1];
			dummyU[3][i][j][0]=-dummyU[3][i][j][1];
			dummyU[4][i][j][0]=dummyU[4][i][j][1];
			dummyU[5][i][j][0]=dummyU[5][i][j][1];
			dummyU[6][i][j][0]=-dummyU[6][i][j][1];
			dummyU[7][i][j][0]=dummyU[7][i][j][1];
			dummyU[8][i][j][0]=dummyU[8][i][j][1];
			
			
		}
	}
	
	
	return 0;
	
}


//============================================================
//	k=kzmax1 fixed boundary condition 
//============================================================


int ZRightReflectedBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j;
	
	
#pragma omp parallel for private(i, j)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			
			dummyU[0][i][j][kzmax1]=dummyU[0][i][j][kzmax1-1];
			dummyU[1][i][j][kzmax1]=dummyU[1][i][j][kzmax1-1];
			dummyU[2][i][j][kzmax1]=dummyU[2][i][j][kzmax1-1];
			dummyU[3][i][j][kzmax1]=-dummyU[3][i][j][kzmax1-1];
			dummyU[4][i][j][kzmax1]=dummyU[4][i][j][kzmax1-1];
			dummyU[5][i][j][kzmax1]=dummyU[5][i][j][kzmax1-1];
			dummyU[6][i][j][kzmax1]=-dummyU[6][i][j][kzmax1-1];
			dummyU[7][i][j][kzmax1]=dummyU[7][i][j][kzmax1-1];
			dummyU[8][i][j][kzmax1]=dummyU[8][i][j][kzmax1-1];
			
			
		}
	}
	
	
	return 0;
	
}





