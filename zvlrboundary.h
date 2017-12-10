

/***********************************************************************
 *
 *	yvlrboundary.h
 *
 *	left & right state boundary condition in z
 *
 *	
 *	i, j : loop index of x, y
 *	m : Vl, Vr vector component index
 *	
 *
 *	2013 Apr. 15 : add ZRightVlrFreeBoundary.
 *	2013 Apr. 15 : add ZLeftVlrVzAntisymBsymBoundary.
 *	2013 Jan. 11 : add fixed boundary.
 *	2012 Dec. 13 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//==================================================
//	k=kzmax1-1 right state free boundary
//==================================================


int ZRightVrFreeBoundary(void)
{		
	
	int i, j;
	int m;
	
	
	//-----Vr free boundary-----
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=0; j<jymax1; j++){
				
				Vr[m][i][j][kzmax1-1]=Vl[m][i][j][kzmax1-1];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}



//======================================================================
//	k=0 left state vz:antisymmetry, bz:symmetry boundary
//======================================================================


int ZLeftVlVzAntisymBzsymBoundary(void)
{		
	
	int i, j;
  
	
	//-----Vl antisymmetry & symmetry boundary-----
	
#pragma omp parallel for private(i, j)
	for(i=0; i<ixmax1; i++){
		for(j=0; j<jymax1; j++){
			
			Vl[0][i][j][0]=Vr[0][i][j][1];
			Vl[1][i][j][0]=Vr[1][i][j][1];
			Vl[2][i][j][0]=Vr[2][i][j][1];
			Vl[3][i][j][0]=-Vr[3][i][j][1];
			Vl[4][i][j][0]=Vr[4][i][j][1];
			Vl[5][i][j][0]=Vr[5][i][j][1];
			Vl[6][i][j][0]=Vr[6][i][j][1];
			Vl[7][i][j][0]=Vr[7][i][j][1];
			Vl[8][i][j][0]=Vr[8][i][j][1];
			
			
		}
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state periodic boundary
//==================================================


int ZVlrPeriodicBoundary(void)
{		
	
	int i, j;
	int m;
	
	
	//-----Vl, Vr periodic boundary in z-----
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=0; j<jymax1; j++){
			
				Vl[m][i][j][0]=Vl[m][i][j][kzmax1-1];
				Vr[m][i][j][kzmax1-1]=Vr[m][i][j][0];
			
        
			}
		}
    
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state free boundary
//==================================================


int ZVlrFreeBoundary(void)
{		
	
	int i, j;
	int m;
	
	
	//-----Vl, Vr free boundary in z-----
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=0; j<jymax1; j++){
			
				Vl[m][i][j][0]=Vr[m][i][j][0];
				Vr[m][i][j][kzmax1-1]=Vl[m][i][j][kzmax1-1];
			
        
			}
		}
    
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state fixed boundary
//==================================================


int ZVlrFixedBoundary(void)
{		
	
	int i, j;
	int m;
	
	
	//-----Vl, Vr fixed boundary in z-----
	
#pragma omp parallel for private(i, j, m)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=0; j<jymax1; j++){
			
				Vl[m][i][j][0]=Vinitial[m][i][j][0];
				Vr[m][i][j][kzmax1-1]=Vinitial[m][i][j][kzmax1];
			
        
			}
		}
    
	}
	
	
	return 0;
	
}




