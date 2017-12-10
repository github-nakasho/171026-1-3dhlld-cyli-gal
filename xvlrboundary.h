

/***********************************************************************
 *
 *	xvlrboundary.h
 *
 *	left & right state boundary condition in x
 *
 *	
 *	j, k : loop index of x, y, z
 *	m : Vl, Vr vector component index
 *	
 *	
 *	2013 Apr. 15 : add XLeftVlAxisBoundary.
 *	2012 Dec. 06 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//==================================================
//	i=ixmax1-1 right state free boundary
//==================================================


int XRightVrFreeBoundary(void)
{	
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr free boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<9; m++){
   
    for(j=0; j<jymax1; j++){
      for(k=0; k<kzmax1; k++){
			
				Vr[m][ixmax1-1][j][k]=Vl[m][ixmax1-1][j][k];
				
				
			}
		}
    
	}
	
	
	return 0;
	
}



//==================================================
//	i=0 left state axis boundary
//==================================================


int XLeftVlAxisBoundary(void)
{		
	
	int j, k;
	
	
	//-----Vl axis boundary-----
	
#pragma omp parallel for private(j, k)
	for(j=0; j<jy/2+1; j++){
		for(k=0; k<kzmax1; k++){
			
			Vl[0][0][j][k]=Vr[0][0][j+jy/2][k];
			Vl[1][0][j][k]=-Vr[1][0][j+jy/2][k];
			Vl[2][0][j][k]=-Vr[2][0][j+jy/2][k];
			Vl[3][0][j][k]=Vr[3][0][j+jy/2][k];
			Vl[4][0][j][k]=-Vr[4][0][j+jy/2][k];
			Vl[5][0][j][k]=-Vr[5][0][j+jy/2][k];
			Vl[6][0][j][k]=Vr[6][0][j+jy/2][k];
			Vl[7][0][j][k]=Vr[7][0][j+jy/2][k];
			Vl[8][0][j][k]=Vr[8][0][j+jy/2][k];
			
			
		}
	}
#pragma omp parallel for private(j, k)
	for(j=jy/2+1; j<jymax1; j++){
		for(k=0; k<kzmax1; k++){
			
			Vl[0][0][j][k]=Vr[0][0][j-jy/2][k];
			Vl[1][0][j][k]=-Vr[1][0][j-jy/2][k];
			Vl[2][0][j][k]=-Vr[2][0][j-jy/2][k];
			Vl[3][0][j][k]=Vr[3][0][j-jy/2][k];
			Vl[4][0][j][k]=-Vr[4][0][j-jy/2][k];
			Vl[5][0][j][k]=-Vr[5][0][j-jy/2][k];
			Vl[6][0][j][k]=Vr[6][0][j-jy/2][k];
			Vl[7][0][j][k]=Vr[7][0][j-jy/2][k];
			Vl[8][0][j][k]=Vr[8][0][j-jy/2][k];
			
			
		}
	}
	
  
  /*
#pragma omp parallel for private(j, k)
  for(j=0; j<jy/2+1; j++){
    for(k=0; k<kzmax1; k++){
      
      Vl[0][0][j][k]=Vr[0][0][j+jy/2][k];
      Vl[1][0][j][k]=Vr[1][0][j+jy/2][k];
      Vl[2][0][j][k]=Vr[2][0][j+jy/2][k];
      Vl[3][0][j][k]=Vr[3][0][j+jy/2][k];
      Vl[4][0][j][k]=-Vr[4][0][j+jy/2][k];
      Vl[5][0][j][k]=-Vr[5][0][j+jy/2][k];
      Vl[6][0][j][k]=Vr[6][0][j+jy/2][k];
      Vl[7][0][j][k]=Vr[7][0][j+jy/2][k];
      Vl[8][0][j][k]=Vr[8][0][j+jy/2][k];
      
      
    }
  }
#pragma omp parallel for private(j, k)
  for(j=jy/2+1; j<jymax1; j++){
    for(k=0; k<kzmax1; k++){
      
      Vl[0][0][j][k]=Vr[0][0][j-jy/2][k];
      Vl[1][0][j][k]=Vr[1][0][j-jy/2][k];
      Vl[2][0][j][k]=Vr[2][0][j-jy/2][k];
      Vl[3][0][j][k]=Vr[3][0][j-jy/2][k];
      Vl[4][0][j][k]=-Vr[4][0][j-jy/2][k];
      Vl[5][0][j][k]=-Vr[5][0][j-jy/2][k];
      Vl[6][0][j][k]=Vr[6][0][j-jy/2][k];
      Vl[7][0][j][k]=Vr[7][0][j-jy/2][k];
      Vl[8][0][j][k]=Vr[8][0][j-jy/2][k];
      
      
    }
  }
*/
	
	return 0;
	
}



//==================================================
//	left & right state periodic boundary
//==================================================


int XVlrPeriodicBoundary(void)
{		
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr periodic boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<9; m++){
   
    for(j=0; j<jymax1; j++){
      for(k=0; k<kzmax1; k++){
				
				Vl[m][0][j][k]=Vl[m][ixmax1-1][j][k];
				Vr[m][ixmax1-1][j][k]=Vr[m][0][j][k];
				
        
			}
		}
    
	}
	
	
	return 0;
	
}



//==================================================
//	left & right state free boundary
//==================================================


int XVlrFreeBoundary(void)
{		
	
	int j, k;
	int m;
	
	
	//-----Vl, Vr free boundary in x-----
	
#pragma omp parallel for private(j, k, m)
	for(m=0; m<9; m++){
   
    for(j=0; j<jymax1; j++){
      for(k=0; k<kzmax1; k++){
			
				Vl[m][0][j][k]=Vr[m][0][j][k];
				Vr[m][ixmax1-1][j][k]=Vl[m][ixmax1-1][j][k];
			
        
			}
		}
    
	}
	
	
	return 0;
	
}







