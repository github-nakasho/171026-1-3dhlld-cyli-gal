

/***********************************************************************
 *
 *	cboundary.h
 *
 *	set central boundary in simulation region.
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : dummyU vector component index
 *	dummyU : U or U1 (conserved variables)
 *
 *
 *  2015 July 20: using absorb_factor array
 *	2013 Apr. 13: coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 


//============================================================
//	central absourbing boundary condition
//============================================================


int CentralAbsorbBoundary(double dummyU[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	int m;
	
	
#pragma omp parallel for private(i, j, k, m)
  for(m=0; m<9; m++){
    for(i=0; i<irin+1; i++){
      for(j=0; j<jymax; j++){
        for(k=0; k<kzmax; k++){
			
					dummyU[m][i][j][k]-=(c_absorb_factor[i][k]
                               *(dummyU[m][i][j][k]-Uinitial[m][i][j][k]));
					
					
				}
			}
		}
	}
	
	
	return 0;
	
}





