

/***********************************************************************
 *
 *	halfstep.h
 *
 *	TVD Runge-Kutta scheme half step time marching function
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2013 Oct. 13 : add CRs effect.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int Halfstep(void)
{
	
	int i, j, k;
	int m;
	
	
	//============================================================
	//	half step GLM source term solver
	//============================================================
	
	GLMSourceSolver(U, V);
	
	
	//============================================================
	//	2nd order slope limiter in x
	//============================================================
	
	MinmodLimiterX();
	//VanLeerLimiterX();
	
	
	//============================================================
	//	HLLD flux in x-direction
	//============================================================
	
	HLLDX();
	
	
	//============================================================
	//	2nd order slope limiter in y
	//============================================================
	
	MinmodLimiterY();
	//VanLeerLimiterY();
	
	
	//============================================================
	//	HLL flux in y-direction
	//============================================================
	
	HLLDY();
	
  
	//============================================================
	//	2nd order slope limiter in z
	//============================================================
	
	MinmodLimiterZ();
	
	
	//============================================================
	//	HLL flux in z-direction
	//============================================================
	
	HLLDZ();
	
	
	//============================================================
	//	calculating source term
	//============================================================
	
	SpiralGravityCooling(V, S);
  //AxiGravity(V, S);
  //AxiGravityCooling(V, S);
  
	
	//-----U^(n+1/2)=U^n-dt/dx(F[i]-F[i-1])-dt/dy(G[j]-G[j-1])-----
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<9; m++){
		for(i=1; i<ixmax1; i++){
			for(j=1; j<jymax1; j++){
				for(k=1; k<kzmax1; k++){
				
          U1[m][i][j][k]=(U[m][i][j][k]
                          -dt*((g2[i][j][k]*F[m][i][j][k]
                                -g1[i][j][k]*F[m][i-1][j][k])
                               +(g4[i][j][k]*G[m][i][j][k]
                                 -g3[i][j][k]*G[m][i][j-1][k])
                               +(g6[i][j][k]*H[m][i][j][k]
                                 -g5[i][j][k]*H[m][i][j][k-1]))
                          +dt*S[m][i][j][k]);
				
					
				}
			}
		}
	}
	
	
	//============================================================
	//	j=0 & j=jymax1 periodic boundary
	//============================================================
	
	YLeftPeriodicBoundary(U1);
	YRightPeriodicBoundary(U1);
	
	
	//==================================================
	//	i=0 axial boundary
	//==================================================
	
	XLeftAxisBoundary(U1);
	
	
	//==================================================
	//	i=ixmax1 free boundary
	//==================================================
	
	XRightFreeBoundary(U1);
	
	
	//============================================================
	//	k=0 vz=antisymmetry, bx by bz=symmetry boundary
	//============================================================
	
	ZLeftVzAntisymBSymBoundary(U1);
	//ZLeftFreeBoundary(U1);
	//ZLeftPeriodicBoundary(U1);
	
	
	//============================================================
	//	k=kzmax1 free boundary
	//============================================================
	
	ZRightFreeBoundary(U1);
	//ZRightPeriodicBoundary(U1);
	
	
	//==================================================
	//	convert from U1 to V.
	//==================================================
	
	UtoV(U1, V);
	
	
	//============================================================
	//	negative density or pressure check
	//============================================================
	
	Check(U1, V);
	
	
	return 0;
  
}












