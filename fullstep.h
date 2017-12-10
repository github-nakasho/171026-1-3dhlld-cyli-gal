

/***********************************************************************
 *
 *	fullstep.h
 *
 *	TVD Runge-Kutta scheme full step time marching function
 *
 *
 *	i, j, k : loop index of x, y, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int Fullstep(void)
{	

	int i, j, k;
	int m;
	
	
	//============================================================
	//	full step GLM source term solver
	//============================================================
	
	GLMSourceSolver(U1, V);
	
	
	//============================================================
	//	2nd order slope limiter in x
	//============================================================

	MinmodLimiterX();
	//VanLeerLimiterX();
	
	
	//============================================================
	//	HLL flux in x-direction
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
  //NoSource(V, S);
	
	
	//-----U^(n+1)[i]=(U^n[i]+U^(n+1/2)[i]
	//							-(dt/dx)(F^(n+1/2)[i]-F^(n+1/2)[i-1])
	//							-(dt/dy)(G^(n+1/2)[j]-G^(n+1/2)[j-1]))/2-----
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<9; m++){
		for(i=1; i<ixmax1; i++){
			for(j=1; j<jymax1; j++){
				for(k=1; k<kzmax1; k++){
				
					U[m][i][j][k]=0.5*(U[m][i][j][k]+U1[m][i][j][k]
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
	
	
	//-----normal boundary conditions-----
	
	//============================================================
	//	j=0 & j=jymax1 periodic boundary
	//============================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	i=0 axial boundary
	//==================================================
	
	XLeftAxisBoundary(U);
	
	
	//==================================================
	//	i=ixmax1 free boundary
	//==================================================
	
	XRightFreeBoundary(U);
	
	
	//============================================================
	//	k=0 vz=antisymmetry, bx by bz=symmetry boundary
	//============================================================
	
	ZLeftVzAntisymBSymBoundary(U);
	//ZLeftFreeBoundary(U);
	//ZLeftPeriodicBoundary(U);
	
	
	//============================================================
	//	k=kzmax1 free boundary
	//============================================================
	
	ZRightFreeBoundary(U);
	//ZRightPeriodicBoundary(U);
	
	
	//-----abnormal boundary conditions-----
	
	//============================================================
	//	absorging boundary near i=ixmax1
	//============================================================
	
	XRightAbsorbBoundary(U);
	
	
	//============================================================
	//	absorging boundary near central region
	//============================================================
	
	CentralAbsorbBoundary(U);
	
	
	//============================================================
	//	absorging boundary near k=kzmax1
	//============================================================
	
	ZRightAbsorbBoundary(U);
	
	
	//==================================================
	//	convert from U to V.
	//==================================================
	
	UtoV(U, V);
	
	
	//============================================================
	//	negative density or pressure check
	//============================================================
	
	Check(U, V);
	
	
	return 0;
  
}















