

/***********************************************************************
 *
 *	firststep.h
 *
 *	TVD Runge-Kutta scheme first step time marching function
 *
 *	
 *	i, j, k : loop index of r, phi, z
 *	m : U, U1, F, G vector component index
 *
 *	
 *	2012 Dec. 26 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int FirstStep(void)
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
	
	//VanLeerLimiterX();
	MinmodLimiterX();
	
	
	//============================================================
	//	HLLD flux in x-direction
	//============================================================
	
	HLLDX();
	//HLLX();
	
	
	//============================================================
	//	2nd order slope limiter in y
	//============================================================
	
	//VanLeerLimiterY();
	MinmodLimiterY();
	
	
	//============================================================
	//	HLLD flux in y-direction
	//============================================================
	
	HLLDY();
	//HLLY();
	
	
	//============================================================
	//	2nd order slope limiter in z
	//============================================================
	
	//VanLeerLimiterZ();
	MinmodLimiterZ();
	
	
	//============================================================
	//	HLL flux in z-direction
	//============================================================
	
	HLLDZ();
	//HLLZ();
	
	
	//============================================================
	//	calculating source term
	//============================================================
	
	SpiralGravityCooling(V, S);
	//AxiGravity(V, S);
	//NoSource(V, S);
	//UniformGz(V, S);
	//Gravity(V, S);
	
	
	//-----U^(n+1/3)=U^n-dt/dx(F[i]-F[i-1])
	//									-dt/dy(G[j]-G[j-1])
	//									-dt/dz(H[k]-H[k-1])
	//									+dt*S[i][j][k]			-----
	
#pragma omp parallel for private(i, j, k, m)
	for(m=0; m<9; m++){
    
    for(i=1; i<ixmax1; i++){
      for(j=1; j<jymax1; j++){
        for(k=1; k<kzmax1; k++){
					
					U1[m][i][j][k]=U[m][i][j][k]
													-dt*((g2[i][j][k]*F[m][i][j][k]
																-g1[i][j][k]*F[m][i-1][j][k])
															 +(g4[i][j][k]*G[m][i][j][k]
																 -g3[i][j][k]*G[m][i][j-1][k])
															 +(g6[i][j][k]*H[m][i][j][k]
																 -g5[i][j][k]*H[m][i][j][k-1]))
													+dt*S[m][i][j][k];
					
					
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
	
	//ZLeftVzAntisymBSymBoundary(U1);
	//ZLeftFreeBoundary(U1);
	ZLeftPeriodicBoundary(U1);
	
	
	//============================================================
	//	k=kzmax1 free boundary
	//============================================================
	
	//ZRightFreeBoundary(U1);
	ZRightPeriodicBoundary(U1);
	
	
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











