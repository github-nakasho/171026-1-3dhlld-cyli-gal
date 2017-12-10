

/***********************************************************************
 *
 *	primitivetof.h 
 *
 *	convert from primitive variables to flux in x
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	
 *	
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int PrimitivetoF(int i, int j, int k, double vdotb, 
								 double rho, double vx, double vy, double vz, 
								 double bx, double by, double bz, 
								 double en, double ptot, 
								 double dummyF[][ixmax1][jymax1][kzmax1])
{
	
	dummyF[0][i][j][k]=rho*vx;
	dummyF[1][i][j][k]=rho*vx*vx+ptot-bx*bx;
	dummyF[2][i][j][k]=xm[i]*(rho*vy*vx-by*bx);
	dummyF[3][i][j][k]=rho*vz*vx-bz*bx;
	dummyF[4][i][j][k]=0.0;
	dummyF[5][i][j][k]=by*vx-vy*bx;
	dummyF[6][i][j][k]=bz*vx-vz*bx;
	dummyF[7][i][j][k]=(en+ptot)*vx-vdotb*bx;
	
	
	return 0;
	
}

