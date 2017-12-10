

/***********************************************************************
 *
 *	primitivetof.h 
 *
 *	convert from primitive variables to flux in y
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	
 *	
 *	2013 Jan. 11 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int PrimitivetoG(int i, int j, int k, double vdotb, 
								 double rho, double vx, double vy, double vz, 
								 double bx, double by, double bz, 
								 double en, double ptot, 
								 double dummyG[][ixmax1][jymax1][kzmax1])
{
	
	dummyG[0][i][j][k]=rho*vy;
	dummyG[1][i][j][k]=rho*vx*vy-bx*by;
	dummyG[2][i][j][k]=x[i]*(rho*vy*vy+ptot-by*by);
	dummyG[3][i][j][k]=rho*vz*vy-bz*by;
	dummyG[4][i][j][k]=bx*vy-vx*by;
	dummyG[5][i][j][k]=0.0;
	dummyG[6][i][j][k]=bz*vy-vz*by;
	dummyG[7][i][j][k]=(en+ptot)*vy-vdotb*by;
	
	
	return 0;
	
}

