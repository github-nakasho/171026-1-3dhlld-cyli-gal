

/***********************************************************************
 *
 *	primitivetoh.h 
 *
 *	convert from primitive variables to flux in z
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	
 *	
 *	2013 Jan. 15 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int PrimitivetoH(int i, int j, int k, double vdotb, 
								 double rho, double vx, double vy, double vz, 
								 double bx, double by, double bz, 
								 double en, double ptot, 
								 double dummyH[][ixmax1][jymax1][kzmax1])
{
	
	dummyH[0][i][j][k]=rho*vz;
	dummyH[1][i][j][k]=rho*vx*vz-bx*bz;
	dummyH[2][i][j][k]=x[i]*(rho*vy*vz-by*bz);
	dummyH[3][i][j][k]=rho*vz*vz+ptot-bz*bz;
	dummyH[4][i][j][k]=bx*vz-vx*bz;
	dummyH[5][i][j][k]=by*vz-vy*bz;
	dummyH[6][i][j][k]=0.0;
	dummyH[7][i][j][k]=(en+ptot)*vz-vdotb*bz;
	
	
	return 0;
	
}

