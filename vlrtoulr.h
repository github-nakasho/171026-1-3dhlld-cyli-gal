

/***********************************************************************
 *
 *	vlrtoulr.h 
 *
 *	convert from Vl to Ul, from Vr to Ur.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	r=x[i] or xm[i] : distance from polar to cell enter, surface
 *	dummyVlr : left or right side V
 *	dummyUlr : left or right side U
 *
 *	
 *	2013 Apr. 13 : add r(=x[i] or xm[i]).
 *	2012 Nov. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int VlrtoUlr(int i, int j, int k, double v2, double b2, double r,
						 double dummyVlr[][ixmax1][jymax1][kzmax1], 
						 double dummyUlr[][ixmax1][jymax1][kzmax1])
{
	
	dummyUlr[0][i][j][k]=dummyVlr[0][i][j][k];
	dummyUlr[1][i][j][k]=dummyVlr[0][i][j][k]*dummyVlr[1][i][j][k];
	dummyUlr[2][i][j][k]=r*dummyVlr[0][i][j][k]*dummyVlr[2][i][j][k];
	dummyUlr[3][i][j][k]=dummyVlr[0][i][j][k]*dummyVlr[3][i][j][k];
	dummyUlr[4][i][j][k]=dummyVlr[4][i][j][k];
	dummyUlr[5][i][j][k]=dummyVlr[5][i][j][k];
	dummyUlr[6][i][j][k]=dummyVlr[6][i][j][k];
	dummyUlr[7][i][j][k]=0.5*(dummyVlr[0][i][j][k]*v2+b2)
												+oneogm1*dummyVlr[7][i][j][k];
	
	
	return 0;
	
}

