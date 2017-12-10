

/***********************************************************************
 *
 *	vlrulrtohlr.h 
 *
 *	convert from Vl & Ul to Hl, from Vr & Ur to Hr.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	b2=bx*bx+by*by+bz*bz : b**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	dummyVlr : left or right side V
 *	dummyUlr : left or right side U
 *	dummyGlr : left or right side H
 *
 *	
 *	2012 Dec. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int VlrUlrtoHlr(int i, int j, int k, double b2, double vdotb, 
						 double dummyVlr[][ixmax1][jymax1][kzmax1], 
						 double dummyUlr[][ixmax1][jymax1][kzmax1], 
						 double dummyHlr[][ixmax1][jymax1][kzmax1])
{
	
	dummyHlr[0][i][j][k]=dummyUlr[3][i][j][k];
	dummyHlr[1][i][j][k]=dummyUlr[1][i][j][k]*dummyVlr[3][i][j][k]
												-dummyVlr[4][i][j][k]*dummyVlr[6][i][j][k];
	dummyHlr[2][i][j][k]=dummyUlr[2][i][j][k]*dummyVlr[3][i][j][k]
												-x[i]
													*dummyVlr[5][i][j][k]*dummyVlr[6][i][j][k];
	dummyHlr[3][i][j][k]=dummyUlr[3][i][j][k]*dummyVlr[3][i][j][k]
												+(dummyVlr[7][i][j][k]+0.5*b2)
												-dummyVlr[6][i][j][k]*dummyVlr[6][i][j][k];
	dummyHlr[4][i][j][k]=dummyVlr[4][i][j][k]*dummyVlr[3][i][j][k]
												-dummyVlr[1][i][j][k]*dummyVlr[6][i][j][k];
	dummyHlr[5][i][j][k]=dummyVlr[5][i][j][k]*dummyVlr[3][i][j][k]
												-dummyVlr[2][i][j][k]*dummyVlr[6][i][j][k];
	dummyHlr[6][i][j][k]=0.0;
	dummyHlr[7][i][j][k]=(dummyUlr[7][i][j][k]
												+dummyVlr[7][i][j][k]+0.5*b2)
												*dummyVlr[3][i][j][k]
												-vdotb*dummyVlr[6][i][j][k];
	
	
	return 0;
	
}

