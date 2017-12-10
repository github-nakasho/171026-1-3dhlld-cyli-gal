

/***********************************************************************
 *
 *	vlrulrtoflr.h 
 *
 *	convert from Vl & Ul to Fl, from Vr & Ur to Fr.
 *
 *
 *	i, j, k : loop index of r, phi, z
 *	b2=br*br+bphi*bphi+bz*bz : b**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	dummyVlr : left or right side V
 *	dummyUlr : left or right side U
 *	dummyFlr : left or right side F
 *
 *	
 *	2012 Nov. 05 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int VlrUlrtoFlr(int i, int j, int k, double b2, double vdotb, 
						 double dummyVlr[][ixmax1][jymax1][kzmax1], 
						 double dummyUlr[][ixmax1][jymax1][kzmax1], 
						 double dummyFlr[][ixmax1][jymax1][kzmax1])
{
	
	dummyFlr[0][i][j][k]=dummyUlr[1][i][j][k];
	dummyFlr[1][i][j][k]=dummyUlr[1][i][j][k]*dummyVlr[1][i][j][k]
												+(dummyVlr[7][i][j][k]+0.5*b2)
												-dummyVlr[4][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[2][i][j][k]=dummyUlr[2][i][j][k]*dummyVlr[1][i][j][k]
												-xm[i]*dummyVlr[5][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[3][i][j][k]=dummyUlr[3][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[6][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[4][i][j][k]=0.0;
	dummyFlr[5][i][j][k]=dummyVlr[5][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[2][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[6][i][j][k]=dummyVlr[6][i][j][k]*dummyVlr[1][i][j][k]
												-dummyVlr[3][i][j][k]*dummyVlr[4][i][j][k];
	dummyFlr[7][i][j][k]=(dummyUlr[7][i][j][k]
												+dummyVlr[7][i][j][k]+0.5*b2)
												*dummyVlr[1][i][j][k]
												-vdotb*dummyVlr[4][i][j][k];
	
	
	return 0;
	
}

