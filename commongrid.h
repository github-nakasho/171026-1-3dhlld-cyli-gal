

/***********************************************************************
 *
 *	commongrid.h
 *	
 *	definition of cell position, cell center interval.
 *
 *
 *	x, y, z : position of cell center
 *	xm[i] : (x[i]+x[i+1])/2 = median of x[i] & x[i+1]
 *	zm[k] : (z[k]+z[k+1])/2 = median of z[k] & z[k+1]
 *	oneox[i] : 1.0/x[i]
 *	oneoxm[i] : 1.0/xm[i];
 *  oneodx[i] : 1.0/dx[i];
 *  oneody : 1.0/dy;
 *	dx : cell center interval of r direction
 *	dy : cell center interval of phi direction
 *	dz : cell center interval of z direction
 *  oneodz[k] : 1.0/dz[k]
 *	dxm[i] : (dx[i]+dx[i+1])/2 = median of dx[i] & dx[i+1]
 *	dzm[k] : (dz[k]+dz[k+1])/2 = median of dz[k] & dz[k+1]
 *	
 *
 *	g1-g6 : S1/dvolume, S2/dvolume, ...S6/dvolume = geometrical factors
 *	
 *	dvolume[i][j][k] : volume of [i][j][k] cell
 *	
 *	minlength : minimum side length of [i][j][k] cell
 *	minminlength : minimum side length of all cell
 *	oneominminlength=1.0/minminlength
 *
 *
 *
 *  2015 July 22: add krin.
 *	2013 Apr. 05: add dx dz array, geometrical factors.
 *	2012 Oct. 03: coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


double x[ixmax], y[jymax], z[kzmax];
double xm[ixmax], zm[kzmax];
double oneox[ixmax], oneoxm[ixmax1], oneodx[ixmax];
double oneody, oneodz[kzmax];
double dx[ixmax], dy, dz[kzmax];
double dxm[ixmax], dzm[kzmax];


double g1[ixmax1][jymax1][kzmax1]={0.0}, 
				g2[ixmax1][jymax1][kzmax1]={0.0}, 
				g3[ixmax1][jymax1][kzmax1]={0.0}, 
				g4[ixmax1][jymax1][kzmax1]={0.0}, 
				g5[ixmax1][jymax1][kzmax1]={0.0}, 
				g6[ixmax1][jymax1][kzmax1]={0.0};

double dvolume[ixmax][jymax][kzmax];

double minlength[ixmax][jymax][kzmax];
double minminlength, oneominminlength;

int krin[irin+1]={0};





