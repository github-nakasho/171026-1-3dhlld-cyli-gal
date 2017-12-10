

/***********************************************************************
 *
 *	input.h
 *
 *	input parameters.
 *
 *
 *	ix, jy, kz : cell number of x, y, z in simulation region
 *	lx, ly, lz : no stretch region.
 *	x0, y0, z0 : r, phi, z origin
 *	rin : central absorbing radius
 *	margin : boundary buffer = 2(2nd order), 4(3rd order)
 *	ixmax, jymax, kzmax : the number of array in x, y, z
 *	ixmax1, jymax1, kzmax1 : ixmax-1, jymax-1, kzmax-1
 *	cfl : CFL condition
 *	dtmin : minimum dt criteria
 *
 *	dx0, dy0, dz0 : cell center interval in x, y, z
 *	dxmax, dzmax : max length of dx, dz
 *
 *	rhoprcr : if rho or pr is negative, set rhoprcr*initial
 *  betacr : if beta is smaller than betacr, inject rho.
 *	nstop : simulation stop n
 *	tint : result output time interval 
 *	tstop : simulation stop time
 *
 *  omegasp: anguler velocity of rotating frame
 *  omegapattern: anguler velocity of pattern (spiral arms, bar)
 *
 *
 *	2015 July 23 : add omegapattern
 *	2015 Feb. 20 : add betacr
 *	2014 July 20 : separated rhoprcr ---> rhocr, prcr
 *	2013 Aug. 22 : add spiral arm scale redius(=scaler) & phasescale
 *	2013 May 10 : add spiral potential parameters.
 *	2013 Apr. 15 : for stretch grid.
 *	2012 Nov. 21 : add margin & generalized (for higher order).
 *	2012 Oct. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


#define ix 240
#define jy 64
#define kz 200
#define lx 2.0
#define ly pi2
#define lz 0.03
#define x0 0.04
#define y0 0.0
#define z0 0.0
#define stretchratio 1.05
#define rin 0.1

#define margin 2
#define cfl 0.4
#define dtmin 1.0e-9

#define ixmax (ix+margin)
#define jymax (jy+margin)
#define kzmax (kz+margin)
#define ixmax1 (ixmax-1)
#define jymax1 (jymax-1)
#define kzmax1 (kzmax-1)

#define  dx0 0.01
#define dy0 (ly/(double)jy)
#define dz0 0.0003

#define irin ((int)(rin/dx0))


double dxmax=50.0*dx0;
double dzmax=50.0*dz0;

#define rhocr 1.0e-4
#define prcr 0.015
#define betacr 0.1
#define oneobetacr (1.0/betacr)
#define nstop 10000000
#define tint 0.05
#define tstop 40.0

double tout=tint;

double epssp0=0.05;
double tausp=200.0*1.0e+6*year/t0;
double scalez=0.03;
double scaler=0.7;
double scalephase=0.1;
int armnum=2;
double omegasp=10.0e+5/omega0*oneokpc;
double omegapattern=10.0e+5/omega0*oneokpc;
double pitch=15.0*M_PI/180.0;


//-----input parameter for restart parameters-----

int restartsw=0;
int nout=0;





