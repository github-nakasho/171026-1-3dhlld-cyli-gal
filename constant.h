

/***********************************************************************
 *
 *	constant.h
 *
 *	definition of physical constant & normalized values.
 *
 *
 *	light : speed of light
 *	kpc : kilo-persec
 *	kb : Boltzmann constant
 *	mu : mean molecular weight
 *	mp : proton mass
 *	g : Gravitational constant
 *	msun : solar mass
 *	year : 365x24x60x60
 *	gm : gamma (specific heat ratio)
 *	gm1 : gm-1
 *	ogm : inverse of gm
 *	ogm1 : inverse of gm1
 *	gmogm1 : gm/gm1
 *	pi : circular constant
 *	pi2, pi4, pi8 : 2xpi, 4xpi, 8xpi
 *	opi2, opi4, opi8 : 1/pi2, 1/pi4, 1/pi8
 *
 *
 *	using GLM-MHD divergence cleaning (Dedner et al. 2002)
 *
 *
 *	cr=cp**2/ch=0.18
 *	onecr=1.0/cr
 *
 *
 *	rho0, l0, t0, b02, v02, te0, lamb0, omega0
 *				: density, length, time, magnetic field, 
 *					velocity**2, temperature, cooling function, 
 *					spiral patern angular velocity normalized value
 *
 *
 *	2012 Nov. 06 add variables using GLM-MHD divergence cleaning. 
 *	2012 Oct. 05 coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


#define light 2.998e+10
#define kpc 3.08568e+21
#define kb 1.3806504e-16
#define year 3.1536e+7
#define mu 0.619047619
#define mp 1.673e-24
#define g 6.67428e-8
#define msun 1.989e+33
#define oneolight (1.0/light)
#define oneokpc (1.0/kpc)
#define oneokb (1.0/kb)
#define oneoyear (1.0/year)
#define mump (mu*mp)
#define oneomump (1.0/mump)

double gm=5.0/3.0;
double gm1=2.0/3.0;
double oneogm=0.6;
double oneogm1=1.5;
double gmogm1=2.5;

double pi=M_PI;
double pi2=2.0*M_PI;
double pi4=4.0*M_PI;
double pi8=8.0*M_PI;
double oneopi2=1.0/(2.0*M_PI);
double oneopi4=1.0/(4.0*M_PI);
double oneopi8=1.0/(8.0*M_PI);

double oneo3=1.0/3.0;

double cr=0.18;
double oneocr=1.0/0.18;


#define rho0 mp
#define l0 (10*kpc)
#define v0 6.558e+6
#define t0 4.7058e+15
#define b0 8.474e-6
#define te0 (mump*oneokb*v0*v0)
#define lamb0 (rho0*v0*v0*v0/l0*(mump/rho0)*(mump/rho0))
#define omega0 (1.0/t0)




