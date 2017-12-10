

/***********************************************************************
 *
 *	commontemp.h
 *
 *	definition of temporary quantities.
 *
 *
 *	I, J, K : negative density or pressure cell point
 *	nstep : number of step
 *	nfile : number of output file
 *	prflag, rhoflag betaflag: pressure, density, beta negative point counter
 *	
 *
 *	Uinitial : conserved variable initial condition
 *
 *	Uinitial[0] : rho 
 *	Uinitial[1], Uinitial[2], Uinitial[3] : rvx, rvy, rvz
 *	Uinitial[4], Uinitial[5], Uinitial[6] : bx, by, bz
 *	Uinitial[7] : energy
 *	Uinitial[8] : psi
 *
 *
 *	Vinitial : primitive variable initial condition
 *
 *	Vinitial[0] : rho 
 *	Vinitial[1], Vinitial[2], Vinitial[3] : vx, vy, vz
 *	Vinitial[4], Uinitial[5], Uinitial[6] : bx, by, bz
 *	Vinitial[7] : energy
 *	Vinitial[8] : psi
 *
 *
 *	rhofloor : minimum rho criteria of each grid
 *	prfloor : minimum pr criteria of each grid
 *
 *	
 *	dt : cfl/cfldt, time step interval
 *	t : simulation time
 *
 *
 *	for MUSCL scheme (interpolate basic quantities)
 *
 *	Vl : left side primitive quantities (between i-1 & i)
 *	Vl[0] : left side rho
 *	Vl[1], Vl[2], Vl[3] : left side vx, vy, vz
 *	Vl[4], Vl[5], Vl[6] : left side bx, by, bz
 *	Vl[7] : left side pr
 *	Vl[8] : left side psi
 *
 *	Vr : right side primitive quantities (between i & i+1)
 *	Vr[0] : right side rho
 *	Vr[1], Vr[2], Vr[3] : right side vx, vy, vz
 *	Vr[4], Vr[5], Vr[6] : right side bx, by, bz
 *	Vr[7] : right side pr
 *	Vr[8] : right side psi
 *
 *
 *	convert from l, r primitive variables to l, r conserved variables
 *
 *	Ul[0] : left side rho
 *	Ul[1], Ul[2], Ul[3] : left side rvx, rvy, rvz
 *	Ul[4], Ul[5], Ul[6] : left side bx, by, bz
 *	Ul[7] : left side en
 *	Ul[8] : left side psi
 *
 *	Ur[0] : right side rho
 *	Ur[1], Ur[2], Ur[3] : right side rvx, rvy, rvz
 *	Ur[4], Ur[5], Ur[6] : right side bx, by, bz
 *	Ur[7] : right side en
 *	Ur[8] : right side psi
 *	
 *
 *	convert from l, r primitive variables to l, r flux in x
 *
 *	Fl[0]=rhol*vxl : left side mass flux 
 *	Fl[1]=rvxl*vxl+ptotl-bxl*bxl : left side x-momentum flux
 *	Fl[2]=rvyl*vxl-byl*bxl : left side y-momentum flux
 *	Fl[3]=rvzl*vxl-bzl*bxl : left side z-momentum flux
 *	Fl[4]=psil : left side bx flux
 *	Fl[5]=byl*vxl-vyl*bxl : left side by flux
 *	Fl[6]=bzl*vxl-vzl*bxl : left side bz flux
 *	Fl[7]=(enl+ptotl)*vxl-vldotbl*bxl : left side en flux
 *	Fl[8]=ch**2*bxl : left side psi flux
 *
 *	Fr[0]=rhor*vxr : right side mass flux 
 *	Fr[1]=rvxr*vxr+ptotr-bxr*bxr : right side x-momentum flux
 *	Fr[2]=rvyr*vxr-byr*bxr : right side y-momentum flux
 *	Fr[3]=rvzr*vxr-bzr*bxr : right side z-momentum flux
 *	Fr[4]=psir : right side bx flux
 *	Fr[5]=byr*vxr-vyr*bxr : right side by flux
 *	Fr[6]=bzr*vxr-vzr*bxr : right side bz flux
 *	Fr[7]=(enr+ptotr)*vxr-vrdotbr*bxr : right side en flux
 *	Fr[8]=ch**2*bxr : righ side psi flux
 *
 *
 *	convert from l, r primitive variables to l, r flux in y
 *
 *	Gl[0]=rhol*vyl : left side mass flux 
 *	Gl[1]=rvxl*vyl-bxl*byl : left side x-momentum flux
 *	Gl[2]=rvyl*vyl+(pgasl+pmagl)-byl*byl : left side y-momentum flux
 *	Gl[3]=rvzl*vyl-bzl*byl : left side z-momentum flux
 *	Gl[4]=bxl*vyl-vxl*byl : left side bx flux
 *	Gl[5]=psil : left side by flux
 *	Gl[6]=bzl*vyl-vzl*byl : left side bz flux
 *	Gl[7]=(enl+(pgasl+pmagl))*vyl-vldotbl*byl : left side en flux
 *	Gl[8]=ch**2*byl : left side psi flux
 *
 *	Gr[0]=rhor*vyr : right side mass flux 
 *	Gr[1]=rvxr*vyr-bxr*byr : right side x-momentum flux
 *	Gr[2]=rvyr*vyr+(pgasr+pmagr)-byr*byr : right side y-momentum flux
 *	Gr[3]=rvzr*vyr-bzr*byr : right side z-momentum flux
 *	Gr[4]=bxr*vyr-vxr*byr : right side bx flux
 *	Gr[5]=psir : right side by flux
 *	Gr[6]=bzr*vyr-vzr*byr : right side bz flux
 *	Gr[7]=(enr+(pgasr+pmagr))*vyr-vrdotbr*byr : right side en flux
 *	Gr[8]=ch**2*byr : right side psi flux
 *
 *
 *	convert from l, r primitive variables to l, r flux in z
 *
 *	Hl[0]=rhol*vzl : left side mass flux 
 *	Hl[1]=rvxl*vzl-bxl*bzl : left side x-momentum flux
 *	Hl[2]=rvyl*vzl-byl*bzl : left side y-momentum flux
 *	Hl[3]=rvzl*vzl+(pgasl+pmagl)-bzl*bzl : left side z-momentum flux
 *	Hl[4]=bxl*vzl-vxl*bzl : left side bx flux
 *	Hl[5]=byl*vzl-vyl*bzl : left side by flux
 *	Hl[6]=psil : left side bz flux
 *	Hl[7]=(enl+(pgasl+pmagl))*vzl-vldotbl*bzl : left side en flux
 *	Hl[8]=ch**2*bzl : left side psi flux
 *
 *	Hr[0]=rhor*vzr : right side mass flux 
 *	Hr[1]=rvxr*vzr-bxr*bzr : right side x-momentum flux
 *	Hr[2]=rvyr*vzr-byr*bzr : right side y-momentum flux
 *	Hr[3]=rvzr*vzr+(pgasr+pmagr)-bzr*bzr : right side z-momentum flux
 *	Hr[4]=bxr*vzr-vxr*bzr : right side bx flux
 *	Hr[5]=byr*vzr-vyr*bzr : right side by flux
 *	Hr[6]=psir : right side bz flux
 *	Hr[7]=(enr+(pgasr+pmagr))*vzr-vrdotbr*bzr : right side en flux
 *	Hr[8]=ch**2*bzr : right side psi flux
 *
 *  oenobeta=1.0/plasmabeta=0.5*B**2/V[7]
 *
 *  x, z, c absorbingfactor=coefficients of absorbing boundary regions
 *
 *	using GLM-MHD source term solver
 *
 *	psicoef : psi diffusion coefficient
 *
 *
 *	2012 Dec. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int I, J, K;
int nstep, nout;
int prflag, rhoflag, oneobetaflag;

double Uinitial[9][ixmax][jymax][kzmax], 
				Vinitial[9][ixmax][jymax][kzmax];

double rhofloor[ixmax][jymax][kzmax], prfloor;
double dt, t;

double Vl[9][ixmax1][jymax1][kzmax1], Vr[9][ixmax1][jymax1][kzmax1];
double Ul[9][ixmax1][jymax1][kzmax1], Ur[9][ixmax1][jymax1][kzmax1];

double Fl[9][ixmax1][jymax1][kzmax1], Fr[9][ixmax1][jymax1][kzmax1];
double Gl[9][ixmax1][jymax1][kzmax1], Gr[9][ixmax1][jymax1][kzmax1];
double Hl[9][ixmax1][jymax1][kzmax1], Hr[9][ixmax1][jymax1][kzmax1];

double oneobeta[ixmax][jymax][kzmax]={0.0};

double x_absorb_factor[11]={0.0};
double z_absorb_factor[11]={0.0};
double c_absorb_factor[irin+1][kzmax]={0.0};

double psicoef;


char filename[64];
char dirname[64];
FILE *fi, *fo;








