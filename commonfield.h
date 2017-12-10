

/***********************************************************************
 *
 *	commonfield.h
 *
 *	definition of field quantities.
 *
 *	
 *	U : conserved variables
 *
 *	U[0]=rho : density
 *	U[1]=rho*vr : momenta in r (=rvr) 
 *	U[2]=r*rho*vphi : angular momenta in z (=amz)
 *	U[3]=r*vz : momenta in z (=rvz)
 *	U[4]=br : magnetic field in r
 *	U[5]=bphi : magnetic field in phi
 *	U[6]=bz : magnetic field in z
 *	U[7]=en : total energy per unit volume
 *	U[8]=psi : GLM-MHD scalor
 *
 *
 *	U1 : conserved variables 
 *						for TVD-RungeKutta half/one-third/two-third step
 *	(notation is the same as U)
 *
 *
 *	V : primitive variables.
 *
 *	V[0]=rho : density
 *	V[1]=vr : velocity in r
 *	V[2]=vphi : velocity in phi
 *	V[3]=vz : velocity in z
 *	V[4]=br : magnetic field in r
 *	V[5]=bphi : magnetic field in phi
 *	V[6]=bz : magnetic field in z
 *	V[7]=pgas : gas pressure
 *	V[8]=psi : GLM-MHD scalor
 *
 *
 *	F : flux in r
 *
 *	F[0]=rvr : mass flux in r
 *	F[1]=rvr*vr+(pgas+pmag)-br*br : r-momenta flux in r
 *	F[2]=amz*vr-r*bphi*br : z-angular momenta flux in r
 *	F[3]=rvz*vr-bz*br : z-momenta flux in r
 *	F[4]=psi : br flux in r
 *	F[5]=bphi*vr-vphi*br : bphi flux in r
 *	F[6]=bz*vr-vz*br : bz flux in r
 *	F[7]=(en+(pgas+pmag))*vr-vdotb*br : energy flux in r
 *	F[8]=ch**2*br : psi flux in r
 *	
 *	
 *	G : flux in phi
 *
 *	G[0]=rho*vphi : mass flux in phi
 *	G[1]=rvr*vphi-bphi*br : r-momenta flux in phi
 *	G[2]=amz*vphi+r*(pgas+pmag)
 *								-r*bphi*bphi : z-angular momenta flux in phi
 *	G[3]=rvz*vphi-bz*bphi : z-momenta flux in phi
 *	G[4]=br*vphi-vr*bphi : br flux in phi
 *	G[5]=psi : bphi flux in phi
 *	G[6]=bz*vphi-vz*bphi : bz flux in phi
 *	G[7]=(en+(pgas+pmag))*vphi-vdotb*bphi : energy flux in phi
 *	G[8]=ch**2*bphi : psi flux in phi
 *
 *	
 *	H : flux in z
 *
 *	H[0]=rvz : mass flux in z
 *	H[1]=rvr*vz-br*bz : r-momenta flux in z
 *	H[2]=amz*vz-r*bphi*bz : z-angular momenta flux in z
 *	H[3]=rvz*vz+(pgas+pmag)-bz*bz : z-momenta flux in z
 *	H[4]=br*vz-vr*bz : br flux in z
 *	H[5]=bphi*vz-vphi*bz : bphi flux in z
 *	H[6]=psi : bz flux in z
 *	H[7]=(en+(pgas+pmag))*vz-vdotb*bz : energy flux in z
 *	H[8]=ch**2*bz : psi flux in z
 *
 *
 *	S : source term
 *
 *	S[0] : mass source
 *	S[1] : r-momenta source
 *	S[2] : z-angular momenta source
 *	S[3] : z-momenta source
 *	S[4] : br source
 *	S[5] : bphi source
 *	S[6] : bz source
 *	S[7] : energy source
 *
 *
 *	Grx : gravity field in r
 *	Gry : gravity field in phi
 *	Grz : gravity field in z
 *	
 *	
 *	Pot : gravitational potential
 *
 *
 *	te : gas temperature
 *	divB : divergence B
 *	
 *	dt : time step
 *	ch : max veocity (using GLM divergence cleaning)
 *
 *
 *	2013 Apr. 05 : notation change to Cylindrical coord.
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


double U[9][ixmax][jymax][kzmax];

double U1[9][ixmax][jymax][kzmax];

double V[9][ixmax][jymax][kzmax];

double F[9][ixmax1][jymax1][kzmax1], 
				G[9][ixmax1][jymax1][kzmax1], 
				H[9][ixmax1][jymax1][kzmax1];

double S[9][ixmax][jymax][kzmax];

double Grx[ixmax][jymax][kzmax],
				Gry[ixmax][jymax][kzmax], 
				Grz[ixmax][jymax][kzmax];

double Grxaxi[ixmax][jymax][kzmax],
				Gryaxi[ixmax][jymax][kzmax],
				Grzaxi[ixmax][jymax][kzmax];

double Pot[ixmax][jymax][kzmax];


double te[ixmax][jymax][kzmax], divB[ix][jy][kz];

double dt;

double ch;




