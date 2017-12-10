


/***********************************************************************
 *
 *	model.h 
 *
 *	set initial model.
 *
 *	
 *	2015 July 11 : add Corotationframe.
 *	2014 Dec. 16 : add AddMagneticFields.
 *	2014 Dec. 02 : add GalacticDisc3.
 *	2014 Apr. 23 : add GalacticDisc2.
 *	2014 Apr. 23 : add DoubleExponentialDisc.
 *	2013 Apr. 13 : add TorusHalo.
 *	2013 Apr. 08 : add SphericalBlastWaveHDCylinder.
 *	2013 Apr. 05 : add MagneticFieldLoopAdvectionCylinder.
 *	2013 Jan. 16 : add Torus.
 *	2013 Jan. 01 : add RTI.
 *	2012 Dec. 24 : add JetMHD.
 *	2012 Dec. 13 : add DW1aShockTube z-directiong (Dai & Woodard 1995)
 *	2012 Dec. 06 : add JetHD.
 *	2012 Nov. 06 : use pointer & reduce useless funcitons.
 *	2012 Nov. 02 : add SphericalBlastWaveHD, MagneticFieldLoopAdvection
 *											SphericalBlastWaveMHD
 *	2012 Oct. 31 : add Orszag-Tang vortex problem (Dai & Woodard 1995)
 *	2012 Oct. 30 : add DW1aShockTube y-direction (Dai & Woodard 1995)
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//============================================================
//	add magnetic fields to equilibrium disc
//============================================================

int Corotationframe(void)
{
  
  int i, j, k;
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        V[2][i][j][k]-=x[i]*omegasp;
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	k=0 impose symmetric boundary
  //==================================================
  
  ZLeftVzAntisymBSymBoundary(U);
  ZLeftFreeBoundary(U);
  
  
  //==================================================
  //	save initial condition.
  //==================================================
  
  SaveInitial();
  
  
  return 0;
  
}



//============================================================
//	add magnetic fields to equilibrium disc
//============================================================

int AddMagneticFields(void)
{

  int i, j, k;
  double discbeta=1.0e+4;
    double r;
  
  
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
          r=sqrt(x[i]*x[i]+z[k]*z[k]);
          
        if(V[0][i][j][k]>0.05 &&
           V[7][i][j][k]/V[0][i][j][k]<0.054 &&
           r>rin){
        //if(V[7][i][j][k]/V[0][i][j][k]<0.755 && r>rin ){
          V[5][i][j][k]=sqrt(2.0*V[7][i][j][k]/discbeta);
          
        
        }
        else{ V[5][i][j][k]=0.0; }
        
      }
    }
  }

  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	k=0 impose symmetric boundary
  //==================================================
  
  ZLeftVzAntisymBSymBoundary(U);
  ZLeftFreeBoundary(U);

  
  //==================================================
  //	save initial condition.
  //==================================================
  
  SaveInitial();
  
  
  return 0;

}


//============================================================
//	realistic disc + halo model 3
//============================================================


int GalacticDisc3(void)
{
  
  int i, j, k;
  
  
  //============================================================
  // a fixed potential parameters
  //
  // vhalo0, dh: the dark logarithmic halo potential parameter
  // a, b, md : parameters of the Miyamoto-Nagai disc
  // db, mb : the spherical Hernquist stellar bulge's parameters
  //
  //============================================================
  
  double vhalo0=131.5e+5/v0, dh=1.2;
  double a=0.65, b=0.026, md=10.0;
  double db=0.07, mb=3.4;
  
  
  //============================================================
  //
  // isothermal gas disc parameters
  //
  // mdisc : parameter of total gas mass
  // rdisc, zdisc : exp(-R/rdisc)sech^2(z/zdisc)
  // discrho0 : gas density @ R=0
  // csdisc2 : speed of sound in disc
  // discbeta : plasma beta @ R=0
  // discbphi0 : toroidal magnetic field strength @ R=0
  //
  //============================================================
  
  double mdisc=0.4*0.75;
  double rdisc=0.7, zdisc=0.03;
  double discrho0=(1.0e+10*msun/(mp*1000*kpc*kpc*kpc)
                   *mdisc/(pi4*rdisc*rdisc*zdisc));
  double csdisc2=1.33e+8*1.0e+4/(v0*v0);
  double discbeta=1.0e+4;
  double discbphi0;
  
  
  //============================================================
  //
  // hydrostatic isothermal halo parameters
  //
  // halorho0 : gas density @ R=0 (input relative value of disc)
  // cshalo2 : speed of sound in halo
  //
  //============================================================
  
  double halorho0=1.0e-4*(1.0e+10*msun/(mp*1000*kpc*kpc*kpc)
                          *mdisc/(pi4*rdisc*rdisc*zdisc));
  double cshalo2=1.33e+8*2.4e+6/(v0*v0);
  
  double x2, z2;
  double potaxi, potaxi0, potaxir;
  double halorho, halopr, discrho, discpr;
  double oneocshalo2=1.0/cshalo2;
  double oneocsdisc2=1.0/csdisc2;
  
  
  discbphi0=sqrt(2.0*csdisc2*discrho0/discbeta);
  //discbphi0=0.0;
  
  
#pragma omp parallel for private(i, j, k, x2, z2, potaxi, potaxi0, potaxir, halorho, halopr, discrho, discpr) firstprivate(vhalo0, dh, a, b, md, db, mb, mdisc, rdisc, zdisc, discrho0, csdisc2, oneocsdisc2, discbeta, halorho0, cshalo2, oneocshalo2, discbphi0)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        x2=x[i]*x[i];
        z2=z[k]*z[k];
        
        
        potaxi=((vhalo0*vhalo0*
                 (2.0*log(l0)+log(x2+z2+dh*dh)))
                -md/sqrt(x2+(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)))
                -mb/(sqrt(x2+z2)+db));
        
        potaxir=((vhalo0*vhalo0*
                 (2.0*log(l0)+log(x2+dh*dh)))
                -md/sqrt(x2+(a+b)*(a+b))
                -mb/(x[i]+db));
        
        potaxi0=((vhalo0*vhalo0*
                  (2.0*log(l0)+log(dh*dh)))
                 -md/(a+b)
                 -mb/db);
        
        
        Grxaxi[i][j][k]=((-2.0*vhalo0*vhalo0/(x2+z2+dh*dh)
                          -md/pow(x2
                                  +(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)), 1.5)
                          -mb/((sqrt(x2+z2)+db)*(sqrt(x2+z2)+db)
                               *sqrt(x2+z2)))
                         *x[i]);
        
        Gryaxi[i][j][k]=0.0;
        
        Grzaxi[i][j][k]=((-2.0*vhalo0*vhalo0/(x2+z2+dh*dh)
                          -md/pow(x2
                                  +(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)), 1.5)
                          *(1.0+a/sqrt(z2+b*b))
                          -mb/((sqrt(x2+z2)+db)*(sqrt(x2+z2)+db)
                               *sqrt(x2+z2)))
                         *z[k]);
        
        
        halorho=halorho0*exp(-oneocshalo2*(potaxi-potaxi0));
        halopr=cshalo2*halorho;
        
        
        discrho=(discrho0
                 *exp(-x[i]/rdisc)
                 *exp(-oneocsdisc2*(potaxi-potaxir)));
        discpr=csdisc2*discrho;
        
        
        //==================================================
        //	disc+halo
        //==================================================
        
        V[0][i][j][k]=discrho+halorho;
        V[1][i][j][k]=0.0;
       
        /*
        V[2][i][j][k]=(sqrt(-x[i]*(Grxaxi[i][j][k]+csdisc2/rdisc))
                       *discrho/V[0][i][j][k]
                       -omegasp*x[i]);
        */
        V[2][i][j][k]=(sqrt(-x[i]*(Grxaxi[i][j][k]+csdisc2/rdisc))
                       *discrho/V[0][i][j][k]
                       -omegasp*x[i]);
        
        V[3][i][j][k]=0.0;
        V[4][i][j][k]=0.0;
        V[5][i][j][k]=(discbphi0
                       *exp(-x[i]/rdisc)
                       /(cosh(z[k]/zdisc)*cosh(z[k]/zdisc)));
        V[6][i][j][k]=0.0;
        V[7][i][j][k]=discpr+halopr;
        V[8][i][j][k]=0.0;
        
        
        
      }
    }
  }
  
  
  //==================================================
  //	convert from V to U.
  //==================================================
  
  VtoU();
  
  
  //==================================================
  //	k=0 impose symmetric boundary
  //==================================================
  
  ZLeftVzAntisymBSymBoundary(U);
  ZLeftFreeBoundary(U);
  //ZLeftPeriodicBoundary(U);
  //ZRightPeriodicBoundary(U);
  
  
  //==================================================
  //	save initial condition.
  //==================================================
  
  SaveInitial();
  
  
  return 0;
  
}



//============================================================
//	realistic disc + halo model 2
//============================================================


int GalacticDisc2(void)
{
	
  int i, j, k;
	
  
  //============================================================
  // a fixed potential parameters
  //
  // vhalo0, dh: the dark logarithmic halo potential parameter
  // a, b, md : parameters of the Miyamoto-Nagai disc
  // db, mb : the spherical Hernquist stellar bulge's parameters
  //
	//============================================================
  
  double vhalo0=131.5e+5/v0, dh=1.2;
  double a=0.65, b=0.026, md=10.0;
  double db=0.07, mb=3.4;
  
 
	//============================================================
  //
  // isothermal gas disc parameters
  //
  // mdisc : total mass of gas disc
  // rdisc, zdisc : exp(-R/rdisc)sech^2(z/zdisc)
  // discrho0 : gas density @ R=0
  // csdisc2 : speed of sound in disc
  // discbeta : plasma beta @ R=0
  // discbphi0 : toroidal magnetic field strength @ R=0
  //
	//============================================================
  
  double mdisc=0.4;
  double rdisc=0.7, zdisc=0.03;
  double discrho0=(1.0e+10*msun/(mp*1000*kpc*kpc*kpc)
									 *mdisc/(pi4*rdisc*rdisc*zdisc));
	double csdisc2=1.33e+8*1.0e+4/(v0*v0);
	double discbeta=1.0e+4;
  double discbphi0;
	
  
  //============================================================
  //
  // hydrostatic isothermal halo parameters
  //
  // halorho0 : gas density @ R=0 (input relative value of disc)
  // cshalo2 : speed of sound in halo
  //
	//============================================================
  
  double halorho0=1.0e-4*(1.0e+10*msun/(mp*1000*kpc*kpc*kpc)
                          *mdisc/(pi4*rdisc*rdisc*zdisc));
  double cshalo2=1.33e+8*2.4e+6/(v0*v0);
	
  double x2, z2;
  double potaxi, potaxi0;
  double halorho, halopr, discrho, discpr;
	double oneocshalo2=1.0/cshalo2;
	
  
  //discbphi0=sqrt(2.0*csdisc2*discrho0/discbeta);
  discbphi0=0.0;
  
  
#pragma omp parallel for private(i, j, k, x2, z2, potaxi, potaxi0, halorho, halopr, discrho, discpr) firstprivate(vhalo0, dh, a, b, md, db, mb, mdisc, rdisc, zdisc, discrho0, csdisc2, discbeta, halorho0, cshalo2, oneocshalo2, discbphi0)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
        x2=x[i]*x[i];
        z2=z[k]*z[k];
      
        
        potaxi=((vhalo0*vhalo0*
                 (2.0*log(l0)+log(x2+z2+dh*dh)))
                -md/sqrt(x2+(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)))
                -mb/(sqrt(x2+z2)+db));
        
        potaxi0=((vhalo0*vhalo0*
                  (2.0*log(l0)+log(dh*dh)))
                 -md/(a+b)
                 -mb/db);
				
				
        Grxaxi[i][j][k]=((-2.0*vhalo0*vhalo0/(x2+z2+dh*dh)
                         -md/pow(x2
                                 +(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)), 1.5)
                         -mb/((sqrt(x2+z2)+db)*(sqrt(x2+z2)+db)
                              *sqrt(x2+z2)))
                         *x[i]);
				
				Gryaxi[i][j][k]=0.0;
				
				Grzaxi[i][j][k]=((-2.0*vhalo0*vhalo0/(x2+z2+dh*dh)
                         -md/pow(x2
                                 +(a+sqrt(z2+b*b))*(a+sqrt(z2+b*b)), 1.5)
                          *(1.0+a/sqrt(z2+b*b))
                         -mb/((sqrt(x2+z2)+db)*(sqrt(x2+z2)+db)
                              *sqrt(x2+z2)))
                         *z[k]);
				
				
				halorho=halorho0*exp(-oneocshalo2*(potaxi-potaxi0));
				halopr=cshalo2*halorho;
				
				
				discrho=(discrho0
								 *exp(-x[i]/rdisc)
								 /(cosh(z[k]/zdisc)*cosh(z[k]/zdisc)));
				discpr=csdisc2*discrho;
				
        
				//==================================================
				//	disc+halo
				//==================================================
        
				V[0][i][j][k]=discrho+halorho;
				//V[0][i][j][k]=discrho;
        V[1][i][j][k]=0.0;
        
        V[2][i][j][k]=(sqrt(-x[i]*(Grxaxi[i][j][k]+csdisc2/rdisc)
                            *discrho/V[0][i][j][k])
                       -omegasp*x[i]);
        
        /*
				V[2][i][j][k]=(sqrt(-x[i]*(Grxaxi[i][j][k]+csdisc2/rdisc))
                       -omegasp*x[i]);
				*/
        
        V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=(discbphi0
											 *exp(-x[i]/rdisc)
                       /(cosh(z[k]/zdisc)*cosh(z[k]/zdisc)));
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=discpr+halopr;
				V[8][i][j][k]=0.0;
				
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	k=0 impose symmetric boundary
	//==================================================
	
	ZLeftVzAntisymBSymBoundary(U);
	ZLeftFreeBoundary(U);
	//ZLeftPeriodicBoundary(U);
	//ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial condition.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	realistic disc + halo model
//============================================================


int DoubleExponentialDisc(void)
{
	
	
	//=================================================================
	// rb : distance from polar to bottom of torus
	// 
	// Miyamoto et al.(1980) galactic potential including DM.
	// m : bulge, disk, halo mass
	// a : bulge, disk, halo distribution parameter in r
	// b : bulge, disk, halo distribution parameter in z
	// 
	// lzdepend : angular momenta velocity dependance(0:const, 1:flat)
	// betab : plasma beta @ bottom of torus
	// rhob : density (g cm^-3) @ bottom of torus (torus component)
	// halorhob : halo density @ bottom of torus (halo component)
	// polytropek : torus temperature
	// polytopeindex : structure of torus parameter
	// cshalo2 : cshalo**2, cshalo is sound of speed in halo
	//=================================================================
	
	
	int i, j, k;
	double rb=1.0;
	double m[3]={1.95, 17.4, 73.5};
	double a[3]={0.0, 0.62, 0.0};
	double b[3]={0.047, 0.015, 3.12};
	double lzdepend=0.95;
	double polyindex=1.05;
	double discbeta=1.0e+6;
	double halorhob=1.0e-4;
	double csdisc2=1.33e+8*1.0e+4/(v0*v0);
	double mdisc=0.2, rdisc=0.4, zdisc=0.03;
	double discrho0=(1.0e+10*msun/(mp*1000*kpc*kpc*kpc)
									 *mdisc/(pi2*zdisc));
	double cshalo2=1.33e+8*2.4e+6/(v0*v0);
	double discbphi0;
	double x2, rb2, oneorb, oneorb2;
	double a0zk2b02, a1zk2b12, a2zk2b22;
	double lzb, potb, psib, dummyrho;
	double discrho, halorho, discpr, halopr;
	double oneobetab, oneocshalo2;
	double oneopoly, polyone, oneopolyone, polyopolyone;
  double potaxi;
	
	
	rb2=rb*rb;
	oneorb=1.0/rb;
	oneorb2=oneorb*oneorb;
	
	oneopoly=1.0/polyindex;
	polyone=polyindex-1.0;
	oneopolyone=1.0/polyone;
	polyopolyone=polyindex/polyone;
	
	
	lzb=-rb2*sqrt(m[0]/pow(rb2+(a[0]+b[0])*(a[0]+b[0]), 1.5)
								+m[1]/pow(rb2+(a[1]+b[1])*(a[1]+b[1]), 1.5)
								+m[2]/pow(rb2+(a[2]+b[2])*(a[2]+b[2]), 1.5));
	
	potb=-(m[0]/sqrt(rb2+(a[0]+b[0])*(a[0]+b[0]))
				 +m[1]/sqrt(rb2+(a[1]+b[1])*(a[1]+b[1]))
				 +m[2]/sqrt(rb2+(a[2]+b[2])*(a[2]+b[2])));
	
	discbphi0=sqrt(2.0*csdisc2*discrho0/discbeta);
	
	oneocshalo2=1.0/cshalo2;
	
	
#pragma omp parallel for private(i, j, k, x2, a0zk2b02, a1zk2b12, a2zk2b22, halorho, halopr, discrho, discpr, potaxi) firstprivate(discbphi0, rb2, oneorb, oneorb2, oneopoly, polyone, oneopolyone, polyopolyone, lzb, potb, oneocshalo2)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				x2=x[i]*x[i];
				
				a0zk2b02=(a[0]+sqrt(z[k]*z[k]+b[0]*b[0]))
				*(a[0]+sqrt(z[k]*z[k]+b[0]*b[0]));
				a1zk2b12=(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]))
				*(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]));
				a2zk2b22=(a[2]+sqrt(z[k]*z[k]+b[2]*b[2]))
				*(a[2]+sqrt(z[k]*z[k]+b[2]*b[2]));
				
				
				potaxi=-(m[0]/sqrt(x2+a0zk2b02)
                  +m[1]/sqrt(x2+a1zk2b12)
                  +m[2]/sqrt(x2+a2zk2b22));
				
				
				Grxaxi[i][j][k]=-(m[0]/pow(x2+a0zk2b02, 1.5)
													+m[1]/pow(x2+a1zk2b12, 1.5)
													+m[2]/pow(x2+a2zk2b22, 1.5))*x[i];
				
				
				Gryaxi[i][j][k]=0.0;
				
				Grzaxi[i][j][k]=-(m[0]/pow(x2+a0zk2b02, 1.5)
													*(1.0+a[0]/sqrt(z[k]*z[k]+b[0]*b[0]))
													+m[1]/pow(x2+a1zk2b12, 1.5)
													*(1.0+a[1]/sqrt(z[k]*z[k]+b[1]*b[1]))
													+m[2]/pow(x2+a2zk2b22, 1.5)
													*(1.0+a[2]/sqrt(z[k]*z[k]+b[2]*b[2])))*z[k];
				
				
				halorho=halorhob*exp(-oneocshalo2*(potaxi-potb));
				halopr=cshalo2*halorho;
				
				
				discrho=(discrho0
								 *exp(-(x[i]-rdisc)/rdisc)
								 *exp(-fabs(z[k])/zdisc));
				discpr=csdisc2*discrho;
				
									
				//==================================================
				//	disc+halo
				//==================================================
					
				V[0][i][j][k]=discrho+halorho;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=((discrho*
												sqrt(-x[i]*(Grxaxi[i][j][k]+csdisc2/rdisc))
												/V[0][i][j][k])
											 -omegasp*x[i]);
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=(discbphi0
											 *exp(-(x[i]-rdisc)*(x[i]-rdisc)/rdisc)
											 *exp(-fabs(z[k])/zdisc));
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=discpr+halopr;
				V[8][i][j][k]=0.0;
				
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	k=0 impose symmetric boundary
	//==================================================
	
	ZLeftVzAntisymBSymBoundary(U);
	ZLeftFreeBoundary(U);
	//ZLeftPeriodicBoundary(U);
	//ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial condition.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	galactic torus+halo model
//============================================================
/*

int TorusHalo(void)
{

	
	//=================================================================
	// rb : distance from polar to bottom of torus
	// 
	// Miyamoto et al.(1980) galactic potential including DM.
	// m : bulge, disk, halo mass
	// a : bulge, disk, halo distribution parameter in r
	// b : bulge, disk, halo distribution parameter in z
	// 
	// lzdepend : angular momenta velocity dependance(0:const, 1:flat)
	// betab : plasma beta @ bottom of torus
	// rhob : density (g cm^-3) @ bottom of torus (torus component)
	// halorhob : halo density @ bottom of torus (halo component)
	// polytropek : torus temperature
	// polytopeindex : structure of torus parameter
	// cshalo2 : cshalo**2, cshalo is sound of speed in halo
	//=================================================================
	
	
	int i, j, k;
	double rb=1.0;
	double m[3]={1.95, 17.4, 73.5};
	double a[3]={0.0, 0.62, 0.0};
	double b[3]={0.047, 0.015, 3.12};
	double lzdepend=0.95;
	double polyindex=1.05;
	double betab=1.0e+4;
	double halorhob=1.0e-5;
	double csdisk2=1.33e+12*1.05/(v0*v0);
	double cshalo2=2.4*1.33e+14/(v0*v0);
	double x2, rb2, oneorb, oneorb2;
	double a0zk2b02, a1zk2b12, a2zk2b22;
	double lzb, potb, psib, dummyrho;
	double diskrho, halorho, diskpr, halopr;
	double oneobetab, oneocshalo2;
	double oneopoly, polyone, oneopolyone, polyopolyone;
	
	
	rb2=rb*rb;
	oneorb=1.0/rb;
	oneorb2=oneorb*oneorb;
	
	oneopoly=1.0/polyindex;
	polyone=polyindex-1.0;
	oneopolyone=1.0/polyone;
	polyopolyone=polyindex/polyone;
	
	
	lzb=-rb2*sqrt(m[0]/pow(rb2+(a[0]+b[0])*(a[0]+b[0]), 1.5)
								+m[1]/pow(rb2+(a[1]+b[1])*(a[1]+b[1]), 1.5)
								+m[2]/pow(rb2+(a[2]+b[2])*(a[2]+b[2]), 1.5));
	
	potb=-(m[0]/sqrt(rb2+(a[0]+b[0])*(a[0]+b[0]))
				 +m[1]/sqrt(rb2+(a[1]+b[1])*(a[1]+b[1]))
				 +m[2]/sqrt(rb2+(a[2]+b[2])*(a[2]+b[2])));
	
	oneobetab=1.0/betab;
	
	oneocshalo2=1.0/cshalo2;


	psib=potb+csdisk2*oneopolyone*(1.0+oneobetab)
				-1.0/(2.0*(lzdepend-1.0))*(lzb*lzb*oneorb2);
	
	
#pragma omp parallel for private(i, j, k, x2, a0zk2b02, a1zk2b12, a2zk2b22, halorho, halopr, dummyrho, diskrho, diskpr) firstprivate(rb2, oneorb, oneorb2, oneopoly, polyone, oneopolyone, polyopolyone, lzb, potb, oneobetab, oneocshalo2, psib)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				x2=x[i]*x[i];
				
				a0zk2b02=(a[0]+sqrt(z[k]*z[k]+b[0]*b[0]))
									*(a[0]+sqrt(z[k]*z[k]+b[0]*b[0]));
				a1zk2b12=(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]))
									*(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]));
				a2zk2b22=(a[2]+sqrt(z[k]*z[k]+b[2]*b[2]))
									*(a[2]+sqrt(z[k]*z[k]+b[2]*b[2]));
				
				
				PotAxi[i][j][k]=-(m[0]/sqrt(x2+a0zk2b02)
													+m[1]/sqrt(x2+a1zk2b12)
													+m[2]/sqrt(x2+a2zk2b22));

				
				GrxAxi[i][j][k]=-(m[0]/pow(x2+a0zk2b02, 1.5)
													+m[1]/pow(x2+a1zk2b12, 1.5)
													+m[2]/pow(x2+a2zk2b22, 1.5))*x[i];

				
				GryAxi[i][j][k]=0.0;				
				
				GrzAxi[i][j][k]=-(m[0]/pow(x2+a0zk2b02, 1.5)
													*(1.0+a[0]/sqrt(z[k]*z[k]+b[0]*b[0]))
													+m[1]/pow(x2+a1zk2b12, 1.5)
													*(1.0+a[1]/sqrt(z[k]*z[k]+b[1]*b[1]))
													+m[2]/pow(x2+a2zk2b22, 1.5)
													*(1.0+a[2]/sqrt(z[k]*z[k]+b[2]*b[2])))*z[k];

				
				halorho=halorhob*exp(-oneocshalo2*(PotAxi[i][j][k]-potb));
				halopr=cshalo2*halorho;
				
				dummyrho=psib-PotAxi[i][j][k]
									+lzb*lzb/(2.0*(lzdepend-1.0)*x2)
										*pow(x2*oneorb2, lzdepend);
				

				if(dummyrho<0.0){
					
					//==================================================
					//	halo
					//==================================================
					
					V[0][i][j][k]=halorho;
					V[2][i][j][k]=0.0;
					V[5][i][j][k]=0.0;
					V[7][i][j][k]=halopr;
					
					
				}
				else{
				
					diskrho=pow(dummyrho/
											(csdisk2*oneopolyone*
											 (1.0+oneobetab*pow(x2*oneorb2, polyone))), 
											oneopolyone);
					
					diskpr=csdisk2*oneopoly*pow(diskrho, polyindex);
					
					
					//==================================================
					//	disk+halo
					//==================================================
					
					V[0][i][j][k]=diskrho+halorho;
					V[7][i][j][k]=diskpr+halopr;
					
					V[2][i][j][k]=lzb*pow(x[i]*oneorb, lzdepend)*oneox[i]
												*diskrho/V[0][i][j][k];
					V[5][i][j][k]=-sqrt(2.0*oneobetab
															*pow(x2*oneorb2, polyone)*diskpr);
					
					
				}
				
				
				V[1][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	k=0 impose symmetric boundary
	//==================================================
	 
	ZLeftVzAntisymBSymBoundary(U);
	//ZLeftFreeBoundary(U);
	
	//==================================================
	//	save initial condition.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}
*/


//============================================================
//	Rayleigh-Taylor instability condition (no magnetic fields)
//============================================================


int RayleighTaylor(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				if(z[k]<0.0)	V[0][i][j][k]=1.0;
				else V[0][i][j][k]=2.0;
				
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.01*(1.0+cos(pi2*x[i]/lx))
											*(1.0+cos(pi2*z[k]/lz))*0.25;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=2.5-V[0][i][j][k]*1.0*z[k];
				V[8][i][j][k]=0.0;

				
			}
		}
	}
	

	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	/*
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	*/
	/*
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	*/
	/*
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	*/
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	 
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	/*
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	 
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	*/
	
	//============================================================
	//	k=0 & k=kzmax1 fixed boundary
	//============================================================
	
	ZLeftFixedBoundary(U);
	ZRightFixedBoundary(U);
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
		
}



//============================================================
//	Jet propagation (including magnetic fields)
//============================================================


int JetMHD(void)
{
	
	int i, j, k;
	double plasmabeta=100.0;
	double oneosqrtbeta;
	double r, oneor;
	
	
	oneosqrtbeta=1.0/sqrt(plasmabeta*0.5);

	
#pragma omp parallel for private(i, j, k, r, oneor) firstprivate(oneosqrtbeta)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
				r=sqrt(x[i]*x[i]+y[j]*y[j]);
				
				
				if(r<=0.05*lx && z[k]<=0.1666*lz){

					oneor=1.0/r;

					
					V[0][i][j][k]=0.1;
					V[3][i][j][k]=6.0*gm;
					V[4][i][j][k]=sqrt(V[7][i][j][k])*oneosqrtbeta*y[j]*oneor;
					V[5][i][j][k]=-sqrt(V[7][i][j][k])*oneosqrtbeta*x[i]*oneor;
					
					
				}
				else{
					
					V[0][i][j][k]=1.0;
					V[3][i][j][k]=0.0;
					V[4][i][j][k]=0.0;
					V[5][i][j][k]=0.0;
				
				
				}
				
				
			}
		}
	}
	

	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	/*
	 //==================================================
	 //	i=0 & i=ixmax1 periodic boundary
	 //==================================================
	 
	 XLeftPeriodicBoundary(U);
	 XRightPeriodicBoundary(U);
	 */
	/*
	 //==================================================
	 //	j=0 & j=jymax1 periodic boundary
	 //==================================================
	 
	 YLeftPeriodicBoundary(U);
	 YRightPeriodicBoundary(U);
	 */
	/*
	 //==================================================
	 //	k=0 & k=kzmax1 periodic boundary
	 //==================================================
	 
	 ZLeftPeriodicBoundary(U);
	 ZRightPeriodicBoundary(U);
	 */
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Jet propagation (no magnetic fields)
//============================================================


int JetHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				

				if(sqrt((x[i]-0.5*lx)*(x[i]-0.5*lx))<=0.05*lx
					 && z[k]<=0.1666*lz){
					V[0][i][j][k]=0.1;
					V[3][i][j][k]=6.0*gm;
				}
				else{
					V[0][i][j][k]=1.0;
					V[3][i][j][k]=0.0;
				}
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	/*
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	*/
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	/*
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	*/
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	/*
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	*/
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Magnetic Field Loop Advection 
//============================================================


int MagneticFieldLoopAdvection(void)
{
	
	int i, j, k;
	double alpha=0.001;
	double loopR=0.3;
	
	
#pragma omp parallel for private(i, j, k) firstprivate(alpha, loopR)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=cos(pi/3.0);
				V[2][i][j][k]=sin(pi/3.0);
				V[3][i][j][k]=0.0;
				
				if(sqrt((x[i]-1.0)*(x[i]-1.0)+(y[j]-0.5)*(y[j]-0.5))<loopR){
					V[4][i][j][k]=-alpha*y[j]
					/sqrt((x[i]-1.0)*(x[i]-1.0)
								+(y[j]-0.5)*(y[j]-0.5));
					V[5][i][j][k]=alpha*x[i]
					/sqrt((x[i]-1.0)*(x[i]-1.0)
								+(y[j]-0.5)*(y[j]-0.5));
				}
				else{
					V[4][i][j][k]=0.0;
					V[5][i][j][k]=0.0;
				}
				
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Magnetic Field Loop Advection in Cylindrical coord.
//============================================================


int MagneticFieldLoopAdvectionCylinder(void)
{
	
	int i, j, k;
	double alpha=0.001;
	double loopR=0.3;
	double cosyj, sinyj;
	double bx, by;
	

#pragma omp parallel for private(i, j, k, cosyj, sinyj, bx, by) firstprivate(alpha, loopR)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				cosyj=cos(y[j]);
				sinyj=sin(y[j]);
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=4.0/3.0*x[i];
				V[3][i][j][k]=0.0;
				
				if(sqrt((x[i]*cosyj-1.5)*(x[i]*cosyj-1.5)
								+x[i]*sinyj*x[i]*sinyj)<loopR){

					bx=-alpha*(x[i]*sinyj
										 /sqrt((x[i]*cosyj-1.5)*(x[i]*cosyj-1.5)
													 +x[i]*sinyj*x[i]*sinyj));
					by=alpha*((x[i]*cosyj-1.5)
										/sqrt((x[i]*cosyj-1.5)*(x[i]*cosyj-1.5)
													+x[i]*sinyj*x[i]*sinyj));
					
					
					V[4][i][j][k]=bx*cosyj+by*sinyj;
					V[5][i][j][k]=-bx*sinyj+by*cosyj;
					
					
				}
				else{
					V[4][i][j][k]=0.0;
					V[5][i][j][k]=0.0;
				}
				
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	//==================================================
	//	i=0 & i=ixmax1 initialize boundary
	//==================================================
	
	XLeftInitializeBoundary(U);
	XRightInitializeBoundary(U);
	
	/*
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	*/
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave 
//		in Cylinder coord. (oblique magnetic fields exist)
//============================================================


int SphericalBlastWaveMHDCylinder(void)
{
	
	int i, j, k;
	double cosyj, sinyj;
	double backb=1.0/sqrt(2.0);
	
	
#pragma omp parallel for private(i, j, k, cosyj, sinyj)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				cosyj=cos(y[j]);
				sinyj=sin(y[j]);
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=backb*(cosyj+sinyj);
				V[5][i][j][k]=backb*(-sinyj+cosyj);
				V[6][i][j][k]=0.0;
				
				if(sqrt((x[i]*cosyj-1.5)*(x[i]*cosyj-1.5)
								+x[i]*sinyj*x[i]*sinyj)<0.1){
					V[7][i][j][k]=10.0;
				}
				else V[7][i][j][k]=0.1;
				
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	//==================================================
	//	i=0 & i=ixmax1 initialize boundary
	//==================================================
	
	XLeftInitializeBoundary(U);
	XRightInitializeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 initialize boundary
	//==================================================
	
	YLeftInitializeBoundary(U);
	YRightInitializeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave (oblique magnetic fields exist)
//============================================================


int SphericalBlastWaveMHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=1.0/sqrt(2.0);
				V[5][i][j][k]=1.0/sqrt(2.0);
				V[6][i][j][k]=0.0;
				
				if(sqrt((x[i]-0.5)*(x[i]-0.5)+(y[j]-0.5)*(y[j]-0.5))<0.1){
					V[7][i][j][k]=10.0;
				}
				else{
					V[7][i][j][k]=0.1;
				}
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave in Cylinder coord(no magnetic fields)
//============================================================


int SphericalBlastWaveHDCylinder(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				
				if(sqrt(pow(x[i]*cos(y[j])-1.5, 2.0)
								+x[i]*sin(y[j])*x[i]*sin(y[j]))<0.1){
					V[7][i][j][k]=10.0;
				}
				else{
					V[7][i][j][k]=0.1;
				}
				
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	/*
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	*/
	/*
	 //==================================================
	 //	k=0 & k=kzmax1 free boundary
	 //==================================================
	 
	 ZLeftFreeBoundary(U);
	 ZRightFreeBoundary(U);
	 */
	/*
	 //==================================================
	 //	i=0 & i=ixmax1 periodic boundary
	 //==================================================
	 
	 XLeftPeriodicBoundary(U);
	 XRightPeriodicBoundary(U);
	 */
	
	 //==================================================
	 //	j=0 & j=jymax1 periodic boundary
	 //==================================================
	 
	 YLeftPeriodicBoundary(U);
	 YRightPeriodicBoundary(U);
	 
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Spherical Blast Wave (no magnetic fields)
//============================================================


int SphericalBlastWaveHD(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				
				if(sqrt(x[i]*x[i]+y[j]*y[j])<0.1){
					V[7][i][j][k]=10.0;
				}
				else{
					V[7][i][j][k]=0.1;
				}
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	/*
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	*/
	/*
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	*/
	/*
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	*/
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Orszag-Tang vortex problem (Orszag & Tang 1998)
//============================================================


int OrszagTangZX(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=gm*gm*oneopi4;
				V[1][i][j][k]=sin(pi2*z[k]);
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=-sin(pi2*x[i]);
				V[4][i][j][k]=sin(pi4*z[k])/sqrt(pi4);
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=-sin(pi2*x[i])/sqrt(pi4);
				V[7][i][j][k]=gm*oneopi4;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Orszag-Tang vortex problem (Orszag & Tang 1998)
//============================================================


int OrszagTangXY(void)
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=gm*gm*oneopi4;
				V[1][i][j][k]=-sin(pi2*y[j]);
				V[2][i][j][k]=sin(pi2*x[i]);
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=-sin(pi2*y[j])/sqrt(pi4);
				V[5][i][j][k]=sin(pi4*x[i])/sqrt(pi4);
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=gm*oneopi4;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 1a MHD shock tube z-direction
//============================================================


int DW1aShockTubeZ(void)
{
	
	int i, j, k;
	
	
	//-----0<=k<kz/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax/2; k++){
				
				
				V[0][i][j][k]=1.08;
				V[1][i][j][k]=0.01;
				V[2][i][j][k]=0.5;
				V[3][i][j][k]=1.2;
				V[4][i][j][k]=3.6/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=0.95;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----kz/2<=k<kz+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=kzmax/2; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=4.0/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	/*
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	*/
	/*
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	*/
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	i=0 & i=ixmax1 periodic boundary
	//==================================================
	
	XLeftPeriodicBoundary(U);
	XRightPeriodicBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 periodic boundary
	//==================================================
	
	YLeftPeriodicBoundary(U);
	YRightPeriodicBoundary(U);
	
	/*
	//==================================================
	//	k=0 & k=kzmax1 periodic boundary
	//==================================================
	
	ZLeftPeriodicBoundary(U);
	ZRightPeriodicBoundary(U);	
	*/
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 1a MHD shock tube y-direction
//============================================================


int DW1aShockTubeY(void)
{
	
	int i, j, k;
	
	
	//-----0<=j<jy/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax/2; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.08;
				V[1][i][j][k]=0.5;
				V[2][i][j][k]=1.2;
				V[3][i][j][k]=0.01;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=3.6/sqrt(pi4);
				V[7][i][j][k]=0.95;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----jy/2<=j<jy+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=jymax/2; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=2.0/sqrt(pi4);
				V[6][i][j][k]=4.0/sqrt(pi4);
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 1a MHD shock tube x-direction
//============================================================


int DW1aShockTubeX(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.08;
				V[1][i][j][k]=1.2;
				V[2][i][j][k]=0.01;
				V[3][i][j][k]=0.5;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=3.6/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=0.95;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ix+1; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=2.0/sqrt(pi4);
				V[5][i][j][k]=4.0/sqrt(pi4);
				V[6][i][j][k]=2.0/sqrt(pi4);
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Dai & Woodward 1995, Table 6 MHD shock tube test
//============================================================


int DW6ShockTube(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=10.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=5.0/sqrt(pi4);
				V[5][i][j][k]=5.0/sqrt(pi4);
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=20.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				V[0][i][j][k]=1.0;
				V[1][i][j][k]=-10.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=5.0/sqrt(pi4);
				V[5][i][j][k]=5.0/sqrt(pi4);
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=1.0;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}



//============================================================
//	Davis (????), HydroDynamic shock tube
//============================================================


int HDShockTube(void)
{
	
	int i, j, k;
	
	
	//-----0<=i<ix/2 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax/2; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
			
				V[0][i][j][k]=0.445;
				V[1][i][j][k]=0.744;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=5.875;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//-----ix/2<=i<ix+1 state-----
	
#pragma omp parallel for private(i, j, k)
	for(i=ixmax/2; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				V[0][i][j][k]=0.5;
				V[1][i][j][k]=0.0;
				V[2][i][j][k]=0.0;
				V[3][i][j][k]=0.0;
				V[4][i][j][k]=0.0;
				V[5][i][j][k]=0.0;
				V[6][i][j][k]=0.0;
				V[7][i][j][k]=0.952;
				V[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	//==================================================
	//	convert from V to U.
	//==================================================
	
	VtoU();
	
	
	//==================================================
	//	i=0 & i=ixmax1 free boundary
	//==================================================
	
	XLeftFreeBoundary(U);
	XRightFreeBoundary(U);
	
	
	//==================================================
	//	j=0 & j=jymax1 free boundary
	//==================================================
	
	YLeftFreeBoundary(U);
	YRightFreeBoundary(U);
	
	
	//==================================================
	//	k=0 & k=kzmax1 free boundary
	//==================================================
	
	ZLeftFreeBoundary(U);
	ZRightFreeBoundary(U);
	
	
	//==================================================
	//	save initial conditions.
	//==================================================
	
	SaveInitial();
	
	
	return 0;
	
}


