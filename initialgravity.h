


/***********************************************************************
 *
 *	initialgravity.h
 *
 *	set initial gravity.
 *
 *	
 *	2014 Dec. 25 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//============================================================
//	initial garvity force
//============================================================


int InitialGravity(void)
{
  
  int i, j, k;
  double x2, z2;
  
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
  
  
  
#pragma omp parallel for private(i, j, k, x2, z2) firstprivate(vhalo0, dh, a, b, md, db, mb)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        
        x2=x[i]*x[i];
        z2=z[k]*z[k];
        
        
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
        
                
        
      }
    }
  }
  
  
  return 0;
  
}



