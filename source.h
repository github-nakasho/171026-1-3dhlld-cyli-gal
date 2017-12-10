

/***********************************************************************
 *
 *	source.h 
 *
 *	set source term.
 *
 *	
 *	2013 Aug. 10 : improve SpiralGravityCooling.
 *	2013 May. 10 : add SpiralGravity.
 *	2013 Apr. 13 : add AxiGravity.
 *	2013 Apr. 05 : change NoSource in Cylindrical coord.
 *	2013 Jan. 16 : add GalacticGravity.
 *	2013 Jan. 01 : coded by Sho Nakamura (Tohoku Univ.).
 *	
 **********************************************************************/



//============================================================
//	existing Spiral Gravity & cooling energy loss
//============================================================


int SpiralGravityCooling(double dummyV[][ixmax][jymax][kzmax], 
												 double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
  
	//-----no-spiral region 0<r<rin-----
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<irin; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				
				//-----GryAxi=0.0-----
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=((dummyV[0][i][j][k]
														 *dummyV[2][i][j][k]*dummyV[2][i][j][k]
														 -dummyV[5][i][j][k]*dummyV[5][i][j][k]
														 +(dummyV[7][i][j][k]
															 +0.5*(dummyV[4][i][j][k]
																		 *dummyV[4][i][j][k]
																		 +dummyV[5][i][j][k]
																		 *dummyV[5][i][j][k]
																		 +dummyV[6][i][j][k]
																		 *dummyV[6][i][j][k])))
														*oneox[i]
														+dummyV[0][i][j][k]*Grxaxi[i][j][k]
														+dummyV[0][i][j][k]*omegasp*omegasp*x[i]
														+(2.0*dummyV[0][i][j][k]
															*dummyV[2][i][j][k]*omegasp));
				dummyS[2][i][j][k]=x[i]*(-2.0*dummyV[0][i][j][k]
                                 *dummyV[1][i][j][k]*omegasp);
				dummyS[3][i][j][k]=dummyV[0][i][j][k]*Grzaxi[i][j][k];
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=(dummyV[5][i][j][k]*dummyV[1][i][j][k]
														-dummyV[2][i][j][k]*dummyV[4][i][j][k])
														*oneox[i];
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=(dummyV[0][i][j][k]
														*(dummyV[1][i][j][k]*Grxaxi[i][j][k]
															+dummyV[3][i][j][k]*Grzaxi[i][j][k])
                                    -(dummyV[0][i][j][k]*dummyV[0][i][j][k]
                                      *Cooling(dummyV[0][i][j][k],
                                               dummyV[7][i][j][k]))
                                    +dummyV[0][i][j][k]
														*x[i]*dummyV[1][i][j][k]*omegasp*omegasp);
        dummyS[8][i][j][k]=0.0;

        
			}
		}
	}
	

	//-----spiral region r>rin-----

#pragma omp parallel for private(i, j, k)
  for(i=irin; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=((dummyV[0][i][j][k]
														 *dummyV[2][i][j][k]*dummyV[2][i][j][k]
														 -dummyV[5][i][j][k]*dummyV[5][i][j][k]
														 +(dummyV[7][i][j][k]
															 +0.5*(dummyV[4][i][j][k]
																		 *dummyV[4][i][j][k]
																		 +dummyV[5][i][j][k]
																		 *dummyV[5][i][j][k]
																		 +dummyV[6][i][j][k]
																		 *dummyV[6][i][j][k])))
														*oneox[i]
														+dummyV[0][i][j][k]*Grx[i][j][k]
														+dummyV[0][i][j][k]*omegasp*omegasp*x[i]
														+(2.0*dummyV[0][i][j][k]
															*dummyV[2][i][j][k]*omegasp));
				dummyS[2][i][j][k]=x[i]*(dummyV[0][i][j][k]*Gry[i][j][k]
                                 -(2.0*dummyV[0][i][j][k]
                                   *dummyV[1][i][j][k]*omegasp));
				dummyS[3][i][j][k]=dummyV[0][i][j][k]*Grz[i][j][k];
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=(dummyV[5][i][j][k]*dummyV[1][i][j][k]
														-dummyV[2][i][j][k]*dummyV[4][i][j][k])
                            *oneox[i];
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=(dummyV[0][i][j][k]
                            *(dummyV[1][i][j][k]*Grx[i][j][k]
                              +dummyV[2][i][j][k]*Gry[i][j][k]
                              +dummyV[3][i][j][k]*Grz[i][j][k])
                            -(dummyV[0][i][j][k]*dummyV[0][i][j][k]
                              *Cooling(dummyV[0][i][j][k], 
                                       dummyV[7][i][j][k]))
														+dummyV[0][i][j][k]
														*x[i]*dummyV[1][i][j][k]*omegasp*omegasp);
        dummyS[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	Axisymmetric Gravity Field with cooling
//============================================================


int AxiGravityCooling(double dummyV[][ixmax][jymax][kzmax],
                      double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
        //-----GryAxi=0.0-----
        
        dummyS[0][i][j][k]=0.0;
        dummyS[1][i][j][k]=((dummyV[0][i][j][k]
                             *dummyV[2][i][j][k]*dummyV[2][i][j][k]
                             -dummyV[5][i][j][k]*dummyV[5][i][j][k]
                             +(dummyV[7][i][j][k]
                               +0.5*(dummyV[4][i][j][k]
                                     *dummyV[4][i][j][k]
                                     +dummyV[5][i][j][k]
                                     *dummyV[5][i][j][k]
                                     +dummyV[6][i][j][k]
                                     *dummyV[6][i][j][k])))
                            *oneox[i]
                            +dummyV[0][i][j][k]*Grxaxi[i][j][k]
                            +dummyV[0][i][j][k]*omegasp*omegasp*x[i]
                            +(2.0*dummyV[0][i][j][k]
                              *dummyV[2][i][j][k]*omegasp));
        dummyS[2][i][j][k]=x[i]*(-2.0*dummyV[0][i][j][k]
                                 *dummyV[1][i][j][k]*omegasp);
        dummyS[3][i][j][k]=dummyV[0][i][j][k]*Grzaxi[i][j][k];
        dummyS[4][i][j][k]=0.0;
        dummyS[5][i][j][k]=((dummyV[5][i][j][k]*dummyV[1][i][j][k]
                            -dummyV[2][i][j][k]*dummyV[4][i][j][k])
                            *oneox[i]);
        dummyS[6][i][j][k]=0.0;
        dummyS[7][i][j][k]=(dummyV[0][i][j][k]
                            *(dummyV[1][i][j][k]*Grxaxi[i][j][k]
                              +dummyV[3][i][j][k]*Grzaxi[i][j][k])
                            -(dummyV[0][i][j][k]*dummyV[0][i][j][k]
                              *Cooling(dummyV[0][i][j][k],
                                       dummyV[7][i][j][k]))
                            +dummyV[0][i][j][k]
                            *x[i]*dummyV[1][i][j][k]*omegasp*omegasp);
        dummyS[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	return 0;
	
}


//============================================================
//	Axisymmetric Gravity Field with cooling
//============================================================


int AxiGravity(double dummyV[][ixmax][jymax][kzmax],
               double dummyS[][ixmax][jymax][kzmax])
{
  
  int i, j, k;
  
  
#pragma omp parallel for private(i, j, k)
  for(i=0; i<ixmax; i++){
    for(j=0; j<jymax; j++){
      for(k=0; k<kzmax; k++){
        
        //-----GryAxi=0.0-----
        
        dummyS[0][i][j][k]=0.0;
        dummyS[1][i][j][k]=((dummyV[0][i][j][k]
                             *dummyV[2][i][j][k]*dummyV[2][i][j][k]
                             -dummyV[5][i][j][k]*dummyV[5][i][j][k]
                             +(dummyV[7][i][j][k]
                               +0.5*(dummyV[4][i][j][k]
                                     *dummyV[4][i][j][k]
                                     +dummyV[5][i][j][k]
                                     *dummyV[5][i][j][k]
                                     +dummyV[6][i][j][k]
                                     *dummyV[6][i][j][k])))
                            *oneox[i]
                            +dummyV[0][i][j][k]*Grxaxi[i][j][k]
                            +dummyV[0][i][j][k]*omegasp*omegasp*x[i]
                            +(2.0*dummyV[0][i][j][k]
                              *dummyV[2][i][j][k]*omegasp));
        dummyS[2][i][j][k]=x[i]*(-2.0*dummyV[0][i][j][k]
                                 *dummyV[1][i][j][k]*omegasp);
        dummyS[3][i][j][k]=dummyV[0][i][j][k]*Grzaxi[i][j][k];
        dummyS[4][i][j][k]=0.0;
        dummyS[5][i][j][k]=((dummyV[5][i][j][k]*dummyV[1][i][j][k]
                             -dummyV[2][i][j][k]*dummyV[4][i][j][k])
                            *oneox[i]);
        dummyS[6][i][j][k]=0.0;
        dummyS[7][i][j][k]=(dummyV[0][i][j][k]
                            *(dummyV[1][i][j][k]*Grxaxi[i][j][k]
                              +dummyV[3][i][j][k]*Grzaxi[i][j][k])
                            +dummyV[0][i][j][k]
                            *x[i]*dummyV[1][i][j][k]*omegasp*omegasp);
        dummyS[8][i][j][k]=0.0;
        
        
      }
    }
  }
  
  
  return 0;
  
}



//============================================================
//	uniform gravitational field in z
//============================================================


int UniformGz(double dummyV[][ixmax][jymax][kzmax], 
							double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=0.0;
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=-dummyV[0][i][j][k]*1.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=0.0;
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=-dummyV[3][i][j][k]*dummyV[0][i][j][k]*1.0;
				dummyS[8][i][j][k]=0.0;
				
				
			}
		}
	}
	
	
	return 0;
	
}



//============================================================
//	no source in cylindrical coord.
//============================================================


int NoSource(double dummyV[][ixmax][jymax][kzmax], 
						 double dummyS[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	
	
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				dummyS[0][i][j][k]=0.0;
				dummyS[1][i][j][k]=(dummyV[0][i][j][k]
														*dummyV[2][i][j][k]*dummyV[2][i][j][k]
														-dummyV[5][i][j][k]*dummyV[5][i][j][k]
														+(dummyV[7][i][j][k]
															+0.5*(dummyV[4][i][j][k]
																		*dummyV[4][i][j][k]
																		+dummyV[5][i][j][k]
																		*dummyV[5][i][j][k]
																		+dummyV[6][i][j][k]
																		*dummyV[6][i][j][k])))
														*oneox[i];
				dummyS[2][i][j][k]=0.0;
				dummyS[3][i][j][k]=0.0;
				dummyS[4][i][j][k]=0.0;
				dummyS[5][i][j][k]=(dummyV[5][i][j][k]*dummyV[1][i][j][k]
														-dummyV[2][i][j][k]*dummyV[4][i][j][k])
														*oneox[i];
				dummyS[6][i][j][k]=0.0;
				dummyS[7][i][j][k]=0.0;
				
				
			}
		}
	}

	return 0;

	
}





