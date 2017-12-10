

/***********************************************************************
 *
 *	spiralgravity.h
 *
 *	spiral arm gavity field.
 *
 *
 *	2015 July 23 : add omegapattern.
 *	2014 July 20 : add epssp, switchfunc.
 *	2013 Aug. 22 : add SpiralGravity model.
 *	2013 Aug. 22 : add SpiralGravity model.
 *	2013 Jan. 16 : coded by Sho Nakamura (Tohoku Univ.)
 *
 **********************************************************************/



int SpiralGravity(void)
{
	
	int i, j, k;
	double m[3]={0.0, 10.0, 0.0};
	double a[3]={0.0, 0.65, 0.0};
	double b[3]={0.0, 0.026, 0.0};
	double oneoscalephase, oneoscaler, cotpitch;
	double x2, xoscalephase, xoscaler, xoscaler2, x2oscaler2;
	double oneoonex2oscaler2, oneoonex2oscaler215;
	double xoscaler2oonex2oscaler215, x2oscaler2oonex2oscaler215;
	double phase, cosphase, sinphase;
	double z2, a1zk2b12, oneosqrtx2a1zk2b12, scalezosqrtz2scalez2;
	double spiralfr, spiralfphi, spiralfz;
  double epssp, switchfunc;
	
	
	oneoscalephase=1.0/scalephase;
	oneoscaler=1.0/scaler;
	cotpitch=cot(pitch);
  
  switchfunc=(t-tausp)/fabs(t-tausp);
  epssp=0.5*epssp0*(t/tausp*(1.0-switchfunc)+(1.0+switchfunc));
	
	
#pragma omp parallel for private(i, j, k, x2, xoscalephase, xoscaler, xoscaler2, x2oscaler2, oneoonex2oscaler2, oneoonex2oscaler215, xoscaler2oonex2oscaler215, x2oscaler2oonex2oscaler215, phase, cosphase, sinphase, z2, a1zk2b12, oneosqrtx2a1zk2b12, scalezosqrtz2scalez2, spiralfr, spiralfphi, spiralfz) firstprivate(oneoscalephase, oneoscaler, cotpitch, epssp)
	for(i=0; i<ixmax; i++){
		
		
		x2=x[i]*x[i];
		xoscalephase=x[i]*oneoscalephase;
		xoscaler=x[i]*oneoscaler;
		xoscaler2=xoscaler*oneoscaler;
		x2oscaler2=xoscaler*xoscaler;
		oneoonex2oscaler2=1.0/(1.0+x2oscaler2);
		oneoonex2oscaler215=1.0*pow(oneoonex2oscaler2, 1.5);
		xoscaler2oonex2oscaler215=xoscaler2*oneoonex2oscaler215;
		x2oscaler2oonex2oscaler215=x2oscaler2*oneoonex2oscaler215;
		
		
		for(j=0; j<jymax; j++){
			
			
			phase=armnum*(y[j]-(omegapattern-omegasp)*t
                    +cotpitch*log(xoscalephase));
			sincos(phase, &sinphase, &cosphase);
			
			
      for(k=0; k<kzmax; k++){
				
				
				z2=z[k]*z[k];
				a1zk2b12=((a[1]+sqrt(z2+b[1]*b[1]))
									*(a[1]+sqrt(z2+b[1]*b[1])));
				
				oneosqrtx2a1zk2b12=1.0/sqrt(x2+a1zk2b12);
				
				scalezosqrtz2scalez2=(scalez/sqrt(z2+scalez*scalez));

				
        
				//-----spiral gravity in r-----
				
				spiralfr=(-epssp*m[1]*oneosqrtx2a1zk2b12
									*scalezosqrtz2scalez2
									*xoscaler2oonex2oscaler215
									*((x[i]*oneosqrtx2a1zk2b12*oneosqrtx2a1zk2b12
										 -((2.0-x2oscaler2)*oneoonex2oscaler2))*cosphase
										+armnum*cotpitch*sinphase));
				
				
				//-----spiral gravity in phi-----
				
				spiralfphi=(epssp*m[1]*oneosqrtx2a1zk2b12
										*scalezosqrtz2scalez2
										*xoscaler2oonex2oscaler215
										*armnum*sinphase);
        
				
				//-----spiral gravity in z-----
				
				spiralfz=(-epssp*m[1]*oneosqrtx2a1zk2b12
									*z[k]*scalezosqrtz2scalez2
									*x2oscaler2oonex2oscaler215
									*((1.0+a[1]/sqrt(z2+b[1]*b[1]))
										*oneosqrtx2a1zk2b12*oneosqrtx2a1zk2b12
										+1.0/(z2+scalez*scalez))*cosphase);
				
				
				//-----bulge, disk, halo & spiral gravity-----
				
				Grx[i][j][k]=Grxaxi[i][j][k]+spiralfr;
				Gry[i][j][k]=Gryaxi[i][j][k]+spiralfphi;
				Grz[i][j][k]=Grzaxi[i][j][k]+spiralfz;
				
        
				
			}
		}
	}

	
  return 0;
	
}


/*
int SpiralGravity(void)
{
	
	int i, j, k;
	double x2;
	double phase, cosphase, sinphase;
	double m[3]={1.95, 17.4, 73.5};
	double a[3]={0.0, 0.62, 0.0};
	double b[3]={0.047, 0.015, 3.12};
	double a1zk2b12, oneosqrtx2a1zk2b12, scalezosqrtz2scalez2;
	double spiralfr, spiralfphi, spiralfz;
	
	
#pragma omp parallel for private(i, j, k, x2, phase, cosphase, sinphase, a1zk2b12, oneosqrtx2a1zk2b12, scalezosqrtz2scalez2, spiralfr, spiralfphi, spiralfz)
	for(i=0; i<ixmax; i++){
		
		x2=x[i]*x[i];
		
		
		for(j=0; j<jymax; j++){
			
			phase=armnum*(-y[j]-omegasp*t+cot(pitch)*log(x[i]/scaler));
			sincos(phase, &sinphase, &cosphase);

			
      for(k=0; k<kzmax; k++){
				
				a1zk2b12=(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]))
									*(a[1]+sqrt(z[k]*z[k]+b[1]*b[1]));
				
				oneosqrtx2a1zk2b12=1.0/sqrt(x2+a1zk2b12);
				
				scalezosqrtz2scalez2=scalez/sqrt(z[k]*z[k]+scalez*scalez);
				
				
				//-----axisymmetric potential + spiral potential-----
				
				Pot[i][j][k]=PotAxi[i][j][k]
											-m[1]*oneosqrtx2a1zk2b12
											*epssp*scalezosqrtz2scalez2*cosphase;
																	
        
				//-----spiral gravity in r-----
				
				spiralfr=-epssp*scalezosqrtz2scalez2
									*m[1]*oneosqrtx2a1zk2b12
									*(x[i]/(x2+a1zk2b12)*cosphase
										+armnum*oneox[i]*cot(pitch)*sinphase);
				
				
				//-----spiral gravity in phi-----
				
				spiralfphi=oneox[i]*epssp
										*m[1]*oneosqrtx2a1zk2b12
										*armnum*scalezosqrtz2scalez2*sinphase;
        
				//-----spiral gravity in z-----
				
				spiralfz=-epssp*m[1]*oneosqrtx2a1zk2b12*cosphase
									*z[k]*scalezosqrtz2scalez2
									*((1.0+a[1]/sqrt(z[k]*z[k]+b[1]*b[1]))
										*oneosqrtx2a1zk2b12*oneosqrtx2a1zk2b12
										+1.0/(z[k]*z[k]+scalez*scalez));
				
				
				//-----bulge, disk, halo & spiral gravity-----
				
				Grx[i][j][k]=GrxAxi[i][j][k]+spiralfr;
				Gry[i][j][k]=GryAxi[i][j][k]+spiralfphi;
				Grz[i][j][k]=GrzAxi[i][j][k]+spiralfz;
				
        
			}
		}
	}
	
	
  return 0;

}
*/




