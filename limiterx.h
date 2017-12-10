

/***********************************************************************
 *
 *	limiterx.h
 *
 *	slope limiter in x
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : V, Vl, Vr vector component index
 *	a : V[i]-V[i-1]
 *	b : V[i+1]-V[i]
 *	c : 0.25*(V[i+1]-V[i-1]) (MC)
 *	grad : primitive variables gradient using interpolation
 *
 *
 *	2012 Nov. 19 : add VanLeerLimiterX, MCLimiterX, SuperbeeLimiterX.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary for Vl[0][][], Vr[ix-1][][]
 *	2012 Oct. 09 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//==================================================
//	minmod slope limiter in x (2nd order)
//==================================================


int MinmodLimiterX(void)
{		

	int i, j, k;
	int m;
	double a, b, grad, curvaturefactor1, curvaturefactor2;
  double b2, gmpr, fast;
  double lambmax, lambmin;
	
	
	//-----minmod limiter (V--->Vl, Vr)-----
	
	
	//-----left & right state of 1<i<ix-1 cells------

	
/*
  i=1;
	curvaturefactor1=dx[i]/(6.0*x[i]);
	
#pragma omp parallel for private(j, k, m, a, b, grad, b2, gmpr, fast, lambmax, lambmin) firstprivate(i, curvaturefactor1)
	for(m=0; m<9; m++){
		for(j=0; j<jymax1; j++){
			for(k=0; k<kzmax1; k++){
        
				
				a=(V[m][i][j][k]-V[m][i-1][j][k])
        /(1.0-dx[i]*dx[i]/(12.0*x[i]*x[i-1]));
				b=(V[m][i+1][j][k]-V[m][i][j][k])
        /(1.0-dx[i]*dx[i]/(12.0*x[i]*x[i+1]));
				
				
				if(a*b<=0.0) grad=0.0;
				else{
					if(a>0.0) grad=min(a, b);
					else grad=max(a, b);
				}
        
        
        b2=(V[4][i][j][k]*V[4][i][j][k]
            +V[5][i][j][k]*V[5][i][j][k]
            +V[6][i][j][k]*V[6][i][j][k]);
        
        gmpr=gm*V[7][i][j][k];
        
        fast=sqrt(((b2+gmpr)
                   +sqrt((b2+gmpr)*(b2+gmpr)
                         -4.0*gmpr*V[4][i][j][k]*V[4][i][j][k]))
                  /(2.0*V[0][i][j][k]));
        
        
        //lambmax=max(0, V[1][i][j][k]+fast);
        //lambmin=min(0, V[1][i][j][k]-fast);
        
				
        lambmax=max(0, V[1][i][j][k]+fast);
        lambmin=-min(0, V[1][i][j][k]-fast);
        
        
				Vl[m][i][j][k]=(V[m][i][j][k]
                        +(0.5-lambmax*oneo3*dt*oneodx[i])
                        *(1.0-curvaturefactor1)*grad);
				Vr[m][i-1][j][k]=(V[m][i][j][k]
                          -(0.5-lambmin*oneo3*dt*oneodx[i-1])
                          *(1.0-curvaturefactor1)*grad);
        
				
			}
		}
	}
	
	
#pragma omp parallel for private(i, j, k, m, a, b, grad, b2, gmpr, fast, lambmax, lambmin, curvaturefactor1, curvaturefactor2)
	for(m=0; m<9; m++){
		for(i=2; i<ixmax1; i++){		
			for(j=0; j<jymax1; j++){
				for(k=0; k<kzmax1; k++){
					
					
					a=1.5*(V[m][i][j][k]-V[m][i-1][j][k])*dxm[i-1]
          /((xm[i]*xm[i]+xm[i]*xm[i-1]+xm[i-1]*xm[i-1])
            /(xm[i]+xm[i-1])
            -(xm[i-1]*xm[i-1]+xm[i-1]*xm[i-2]+xm[i-2]*xm[i-2])
            /(xm[i-1]+xm[i-2]));
					b=1.5*(V[m][i+1][j][k]-V[m][i][j][k])*dxm[i-1]
          /((xm[i+1]*xm[i+1]+xm[i+1]*xm[i]+xm[i]*xm[i])
            /(xm[i+1]+xm[i])
            -(xm[i]*xm[i]+xm[i]*xm[i-1]+xm[i-1]*xm[i-1])
            /(xm[i]+xm[i-1]));
          
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(a>0.0) grad=min(a, b);
						else grad=max(a, b);
					}
					
					
					curvaturefactor1=((x[i]+oneo3*(0.5*dx[i]-dx[i-1]))
														/(x[i]+0.25*(dx[i]-dx[i-1])));
					curvaturefactor2=((x[i]+oneo3*(dx[i]-0.5*dx[i-1]))
														/(x[i]+0.25*(dx[i]-dx[i-1])));;
          
          
          b2=(V[4][i][j][k]*V[4][i][j][k]
              +V[5][i][j][k]*V[5][i][j][k]
              +V[6][i][j][k]*V[6][i][j][k]);
          
          gmpr=gm*V[7][i][j][k];
          
          fast=sqrt(((b2+gmpr)
                     +sqrt((b2+gmpr)*(b2+gmpr)
                           -4.0*gmpr*V[4][i][j][k]*V[4][i][j][k]))
                    /(2.0*V[0][i][j][k]));
          
          
          lambmax=max(0, V[1][i][j][k]+fast);
          lambmin=min(0, V[1][i][j][k]-fast);
          
					
					Vl[m][i][j][k]=(V[m][i][j][k]
                          +(0.5-lambmax*oneo3*dt*oneodx[i])
                          *grad
                          *curvaturefactor1);
					Vr[m][i-1][j][k]=(V[m][i][j][k]
                            -(0.5-lambmin*oneo3*dt*oneodx[i-1])
                            *grad
                            *curvaturefactor2);
					
				}
				
				
			}
		}
	}
 
*/
	
	
#pragma omp parallel for private(i, j, k, m, a, b, grad, b2, gmpr, fast, lambmax, lambmin)
	for(m=0; m<9; m++){
		for(i=1; i<ixmax1; i++){		
			for(j=0; j<jymax1; j++){
				for(k=0; k<kzmax1; k++){
					
					
					a=V[m][i][j][k]-V[m][i-1][j][k];
					b=V[m][i+1][j][k]-V[m][i][j][k];
          
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(a>0.0) grad=min(a, b);
						else grad=max(a, b);
					}
					
          
          b2=(V[4][i][j][k]*V[4][i][j][k]
              +V[5][i][j][k]*V[5][i][j][k]
              +V[6][i][j][k]*V[6][i][j][k]);
          
          gmpr=gm*V[7][i][j][k];
          
          fast=sqrt(((b2+gmpr)
                     +sqrt((b2+gmpr)*(b2+gmpr)
                           -4.0*gmpr*V[4][i][j][k]*V[4][i][j][k]))
                    /(2.0*V[0][i][j][k]));
          
          
          lambmax=max(0, V[1][i][j][k]+fast);
          lambmin=min(0, V[1][i][j][k]-fast);
          
					
					Vl[m][i][j][k]=(V[m][i][j][k]
                          +0.5*(1.0-lambmax*dt*oneodx[i])
                          *grad);
					Vr[m][i-1][j][k]=(V[m][i][j][k]
                            -0.5*(1.0+lambmin*dt*oneodx[i-1])
                            *grad);
					
				}
				
				
			}
		}
	}	
	
	/*
	//==================================================
	//	i=0 Vl axis boundary
	//==================================================
	
	XLeftVlAxisBoundary();
	*/
	/*
	//==================================================
	//	i=ixmax1-1 Vr free boundary
	//==================================================
	
	XRightVrFreeBoundary();
	*/
	/*
	//==================================================
	//	Vl, Vr initialize boundary in x
	//==================================================
	 
	XVlrInitializeBoundary();
	*/
	/*
	//==================================================
	//	Vl, Vr free boundary in x
	//==================================================
	
	XVlrFreeBoundary();
	*/
	
	//==================================================
	//	Vl, Vr periodic boundary in x
	//==================================================
	
	XVlrPeriodicBoundary();
	
	
	return 0;

}



//==================================================
//	vanLeer slope limiter in x (2nd order)
//==================================================


int VanLeerLimiterX(void)
{		
	
	int i, j, k;
	int m;
	double a, b, grad;
	
	
	//-----vanLeer limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<i<ix-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, grad)
	for(m=0; m<9; m++){
    
    for(i=1; i<ixmax1; i++){		
      for(j=0; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i-1][j][k];
					b=V[m][i+1][j][k]-V[m][i][j][k];
					
					
					if(a*b<=0.0)	grad=0.0;
					else	grad=a*b/(a+b);
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+grad;
					Vr[m][i-1][j][k]=V[m][i][j][k]-grad;
				
          
				}
			}
		}
    
	}
	
	/*
	//==================================================
	//	Vl, Vr initialize boundary in x
	//==================================================
	
	XVlrInitializeBoundary();
	*/
	
	//==================================================
	//	Vl, Vr free boundary in x
	//==================================================
	
	XVlrFreeBoundary();
	
	/*
	//==================================================
	//	Vl, Vr periodic boundary in x
	//==================================================
	
	XVlrPeriodicBoundary();
	*/
	
	return 0;
	
}



//============================================================
//	MC (Monotonized Central)slope limiter in x (2nd order)
//============================================================


int MCLimiterX(void)
{	
	
	int i, j, k;
	int m;
	double a, b, c, grad;
	
	
	//-----MC limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<i<ix-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, c, grad)
	for(m=0; m<9; m++){
    
    for(i=1; i<ixmax1; i++){		
      for(j=0; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i-1][j][k];
					b=V[m][i+1][j][k]-V[m][i][j][k];
					c=0.25*(V[m][i+1][j][k]-V[m][i-1][j][k]);
					
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(a>0.0) grad=min(a, min(b, c));
						else grad=max(a, max(b, c));
					}
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+grad;
					Vr[m][i-1][j][k]=V[m][i][j][k]-grad;
				
          
				}
			}
		}
    
	}
	
	
	//==================================================
	//	Vl, Vr periodic boundary in x
	//==================================================
	
	XVlrPeriodicBoundary();

	
	return 0;
	
}



//============================================================
//	superbee slope limiter in x (2nd order)
//============================================================


int SuperbeeLimiterX(void)
{	
	
	int i, j, k;
	int m;
	double a, b, grad;
	
	
	//-----superbee limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<i<ix-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, grad)
	for(m=0; m<9; m++){
    
    for(i=1; i<ixmax1; i++){		
      for(j=0; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i-1][j][k];
					b=V[m][i+1][j][k]-V[m][i][j][k];
					
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(fabs(a)>=fabs(b)){
							if(a>0.0)	grad=min(a, 2.0*b);
							else	grad=max(a, 2.0*b);
						}
						else{
							if(a>0.0) grad=min(2.0*a, b);
							else grad=max(2.0*a, b);
						}
					}
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+0.5*grad;
					Vr[m][i-1][j][k]=V[m][i][j][k]-0.5*grad;
				
          
				}
			}
		}
    
	}
	
	
	//==================================================
	//	Vl, Vr periodic boundary in x
	//==================================================
	
	XVlrPeriodicBoundary();

	
	return 0;
	
}






