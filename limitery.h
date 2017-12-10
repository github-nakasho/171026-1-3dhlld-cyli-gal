

/***********************************************************************
 *
 *	limitery.h
 *
 *	slope limiter in y
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
 *	2012 Dec. 13 : add Vl, Vr boundary condition
 *	2012 Nov. 19 : add VanLeerLimiterX, MCLimiterX, SuperbeeLimiterX.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 31 : add periodic boundary for Vl[][0][], Vr[][jy-1][].
 *	2012 Oct. 30 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



//==================================================
//	minmod slope limiter in y (2nd order)
//==================================================


int MinmodLimiterY(void)
{		

	int i, j, k;
	int m;
	double a, b, grad;
  double b2, gmpr, fast;
  double lambmax, lambmin;
	
	
	//-----minmod limiter (V--->Vl, Vr)-----	
	//-----left & right state of 1<j<jy-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, grad, b2, gmpr, fast, lambmax, lambmin)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){		
      for(j=1; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i][j-1][k];
					b=V[m][i][j+1][k]-V[m][i][j][k];
				
					
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
                           -4.0*gmpr*V[5][i][j][k]*V[5][i][j][k]))
                    /(2.0*V[0][i][j][k]));
          
          /*
           lambmax=max(0, V[1][i][j][k]+fast);
           lambmin=min(0, V[1][i][j][k]-fast);
           */
          
          lambmax=max(0, V[2][i][j][k]+fast);
          lambmin=min(0, V[2][i][j][k]-fast);
          
          
					Vl[m][i][j][k]=(V[m][i][j][k]
                          +0.5*(1.0-lambmax*dt*oneox[i]*oneody)*grad);
					Vr[m][i][j-1][k]=(V[m][i][j][k]
                            -0.5*(1.0+lambmin*dt*oneox[i]*oneody)*grad);
					
          
				}
			}
		}

	}
	
	/*
	//==================================================
	//	Vl, Vr free boundary in y
	//==================================================
	
	YVlrFreeBoundary();
	*/
	
	//==================================================
	//	Vl, Vr periodic boundary in y
	//==================================================
	
	YVlrPeriodicBoundary();
	
	/*
	//==================================================
	//	Vl, Vr initialize boundary in y
	//==================================================
	
	YVlrInitializeBoundary();
	*/
	
	return 0;

}



//==================================================
//	vanLeer slope limiter in y (2nd order)
//==================================================


int VanLeerLimiterY(void)
{
	
	int i, j, k;
	int m;
	double a, b, grad;
	
	
	//-----vanLeer limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<j<jy-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, grad)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=1; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i][j-1][k];
					b=V[m][i][j+1][k]-V[m][i][j][k];
					
					
					if(a*b<=0.0) grad=0.0;
					else grad=a*b/(a+b);
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+grad;
					Vr[m][i][j-1][k]=V[m][i][j][k]-grad;
					
				}
			}
		}
    
	}
	
	/*
	//==================================================
	//	Vl, Vr free boundary in y
	//==================================================
	
	YVlrFreeBoundary();
	*/
	
	//==================================================
	//	Vl, Vr periodic boundary in y
	//==================================================
	
	YVlrPeriodicBoundary();
	
	
	//==================================================
	//	Vl, Vr initialize boundary in y
	//==================================================
	
	YVlrInitializeBoundary();
	
	
	return 0;
	
}



//============================================================
//	MC (Monotonized Central) slope limiter in y (2nd order)
//============================================================


int MCLimiterY(void)
{
	
	int i, j, k;
	int m;
	double a, b, c, grad;
	
	
	//-----MC limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<j<jy-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, c, grad)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=1; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
				
					a=V[m][i][j][k]-V[m][i][j-1][k];
					b=V[m][i][j+1][k]-V[m][i][j][k];
					c=0.25*(V[m][i][j+1][k]-V[m][i][j-1][k]);
					
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(a>0.0) grad=min(a, min(b, c));
						else grad=max(a, max(b, c));
					}
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+grad;
					Vr[m][i][j-1][k]=V[m][i][j][k]-grad;
					
          
				}
			}
		}
    
	}
	
	
	//==================================================
	//	Vl, Vr periodic boundary in y
	//==================================================
	
	YVlrPeriodicBoundary();
	
	
	return 0;
	
}



//==================================================
//	superbee slope limiter in y (2nd order)
//==================================================


int SuperbeeLimiterY(void)
{
	
	int i, j, k;
	int m;
	double a, b, grad;
	
	
	//-----superbee limiter (V--->Vl, Vr)-----
	//-----left & right state of 1<j<jy-1 cells------
	
#pragma omp parallel for private(i, j, k, m, a, b, grad)
	for(m=0; m<9; m++){
    
    for(i=0; i<ixmax1; i++){
      for(j=1; j<jymax1; j++){
        for(k=0; k<kzmax1; k++){
					
					a=V[m][i][j][k]-V[m][i][j-1][k];
					b=V[m][i][j+1][k]-V[m][i][j][k];
					
					
					if(a*b<=0.0) grad=0.0;
					else{
						if(fabs(a)>=fabs(b)){
							if(a>0.0) grad=min(a, 2.0*b);
							else grad=max(a, 2.0*b);
						}
						else{
							if(a>0.0) grad=min(2.0*a, b);
							else grad=max(2.0*a, b);
						}
					}
					
					
					Vl[m][i][j][k]=V[m][i][j][k]+0.5*grad;
					Vr[m][i][j-1][k]=V[m][i][j][k]-0.5*grad;
					
          
				}
			}
		}
    
	}
	
	
	//==================================================
	//	Vl, Vr periodic boundary in y
	//==================================================
	
	YVlrPeriodicBoundary();
	
	
	return 0;
	
}



