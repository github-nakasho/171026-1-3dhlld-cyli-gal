

/***********************************************************************
 *
 *	cylindricalgrid.h
 *	
 *	making cylindrical coordinate grid.
 *
 *
 *	i, j, k : loop index of r, phi, z
 *	dxpdxo2, dxmdxo4 : 0.5*(dx[i]+dx[i-1]), 0.25*(dx[i]-dx[i-1])
 *	dxmdxo2 : 0.5*(dz[k]+dz[k-1])
 *	oneminlength=1.0/minlength
 *  r: distance from origin.
 *
 *	
 *	2015 July 22: add krin.
 *	2015 July 20: add absorb_factor calcuration.
 *	2013 Apr. 13: add stretch.
 *	2013 Apr. 05: coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int CylindricalGrid(void)
{

	int i, j, k;
	double oneodvol;
	double oneominlength;
  double r;
	
	
	//-----r, stretch length-----
	
	dx[0]=dx0;
	x[0]=x0-dx0;
  
	for(i=0; i<ixmax1; i++){
		
		
		x[i+1]=x[i]+dx[i];
    dx[i+1]=stretchratio*dx[i];
		
    if(x[i]>lx) dx[i+1]=min(stretchratio*dx[i], dxmax);
    //if(x[i]>lx) dx[i+1]=stretchratio*dx[i];
    else dx[i+1]=dx[i];
		
		
		xm[i]=(x[i]+x[i+1])*0.5;
		dxm[i]=(dx[i]+dx[i+1])*0.5;
		oneox[i]=1.0/x[i];
		oneoxm[i]=1.0/xm[i];
    oneodx[i]=1.0/dx[i];
		
		
	}

	
	xm[ixmax1]=x[ixmax1]+0.5*dx[ixmax1];
	dxm[ixmax1]=dx[ixmax1];
	
	
	oneox[ixmax1]=1.0/x[ixmax1];
  oneodx[ixmax1]=1.0/dx[ixmax1];
	

	//-----phi, no stretch length-----
	
	dy=dy0;
	y[0]=y0-dy;
  oneody=1.0/dy;
	
	for(j=0; j<jymax1; j++) y[j+1]=y[j]+dy;
	

	//-----z, stretch length-----
	
	dz[0]=dz0;
	z[0]=z0-dz0;
	
	for(k=0; k<kzmax1; k++){
		
		z[k+1]=z[k]+dz[k];
		
    
    if(z[k]>lz) dz[k+1]=min(stretchratio*dz[k], dzmax);
    //if(z[k]>lz) dz[k+1]=stretchratio*dz[k];
    else dz[k+1]=dz[k];
		
    
    
		zm[k]=(z[k]+z[k+1])*0.5;
		dzm[k]=(dz[k]+dz[k+1])*0.5;
    oneodz[k]=1.0/dz[k];
		
		
	}
	
  
  zm[kzmax1]=z[kzmax1]+0.5*dz[kzmax1];
  oneodz[kzmax1]=1.0/dz[kzmax1];

  
	//-----set [i][j][k] volume, geometrical factors and so on-----
	
#pragma omp parallel for private(i, j, k, oneodvol)
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			for(k=1; k<kzmax1; k++){
				
				//-----set cell volume of i, j, k & inverse-----
				
				dvolume[i][j][k]=0.5*(xm[i]*xm[i]-xm[i-1]*xm[i-1])*dy*dzm[k-1];
				oneodvol=1.0/dvolume[i][j][k];
				
				
				//-----set geometrical factor-----
				
				g1[i][j][k]=xm[i-1]*dy*dzm[k-1]*oneodvol;
				g2[i][j][k]=xm[i]*dy*dzm[k-1]*oneodvol;
				g3[i][j][k]=dxm[i-1]*dzm[k-1]*oneodvol;
				g4[i][j][k]=g3[i][j][k];
				g5[i][j][k]=0.5*(xm[i]*xm[i]-xm[i-1]*xm[i-1])*dy*oneodvol;
				g6[i][j][k]=g5[i][j][k];
				
				
			}
		}
	}
	
	
  
  //-----searching minimum length of grids-----
  
	oneominminlength=0.0;
	
	
#pragma omp parallel for private(i, j, k, oneominlength) shared(oneominminlength)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
				//-----searching minimum length of [i][j][k] cell-----
				
				minlength[i][j][k]=min(dx[i], min(x[i]*dy, dz[k]));
				
				
				//-----searching minimum length of all cell-----
				
				oneominlength=1.0/minlength[i][j][k];
				
				
#pragma omp flush(oneominminlength)
				if(oneominlength>oneominminlength)	oneominminlength=oneominlength;
				
				
			}
		}
	}

	
	minminlength=1.0/oneominminlength;
  
  
  //-----setting absorbing coefficients-----
	
#pragma omp parallel for private(i)
  for(i=ixmax1-10; i<ixmax; i++){
    
    x_absorb_factor[i-(ixmax1-10)]=(0.5*
                                    (1.0-tanh((x[ixmax1]-x[i]-5.0*dx[ixmax1])
                                              /(2.0*dx[ixmax1]))));
    
    //printf("%d\t%e\n", i, x_absorb_factor[i-(ixmax1-10)]);
  }
	

#pragma omp parallel for private(k)
  for(k=kzmax1-10; k<kzmax; k++){
          
    z_absorb_factor[k-(kzmax1-10)]=(0.5
                                    *(1.0-tanh((z[kzmax1]-z[k]-5.0*dz[kzmax1])
                                               /(2.0*dz[kzmax1]))));
    
    //printf("%d\t%e\n", k, z_absorb_factor[k-(kzmax1-10)]);
  }
  
  //printf("%d\n", irin);
#pragma omp parallel for private(i, k, r)
  for(i=0; i<irin+1; i++){
    for(k=0; k<kzmax; k++){
          
      r=sqrt(x[i]*x[i]+z[k]*z[k]);
          
      c_absorb_factor[i][k]=(0.5
                             *(1.0-tanh((r-rin+5.0*dx0)/(2.0*dx0))));
          
          
    }
  }
  
  
  //-----setting krin-----
#pragma omp parallel for private(i, k, r)
  for(i=0; i<irin+1; i++){
    for(k=0; k<kzmax; k++){
      
      r=sqrt(x[i]*x[i]+z[k]*z[k]);
      
      
      if(r>rin){
      
        krin[i]=k;
      
        
        break;
        
      }
      
    }
  }
  

	return 0;

}



