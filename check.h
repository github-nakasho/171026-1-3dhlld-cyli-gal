

/***********************************************************************
 *
 *	check.h
 *
 *	checking negative density or pressure.
 *
 *
 *	i, j, k : loop index of x, y, z
 *	v2=vx*vx+vy*vy+vz*vz
 *	b2=bx*bx+by*by+bz*bz
 *	
 *
 *	2015 Feb. 20 : add beta check.
 *	2012 Dec. 19 : using dummyU, dummyV
 *	2012 Nov. 06 : improved.
 *	2012 Oct. 08 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 

int Check(double dummyU[][ixmax][jymax][kzmax], 
					double dummyV[][ixmax][jymax][kzmax])
{
	
	int i, j, k;
	double v2, b2;
	
	
	rhoflag=0;
	prflag=0;
    oneobetaflag=0;
	

#pragma omp parallel for private(i, j, k, v2, b2, prfloor)
	for(i=0; i<ixmax; i++){
		for(j=0; j<jymax; j++){
			for(k=0; k<kzmax; k++){
				
                
                v2=dummyV[2][i][j][k]*dummyV[2][i][j][k];
                
                b2=(dummyV[4][i][j][k]*dummyV[4][i][j][k]
                    +dummyV[5][i][j][k]*dummyV[5][i][j][k]
                    +dummyV[6][i][j][k]*dummyV[6][i][j][k]);
                
                
				//-----checking negative density-----
				
				if(dummyV[0][i][j][k]<rhofloor[i][j][k]){
					
					rhoflag++;
				
					/*
					printf("negative density at %d %d %d "
								 "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", 
								 i, j, k, 
								 dummyU[0][i][j][k],
								 dummyU[1][i][j][k], 
								 dummyU[2][i][j][k], 
								 dummyU[3][i][j][k], 
								 dummyV[0][i][j][k], 
								 dummyV[1][i][j][k],
								 dummyV[2][i][j][k],
								 dummyV[3][i][j][k]);
					*/
					
					
					dummyV[0][i][j][k]=rhofloor[i][j][k];
					dummyV[1][i][j][k]=0.0;
					dummyV[3][i][j][k]=0.0;
					
					
					//-----V--->U-----
     

					dummyU[0][i][j][k]=dummyV[0][i][j][k];
					dummyU[1][i][j][k]=0.0;
					dummyU[2][i][j][k]=x[i]*dummyV[0][i][j][k]*dummyV[2][i][j][k];
					dummyU[3][i][j][k]=0.0;
					dummyU[7][i][j][k]=0.5*(dummyV[0][i][j][k]*v2+b2)
															+oneogm1*dummyV[7][i][j][k];
					
					
				}
                
                
                //-----checking negative pressure-----
                
                prfloor=dummyV[0][i][j][k]*prcr;
                
                
                if(dummyV[7][i][j][k]<prfloor){
                    
                    prflag++;
                    
                    /*
                     printf("negative pressure at %d %d %d "
                     "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
                     i, j, k,
                     dummyU[1][i][j][k],
                     dummyU[2][i][j][k],
                     dummyU[3][i][j][k],
                     dummyU[7][i][j][k],
                     dummyV[1][i][j][k],
                     dummyV[2][i][j][k],
                     dummyV[3][i][j][k],
                     dummyV[7][i][j][k]);
                     */
                    
                    //printf("negative pressure at %d %d %d\n", i, j, k);
                    
                    
                    dummyV[1][i][j][k]=0.0;
                    dummyV[3][i][j][k]=0.0;
                    dummyV[7][i][j][k]=prfloor;
                    
                    
                    //-----V--->U-----
                    
                    
                    dummyU[1][i][j][k]=0.0;
                    dummyU[3][i][j][k]=0.0;
                    dummyU[7][i][j][k]=(0.5*(dummyV[0][i][j][k]*v2+b2)
                                        +oneogm1*dummyV[7][i][j][k]);
                    
                    
                }
                
                
                //-----checking beta criteria-----
                
                oneobeta[i][j][k]=0.5*b2/dummyV[7][i][j][k];
                
                if(oneobeta[i][j][k]>oneobetacr){
                    
                    
                    oneobetaflag++;
                    
                    
                    /*
                    printf("lower beta at %d %d %d "
                    "%.3e %.3e %.3e %.3e %.3e\n",
                    i, j, k,
                    dummyU[1][i][j][k],
                    dummyU[4][i][j][k],
                    dummyU[5][i][j][k],
                    dummyU[6][i][j][k],
                    dummyV[7][i][j][k]);
                    */
                    
                    
                    dummyV[0][i][j][k]=betacr*oneobeta[i][j][k]*dummyV[0][i][j][k];
                    dummyV[1][i][j][k]=0.0;
                    dummyV[3][i][j][k]=0.0;
                    dummyV[7][i][j][k]=betacr*oneobeta[i][j][k]*dummyV[7][i][j][k];
                    
                    
                    //-----V--->U-----
                    
                    dummyU[0][i][j][k]=dummyV[0][i][j][k];
                    dummyU[1][i][j][k]=0.0;
                    dummyU[2][i][j][k]=x[i]*dummyV[0][i][j][k]*dummyV[2][i][j][k];
                    dummyU[3][i][j][k]=0.0;
                    dummyU[7][i][j][k]=(0.5*(dummyV[0][i][j][k]*v2+b2)
                                        +oneogm1*dummyV[7][i][j][k]);
                    
                    
                    oneobeta[i][j][k]=oneobetacr;
                    
                    
                }
                
			}
		}
	}
	
	
	printf("negative density point : %d\n"
				 "negative pressure point : %d\n"
                    "critical beta point : %d\n",
				 rhoflag, prflag, oneobetaflag);
	
	
	return 0;
	
}


