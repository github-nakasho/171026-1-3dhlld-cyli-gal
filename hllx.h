

/***********************************************************************
 *
 *	hllx.h
 *
 *	HLL flux in x-direction
 *
 *	
 *	i, j, k : loop index of x, y, z
 *	m : Vl, Ul, Fl, Vr, Ur, Fr vector component index
 *	v2=vx*vx+vy*vy+vz*vz : v**2
 *	b2=bx*bx+by*by+bz*bz : B**2
 *	vdotb=vx*bx+vy*by+vz*bz : dot_product(v, B)
 *	gmpr=gamma*pr
 *	vfl, vfr : left & right side fast mode phase speed
 *	sl, sr : left & right side Riemann fan speed
 *
 *	
 *	2014 Apr. 16 : adapting for cylindrical coordinate.
 *	2012 Dec. 28 : changing flux decision.
 *	2012 Nov. 06 : using pointer & add function.
 *	2012 Nov. 03 : add GLM-MHD divergence cleaning.
 *	2012 Oct. 09 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int HLLX(void)
{	
	
	int i, j, k;
	int m;
	double v2, b2, vdotb;
	double gmpr;
	double vfl, vfr;
	double sl, sr;
	

	//-----HLL flux calculation in x-direction-----	
	
#pragma omp parallel for private(i, j, k, m, v2, b2, vdotb, gmpr, vfl, vfr, sl, sr)
	for(i=0; i<ixmax1; i++){
		for(j=0; j<jymax1; j++){
			for(k=0; k<kzmax1; k++){
				
				
				//-----left side-----
				
				v2=Vl[1][i][j][k]*Vl[1][i][j][k]
						+Vl[2][i][j][k]*Vl[2][i][j][k]
						+Vl[3][i][j][k]*Vl[3][i][j][k];
				
				b2=Vl[4][i][j][k]*Vl[4][i][j][k]
						+Vl[5][i][j][k]*Vl[5][i][j][k]
						+Vl[6][i][j][k]*Vl[6][i][j][k];
				
				vdotb=Vl[1][i][j][k]*Vl[4][i][j][k]
							+Vl[2][i][j][k]*Vl[5][i][j][k]
							+Vl[3][i][j][k]*Vl[6][i][j][k];
				
				
				//-----Vl--->Ul-----
				
				VlrtoUlr(i, j, k, v2, b2, xm[i], Vl, Ul);
				
				
				//-----Vl, Ul--->Fl-----
				
				VlrUlrtoFlr(i, j, k, b2, vdotb, Vl, Ul, Fl);
				
				
				//-----left side fast mode phase speed-----
				
				gmpr=gm*Vl[7][i][j][k];
				
				vfl=sqrt(((b2+gmpr)
									+sqrt((b2+gmpr)*(b2+gmpr)
												-4.0*gmpr*Vl[4][i][j][k]*Vl[4][i][j][k]))
								 /(2.0*Vl[0][i][j][k]));
				
				
				//-----right side-----
				
				v2=Vr[1][i][j][k]*Vr[1][i][j][k]
						+Vr[2][i][j][k]*Vr[2][i][j][k]
						+Vr[3][i][j][k]*Vr[3][i][j][k];
				
				b2=Vr[4][i][j][k]*Vr[4][i][j][k]
						+Vr[5][i][j][k]*Vr[5][i][j][k]
						+Vr[6][i][j][k]*Vr[6][i][j][k];
				
				vdotb=Vr[1][i][j][k]*Vr[4][i][j][k]
								+Vr[2][i][j][k]*Vr[5][i][j][k]
								+Vr[3][i][j][k]*Vr[6][i][j][k];
				
				
				//-----Vr--->Ur-----
				
				VlrtoUlr(i, j, k, v2, b2, xm[i], Vr, Ur);
				
				
				//-----Vr, Ur--->Fr-----
				
				VlrUlrtoFlr(i, j, k, b2, vdotb, Vr, Ur, Fr);
				
				
				//-----right side fast mode phase speed-----
				
				gmpr=gm*Vr[7][i][j][k];
				
				vfr=sqrt(((b2+gmpr)
									+sqrt((b2+gmpr)*(b2+gmpr)
												-4.0*gmpr*Vr[4][i][j][k]*Vr[4][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[1][i][j][k], Vr[1][i][j][k])-max(vfl, vfr);
				sr=max(Vl[1][i][j][k], Vr[1][i][j][k])+max(vfl, vfr);
				
				
				if(sl>=0.0){
					
					//-----F=Fl-----
					
					for(m=0; m<8; m++)	F[m][i][j][k]=Fl[m][i][j][k];
					
					
				}
				else{
					
					if(sr>=0.0){
					
						//-----F=Fhll-----
						
						for(m=0; m<8; m++){
							
							F[m][i][j][k]=(sr*Fl[m][i][j][k]-sl*Fr[m][i][j][k]
														 +sr*sl*(Ur[m][i][j][k]-Ul[m][i][j][k]))
														/(sr-sl);
							
							
						}
						
						
					}
					else{
						
						//-----F=Fr-----
					
						for(m=0; m<8; m++)	F[m][i][j][k]=Fr[m][i][j][k];
					
					
					}
				}
				
				
			}
		}
	}
	
	
	//-----GLM-MHD flux-----
		
#pragma omp parallel for private(i, j, k)
	for(i=0; i<ixmax1; i++){
		for(j=0; j<jymax1; j++){
			for(k=0; k<kzmax1; k++){
				
				F[4][i][j][k]+=0.5*(Vl[8][i][j][k]+Vr[8][i][j][k]
														+ch*(Vl[4][i][j][k]-Vr[4][i][j][k]));
				F[8][i][j][k]=0.5*ch*(ch*(Vl[4][i][j][k]+Vr[4][i][j][k])
															+(Vl[8][i][j][k]-Vr[8][i][j][k]));
			
			}
		}
	}

	
	return 0;

}























