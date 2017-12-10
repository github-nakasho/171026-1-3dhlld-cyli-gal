

/***********************************************************************
 *
 *	hllz.h
 *
 *	HLL flux in z-direction
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
 *	2012 Dec. 28 : change flux decision.
 *	2012 Dec. 03 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int HLLZ(void)
{	
	
	int i, j, k;
	int m;
	double v2, b2, vdotb;
	double gmpr;
	double vfl, vfr;
	double sl, sr;
	
	
	//-----HLL flux calculation in z-direction-----	
	
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
				
				VlrtoUlr(i, j, k, v2, b2, x[i], Vl, Ul);
				
				
				//-----Vl, Ul--->Hl-----
				
				VlrUlrtoHlr(i, j, k, b2, vdotb, Vl, Ul, Hl);
				
			
				//-----left side fast mode phase speed-----
				
				gmpr=gm*Vl[7][i][j][k];
				
				vfl=sqrt(((b2+gmpr)
									+sqrt((b2+gmpr)*(b2+gmpr)
												-4.0*gmpr*Vl[6][i][j][k]*Vl[6][i][j][k]))
								 /(2.0*Vl[0][i][j][k]));
				
				
				//-----right side----
				
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
				
				VlrtoUlr(i, j, k, v2, b2, x[i], Vr, Ur);
				
				
				//-----Vr, Ur--->Hr-----
				
				VlrUlrtoHlr(i, j, k, b2, vdotb, Vr, Ur, Hr);
				
				
				//-----right side fast mode phase speed-----
				
				gmpr=gm*Vr[7][i][j][k];
				
				vfr=sqrt(((b2+gmpr)
									+sqrt((b2+gmpr)*(b2+gmpr)
												-4.0*gmpr*Vr[6][i][j][k]*Vr[6][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[3][i][j][k], Vr[3][i][j][k])-max(vfl, vfr);
				sr=max(Vl[3][i][j][k], Vr[3][i][j][k])+max(vfl, vfr);
				
				
				if(sl>=0.0){
					
					//-----H=Hl-----
					
					for(m=0; m<8; m++)	H[m][i][j][k]=Hl[m][i][j][k];
				
				
				}
				else{
					
					if(sr>=0.0){
					
						//-----H=Hhll-----
						
						for(m=0; m<8; m++){
							
							H[m][i][j][k]=(sr*Hl[m][i][j][k]-sl*Hr[m][i][j][k]
														 +sr*sl*(Ur[m][i][j][k]-Ul[m][i][j][k]))
														/(sr-sl);
							
							
						}
						
						
					}
					else{
						
						//-----H=Hr-----
					
						for(m=0; m<8; m++)	H[m][i][j][k]=Hr[m][i][j][k];
				
					
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
				
				H[6][i][j][k]+=0.5*(Vl[8][i][j][k]+Vr[8][i][j][k]
														+ch*(Vl[6][i][j][k]-Vr[6][i][j][k]));
				H[8][i][j][k]=0.5*ch*(ch*(Vl[6][i][j][k]+Vr[6][i][j][k])
															+(Vl[8][i][j][k]-Vr[8][i][j][k]));
				
			}
		}
	}

	
	return 0;

}












