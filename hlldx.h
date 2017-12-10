

/***********************************************************************
 *
 *	hlldx.h
 *
 *	HLLD flux in x-direction
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
 *	ptl, ptr : left & right side total pressure
 *	ptint : intermediate total pressure
 *
 *
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int HLLDX(void)
{	

	int i, j, k, m;
	double vl2, bl2, vldotbl;
	double vr2, br2, vrdotbr;
	double vlldotbll, vrrdotbrr;
	double vintdotbint;
	double gmpr;
	double vfl, vfr;
	double slvl, srvr, oneoslsm, oneosrsm;
	double sl, sll, sm, srr, sr;
	double ptl, ptr, ptint;
	double rholl, rhorr;
	double bxint;
	double vyll, vyrr, vzll, vzrr, vyint, vzint;
	double byll, byrr, bzll, bzrr, byint, bzint;
	double enll, enrr, enlll, enrrr;
	double sqrtrholl, sqrtrhorr;
	double sgn, oneosqrtrho;
	double denominator, numerator;
	
	
	//-----HLLD flux calculation in x-direction-----	
	
#pragma omp parallel for private(i, j, k, m, rholl, rhorr, enll, enrr, enlll, enrrr, vyll, vyrr, byll, byrr, vzll, vzrr, bzll, bzrr, bxint, vyint, byint, vzint, bzint, vl2, bl2, vr2, br2, vrdotbr, vldotbl, vlldotbll, vrrdotbrr, vintdotbint, gmpr, ptl, ptr, vfl, vfr, sl, sll, sr, srr, slvl, srvr, sm, oneoslsm, oneosrsm, ptint, sgn, sqrtrholl, sqrtrhorr, oneosqrtrho, denominator, numerator)
	for(i=0; i<ixmax1; i++){
		for(j=0; j<jymax1; j++){
			for(k=0; k<kzmax1; k++){
	
				//------left side v**2, b**2, scalor_product(v, b)-----
				
				vl2=Vl[1][i][j][k]*Vl[1][i][j][k]
						+Vl[2][i][j][k]*Vl[2][i][j][k]
						+Vl[3][i][j][k]*Vl[3][i][j][k];
				
				bl2=Vl[4][i][j][k]*Vl[4][i][j][k]
						+Vl[5][i][j][k]*Vl[5][i][j][k]
						+Vl[6][i][j][k]*Vl[6][i][j][k];
				
				vldotbl=Vl[1][i][j][k]*Vl[4][i][j][k]
								+Vl[2][i][j][k]*Vl[5][i][j][k]
								+Vl[3][i][j][k]*Vl[6][i][j][k];
				
			
				//-----left side gamma*pressure-----
				
				gmpr=gm*Vl[7][i][j][k];
				
				
				//----left side total pressure-----
				
				ptl=Vl[7][i][j][k]+0.5*bl2;
				
				
				//-----left side fast mode phase speed-----
				
				vfl=sqrt(((bl2+gmpr)
									+sqrt((bl2+gmpr)*(bl2+gmpr)
												-4.0*gmpr*Vl[4][i][j][k]*Vl[4][i][j][k]))
								 /(2.0*Vl[0][i][j][k]));

				
				//-----right side-----
				
				vr2=Vr[1][i][j][k]*Vr[1][i][j][k]
						+Vr[2][i][j][k]*Vr[2][i][j][k]
						+Vr[3][i][j][k]*Vr[3][i][j][k];
				
				br2=Vr[4][i][j][k]*Vr[4][i][j][k]
						+Vr[5][i][j][k]*Vr[5][i][j][k]
						+Vr[6][i][j][k]*Vr[6][i][j][k];
				
				vrdotbr=Vr[1][i][j][k]*Vr[4][i][j][k]
								+Vr[2][i][j][k]*Vr[5][i][j][k]
								+Vr[3][i][j][k]*Vr[6][i][j][k];
				
				
				//-----right side gamma*pressure-----
				
				gmpr=gm*Vr[7][i][j][k];
				
				
				//----right side total pressure-----
				
				ptr=Vr[7][i][j][k]+0.5*br2;
				
				
				//-----right side fast mode phase speed-----
				
				vfr=sqrt(((br2+gmpr)
									+sqrt((br2+gmpr)*(br2+gmpr)
												-4.0*gmpr*Vr[4][i][j][k]*Vr[4][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[1][i][j][k], Vr[1][i][j][k])-max(vfl, vfr);
				sr=max(Vl[1][i][j][k], Vr[1][i][j][k])+max(vfl, vfr);
				
				
				if(sl>=0.0){
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, xm[i], Vl, Ul);
					
					
					//-----Vl, Ul--->F=Fhlll-----
					
					VlrUlrtoFlr(i, j, k, bl2, vldotbl, Vl, Ul, F);
					
					
				}
				else if(sr<=0.0){
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, xm[i], Vr, Ur);
					
					
					//-----Vr, Ur--->F=Fhllr-----
					
					VlrUlrtoFlr(i, j, k, br2, vrdotbr, Vr, Ur, F);
					
					
				}
				else{
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, xm[i], Vl, Ul);
					
					
					//-----Vl, Ul--->Fl-----
					
					//VlrUlrtoFlr(i, j, k, bl2, vldotbl, Vl, Ul, Fl);
					
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, xm[i], Vr, Ur);
					
					
					//-----Vr, Ur--->Fr-----
					
					//VlrUlrtoFlr(i, j, k, br2, vrdotbr, Vr, Ur, Fr);
					
					
					//-----entropy speed, vxll=vxlll=vxrrr=vxrr=sm-----
					
					slvl=sl-Vl[1][i][j][k];
					srvr=sr-Vr[1][i][j][k];
					
					sm=(srvr*Ur[1][i][j][k]-slvl*Ul[1][i][j][k]-ptr+ptl)
							/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);

					oneoslsm=1.0/(sl-sm);
					oneosrsm=1.0/(sr-sm);
					
					
					//-----intermediate total pressure, 
					//											ptll=ptlll=ptrrr=ptrr=ptint----- 
					
					ptint=(srvr*Ur[0][i][j][k]*ptl-slvl*Ul[0][i][j][k]*ptr
								 +Ur[0][i][j][k]*Ul[0][i][j][k]*srvr*slvl
								 *(Vr[1][i][j][k]-Vl[1][i][j][k]))
								/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);
					
					
					//-----intermediate bx=bx_hll(Fr[4]=0, Fl[4]=0)-----
					/*
					bxint=(sr*Ur[4][i][j][k]-sl*Ul[4][i][j][k]
								 -Fr[4][i][j][k]+Fl[4][i][j][k])
								/(sr-sl);
					*/
					
					bxint=(sr*Ur[4][i][j][k]-sl*Ul[4][i][j][k])/(sr-sl);
					

					//-----left left left & right right right side density-----
					//-----rholl=rholll, rhorr=rhorrr-----
					
					rholl=Ul[0][i][j][k]*slvl*oneoslsm;
					rhorr=Ur[0][i][j][k]*srvr*oneosrsm;
					
					
					//-----left left side vy, vz, by, bz, vlldotbll-----
          
          //-----caution! if Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vyll=vyl, vzll=vzl, byll=byl, bzll=bzl-----
          					
           sgn=(fabs(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint)
                     -1.0e-8*bxint*bxint)
                /(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint)
                  -1.0e-8*bxint*bxint));
           
           denominator=((1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint)
                         +1.0)
                        +(1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint)
                          -1.0)*sgn)*0.5;
           
           numerator=((((Ul[0][i][j][k]*slvl*slvl-bxint*bxint)+1.0)
                       +((Ul[0][i][j][k]*slvl*slvl-bxint*bxint)-1.0)*sgn)
                      *0.5);
           
           
           vyll=Vl[2][i][j][k]-(bxint*Ul[5][i][j][k]
                                *(sm-Vl[1][i][j][k])
                                *((denominator+denominator*sgn)*0.5));
           
           vzll=Vl[3][i][j][k]-(bxint*Ul[6][i][j][k]
                                *(sm-Vl[1][i][j][k])
                                *((denominator+denominator*sgn)*0.5));
           
           byll=Ul[5][i][j][k]*numerator*denominator;
           
           bzll=Ul[6][i][j][k]*numerator*denominator;
           				

/*          
					//caution! if denominator==0.0, degenerate HLLD & HLL
					
					if(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint)
						 >1e-8*bxint*bxint){
					
						denominator=1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint);
						numerator=Ul[0][i][j][k]*slvl*slvl-bxint*bxint;

						
						vyll=Vl[2][i][j][k]-(bxint*Ul[5][i][j][k]
																 *(sm-Vl[1][i][j][k])*denominator);
						
						vzll=Vl[3][i][j][k]-(bxint*Ul[6][i][j][k]
																 *(sm-Vl[1][i][j][k])*denominator);
						
						byll=Ul[5][i][j][k]*numerator*denominator;
						
						bzll=Ul[6][i][j][k]*numerator*denominator;
					
						
					}
					else{
						
						vyll=Vl[2][i][j][k];
						
						vzll=Vl[3][i][j][k];
						
						byll=Ul[5][i][j][k];
						
						bzll=Ul[6][i][j][k];
						
						
					}
*/
					
					vlldotbll=sm*bxint+vyll*byll+vzll*bzll;
					

          //-----right right side vy, vz, by, bz, vrrdotbrr-----
					
					//-----caution! if Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vyrr=vyr, vzrr=vzr, byrr=byr, bzrr=bzr-----
          					
           sgn=(fabs(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint)
           -1.0e-8*bxint*bxint)
           /(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint)
           -1.0e-8*bxint*bxint));
           
           denominator=((1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint)
           +1.0)
           +(1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint)
           -1.0)*sgn)*0.5;
           
           numerator=((((Ur[0][i][j][k]*srvr*srvr-bxint*bxint)+1.0)
           +((Ur[0][i][j][k]*srvr*srvr-bxint*bxint)-1.0)*sgn)
           *0.5);
           
           
           vyrr=Vr[2][i][j][k]-(bxint*Ur[5][i][j][k]
           *(sm-Vr[1][i][j][k])
           *((denominator+denominator*sgn)*0.5));
           
           vzrr=Vr[3][i][j][k]-(bxint*Ur[6][i][j][k]
           *(sm-Vr[1][i][j][k])
           *((denominator+denominator*sgn)*0.5));
           
           byrr=Ur[5][i][j][k]*numerator*denominator;
           
           bzrr=Ur[6][i][j][k]*numerator*denominator;
           		
          /*
					//-----right right side vy, vz, by, bz, vrrdotbrr-----
					
					//caution! if denominator==0.0, degenerate HLLD & HLL
					
					if(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint)
						 >1.0e-8*bxint*bxint){
					
						denominator=1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bxint*bxint);
						numerator=Ur[0][i][j][k]*srvr*srvr-bxint*bxint;

						
						vyrr=Vr[2][i][j][k]-(bxint*Vr[5][i][j][k]
																 *(sm-Vr[1][i][j][k])*denominator);
						
						vzrr=Vr[3][i][j][k]-(bxint*Vr[6][i][j][k]
																 *(sm-Vr[1][i][j][k])*denominator);
						
						byrr=Ur[5][i][j][k]*numerator*denominator;
						
						bzrr=Ur[6][i][j][k]*numerator*denominator;
						
						
					}
					else{
					
						vyrr=Vr[2][i][j][k];
						
						vzrr=Vr[3][i][j][k];
						
						byrr=Ur[5][i][j][k];
						
						bzrr=Ur[6][i][j][k];
						
						
					}
					*/
					
					vrrdotbrr=sm*bxint+vyrr*byrr+vzrr*bzrr;
					
					
					//-----left left side en-----
					
					enll=(slvl*Ul[7][i][j][k]-ptl*Vl[1][i][j][k]+ptint*sm
								+bxint*(vldotbl-vlldotbll))*oneoslsm;

					
					//-----right right side en-----
					
					enrr=(srvr*Ur[7][i][j][k]-ptr*Vr[1][i][j][k]+ptint*sm
								+bxint*(vrdotbr-vrrdotbrr))*oneosrsm;
					
					
					//-----left & right direction alfven speed-----
					
					sqrtrholl=sqrt(rholl);
					sqrtrhorr=sqrt(rhorr);
					
					sll=sm-fabs(bxint)/sqrtrholl;
					srr=sm+fabs(bxint)/sqrtrhorr;
					

					if(sm>=0.0){
						
						if(sll>=0.0){
						
							//-----F=Fhll_l-----
						
							PrimitivetoF(i, j, k, vlldotbll, rholl, sm, vyll, vzll, 
													 bxint, byll, bzll, enll, ptint, F);
						
							
						}
						else{
							
							sgn=fabs(bxint)/bxint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
							
							//-----F=Fhlldl-----
							
							vyint=(sqrtrholl*vyll+sqrtrhorr*vyrr+(byrr-byll)*sgn)
										*oneosqrtrho;
							vzint=(sqrtrholl*vzll+sqrtrhorr*vzrr+(bzrr-bzll)*sgn)
										*oneosqrtrho;
							byint=(sqrtrholl*byrr+sqrtrhorr*byll
										 +sqrtrholl*sqrtrhorr*(vyrr-vyll)*sgn)
										*oneosqrtrho;
							bzint=(sqrtrholl*bzrr+sqrtrhorr*bzll
										 +sqrtrholl*sqrtrhorr*(vzrr-vzll)*sgn)
										*oneosqrtrho;
							
							
							vintdotbint=sm*bxint+vyint*byint+vzint*bzint;
							
							
							enlll=enll-sqrtrholl*(vlldotbll-vintdotbint)*sgn;
							
							
							PrimitivetoF(i, j, k, vintdotbint, rholl, sm, vyint, vzint, 
													 bxint, byint, bzint, enlll, ptint, F);
							
							
						}
						
					}
					else{
						
						if(srr<=0.0){
							
							//-----F=Fhll_r----
							
							PrimitivetoF(i, j, k, vrrdotbrr, rhorr, sm, vyrr, vzrr, 
													 bxint, byrr, bzrr, enrr, ptint, F);
						
							
						}
						else{
						
							sgn=fabs(bxint)/bxint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
						
							
							//-----F=Fhlldr-----
								
							vyint=(sqrtrholl*vyll+sqrtrhorr*vyrr+(byrr-byll)*sgn)
										*oneosqrtrho;
							vzint=(sqrtrholl*vzll+sqrtrhorr*vzrr+(bzrr-bzll)*sgn)
										*oneosqrtrho;
							byint=(sqrtrholl*byrr+sqrtrhorr*byll
											+sqrtrholl*sqrtrhorr*(vyrr-vyll)*sgn)
										*oneosqrtrho;
							bzint=(sqrtrholl*bzrr+sqrtrhorr*bzll
											+sqrtrholl*sqrtrhorr*(vzrr-vzll)*sgn)
										*oneosqrtrho;
								
								
							vintdotbint=sm*bxint+vyint*byint+vzint*bzint;
								
								
							enrrr=enrr+sqrtrhorr*(vrrdotbrr-vintdotbint)*sgn;
								
								
							PrimitivetoF(i, j, k, vintdotbint, rhorr, sm, vyint, vzint, 
													 bxint, byint, bzint, enrrr, ptint, F);
								
							
						}
							
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









