

/***********************************************************************
 *
 *	hlldy.h
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
 *	2013 Apr. 13 : HLLD flux for phi direction in Cylindrical cood.
 *	2013 Mar. 25 : improved HLLD flux
 *	2012 Dec. 20 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/



int HLLDY(void)
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
	double byint;
	double vzll, vzrr, vxll, vxrr, vzint, vxint;
	double bzll, bzrr, bxll, bxrr, bzint, bxint;
	double enll, enrr, enlll, enrrr;
	double sqrtrholl, sqrtrhorr;
	double sgn, oneosqrtrho;
	double denominator, numerator;
	
	
	//-----HLLD flux calculation in y-direction-----	
	
#pragma omp parallel for private(i, j, k, m, rholl, rhorr, enll, enrr, enlll, enrrr, vzll, vzrr, bzll, bzrr, vxll, vxrr, bxll, bxrr, byint, vzint, bzint, vxint, bxint, vl2, bl2, vr2, br2, vrdotbr, vldotbl, vlldotbll, vrrdotbrr, vintdotbint, gmpr, ptl, ptr, vfl, vfr, sl, sll, sr, srr, slvl, srvr, sm, oneoslsm, oneosrsm, ptint, sgn, sqrtrholl, sqrtrhorr, oneosqrtrho, denominator, numerator)
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
												-4.0*gmpr*Vl[5][i][j][k]*Vl[5][i][j][k]))
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
												-4.0*gmpr*Vr[5][i][j][k]*Vr[5][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[2][i][j][k], Vr[2][i][j][k])-max(vfl, vfr);
				sr=max(Vl[2][i][j][k], Vr[2][i][j][k])+max(vfl, vfr);
				
				
				if(sl>=0.0){
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, x[i], Vl, Ul);
					
					
					//-----Vl, Ul--->G=Gl-----
					
					VlrUlrtoGlr(i, j, k, bl2, vldotbl, Vl, Ul, G);
					
					
				}
				else if(sr<=0.0){
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, x[i], Vr, Ur);
					
					
					//-----Vr, Ur--->G=Gr-----
					
					VlrUlrtoGlr(i, j, k, br2, vrdotbr, Vr, Ur, G);
				
				
				}
				else{
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, x[i], Vl, Ul);
					
					
					//-----Vl, Ul--->Gl-----
					
					//VlrUlrtoGlr(i, j, k, bl2, vldotbl, Vl, Ul, Gl);
					
					
					//-----Vr--->Ur-----
						
					VlrtoUlr(i, j, k, vr2, br2, x[i], Vr, Ur);
					
					
					//-----Vr, Ur--->Gr-----
						
					//VlrUlrtoGlr(i, j, k, br2, vrdotbr, Vr, Ur, Gr);
					
					
					//-----entropy speed, vxll=vxlll=vxrrr=vxrr=sm-----
					
					slvl=sl-Vl[2][i][j][k];
					srvr=sr-Vr[2][i][j][k];
					
					sm=((srvr*Vr[0][i][j][k]*Vr[2][i][j][k]
							 -slvl*Vl[0][i][j][k]*Vl[2][i][j][k])
							-ptr+ptl)
							/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);

					oneoslsm=1.0/(sl-sm);
					oneosrsm=1.0/(sr-sm);
					
						
					//-----intermediate total pressure, 
					//											ptll=ptlll=ptrrr=ptrr=ptint----- 
					
					ptint=(srvr*Ur[0][i][j][k]*ptl-slvl*Ul[0][i][j][k]*ptr
								 +Ur[0][i][j][k]*Ul[0][i][j][k]*srvr*slvl
								 *(Vr[2][i][j][k]-Vl[2][i][j][k]))
								/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);
					
					
					//-----intermediate by=by_hll(Gr[5]=0, Gl[5]=0)-----
					/*
					byint=(sr*Ur[5][i][j][k]-sl*Ul[5][i][j][k]
								 -Gr[5][i][j][k]+Gl[5][i][j][k])
								/(sr-sl);
					*/
					
					byint=(sr*Ur[5][i][j][k]-sl*Ul[5][i][j][k])/(sr-sl);
					
					
					//-----left left left & right right right side density-----
					//-----rholl=rholll, rhorr=rhorrr-----
					
					rholl=Ul[0][i][j][k]*slvl*oneoslsm;
					rhorr=Ur[0][i][j][k]*srvr*oneosrsm;
					
					
          //-----left left side vz, vx, bz, bx, vlldotbll-----
					
					//-----caution! if Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vzll=vzl, vxll=vxl, bzll=bzl, bxll=bxl-----
          					
          sgn=(fabs(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint)
                    -1.0e-8*byint*byint)
               /(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint)
                 -1.0e-8*byint*byint));
           
          denominator=((1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint)
                        +1.0)
                       +(1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint)
                         -1.0)*sgn)*0.5;
           
          numerator=((((Ul[0][i][j][k]*slvl*slvl-byint*byint)+1.0)
                      +((Ul[0][i][j][k]*slvl*slvl-byint*byint)-1.0)*sgn)
                     *0.5);
           
           
          vzll=Vl[3][i][j][k]-(byint*Ul[6][i][j][k]
                               *(sm-Vl[2][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          vxll=Vl[1][i][j][k]-(byint*Ul[4][i][j][k]
                               *(sm-Vl[2][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          bzll=Ul[6][i][j][k]*numerator*denominator;
           
          bxll=Ul[4][i][j][k]*numerator*denominator;
           				
          
					
          /*
					//-----left left side vy, vz, by, bz, vlldotbll-----
					
					//caution! if denominator==0.0, degenerate HLLD & HLL
					
					if(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint)
						 >1.0e-8*byint*byint){
						
						denominator=1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-byint*byint);
						numerator=Ul[0][i][j][k]*slvl*slvl-byint*byint;

						
						vzll=Vl[3][i][j][k]-(byint*Ul[6][i][j][k]
																 *(sm-Vl[2][i][j][k])*denominator);
						
						vxll=Vl[1][i][j][k]-(byint*Ul[4][i][j][k]
																 *(sm-Vl[2][i][j][k])*denominator);
						
						bzll=Ul[6][i][j][k]*numerator*denominator;
						
						bxll=Ul[4][i][j][k]*numerator*denominator;
						
						
					}
					else{
						
						vzll=Vl[3][i][j][k];
						
						vxll=Vl[1][i][j][k];
						
						bzll=Ul[6][i][j][k];
						
						bxll=Ul[4][i][j][k];
						
						
					}
					*/
					
          
					vlldotbll=sm*byint+vzll*bzll+vxll*bxll;
					
          
					//-----right right side vz, vx, bz, bx, vrrdotbrr-----
					
					//-----caution! if Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vzrr=vzr, vxrr=vxr, bzrr=bzr, bxrr=bxr-----
          				
           sgn=(fabs(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint)
           -1.0e-8*byint*byint)
           /(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint)
           -1.0e-8*byint*byint));
           
           denominator=((1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint)
           +1.0)
           +(1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint)
           -1.0)*sgn)*0.5;
           
           numerator=((((Ur[0][i][j][k]*srvr*srvr-byint*byint)+1.0)
           +((Ur[0][i][j][k]*srvr*srvr-byint*byint)-1.0)*sgn)
           *0.5);
           
           
           vzrr=Vr[3][i][j][k]-(byint*Ur[6][i][j][k]
           *(sm-Vr[2][i][j][k])
           *((denominator+denominator*sgn)*0.5));
           
           vxrr=Vr[1][i][j][k]-(byint*Ur[4][i][j][k]
           *(sm-Vr[2][i][j][k])
           *((denominator+denominator*sgn)*0.5));
           
           bzrr=Ur[6][i][j][k]*numerator*denominator;
           
           bxrr=Ur[4][i][j][k]*numerator*denominator;
                
					
/*
					//-----right right side vy, vz, by, bz, vrrdotbrr-----
					
					//caution! if denominator==00, degenerate HLLD & HLL
					
					if(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint)
						 >1.0e-8*byint*byint){
						
						denominator=1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-byint*byint);
						numerator=Ur[0][i][j][k]*srvr*srvr-byint*byint;
					
						
						vzrr=Vr[3][i][j][k]-(byint*Vr[6][i][j][k]
																 *(sm-Vr[2][i][j][k])*denominator);
						
						vxrr=Vr[1][i][j][k]-(byint*Vr[4][i][j][k]
																 *(sm-Vr[2][i][j][k])*denominator);
						
						bzrr=Ur[6][i][j][k]*numerator*denominator;
						
						bxrr=Ur[4][i][j][k]*numerator*denominator;
						
						
					}
					else{	
						
						vzrr=Vr[3][i][j][k];
						
						vxrr=Vr[1][i][j][k];
						
						bzrr=Ur[6][i][j][k];
						
						bxrr=Ur[4][i][j][k];
						
						
					}
*/					
					
					vrrdotbrr=sm*byint+vzrr*bzrr+vxrr*bxrr;
					
					
					//-----left left side en-----
					
					enll=(slvl*Ul[7][i][j][k]-ptl*Vl[2][i][j][k]+ptint*sm
								+byint*(vldotbl-vlldotbll))*oneoslsm;

					
					//-----right right side en-----
					
					enrr=(srvr*Ur[7][i][j][k]-ptr*Vr[2][i][j][k]+ptint*sm
								+byint*(vrdotbr-vrrdotbrr))*oneosrsm;
					
					
					//-----left & right direction alfven speed-----
					
					sqrtrholl=sqrt(rholl);
					sqrtrhorr=sqrt(rhorr);
					
					sll=sm-fabs(byint)/sqrtrholl;
					srr=sm+fabs(byint)/sqrtrhorr;
					

					if(sm>=0.0){
						
						if(sll>=0.0){
						
							//-----G=Ghlll-----
						
							PrimitivetoG(i, j, k, vlldotbll, rholl, vxll, sm, vzll, 
													 bxll, byint, bzll, enll, ptint, G);
						
							
						}
						else{
							
							sgn=fabs(byint)/byint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
							
							//-----G=Ghlldl-----
							
							vzint=(sqrtrholl*vzll+sqrtrhorr*vzrr+(bzrr-bzll)*sgn)
										*oneosqrtrho;
							vxint=(sqrtrholl*vxll+sqrtrhorr*vxrr+(bxrr-bxll)*sgn)
										*oneosqrtrho;
							bzint=(sqrtrholl*bzrr+sqrtrhorr*bzll
										 +sqrtrholl*sqrtrhorr*(vzrr-vzll)*sgn)
										*oneosqrtrho;
							bxint=(sqrtrholl*bxrr+sqrtrhorr*bxll
										 +sqrtrholl*sqrtrhorr*(vxrr-vxll)*sgn)
										*oneosqrtrho;
							
							
							vintdotbint=sm*byint+vzint*bzint+vxint*bxint;
							
							
							enlll=enll-sqrtrholl*(vlldotbll-vintdotbint)*sgn;
							
							
							PrimitivetoG(i, j, k, vintdotbint, rholl, vxint, sm, vzint, 
													 bxint, byint, bzint, enlll, ptint, G);
							
							
						}
					}
					else{
						
						if(srr<=0.0){
							
							//-----G=Ghllr----
							
							PrimitivetoG(i, j, k, vrrdotbrr, rhorr, vxrr, sm, vzrr, 
													 bxrr, byint, bzrr, enrr, ptint, G);
						
							
						}
						else{
						
							sgn=fabs(byint)/byint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
						
							
							//-----G=Ghlldr-----
								
							vzint=(sqrtrholl*vzll+sqrtrhorr*vzrr+(bzrr-bzll)*sgn)
										*oneosqrtrho;
							vxint=(sqrtrholl*vxll+sqrtrhorr*vxrr+(bxrr-bxll)*sgn)
										*oneosqrtrho;
							bzint=(sqrtrholl*bzrr+sqrtrhorr*bzll
										 +sqrtrholl*sqrtrhorr*(vzrr-vzll)*sgn)
										*oneosqrtrho;
							bxint=(sqrtrholl*bxrr+sqrtrhorr*bxll
										 +sqrtrholl*sqrtrhorr*(vxrr-vxll)*sgn)
										*oneosqrtrho;
								
								
							vintdotbint=sm*byint+vzint*bzint+vxint*bxint;
								
								
							enrrr=enrr+sqrtrhorr*(vrrdotbrr-vintdotbint)*sgn;
								
								
							PrimitivetoG(i, j, k, vintdotbint, rhorr, vxint, sm, vzint, 
													 bxint, byint, bzint, enrrr, ptint, G);
								
							
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
				
				G[5][i][j][k]+=0.5*(Vl[8][i][j][k]+Vr[8][i][j][k]
														+ch*(Vl[5][i][j][k]-Vr[5][i][j][k]));
				G[8][i][j][k]=0.5*ch*(ch*(Vl[5][i][j][k]+Vr[5][i][j][k])
															+(Vl[8][i][j][k]-Vr[8][i][j][k]));
				
				
			}
		}
	}
	
	
	
	return 0;

}









