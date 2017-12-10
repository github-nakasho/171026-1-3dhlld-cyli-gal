

/***********************************************************************
 *
 *	hlldx.h
 *
 *	HLLD flux in z-direction
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



int HLLDZ(void)
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
	double bzint;
	double vxll, vxrr, vyll, vyrr, vxint, vyint;
	double bxll, bxrr, byll, byrr, bxint, byint;
	double enll, enrr, enlll, enrrr;
	double sqrtrholl, sqrtrhorr;
	double sgn, oneosqrtrho;
	double denominator, numerator;
	
	
	//-----HLLD flux calculation in y-direction-----	
	
#pragma omp parallel for private(i, j, k, m, rholl, rhorr, enll, enrr, enlll, enrrr, vxll, vxrr, bxll, bxrr, vyll, vyrr, byll, byrr, bzint, vxint, bxint, vyint, byint, vl2, bl2, vr2, br2, vrdotbr, vldotbl, vlldotbll, vrrdotbrr, vintdotbint, gmpr, ptl, ptr, vfl, vfr, sl, sll, sr, srr, slvl, srvr, sm, oneoslsm, oneosrsm, ptint, sgn, sqrtrholl, sqrtrhorr, oneosqrtrho, denominator, numerator)
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
												-4.0*gmpr*Vl[6][i][j][k]*Vl[6][i][j][k]))
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
												-4.0*gmpr*Vr[6][i][j][k]*Vr[6][i][j][k]))
								 /(2.0*Vr[0][i][j][k]));

				
				//-----Riemann fan speed-----
				
				sl=min(Vl[3][i][j][k], Vr[3][i][j][k])-max(vfl, vfr);
				sr=max(Vl[3][i][j][k], Vr[3][i][j][k])+max(vfl, vfr);
				
				
				if(sl>=0.0){
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, x[i], Vl, Ul);
					
					
					//-----Vl, Ul--->H=Hl-----
					
					VlrUlrtoHlr(i, j, k, bl2, vldotbl, Vl, Ul, H);
					
					
				}
				else if(sr<=0.0){
					
					//-----Vr--->Ur-----
					
					VlrtoUlr(i, j, k, vr2, br2, x[i], Vr, Ur);
					
					
					//-----Vr, Ur--->H=Hr-----
					
					VlrUlrtoHlr(i, j, k, br2, vrdotbr, Vr, Ur, H);
				
				
				}
				else{
					
					//-----Vl--->Ul-----
					
					VlrtoUlr(i, j, k, vl2, bl2, x[i], Vl, Ul);
					
					
					//-----Vl, Ul--->Hl-----
					
					//VlrUlrtoHlr(i, j, k, bl2, vldotbl, Vl, Ul, Hl);
					
					
					//-----Vr--->Ur-----
						
					VlrtoUlr(i, j, k, vr2, br2, x[i], Vr, Ur);
					
					
					//-----Vr, Ur--->Hr-----
						
					//VlrUlrtoHlr(i, j, k, br2, vrdotbr, Vr, Ur, Hr);
					
					
					//-----entropy speed, vxll=vxlll=vxrrr=vxrr=sm-----
					
					slvl=sl-Vl[3][i][j][k];
					srvr=sr-Vr[3][i][j][k];
					
					sm=(srvr*Ur[3][i][j][k]-slvl*Ul[3][i][j][k]-ptr+ptl)
							/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);

					oneoslsm=1.0/(sl-sm);
					oneosrsm=1.0/(sr-sm);
					
						
					//-----intermediate total pressure, 
					//											ptll=ptlll=ptrrr=ptrr=ptint----- 
					
					ptint=(srvr*Ur[0][i][j][k]*ptl-slvl*Ul[0][i][j][k]*ptr
								 +Ur[0][i][j][k]*Ul[0][i][j][k]*srvr*slvl
								 *(Vr[3][i][j][k]-Vl[3][i][j][k]))
								/(srvr*Ur[0][i][j][k]-slvl*Ul[0][i][j][k]);
					
					
					//-----intermediate bz=bz_hll(Hr[6]=0, Hl[6]=0)-----
					/*
					bzint=(sr*Ur[6][i][j][k]-sl*Ul[6][i][j][k]
								 -Hr[6][i][j][k]+Hl[6][i][j][k])
								/(sr-sl);
					*/
					
					bzint=(sr*Ur[6][i][j][k]-sl*Ul[6][i][j][k])/(sr-sl);
					
					
					//-----left left left & right right right side density-----
					//-----rholl=rholll, rhorr=rhorrr-----
					
					rholl=Ul[0][i][j][k]*slvl*oneoslsm;
					rhorr=Ur[0][i][j][k]*srvr*oneosrsm;
					
					
          //-----left left side vx, vy, bx, by, vlldotbll-----
					
					//-----caution! if Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vxll=vxl, vyll=vyl, bxll=bxl, byll=byl-----
          					
          sgn=(fabs(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint)
                    -1.0e-8*bzint*bzint)
               /(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint)
                 -1.0e-8*bzint*bzint));
           
          denominator=((1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint)
                        +1.0)
                       +(1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint)
                         -1.0)*sgn)*0.5;
           
          numerator=((((Ul[0][i][j][k]*slvl*slvl-bzint*bzint)+1.0)
                      +((Ul[0][i][j][k]*slvl*slvl-bzint*bzint)-1.0)*sgn)
                     *0.5);
           
           
          vxll=Vl[1][i][j][k]-(bzint*Ul[4][i][j][k]
                               *(sm-Vl[3][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          vyll=Vl[2][i][j][k]-(bzint*Ul[5][i][j][k]
                               *(sm-Vl[3][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          bxll=Ul[4][i][j][k]*numerator*denominator;
           
          byll=Ul[5][i][j][k]*numerator*denominator;
           					
					
/*          
					//-----left left side vy, vz, by, bz, vlldotbll-----
					
					//caution! if denominator==0.0, degenerate HLLD & HLL
					
					if(fabs(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint)
						 >1.0e-8*bzint*bzint){
						
						denominator=1.0/(Ul[0][i][j][k]*slvl*(sl-sm)-bzint*bzint);
						numerator=Ul[0][i][j][k]*slvl*slvl-bzint*bzint;

						
						vxll=Vl[1][i][j][k]-(bzint*Ul[4][i][j][k]
																 *(sm-Vl[3][i][j][k])*denominator);
						
						vyll=Vl[2][i][j][k]-(bzint*Ul[5][i][j][k]
																 *(sm-Vl[3][i][j][k])*denominator);
						
						bxll=Ul[4][i][j][k]*numerator*denominator;
						
						byll=Ul[5][i][j][k]*numerator*denominator;
						
						
					}
					else{
						
						vxll=Vl[1][i][j][k];
						
						vyll=Vl[2][i][j][k];
						
						bxll=Ul[4][i][j][k];
						
						byll=Ul[5][i][j][k];
						
						
					}
*/
					
					
					vlldotbll=sm*bzint+vxll*bxll+vyll*byll;
					

          //-----right right side vx, vy, bx, by, vrrdotbrr-----
					
					//-----caution! if Ul[0][i][j][k]*slvl*(sl-sm)-bxint*bxint=0,
					//		  degenarate HLLD & HLL, so denominator=numerator=1, 
					//			vxrr=vxr, vyrr=vyr, bxrr=bxr, byrr=byr-----
          					
          sgn=(fabs(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint)
                    -1.0e-8*bzint*bzint)
               /(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint)
                 -1.0e-8*bzint*bzint));
           
          denominator=((1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint)
                        +1.0)
                       +(1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint)
                         -1.0)*sgn)*0.5;
           
          numerator=((((Ur[0][i][j][k]*srvr*srvr-bzint*bzint)+1.0)
                      +((Ur[0][i][j][k]*srvr*srvr-bzint*bzint)-1.0)*sgn)
                     *0.5);
           
           
          vxrr=Vr[1][i][j][k]-(bzint*Ur[4][i][j][k]
                               *(sm-Vr[3][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          vyrr=Vr[2][i][j][k]-(bzint*Ur[5][i][j][k]
                               *(sm-Vr[3][i][j][k])
                               *((denominator+denominator*sgn)*0.5));
           
          bxrr=Ur[4][i][j][k]*numerator*denominator;
           
          byrr=Ur[5][i][j][k]*numerator*denominator;
           				

/*          
					//-----right right side vy, vz, by, bz, vrrdotbrr-----
					
					//caution! if denominator==0.0, degenerate HLLD & HLL
					
					if(fabs(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint)
						 >1.0e-8*bzint*bzint){
						
						denominator=1.0/(Ur[0][i][j][k]*srvr*(sr-sm)-bzint*bzint);
						numerator=Ur[0][i][j][k]*srvr*srvr-bzint*bzint;
						
						
						vxrr=Vr[1][i][j][k]-(bzint*Ur[4][i][j][k]
																 *(sm-Vr[3][i][j][k])*denominator);
						
						vyrr=Vr[2][i][j][k]-(bzint*Ur[5][i][j][k]
																 *(sm-Vr[3][i][j][k])*denominator);
						
						bxrr=Ur[4][i][j][k]*numerator*denominator;
						
						byrr=Ur[5][i][j][k]*numerator*denominator;
						
						
					}
					else{	
						
						vxrr=Vr[1][i][j][k];
						
						vyrr=Vr[2][i][j][k];
						
						bxrr=Ur[4][i][j][k];
						
						byrr=Ur[5][i][j][k];
						
						
					}
*/
          
					
					vrrdotbrr=sm*bzint+vxrr*bxrr+vyrr*byrr;
					
					
					//-----left left side en-----
					
					enll=(slvl*Ul[7][i][j][k]-ptl*Vl[3][i][j][k]+ptint*sm
								+bzint*(vldotbl-vlldotbll))*oneoslsm;

					
					//-----right right side en-----
					
					enrr=(srvr*Ur[7][i][j][k]-ptr*Vr[3][i][j][k]+ptint*sm
								+bzint*(vrdotbr-vrrdotbrr))*oneosrsm;
					
					
					//-----left & right direction alfven speed-----
					
					sqrtrholl=sqrt(rholl);
					sqrtrhorr=sqrt(rhorr);
					
					sll=sm-fabs(bzint)/sqrtrholl;
					srr=sm+fabs(bzint)/sqrtrhorr;
					

					if(sm>=0.0){
						
						if(sll>=0.0){
						
							//-----H=Hhlll-----
						
							PrimitivetoH(i, j, k, vlldotbll, rholl, vxll, vyll, sm, 
													 bxll, byll, bzint, enll, ptint, H);
						
							
						}
						else{
							
							sgn=fabs(bzint)/bzint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
							
							//-----H=Hhlldl-----
							
							vxint=(sqrtrholl*vxll+sqrtrhorr*vxrr+(bxrr-bxll)*sgn)
										*oneosqrtrho;
							vyint=(sqrtrholl*vyll+sqrtrhorr*vyrr+(byrr-byll)*sgn)
										*oneosqrtrho;
							bxint=(sqrtrholl*bxrr+sqrtrhorr*bxll
										 +sqrtrholl*sqrtrhorr*(vxrr-vxll)*sgn)
										*oneosqrtrho;
							byint=(sqrtrholl*byrr+sqrtrhorr*byll
										 +sqrtrholl*sqrtrhorr*(vyrr-vyll)*sgn)
										*oneosqrtrho;
							
							
							vintdotbint=sm*bzint+vxint*bxint+vyint*byint;
							
							
							enlll=enll-sqrtrholl*(vlldotbll-vintdotbint)*sgn;
							
							
							PrimitivetoH(i, j, k, vintdotbint, rholl, vxint, vyint, sm, 
													 bxint, byint, bzint, enlll, ptint, H);
							
							
						}
					}
					else{
						
						if(srr<=0.0){
							
							//-----H=Hhllr----
							
							PrimitivetoH(i, j, k, vrrdotbrr, rhorr, vxrr, vyrr, sm, 
													 bxrr, byrr, bzint, enrr, ptint, H);
						
							
						}
						else{
						
							sgn=fabs(bzint)/bzint;
							oneosqrtrho=1.0/(sqrtrholl+sqrtrhorr);
						
							
							//-----H=Hhlldr-----
								
							vxint=(sqrtrholl*vxll+sqrtrhorr*vxrr+(bxrr-bxll)*sgn)
										*oneosqrtrho;
							vyint=(sqrtrholl*vyll+sqrtrhorr*vyrr+(byrr-byll)*sgn)
										*oneosqrtrho;
							bxint=(sqrtrholl*bxrr+sqrtrhorr*bxll
										 +sqrtrholl*sqrtrhorr*(vxrr-vxll)*sgn)
										*oneosqrtrho;
							byint=(sqrtrholl*byrr+sqrtrhorr*byll
										 +sqrtrholl*sqrtrhorr*(vyrr-vyll)*sgn)
										*oneosqrtrho;
								
								
							vintdotbint=sm*bzint+vxint*bxint+vyint*byint;
								
								
							enrrr=enrr+sqrtrhorr*(vrrdotbrr-vintdotbint)*sgn;
								
								
							PrimitivetoH(i, j, k, vintdotbint, rhorr, vxint, vyint, sm, 
													 bxint, byint, bzint, enrrr, ptint, H);
								
							
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
				
				H[6][i][j][k]+=0.5*(Vl[8][i][j][k]+Vr[8][i][j][k]
														+ch*(Vl[6][i][j][k]-Vr[6][i][j][k]));
				H[8][i][j][k]=0.5*ch*(ch*(Vl[6][i][j][k]+Vr[6][i][j][k])
															+(Vl[8][i][j][k]-Vr[8][i][j][k]));
				
			}
		}
	}
	
	
	
	return 0;

}









