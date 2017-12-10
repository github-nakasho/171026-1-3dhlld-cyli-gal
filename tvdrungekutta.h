

/***********************************************************************
 *
 *	tvdrungekutta.h
 *
 *	TVD 3rd order Runge-Kutta time marching scheme
 *
 *
 *	2013 Apr. 05 : changed to Cylindrical coord.
 *	2012 Dec. 26 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


int TVDRK(void)
{
		
	//-----GLM-MHD source term diffusion coefficient-----
	
	//psicoef=exp(-dt*ch*oneocr*oneo3);
	psicoef=exp(-dt*ch*oneocr*0.5);
  
	
  //============================================================
	//	U^(n+1/2)[i]=U^n[i]-(dt/dx)(F^n[i]-F^n[i-1])
	//	F[i-1] : left side flux (between i-1 & i)
	//	F[i] : right side flux (between i & i+1)
	//============================================================
	
	Halfstep();
	
	
	//============================================================
	//	U^(n+1)[i]=(U^n[i]+U^(n+1/2)[i]
	//							-(dt/dx)(F^(n+1/2)[i]-F^(n+1/2)[i-1]))/2
	//	F[i-1] : left side flux (between i-1 & i)
	//	F[i] : right side flux (between i & i+1)
	//============================================================
	
	Fullstep();

  
  /*
	//============================================================
	//	U^(n+1/3)[i]=U^n[i]-(dt/dx)(F^n[i]-F^n[i-1])+dt*S^n[i]
	//	F[i-1] : left side flux (between i-1 & i)
	//	F[i] : right side flux (between i & i+1)
	//============================================================
	
	FirstStep();
	
	
	//============================================================
	//	U^(n+2/3)[i]=3/4*U^n[i]+(U^(n+1/3)[i]
	//							-(dt/dx)(F^(n+1/3)[i]-F^(n+1/3)[i-1]))/4
	//	F[i-1] : left side flux (between i-1 & i)
	//	F[i] : right side flux (between i & i+1)
	//============================================================
	
	SecondStep();
	
	
	//============================================================
	//	U^(n+1)[i]=1/3*U^n[i]+(U^(n+2/3)[i]
	//							-2/3*(dt/dx)(F^(n+2/3)[i]-F^(n+2/3)[i-1]))
	//	F[i-1] : left side flux (between i-1 & i)
	//	F[i] : right side flux (between i & i+1)
	//============================================================
	
	ThirdStep();
	*/
  
  
	
	return 0;
	
}




