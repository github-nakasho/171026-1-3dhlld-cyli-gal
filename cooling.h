

/***********************************************************************
 *
 *	cooling.h 
 *
 *	setting cooling energy loss term.
 *
 *	
 *	2015 May 20 : improved Cooling function.
 *	2013 May 10 : coded by Sho Nakamura (Tohoku Univ.).
 *	
 **********************************************************************/



//============================================================
//	cooling energy loss (Raymond 1978)
//============================================================


int Cooling(double dummyrho, double dummypr)
{

	double temperature=dummypr/dummyrho;
	double tempma, fabstempma, tempmb, fabstempmb;
	double tempmc, fabstempmc, tempmd, fabstempmd;
	double signa, signb, signc, signd;
	double A, B, C;

	//-----a: 10^4K, b: 10^4.24K, c: 10^4.4K, d: 10^4.88K-----

	double a=0.031, b=0.054, c=0.078, d=0.235;
	
	tempma=temperature-a;
	fabstempma=fabs(tempma);
	tempmb=temperature-b;
	fabstempmb=fabs(tempmb);
	tempmc=temperature-c;
	fabstempmc=fabs(tempmc);
	tempmd=temperature-d;
	fabstempmd=fabs(tempmd);
	

	signa=tempma/fabstempma;
	signb=tempmb/fabstempmb;
	signc=tempmc/fabstempmc;
	signd=tempmd/fabstempmd;

	A=0.108e-37*pow(te0*temperature, 10.26);
	//B=0.325e+14*pow(te0*temperature, -2.036);
	//C=0.094e-14*pow(te0*temperature, 2.165);
    B=0.0;
    C=0.0;

    
	return (0.5*(A*(1.0-signa*signb)
               +B*(1.0-signb*signc)
               +C*(1.0-signc*signd)));


	//=====a new cooling function(10^4K<T<10^5K)


	/*
	if(temperature>=0.775) return 0.0;
	else if(temperature>=0.245 && temperature<0.775) return 1.08e+5;
	else if(temperature>=0.124 && temperature<0.245){
		
		return 1.78e+6*temperature*temperature;
		
		
	}
	else if(temperature>=0.062 && temperature<0.124) return 2.41e+4;
	else return 0.0;
	*/
  
}







