

/***********************************************************************
 *
 *	3-dimensional MHD simulation code in Cylindrical coord.
 *
 *	numerical flux : HLLD
 *									(Miyosho & Kusano, 2005)
 *
 *	flux limiter : 2nd order limiters
 *
 *	time marching : 2nd order TVD Runge-Kutta
 *
 *	
 *	2013 Apr. 05 : making cylindrical grid.
 *	2012 Dec. 26 : time 2nd order ---> 3rd order.
 *	2012 Nov. 21 : using margin & generalized (for higher order).
 *	2012 Nov. 19 : add vanLeer, MC, superbee limiter(2nd order).
 *	2012 Nov. 06 : using pointer & function improvement.
 *	2012 Nov. 03 : add GLM divergence cleaning.
 *	2012 Oct. 10 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <omp.h>


//============================================================
//	definition of max(a, b), min(a, b) functions.
//============================================================

#include "maxmin.h"


//============================================================
//	definition of physical constant
//============================================================

#include "constant.h"


//============================================================
//	input parameters
//============================================================

#include "input.h"


//============================================================
//	definition of quantities
//============================================================

#include "commonfield.h"
#include "commontemp.h"
#include "commongrid.h"


//============================================================
//	make cell center
//============================================================

#include "cartesiangrid.h"
#include "cylindricalgrid.h"


//============================================================
//	calculate dt from CFL condition
//============================================================

#include "setdt.h"


//============================================================
//	set boundary conditions
//============================================================

#include "xboundary.h"
#include "yboundary.h"
#include "zboundary.h"
#include "centralboundary.h"
#include "xvlrboundary.h"
#include "yvlrboundary.h"
#include "zvlrboundary.h"


//============================================================
//	convert from V to U, from U to V, from U1 to V.
//============================================================

#include "vtou.h"
#include "utov.h"


//============================================================
//	save initial conditions.
//============================================================

#include "saveinitial.h"


//============================================================
//	set initial model
//============================================================

#include "model.h"


//============================================================
//	set initial gravity
//============================================================

#include "initialgravity.h"


//============================================================
//	set spiral arm gravity
//============================================================

#include "spiralgravity.h"


//============================================================
//	set source term & cooling energy loss
//============================================================

#include "cooling.h"
#include "source.h"


//============================================================
//	convert from Vl to Ul, from Vr to Ur
//============================================================

#include "vlrtoulr.h"


//============================================================
//	convert from Vl & Ul to Fl, from Vr & Ur to Fr
//============================================================

#include "vlrulrtoflr.h"


//============================================================
//	convert from Vl & Ul to Gl, from Vr & Ur to Gr
//============================================================

#include "vlrulrtoglr.h"


//============================================================
//	convert from Vl & Ul to Hl, from Vr & Ur to Hr
//============================================================

#include "vlrulrtohlr.h"


//============================================================
//	convert from primitive variables to F, G, H
//============================================================

#include "primitivetof.h"
#include "primitivetog.h"
#include "primitivetoh.h"


//============================================================
//	main solvers
//============================================================

#include "glmsourcesolver.h"
#include "limiterx.h"
#include "limitery.h"
#include "limiterz.h"
#include "hlldx.h"
#include "hlldy.h"
#include "hlldz.h"
#include "hllx.h"
#include "hlly.h"
#include "hllz.h"
#include "halfstep.h"
#include "fullstep.h"
#include "firststep.h"
#include "secondstep.h"
#include "thirdstep.h"
#include "tvdrungekutta.h"


//============================================================
//	negative density & pressure checker.
//============================================================

#include "check.h"


//============================================================
//	define output function.
//============================================================

#include "output.h"


//============================================================
//	initial perturbation.
//============================================================

#include "perturbation.h"


//============================================================
//	read binary files for restart simulation.
//============================================================

#include "restartread.h"



int main(void)
{
	
	//============================================================
	//	set cell center
	//============================================================
	
	//CartesianGrid();
	CylindricalGrid();
	
	
	printf("complete grid setting.\n");
  
  
  //-----set t=0, nstep=0, outputfile number=0-----
  
  nstep=0;
  t=0.0;
  
  
  if(restartsw){
   
    //==================================================
    //  inputting initial rho & pr
    //==================================================
      
    RestartRead(0);
    SaveInitial();
      
    
    //==================================================
    //  read pre-simulation binary data
    //==================================================
    
    RestartRead(nout);
    
    
    //=======================================================
    //	add toroidal B-fields & change to co-rotation frame
    //=======================================================
    
    //AddMagneticFields();
    //Corotationframe();
    
    
    t=nout*tint;
    tout=(nout+1)*tint;
    
    
    //==================================================
    //	set axisymmetric gravity
    //==================================================
    
    InitialGravity();
    
    
    //GNUOutputCylinder();
    //GDLOutputCylinder();
    
    
    //BinaryOutput();
      
    //exit(1);
      
  }
  else{
    
	
    //============================================================
    //	set initial condition
    //============================================================
	
    //TorusHalo();
    //DoubleExponentialDisc();
    //GalacticDisc2();
    GalacticDisc3();
    
    
    //GNUOutputCylinder();
    //GDLOutputCylinder();
    BinaryOutput();
    
    
    //exit(1);
  
  }
  
  
	printf("complete model setting.\n");
  
	
	
	while(nstep<nstop+1){
    
    
		nstep++;
		
		printf("nstep = %d\n", nstep);
		
		
		Setdt();
				
		printf("dt = %.3e\n", dt);
	
		
		//-----if dt=inf or dt<dtmin, program is force-quite-----
		
		if(isinf(dt) || dt<dtmin){
			printf("stop due to small dt or inf. program is done.\n");
			exit(1);
		}

		
    //============================================================
    //	set spiral arm gravity
    //============================================================
    
    SpiralGravity();
    
    
		//============================================================
		//	3rd order TVD Runge-Kutta time marching
		//============================================================
		
		TVDRK();
		
		
		t+=dt;
		
        
    
		
		//-----if Time>outputtime, data set output-----
		
		if(t>tout){
        	
            printf("time=%e.\n output complete.\n", t);
			
			
			//============================================================
			//	output for analysis & visualization
			//============================================================
			
			tout+=tint;
			nout++;
			
			
			//GNUOutputCylinder();
            //GDLOutputCylinder();
			//AVSOutput();
			BinaryOutput();
			
        	//exit(0);
            //return 0;
		}
		
		
		//-----if Time>simulationtime, program is done-----
		
		if(t>tstop){
			
			printf("time is over %.3e.\n", tstop);
			
			
			break;
			
		}
		
		
	}

	
	return 0;

}


