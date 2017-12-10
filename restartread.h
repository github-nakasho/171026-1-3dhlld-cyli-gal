/*
 restartread.h
 
 read binary datas.
 
 11/Feb/2012 coded by Sho Nakamura.
*/

int RestartRead(int num)
{
	
	//========================================
	// read field.
	//========================================
	
	sprintf(filename, "%d/rho", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[0], sizeof V[0], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/vr", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[1], sizeof V[1], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/vphi", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[2], sizeof V[2], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/vz", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[3], sizeof V[3], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/br", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[4], sizeof V[4], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/bphi", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[5], sizeof V[5], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/bz", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[6], sizeof V[6], 1, fi);
	
	fclose(fi);
	
	
	
	sprintf(filename, "%d/pr", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[7], sizeof V[7], 1, fi);
	
	fclose(fi);
	
	
	sprintf(filename, "%d/psi", num);
	fi=fopen(filename, "rb");
	if(fi==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fread(V[8], sizeof V[8], 1, fi);
	
	fclose(fi);
	
	
	printf("%d binary read complete.\n", num);
    
    
    //============================================================
    //  converting from primitive to conservative variables
    //============================================================
    
    VtoU();
    
    
    //========================================
    //  checking rho, pr & plasma beta
    //========================================
    
    Check(U, V);
	
	
	return 0;
}



