

/***********************************************************************
 *
 *	output.h
 *
 *	data set output for analitics, visualization.
 *
 *
 *	gnuplot, IDL(GDL), AVS, matplotlib
 *
 *
 *	2013 Apr. 08 : add data output in Cylindrical coord.
 *	2012 Nov. 02 : add divB output for GDL 2-d visualization.
 *	2012 Oct. 10 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/
 

//==================================================
// data set output for gnuplot in Cartesian
//==================================================


int GNUOutputCylinder(void)
{
	
	int i, j, k;
	
	
	//-----0<i<ixmax1, j=jy/2, k=kz/2-----
	
	sprintf(filename, "./1d-r-%d.dat", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open .%s.\n", filename);
		exit(1);
	}
	
	
	j=1;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, j, k);
	fprintf(fo, "#x rho vx vy vz bx by bz en pr\n");
	
	
	for(i=1; i<ixmax1; i++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						x[i], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	//-----i=ix/2, 0<j<jymax1, k=kz/2-----
	
	sprintf(filename, "./1d-phi-%d.dat", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	
	i=160;
	k=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, i, k);
	fprintf(fo, "#y rho vx vy vz bx by bz en pr\n");
	
	
	for(j=1; j<jymax1; j++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						y[j], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	//-----i=ix/2, j=jy/2, 0<k<kzmax1-----
	
	sprintf(filename, "./1d-z-%d.dat", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	
	i=160;
	j=1;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, i, j);
	fprintf(fo, "#z rho vx vy vz bx by bz en pr\n");
	
	
	for(k=1; k<kzmax1; k++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						z[k], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	printf("\n output %d %.3e\n", nout, t);
	
	
	return 0;
	
}



//==================================================
// data set output for gnuplot in Cartesian
//==================================================


int GNUOutputCartesian(void)
{
	
	int i, j, k;
	
	
//-----0<i<ixmax1, j=jy/2, k=kz/2-----
	
	sprintf(filename, "./1d-x-%d.dat", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}

	
	j=jy/2;
	k=kz/2;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, j, k);
	fprintf(fo, "#x rho vx vy vz bx by bz en pr\n");
	
	
	for(i=1; i<ixmax1; i++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						x[i], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	//-----i=ix/2, 0<j<jymax1, k=kz/2-----
	
	sprintf(filename, "./1d-y-%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	i=ix/2;
	k=kz/2;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, i, k);
	fprintf(fo, "#y rho vx vy vz bx by bz en pr\n");
	
	
	for(j=1; j<jymax1; j++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						y[j], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	//-----i=ix/2, j=jy/2, 0<k<kz+1-----
	
	sprintf(filename, "./1d-z-%d.txt", nout);
	
	fo=fopen(filename, "w");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	
	i=ix/2;
	j=jy/2;
	
	fprintf(fo, "#\t nstep=%d,\ttime=%.3e,\tj=%d,\tk=%d \n", 
					nstep, t, i, j);
	fprintf(fo, "#z rho vx vy vz bx by bz en pr\n");
	
	
	for(k=1; k<kzmax1; k++){
		
		fprintf(fo, "%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\t%.3e\t%.3e\t%.3e\t"
						"%.3e\n", 
						z[k], 
						V[0][i][j][k], V[1][i][j][k], V[2][i][j][k], V[3][i][j][k], 
						V[4][i][j][k], V[5][i][j][k], V[6][i][j][k], U[7][i][j][k],
						V[7][i][j][k]);
		
		
	}
	
	fclose(fo);
	
	
	printf("\n output %d %.3e\n", nout, t);
	
	
	return 0;

}



//==================================================
// output for IDL, GDL in Cylindrical coord.
//==================================================


int GDLOutputCylinder(void)
{
	
	int i, j, k;
	
	
	//-----output x-----
	
	fo=fopen("r.dat", "w");	
	for(i=1; i<ixmax1; i++){
		fprintf(fo, "%e\n", x[i]);
	}
	
	fclose(fo);
	
	
	//-----output y-----
	
	fo=fopen("phi.dat", "w");
	for(j=1; j<jymax1; j++){
		fprintf(fo, "%e\n", y[j]);
	}
	
	fclose(fo);
	
	
	//-----output z-----
	
	fo=fopen("z.dat", "w");
	for(k=1; k<kzmax1; k++){
		fprintf(fo, "%e\n", z[k]);
	}	
	
	fclose(fo);
	
	
	//-----output k=1 density-----
	
	sprintf(filename, "./2d-density-rphi-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=1;
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output k=1 pressure-----
	
	sprintf(filename, "./2d-pressure-rphi-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=1;
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			fprintf(fo, "%e\n", V[7][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output k=1 benergy-----
	
	sprintf(filename, "./2d-benergy-rphi-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=1;
	for(i=1; i<ixmax1; i++){
		for(j=1; j<jymax1; j++){
			fprintf(fo, "%e\n", 0.5*(V[4][i][j][k]*V[4][i][j][k]
															 +V[5][i][j][k]*V[5][i][j][k]
															 +V[6][i][j][k]*V[6][i][j][k]));
		}
	}
	
	fclose(fo);	
	
	
	//-----output j=0 density-----
	
	sprintf(filename, "./2d-density-rz-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	j=0;
	for(k=1; k<kzmax1; k++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	return 0;
	
}



//==================================================
// data set output for IDL, GDL in Cartesian coord.
//==================================================


int GDLOutputCartesian(void)
{

	int i, j, k;
	
	
	//-----output x-----
	
	fo=fopen("x.dat", "w");	
	for(i=1; i<ixmax1; i++){
		fprintf(fo, "%e\n", x[i]);
	}
	
	fclose(fo);
	
	
	//-----output y-----
	
	fo=fopen("y.dat", "w");
	for(j=1; j<jymax1; j++){
		fprintf(fo, "%e\n", y[j]);
	}
	
	fclose(fo);
	
	
	//-----output z-----
	
	fo=fopen("z.dat", "w");
	for(k=1; k<kzmax1; k++){
		fprintf(fo, "%e\n", z[k]);
	}	
	
	fclose(fo);
	
	
	//-----output density-----
	
	sprintf(filename, "./2d-density-xy-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=kz/2;
	for(j=1; j<jymax1; j++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output density-----
	
	sprintf(filename, "./2d-density-zx-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	j=jy/2;
	for(k=1; k<kzmax1; k++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output pressure-----
	
	sprintf(filename, "./2d-pressure-xy-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=kz/2;
	for(j=1; j<jymax1; j++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[7][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output pressure-----
	
	sprintf(filename, "./2d-pressure-zx-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	j=jy/2;
	for(k=1; k<kzmax1; k++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[7][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output benergy-----
	
	sprintf(filename, "./2d-benergy-xy-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	k=kz/2;
	for(j=1; j<jymax1; j++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", 0.5*(V[4][i][j][k]*V[4][i][j][k]
															 +V[5][i][j][k]*V[5][i][j][k]
															 +V[6][i][j][k]*V[6][i][j][k]));
		}
	}
	
	fclose(fo);
	
	
	//-----output by-----
	
	sprintf(filename, "./2d-by-zx-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	j=jy/2;
	for(k=1; k<kzmax1; k++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[5][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output bz-----
	
	sprintf(filename, "./2d-bz-zx-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	j=jy/2;
	for(k=1; k<kzmax1; k++){
		for(i=1; i<ixmax1; i++){
			fprintf(fo, "%e\n", V[6][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----output density-----
	
	sprintf(filename, "./2d-density-zy-%d.dat", nout);
	
	fo=fopen(filename, "w");
	
	i=ix/2;
	for(k=1; k<kzmax1; k++){
		for(j=1; j<jymax1; j++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	
	/*
	//-----output density-----
	
	sprintf(filename, "./2d-density-zx-%d.txt", nout);
	
	fo=fopen(filename, "w");
	
	j=1;
	for(i=1; i<ix+1; i++){
		for(k=1; k<kz+1; k++){
			fprintf(fo, "%e\n", V[0][i][j][k]);
		}
	}
	
	fclose(fo);
	*/
	
	/*
	//-----output pressure-----
	
	sprintf(filename, "./2d-pressure%d.txt", nout);
	
	fo=fopen(filename, "w");
	
	k=kz/2;
	for(j=1; j<jy+1; j++){
		for(i=1; i<ix+1; i++){
			fprintf(fo, "%e\n", V[7][i][j][k]);
		}
	}
	
	fclose(fo);
	
	
	//-----calculate & output divB-----
	
	MinmodLimiterX();
	
	
	for(i=0; i<ix; i++){
		for(j=0; j<jy; j++){
			for(k=0; k<kz; k++){
				
				divB[i][j][k]=(Vr[4][i+1][j][k]-Vl[4][i][j][k])/dx;
			
			
			}
		}
	}
	
	
	MinmodLimiterY();
	
	
	for(i=0; i<ix; i++){
		for(j=0; j<jy; j++){
			for(k=0; k<kz; k++){
				
				divB[i][j][k]+=(Vr[5][i][j+1][k]-Vl[5][i][j][k])/dy;
				
				
			}
		}
	}
	
	
	sprintf(filename, "./2d-divb%d.txt", nout);
	
	fo=fopen(filename, "w");
	
	k=kz/2;
	for(j=0; j<jy; j++){
		for(i=0; i<ix; i++){
			fprintf(fo, "%e\n", divB[i][j][k]);
		}
	}
	
	fclose(fo);
	*/
	
	return 0;

}


/*
// for AVS output

int AVSOutput(void)
{
	
//============================================================
//	convert double ---> float & thin out
//============================================================

#pragma omp parallel for private(i, j, k, I, J, K, onerho)
	for(i=0; i<avsix+1; i++){
		for(j=0; j<avsjy+1; j++){
			for(k=0; k<avskz+1; k++){
				
				I=2*i;
				J=j;
				K=2*k;
				
				avsx[i][j][k]=x[I]*cos(y[J])*onekpc;
				avsy[i][j][k]=x[I]*sin(y[J])*onekpc;
				avsz[i][j][k]=z[K]*onekpc;
				
				onerho=1.0/rho[I][J][K];
				
				avsrho[i][j][k]=rho[I][J][K];
				
				
				avsvx[i][j][k]=(rvx[I][J][K]*cos(y[J])-rvy[I][J][K]*sin(y[J]))*onerho;
				avsvy[i][j][k]=(rvx[I][J][K]*sin(y[J])+rvy[I][J][K]*cos(y[J]))*onerho;
				avsvz[i][j][k]=rvz[I][J][K]*onerho;
				avsbx[i][j][k]=bx[I][J][K]*cos(y[J])-by[I][J][K]*sin(y[J]);
				avsby[i][j][k]=bx[I][J][K]*sin(y[J])+by[I][J][K]*cos(y[J]);
				avsbz[i][j][k]=bz[I][J][K];
				avsee[i][j][k]=ee[I][J][K];
				avspr[i][j][k]=pr[I][J][K];
				avste[i][j][k]=pr[I][J][K]*onerho/(onemump*kb);
				avsbeta[i][j][k]=pr[I][J][K]
													/(onepi8*(bx[I][J][K]*bx[I][J][K]
																		+by[I][J][K]*by[I][J][K]
																		+bz[I][J][K]*bz[I][J][K]));
				avsBenergy[i][j][k]=onepi8*(bx[I][J][K]*bx[I][J][K]
																		+by[I][J][K]*by[I][J][K]
																		+bz[I][J][K]*bz[I][J][K]);
				avspot[i][j][k]=pot[I][J][K];
				avsspiral[i][j][k]=diskspiral[I][J][K];
				
			}
		}
	}
	
	
//============================================================
//	make output step directory
//============================================================
	
	
	sprintf(dirname, "./avs%d", noutput);
	mkdir(dirname, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP |
				S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH | S_IXOTH);
	
	
//============================================================
//	output files in maked directory. 
//============================================================
	
	
	//========================================
	//	output files in maked directory. 
	//========================================
	
	
	//x
	sprintf(filename, "%s/coord1", dirname);
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				fwrite(&avsx[i][j][k], sizeof avsx[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	

	//y
	sprintf(filename, "%s/coord2", dirname);
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}

	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				fwrite(&avsy[i][j][k], sizeof avsy[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}

	
	fclose(fo);

	
	//z
	sprintf(filename, "%s/coord3", dirname);
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}

	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				fwrite(&avsz[i][j][k], sizeof avsz[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//rho
	sprintf(filename, "%s/rho", dirname);
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsrho[i][j][k], sizeof avsrho[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	

	//vx
	sprintf(filename, "%s/vx", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsvx[i][j][k], sizeof avsvx[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//vy
	sprintf(filename, "%s/vy", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsvy[i][j][k], sizeof avsvy[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//vz
	sprintf(filename, "%s/vz", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsvz[i][j][k], sizeof avsvz[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//bx
	sprintf(filename, "%s/bx", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsbx[i][j][k], sizeof avsbx[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	//by
	sprintf(filename, "%s/by", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsby[i][j][k], sizeof avsby[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//bz
	sprintf(filename, "%s/bz", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsbz[i][j][k], sizeof avsbz[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//ee
	sprintf(filename, "%s/ee", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsee[i][j][k], sizeof avsee[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//pr
	sprintf(filename, "%s/pr", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avspr[i][j][k], sizeof avspr[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//te
	sprintf(filename, "%s/te", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avste[i][j][k], sizeof avste[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//beta
	sprintf(filename, "%s/beta", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsbeta[i][j][k], sizeof avsbeta[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	//Benergy
	sprintf(filename, "%s/Benergy", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsBenergy[i][j][k], sizeof avsBenergy[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	//pot
	sprintf(filename, "%s/pot", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avspot[i][j][k], sizeof avspot[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	//spiral
	sprintf(filename, "%s/spiral", dirname);
	fo=fopen(filename, "wb");
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	k=0;
	while(k<avskz+1){
		j=0;
		while(j<avsjy+1){
			i=0;
			while(i<avsix+1){
				
				fwrite(&avsspiral[i][j][k], sizeof avsspiral[i][j][k], 1, fo);
				
				i++;
			}
			j++;
		}
		k++;
	}
	
	fclose(fo);
	
	
	return 0;
}

/*
// for IDL output

int IDLOutput(void)
{

	return 0;
}


// for OpenDX output

int ODXOutput(void)
{
	
	return 0;
}
*/




//all data binary output for analysis & restart 

int BinaryOutput(void)
{
	//========================================
	//	make directory (output step)
	//========================================
	
	sprintf(dirname, "./%d", nout);
	mkdir(dirname, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP |
				S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH | S_IXOTH);
	
	//========================================
	//	output files in maked directory. 
	//========================================
	
	//grid number, Time
	sprintf(filename, "%s/GridTime.dat", dirname);
	
	fo=fopen(filename, "w");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fprintf(fo, "%d\t%d\t%d\t%.3e", ixmax, jymax, kzmax, t);
	
	fclose(fo);
	
	
	//r
	sprintf(filename, "%s/r", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(x, sizeof x, 1, fo);
	
	fclose(fo);
	
	
	//phi
	sprintf(filename, "%s/phi", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(y, sizeof y, 1, fo);
	
	fclose(fo);
	
	
	//z
	sprintf(filename, "%s/z", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(z, sizeof z, 1, fo);
	
	fclose(fo);
	
	
	//rho
	sprintf(filename, "%s/rho", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[0], sizeof V[0], 1, fo);
	
	fclose(fo);
	
	
	//vr
	sprintf(filename, "%s/vr", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[1], sizeof V[1], 1, fo);
	
	fclose(fo);
	
	
	//vphi
	sprintf(filename, "%s/vphi", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[2], sizeof V[1], 1, fo);
	
	fclose(fo);
	
	
	//vz
	sprintf(filename, "%s/vz", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[3], sizeof V[3], 1, fo);
	
	fclose(fo);
	
	
	//bx
	sprintf(filename, "%s/br", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[4], sizeof V[4], 1, fo);
	
	fclose(fo);
	
	
	//by
	sprintf(filename, "%s/bphi", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[5], sizeof V[5], 1, fo);
	
	fclose(fo);
	
	
	//bz
	sprintf(filename, "%s/bz", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[6], sizeof V[6], 1, fo);
	
	fclose(fo);
	
	
	//pr
	sprintf(filename, "%s/pr", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[7], sizeof V[7], 1, fo);
	
	fclose(fo);
	
	
	//psi
	sprintf(filename, "%s/psi", dirname);
	
	fo=fopen(filename, "wb");
	
	if(fo==NULL){
		printf("I can't open %s.\n", filename);
		exit(1);
	}
	
	fwrite(V[8], sizeof V[8], 1, fo);
	
	fclose(fo);
	
	
	printf("%d step binary output complete.\n", nout);
	
	return 0;
}




