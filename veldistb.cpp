/* Program to calculate velocities of the particles according to */
/* Maxwell-Boltzmann distributions */
void veldismb(int np, double **r, double te,
		int dim, long int *seed, int state)
{
	int i,j;
	double rtmp,ke,tinst,rscal,sumpx,
				 sumpy,sumpz,totmass,m;
	int is,x;

	// Mass of the particles
	m = 1.0;

	sumpx=sumpy=sumpz=0.0;

	for(i=0;i<np;i++){

		rtmp = sqrt(te/m);

		// Generate Gaussian distributed velocities if it is first run
		if(state != 1) for(j=dim;j<(2*dim);j++) r[i][j] = rtmp*gaus(seed);

		sumpx += (m*r[i][3]);
		sumpy += (m*r[i][4]);
		sumpz += (m*r[i][5]);

	}

	totmass = (double)np*m;

	sumpx = sumpx/totmass;
	sumpy = sumpy/totmass;
	sumpz = sumpz/totmass;

	for(i=0;i<np;i++){
		r[i][3] -= sumpx;
		r[i][4] -= sumpy;
		r[i][5] -= sumpz;
	}

	/* ############## checking the temperature and scaling ############### */

	ke = 0.0;
	for(i=0;i<np;i++) ke += m*(r[i][3]*r[i][3] + r[i][4]*r[i][4] + r[i][5]*r[i][5]);


	ke = ke*0.5;
	tinst = 2.0*ke/(3.0*((double)np - 1.0));
	cout << "veldistb.c:Instantaneous temperature before velocity scaling = " << tinst << endl;

	rscal = sqrt(te/tinst);

	for(i=0;i<np;i++)for(j=dim;j<(2*dim);j++) r[i][j] *= rscal;

	ke = 0.0;
	for(i=0;i<np;i++) ke += m*(r[i][3]*r[i][3]+r[i][4]*r[i][4]+r[i][5]*r[i][5]);


	ke = ke*0.5;
	tinst = 2.0*ke/(3.0*((double)np - 1.0));
	cout << "veldistb.c:Instantaneous temperature after velocity scaling = " << tinst << endl;

}
