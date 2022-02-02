/* Velocity-verlet integrator Langevin bath (Langevin Equation) */
//
double move(int np, double **r, double **zet, double **tht, double dt,
		double te, double Lx, double Ly, double Lz, int npc, int nc,
		double *rcut2, double *sig2, long int *seed, double gama, 
		double kb, int dim, double *eps, int nspacBeads, 
		double *Diameter, double celLenList)
{
	/*
	 * nn - positions
	 * crd - no. Coordinate of each particle
	 * np - number of particles
	 * i,j  - indices
	 */
	int i,j;
	static int k=0;
	double pe;
	static double dtb2,dt2b2,pf,dtsq,
								dt3b2,sigma,**ct;

	/* position r[..][0-1], velocity r[..][2-3], force r[..][4-5], 
	 * random dis. r[..][6-7]...
	 * rad - sqrt(2dt)
	 * dt - time step 
	 * te - temperature
	 */

	if(k == 0){
		
		pf = 1.0/2.0/sqrt(3.0);
		dtsq = sqrt(dt);
		
		dtb2 = dt/2.0;
		dt2b2 = dt*dtb2;
		dt3b2 = pow(dt,3.0/2.0);
		
		sigma = sqrt(2.0*te*gama); // k_B and mass is unity

		ct = new double*[np];
		
		for(i=0;i<np;i++) ct[i] = new double[dim];
		k++;
	}

#pragma omp parallel for simd firstprivate(dt2b2,gama,sigma,dt3b2,pf,dt,np) schedule(dynamic)

	for(i=0;i<np;i++){

		ct[i][0] = dt2b2*(r[i][6] - gama*r[i][3]) + sigma*dt3b2*(0.50*zet[i][0] + pf*tht[i][0]);
		ct[i][1] = dt2b2*(r[i][7] - gama*r[i][4]) + sigma*dt3b2*(0.50*zet[i][1] + pf*tht[i][1]);
		ct[i][2] = dt2b2*(r[i][8] - gama*r[i][5]) + sigma*dt3b2*(0.50*zet[i][2] + pf*tht[i][2]);

		r[i][0] = r[i][0] + r[i][3]*dt + ct[i][0];
		r[i][1] = r[i][1] + r[i][4]*dt + ct[i][1];
		r[i][2] = r[i][2] + r[i][5]*dt + ct[i][2];
		
		// Assign forces at updated step to the previous step 
		r[i][9] = r[i][6];
		r[i][10] = r[i][7];
		r[i][11] = r[i][8];
	}

	// Calculate forces using cell list at updated positions 
	pe = forcecell(np,Lx,Ly,Lz,r,npc,nc,rcut2,sig2,kb,dim,eps,nspacBeads,Diameter,celLenList);

#pragma omp parallel for simd firstprivate(dtb2,gama,sigma,dtsq,dt,np) schedule(dynamic)

	for(i=0;i<np;i++){
		r[i][3] = r[i][3] + dtb2*(r[i][6] + r[i][9]) - dt*gama*r[i][3] + sigma*dtsq*zet[i][0] - gama*ct[i][0]; 
		r[i][4] = r[i][4] + dtb2*(r[i][7] + r[i][10]) - dt*gama*r[i][4] + sigma*dtsq*zet[i][1] - gama*ct[i][1]; 
		r[i][5] = r[i][5] + dtb2*(r[i][8] + r[i][11]) - dt*gama*r[i][5] + sigma*dtsq*zet[i][2] - gama*ct[i][2]; 

	}

	return pe/(double)np;
}
