// Compute forces by forming cell list 
// for non-bonding interactions
//
double forcecell(int np, double Lx, double Ly, double Lz, double **r,
		int npc, int nc, double *rcut2, double *sig2, double kb, int dim,
		double *eps, int nspacBeads, double *sigma, double celLenList)
{
	static double clx,cly,clz,***g,**rcut2_OMP,**sig2_OMP,**eps_OMP;

	double vLJ,dx,dy,dz,lde,pe,peBond;

	int dkm,cnm,cnk,i,j,k,l,m,n,o,ix,iy,iz,tn;

	static int ncx,ncy,ncz,mpc,mci,**cn,p=0,**cli,mcl,lm;

	Distances distances;
	Potentials potentials;
	Force force;

	// cnm - Chain number of m'th bead
	// cnk - Chain number of k'th bead
	// cl - cell length for listing particles to
	// calculate non-bonding forces
	// nc - number of cells along an axis
	// mci - maximum number of cells to be created
	// mpc - maximum number of particles expected in a cell
	// cn - array to store neighbours of each cell
	// k - index for number of times this function is called
	// cli - cell index: stores indices of particles in each cell

	if(!p){

		ncx = (int)(Lx/celLenList); //no cells on one side
		ncy = (int)(Ly/celLenList); //no cells on one side
		ncz = (int)(Lz/celLenList); //no cells on one side

		// Recalculate exact size of cell

		clx = Lx/(double)ncx;
		cly = Ly/(double)ncy;
		clz = Lz/(double)ncz;

		// Check minimum size of cell list
		if(ncx < 4 || ncy < 4 || ncz < 4){
			cout << "Either Box length(s) are small OR cell length(s) are large for non-Bonding force cell list(s)" << endl;
			exit(EXIT_FAILURE);
		}

		//mci - total no of cells is equal to cube of 
		// number of cells in each direction
		mci = ncx*ncy*ncz;

		//Allocate space for assigning neighbours of cells
		//first dimension has mci number of cells
		//second dimension has (3^3-1)/2=13 cells

		cn = new int*[mci];
		for(i=0;i<mci;i++) cn[i] = new int[13];

		// Modified code for cuboidal box
		gencell(ncx,ncy,ncz,cn);

		// mcl - maximum number of particles that can be 
		// placed inside a cell with maximum permissible density
		//as 3 times actual density

		lde = 2.0; // Consider reduced upper density as 2.00

		mcl = (int)(lde*clx*cly*clz); // density * area

		//Allocate space for storing index of particles
		//in each cell
		cli = new int*[mci];
		for(i=0;i<mci;i++) cli[i] = new int[mcl+1];

		// Note:cli[mcl] stores number of particle in a 
		//particular cell, therefore size is mcl+1
		// Arrays to be defined for parallel execution
		lm = omp_get_max_threads();

		g = new double**[lm];
		sig2_OMP = new double*[lm];
		rcut2_OMP = new double*[lm];
		eps_OMP = new double*[lm];
		for(i=0;i<lm;i++){
			g[i] = new double*[np];
			sig2_OMP[i] = new double[3];
			rcut2_OMP[i] = new double[3];
			eps_OMP[i] = new double[3];
			for(j=0;j<np;j++) g[i][j] = new double[2*dim];
		}

		p++;

		// Store diameters, cut-offs, and interaction 
		// energies per thread at once locally 
		for(i=0;i<lm;i++){
			for(j=0;j<3;j++){
				sig2_OMP[i][j] = sig2[j];
				rcut2_OMP[i][j] = rcut2[j];
				eps_OMP[i][j] = eps[j];
			}
		}
	}

	// Clear first list info

	for(i = 0; i < mci; i++) cli[i][mcl] = 0;

	//Sort particles into the cells
	
	for(i = 0; i < np; i++){

		// First apply periodic boundary condition to 
		// bring particles to box located at the +ve quadrant
		
		dx = (r[i][0] - Lx*(rint(r[i][0]/Lx) - 0.5))/clx;
		dy = (r[i][1] - Ly*(rint(r[i][1]/Ly) - 0.5))/cly;
		dz = (r[i][2] - Lz*(rint(r[i][2]/Lz) - 0.5))/clz;

		ix = (int)dx;
		iy = (int)dy;
		iz = (int)dz;

		//printf("x-index = %d y-index = %d z-index = %d\n",ix,iy,iz);

		//Index of cell which need to be updated
		j = icell(ix,iy,iz,ncx,ncy,ncz);

		//Increase the number of particle in j'th cell
		cli[j][mcl]++;

		//Check for dimension beyond the pre-defined index
		if(cli[j][mcl] >= mcl){
			cout << "forcecell.c: cell size crossed dimension of the array" << endl;
			exit(EXIT_FAILURE);
		}

		//Array index position in which i th particle is to be
		//saved - for eg first particle on 0'th position
		k = cli[j][mcl] - 1;
		//Store i'th particle in ix'th cell and iy'th position of cell index array
		cli[j][k] = i;

	}

	//Initialize force array
	
	for(i = 0; i < np; i++)for(j = (2*dim); j < (3*dim); j++) r[i][j] = 0.0;

	//Private read data for each thread

	for(i = 0; i < lm; i++){
		for(j = 0; j < np; j++){
			g[i][j][0] = r[j][0];
			g[i][j][1] = r[j][1];
			g[i][j][2] = r[j][2];
			g[i][j][3] = 0.0;
			g[i][j][4] = 0.0;
			g[i][j][5] = 0.0;
		}
	}

	// Calculate non-bonding forces from the cell list

	// Loop over all cells
	vLJ = 0.0 ; // Initialize potential energy for LJ pairs

	// Loop over all cells

#pragma omp parallel for simd firstprivate(npc,nspacBeads,Lx,Ly,Lz,mci,mcl) \
	private(j,k,tn,l,m,n,o,cnm,cnk,dkm) reduction(+:vLJ) schedule(dynamic)

	for(i = 0; i < mci; i++){
		tn = omp_get_thread_num();

		for(j = 0; j < cli[i][mcl]; j++){
			// k - index of the first particle
			k = cli[i][j];
			cnk = k/npc; // Chain number
			//Pairs from the self count
			//Loop for remaining particles in the same cell
			//avoiding self-count
			for(l = j+1; l < cli[i][mcl]; l++){
				//m- index of the second particle
				m = cli[i][l];
				cnm = m/npc; // Chain number
				dkm = abs(m-k); // Difference in particles' index

				// Find the distance and L-J force
				if(dkm != 1 || cnm != cnk) vLJ += ljf(npc,g[tn],Lx,Ly,Lz,m,k,nspacBeads,rcut2_OMP[tn],sig2_OMP[tn],eps_OMP[tn]); 

			}

			//Pairs from boxes on the neighborhood
			for(n = 0; n < 13; n++){
				//o index of neighboring boxes
				o = cn[i][n];
				for(l = 0; l < cli[o][mcl]; l++){
					// m - index of the second particles
					m = cli[o][l];
					cnm = m/npc; // Chain number
					dkm = abs(m-k); // Difference in particles' index

					//Find the distance and L-J force
					if(dkm != 1 || cnm != cnk) vLJ += ljf(npc,g[tn],Lx,Ly,Lz,m,k,nspacBeads,rcut2_OMP[tn],sig2_OMP[tn],eps_OMP[tn]);

				}//End of loop over neighbour boxes
			}//End of loop over particles - second
		}//End of loop over particle - first
	}//End of loop over all boxes

	// Bonding forces - Harmonic
	peBond = harmonic(np,g,kb,sigma,nspacBeads,npc);

	// Collecting LJ forces from all threads
	for(i = 0; i < np; i++)for(j = 0; j < lm; j++)for(k = dim; k < (2*dim); k++) r[i][k+dim] += g[j][i][k];

	pe = vLJ + peBond; // Total potential energy 

	return pe;
}
