#include "interactions.h"
/*
 * Radial distribution function
 */
void rdf(double Lx, double Ly, double Lz, int np, double **r,
		int dim, int nspacBeads, int nStikers, int nSpacers, int npc)
{

	Distances distances;
	
	/*
	 * i,j,k,jj indices
	 * mbil - max. bin length
	 * np - number of particles
	 */

	fstream out;
	int i,j,bn,tn,l,m,n;
	static int mbil,csize,mt;
	static unsigned long int k = 0;

	/*
	 * cvn - nor. constant
	 * bl - box length
	 * dr - size of a bin
	 * dis1,dis2,dis3,dis4 - temp var.s
	 * cbn - norm. for the bins
	 * gr - rad. dis. fn.
	 * ds2 - ds*ds
	 * ds  - distance
	 */

	double dis1,dis2,dis3,dis4;
	double cvn,ds2,ds,dx,dy,dz;
	static double ***g,*cbn,**gr,***lgr,
								dr = 0.01,pi,V;

	// Normalization constants for bins
	if(!k){

		mt = omp_get_max_threads();

		V = Lx*Ly*Lz;

		csize = 1;

		pi = 4.0*atan(1.0);

		mbil = (int)(Lz/2.0/dr);

		cbn = new double[mbil];
		gr = new double*[mbil];
		for(i = 0; i < mbil; i++) gr[i] = new double[3];
		
		lgr = new double**[mt];
		g = new double**[mt];

		for(i = 0; i < mt; i++){
			lgr[i] = new double*[mbil];
			g[i] = new double*[np];
			for(j = 0; j < mbil; j++) lgr[i][j] = new double[3];
			for(j = 0; j < np; j++) g[i][j] = new double[dim];
		}

		for(i = 0; i < mbil; i++){
			cbn[i] = 0.0;
			for( j = 0; j < 3; j++) gr[i][j] = 0.0;
		}

		cvn = 4.0*pi/3.0;
		dis1 = dr;
		dis2 = 0.0;

		for(j = 0; j < mbil; j++){
			dis3 = dis1*dis1;
			dis4 = dis2*dis2;
			cbn[j] = cvn*(dis1*dis3 - dis2*dis4);
			dis1 += dr;
			dis2 += dr;
		}
	}

	k++;

	// Allocate co-ordinates to all the threads
	for(j = 0; j < mt; j++){
		for(i = 0; i < np; i++)for(l = 0; l < dim; l++) g[j][i][l] = r[i][l];
		for(i = 0; i < mbil; i++)for(l = 0; l < 3; l++) lgr[j][i][l] = 0.0;
	}

#pragma omp parallel for firstprivate(np,npc,dr,Lx,Ly,Lz,mbil,nspacBeads)\
	private(j,dx,dy,dz,ds2,ds,bn,tn,m,n) schedule(dynamic,csize)
	for(i = 0; i < np; i++){

		m = beadIndex(i,npc);

		tn = omp_get_thread_num();

		for(j = i + 1; j < np; j++){
		
			n = beadIndex(j,npc);

			dx = g[tn][i][0] - g[tn][j][0];
			dy = g[tn][i][1] - g[tn][j][1];
			dz = g[tn][i][2] - g[tn][j][2];

	    dx = distances.minimumImage(Lx,dx);
	    dy = distances.minimumImage(Ly,dy);
	    dz = distances.minimumImage(Lz,dz);

			ds2 = distances.squareDistance(dx,dy,dz);
			ds = sqrt(ds2);

			if(ds < 0.25){
				cout << "grs.c: " << ds <<" "<< dx <<" "<< dy <<" "<< dz << endl;
				cout << "grs.c: " << i <<" "<< j << endl;
				exit(EXIT_FAILURE);
			}

			bn = (int)(ds/dr);

			if(bn < mbil){
			 if(m%nspacBeads != 0 && n%nspacBeads != 0) lgr[tn][bn][0] += 2.0;
			 else if(!(m%nspacBeads) && !(n%nspacBeads)) lgr[tn][bn][1] += 2.0;
			 else lgr[tn][bn][2] += 1.0;
			}
		}
	}

	// Collect data from all threads
	for(i = 0; i < mbil; i++)for(j = 0; j < mt; j++)for(l = 0; l < 3; l++) gr[i][l] += lgr[j][i][l];

	if(!(k%100)){
		
		out.open("rdf.dat", ios::out);
		
		for(i = 0; i < mbil; i++){
			
			dis1 = ((double)i + 0.5)*dr;
		
			dis2 = (V*gr[i][0])/cbn[i]/((double)(k*nStikers*(nStikers - 1)));
			dis3 = (V*gr[i][1])/cbn[i]/((double)(k*nSpacers*(nSpacers - 1)));
			dis4 = (V*gr[i][2])/cbn[i]/((double)(k*nStikers*nSpacers));
			
			out << dis1 <<"  "<< dis2 <<"  "<< dis3 <<"  "<< dis4 << endl;
		}

		out.close();
	}

}
