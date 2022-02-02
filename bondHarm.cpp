#include "interactions.h"
/* Harmonic Bonding Force */
double harmonic(int np, double ***r, double kb, double *sig,
		int nspacBeads, int npc)
{

	Distances distances;
	Potentials potentials;
	Force force;

	int i,k,tn,m,n;
	double rx,ry,rz,r0br,ds2,ds,fij,fijx,fijy,fijz,pe;
	double sig_sticker,sig_spacer,sig_stikSpac,
				 sigma,sigbr;

	sig_sticker = sig[0]; // Bond length between two sticker beads
	sig_spacer = sig[1]; // Bond between two spacer beads
	sig_stikSpac = sig[2]; // Bond between sticker and spacer beads

	pe = 0.0;

#pragma omp parallel for firstprivate(nspacBeads,sig_sticker,sig_spacer,sig_stikSpac,np,npc,kb)\
	private(k,m,n,rx,ry,rz,ds2,ds,sigma,sigbr,fij,fijx,fijy,fijz) reduction(+:pe) schedule(dynamic)

	for(i = 0; i < np-1; i++){

		tn = omp_get_thread_num();

		k = i + 1;

		if(k%npc != 0){

			m = beadIndex(i,npc);
			n = beadIndex(k,npc);

			// Identify the covalent bond between particles' type
			if(m%nspacBeads != 0 && n%nspacBeads != 0) sigma = sig_sticker;
			else if (!(m%nspacBeads) && !(n%nspacBeads)) sigma = sig_spacer;
			else sigma = sig_stikSpac;

			rx = r[tn][k][0] - r[tn][i][0];
			ry = r[tn][k][1] - r[tn][i][1];
			rz = r[tn][k][2] - r[tn][i][2];

			ds2 = distances.squareDistance(rx,ry,rz);
			// Radial distance between two covalently bonded beads
			ds = sqrt(ds2); 
			sigbr = sigma/ds;

			pe += potentials.peHarmonic(kb,ds,sigma); 
			
			fijx = force.harmonicBond(kb,sigbr)*rx;
			fijy = force.harmonicBond(kb,sigbr)*ry;
			fijz = force.harmonicBond(kb,sigbr)*rz;

			r[tn][i][3] -= fijx;
			r[tn][i][4] -= fijy;
			r[tn][i][5] -= fijz;

			r[tn][k][3] += fijx;
			r[tn][k][4] += fijy;
			r[tn][k][5] += fijz;

		}
	}

	return (pe*0.50);
}
