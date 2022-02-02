#include "interactions.h"
// Lennard-Jones and Weeks-Chandler-Anderson Force Module

double ljf(int npc, double **r, double Lx, double Ly, double Lz, 
		int i, int j, int nspacBeads, double *rcut2, double *sig2, 
		double *eps)
{

	Distances distances;
	Potentials potentials;
	Force force;

	double ds2,ds2s,dsi6,vij,fij,
				 fijx,fijy,fijz,dx,dy,dz;
	int m,n;

	/*
	 * rx,ry - distaces
	 * ds - distance
	 * ds2 - sqare of the distance
	 * ds2 - sqare of the distance
	 * dsi6 - 1/ds^6
	 * rc - square of potential cut off
	 * vij - pair potential
	 * fij - pair force
	 * fijx, fijy - force components
	 */

	dx = r[j][0] - r[i][0];
	dy = r[j][1] - r[i][1];
	dz = r[j][2] - r[i][2];

	dx = distances.minimumImage(Lx,dx);
	dy = distances.minimumImage(Ly,dy);
	dz = distances.minimumImage(Lz,dz);

	ds2 = distances.squareDistance(dx,dy,dz);

	// Identify the bead type whether sticker OR spacer 
	// along a chain 

	m = beadIndex(i,npc);
	n = beadIndex(j,npc);

	// Interactions between sticker beads
	if(m%nspacBeads != 0 && n%nspacBeads != 0 ){ 

		if(ds2 <= rcut2[0]){

			ds2s = distances.inverseSquare(sig2[0],ds2);
			dsi6 = distances.inverseSquareCube(ds2s);

			vij = potentials.peLJ(eps[0],dsi6);
			fij = force.nonBonding(eps[0],ds2,dsi6);

			if(fabs(ds2) < 0.25){
				cout << "flj.cpp: vij =" << vij << " fij = " << fij << " ds2 = " << ds2 << " ds2s = " << ds2s << endl;
				cout << "flj.cpp: dx = " << dx << " dy = " << dy << " dz = " << dz << "sig2LJ = " << sig2[0] << endl;
				cout << "flj:i = " << i << " j = " << j << endl;
				exit(EXIT_FAILURE);
			}

			fijx = fij*dx;
			fijy = fij*dy;
			fijz = fij*dz;

			r[i][3] -= fijx;
			r[i][4] -= fijy;
			r[i][5] -= fijz;

			r[j][3] += fijx;
			r[j][4] += fijy;
			r[j][5] += fijz;

			return vij;

		}

		else return 0.0;

	}else if(!(m%nspacBeads) && !(n%nspacBeads)){ 

		// Interactions between spacer equal sizes beads
		if(ds2 <= rcut2[1]){

			ds2s = distances.inverseSquare(sig2[1],ds2);
			dsi6 = distances.inverseSquareCube(ds2s);

			vij = potentials.peWCA(eps[1],dsi6);
			fij = force.nonBonding(eps[1],ds2,dsi6);
			
			if(fabs(ds2) < 0.25){
				cout << "flj:fij = " << fij <<" vij = " << vij << " ds2 = " << ds2 << " ds2s = " << ds2s << endl;
				cout << "flj:dx = " << dx << " dy = " << dy << " dz = " << dz << "sig2LJ = " << sig2[0] << endl;
				cout << "flj:i = " << i << " j = " << j << endl;
				exit(EXIT_FAILURE);
			}

			fijx = fij*dx;
			fijy = fij*dy;
			fijz = fij*dz;

			r[i][3] -= fijx;
			r[i][4] -= fijy;
			r[i][5] -= fijz;

			r[j][3] += fijx;
			r[j][4] += fijy;
			r[j][5] += fijz;

			return vij;

		}

		else return 0.0;

	}else{ 

		// Interactions between spacer Unequal sizes beads
		if(ds2 <= rcut2[2]){

			ds2s = distances.inverseSquare(sig2[2],ds2);
			dsi6 = distances.inverseSquareCube(ds2s);

			vij = potentials.peWCA(eps[2],dsi6);
			fij = force.nonBonding(eps[2],ds2,dsi6);

			if(fabs(ds2) < 0.25){
				cout << "flj:fij = " << fij <<" vij = " << vij << " ds2 = " << ds2 << " ds2s = " << ds2s << endl;
				cout << "flj:dx = " << dx << " dy = " << dy << " dz = " << dz << "sig2LJ = " << sig2[0] << endl;
				cout << "flj:i = " << i << " j = " << j << endl;
				exit(EXIT_FAILURE);
			}

			fijx = fij*dx;
			fijy = fij*dy;
			fijz = fij*dz;

			r[i][3] -= fijx;
			r[i][4] -= fijy;
			r[i][5] -= fijz;

			r[j][3] += fijx;
			r[j][4] += fijy;
			r[j][5] += fijz;

			return vij;

		}

		else return 0.0;
	}

}
