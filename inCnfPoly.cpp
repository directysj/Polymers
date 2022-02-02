#include "interactions.h"
/* 
 * Program to generate initial configuration of system
 * of particles (e.g., linear polymers) based on random 
 * packing interpartice separation criteria.
 */
void genCnfPoly(int np, int npc, double blx, double bly,
		double blz, double **r, long *seed, long *seed1)
{
	
	Distances distances;

	int i,j,k,flg;
	double dr2,a,b,c,x1,x2,x3,sc,rx,ry,rz,ds2,
				 ran_x,ran_y,ran_z;
	long int idum,idum_new;

	idum = *seed;
	idum_new = *seed1;

	for(i = 0; i < np; i++){

		dr2 = 0.80;
		k = 0;
		do{
			flg = 1;

			if(i%npc == 0){

				a = ran3(&idum)*blx;
				b = ran3(&idum)*bly;
				c = ran3(&idum)*blz;

			}
			else
			{

				do{
					x1 = 2.0*ran3(&idum) - 1.0;
					x2 = 2.0*ran3(&idum) - 1.0;
					x3 = x1*x1 + x2*x2;

				}while(x3 > 1.0);

				c = 1.0 - 2.0*x3;
				x3 = sqrt(1.0 - x3);
				a = 2.0*x1*x3;
				b = 2.0*x2*x3;

				sc = 0.8 + 0.50*(2.0*ran3(&idum) - 1.0);  // offset + random number

        // Randomize the new bead position by adding/subtracting
				// from the previous bead position
				if(rint(ran3(&idum_new)) != 0) ran_x = -1.0;
				else ran_x = 1.0;

				if(rint(ran3(&idum_new)) != 0) ran_y = -1.0;
				else ran_y = 1.0;

				if(rint(ran3(&idum_new)) != 0) ran_z = -1.0;
				else ran_z = 1.0;

				a = r[i-1][0] + sc*a*ran_x;
				b = r[i-1][1] + sc*b*ran_y;
				c = r[i-1][2] + sc*c*ran_z;

			}

			for(j = i - 1; j > -1; j--){

				rx = r[j][0] - a;
				ry = r[j][1] - b;
				rz = r[j][2] - c;

        rx = distances.minimumImage(blx,rx);
        ry = distances.minimumImage(bly,ry);
        rz = distances.minimumImage(blz,rz);

				ds2 = distances.squareDistance(rx,ry,rz);

				flg = flg && ds2 > dr2;
			}

			// Change the flg if a bead is placed out of the central box
			if(a >= blx || b >= bly || c >= blz) flg = 0;
			if(a < 0.0 || b < 0.0 || c < 0.0) flg = 0;

			k++;

			if(k > 1000){
				dr2 *= dr2*0.95;
				k = 1;
			}
		}while(!flg);
		k = 1;

		r[i][0] = a;
		r[i][1] = b;
		r[i][2] = c;

		//	cout << i << endl; 
	}

	cout << "Succesfully generated initial configuration of linear polymers" << endl;

	fstream out;

	out.open("icrd.dat", ios::out);
	for(i = 0; i < np; i++) out << r[i][0] <<"  "<< r[i][1] <<"  "<< r[i][2] << endl;
	out.close();

}	
