/* 
 * ***************************************************************
 * Underdamped Langevin dynamics simulation of linear homopolymers
 * ***************************************************************
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <omp.h>

#define icell(i,j,k,m,n,o) (i+m)%m+((j+n)%n)*m+((k+o)%o)*n*m

using namespace std;

#include "simparin.cpp"
#include "ran3.cpp"
#include "gencell.cpp"
#include "gauss.cpp"
#include "iseed.cpp"
#include "inCnfPoly.cpp"
#include "pbcPolyCheck.cpp"
#include "volFracBeads.cpp"
#include "veldistb.cpp"
#include "beadIndex.cpp"
#include "bondHarm.cpp"
#include "fljp.cpp"
#include "forcecellmp.cpp"
#include "gausranlang.cpp"
#include "mv_vvlt_ld.cpp"
#include "ke.cpp"
#include "removeRecreateTrj.cpp"
#include "resEqSim.cpp"
#include "trjrst.cpp"
#include "grs.cpp"

int main()
{
	/*
	 * ns - number of steps
	 * tf - Final time of the simulation
	 * nsi - number of sampling interval
	 * si - sampling interval in time
	 * neq - equilibration steps
	 * teq - equilibration time
	 * np - number of particles
	 * fsize - Size of one trajectory file in MB
	 */ 
	int j,n,flg,nsi,np,nc,npc,neq,
			nav,dim,nap,*api,res,mt,nspacBeads,
			state,fsize,grCal,nStikers,nSpacers;

	long int i,ns,nsum_runvar,seed,seed1;

	/*
	 * te - temperature
	 * dt - time step
	 * Lx, Ly, Lz - box lengths 
	 * rcut - P.E. cut off
	 * rcld - Cutoff radius for the active particles
	 */

	double si,siRDF,te,dt,tf,teq,af,pe,ke,etot,tinst; 
	double pesum,kesum,etotsum,tmpsum,ptaf,
				 pesqsum,kesqsum,etotsqsum,tmpsqsum;
	double Lx,Ly,Lz,volume,phi_st,phi_sp;
	double peav,keav,etotav,tmpav,pesqav,kesqav,
				 etotsqav,tmpsqav,flpe,flke,fletot,
				 fltmp,napsum,napsqsum,napav,napsqav,
				 flnap;
	/*
	 * nn - postions 
	 * crd - no. coordinate of each particle
	 * position r[..,0-2],velocity r[..,3-5],
	 * force r[..,6-8]
	 * two uncorrelated random distributions
	 * zet[i][0-2] and tht[i][0-2]
	 * eps - Strength of Lennard-Jones interactions
	 * celLenList - Cell length for force calculation
	 * Index position of spacer beads along a polymer : nspacBeads
	 */
	double **r,tc,*sig2,*sigma,*rcut2,*eps,**fa,sigLJ,sigWCA;
	double epsLJ,epsWCA,epsCI,**zet,**tht,rcLJ,rcWCA,gama,kb,celLenList;

	fstream inp,sys,pos,out,stat;

	/* Read input parameters for simulation(s) */
	inp.open("inp.dat",ios::in);
	sys.open("system.dat",ios::in);
	if(inp){
		inp >> nav >> epsWCA >> epsLJ >> si >> siRDF >> tf >> teq >> te >> dt >> gama >> rcLJ >> kb;
		sys >> nc >> npc >> dim >> nspacBeads >> Lx >> Ly >> Lz >> sigLJ >> sigWCA >> celLenList;
		inp.close();
		sys.close();
	}else{
		cout << "inp.dat - system.dat:: ERROR in reading input simulation parameters" << endl;
		exit(EXIT_FAILURE);
	}

	// Read state of simulation
	stat.open("state.sim",ios::in); 
	if(stat){
		stat >> state >> fsize;
		stat.close();
	}else{
		cout << "state.sim:: simulation state could not be read" << endl;
		exit(EXIT_FAILURE);
	}

	sigma = new double[3];
	sig2 = new double[3];
	rcut2 = new double[3];
	eps = new double[3];

	// Total number of particles in the system to be simulated
	np = nc*npc;

	// Volume of the central simulation box
	volume = Lx*Ly*Lz;

	// Calculate volume fraction of sticker and spacer beads
	phi_st = fracBeads(nc,npc,nspacBeads,volume,sigLJ,sigWCA,&phi_sp,&nStikers,&nSpacers);

	rcWCA = pow(2.0,1.0/6.0); // Cut-off for WCA potential

	// Get and assemble simulation parameters
	simparinit(sigLJ,sigWCA,sig2,epsWCA,epsLJ,epsCI,rcLJ,rcWCA,rcut2,eps,sigma);

	if(celLenList < sqrt(rcut2[0])){
		cout << "Cell length for non-Bonding force calculation is smaller that RCUT" << endl;
		cout << "Rerun after incresing the cell length in file system.dat" << endl;
		exit(EXIT_FAILURE);
	} 

	//api = new int[npc];
	//nnb = new int[npc];
	r = new double*[np];
	//fa = new double*[npc];
	zet = new double*[np];
	tht = new double*[np];
	for(i=0;i<np;i++){
		r[i] = new double[4*dim];
		//fa[i] = new double[dim];
		zet[i] = new double[dim];
		tht[i] = new double[dim];
	}

	inp.open("checkParSim.dat",ios::out);
	inp << "Temperature = " << te << endl;
	inp << "chains = " << nc << " beads in chain = " << npc << " total particles = " << np << endl;
	inp << "Lx = " << Lx <<  " Ly = " << Ly << " Lz = " << Lz << endl;
	inp << "Time step = " << dt << " Final time = " << tf << " sampling interval = " << si << endl;
	inp << "LJ CutOff = " << rcLJ << " WCA CutOff = " << rcWCA << endl;
	inp << "Square of LJ CutOff = " << rcut2[0] << " WCA CutOff = " << rcut2[1] << endl;
	inp << "volume fraction of stickers = " << phi_st << " and spacers = " << phi_sp << endl;
	inp << "sig_LJ = " << sigma[0] << " sig_WCA = " << sigma[1] << " sig_LB_MIX = " << sigma[2] << endl;
	inp << "eps_LJ = " << eps[0] << " eps_WCA = " << eps[1] <<  " eps_LB_MIX = " << eps[2] << endl;
	inp.close();

	// Generate two seeds for uncorrelated Gaussian random variables
	seed = seedrv();
	seed1 = seedrv();

	// Generate/Load initial configuration of the system 

	if(!state) genCnfPoly(np,npc,Lx,Ly,Lz,r,&seed,&seed1);
	else if(state == 1){
		pos.open("crdOutLast.dat",ios::in);
		if(pos){
			for(i=0;i<np;i++) pos >> r[i][0] >> r[i][1] >> r[i][2] >> r[i][3] >> r[i][4] >> r[i][5];
			pos.close();
			cout << "Initial system configuration loaded successfully" << endl;
		}else{
			cout << "Initial system configuration could not be loaded/generated" << endl;
			exit(EXIT_FAILURE);
		}
	}else{
		pos.open("icrd.dat", ios::in);
		if(pos){
			for(i=0;i<np;i++) pos >> r[i][0] >> r[i][1] >> r[i][2];
			pos.close();
			cout << "Initial system configuration loaded successfully" << endl;
		}else{
			cout << "Initial system configuration could not be loaded/generated" << endl;
			exit(EXIT_FAILURE);
		}
	}

	// Random Gaussian distribution of velocities
	veldismb(np,r,te,dim,&seed,state);

	nsi = (int)(si/dt); // No of sampling intervals
	grCal = (int)(siRDF/dt); // Sampling intervals for g(r) calculations
	neq = (int)(teq/dt); // No of equilibration time
	tf += teq;

	nsum_runvar = 0; // Running variable summation number
	tc = 0.0; // Current time initialization
	i = 0; // Step number initialization

	pesum = kesum = etotsum = tmpsum = napsum = 0.0;
	pesqsum = kesqsum = etotsqsum = tmpsqsum = napsqsum = 0.0;

	// Restart the code if required that is decided by 0/1 and state of simulation
	// Here state -> 2 means if system fails while equilibration runs are going on
	// ELSE state -> 0 & 1 applied 
	if(state == 2) saveResInfoEq(&tc,&i,np,fsize,r,dt);
	else traj(&tc,&i,np,r,0,fsize,dt);

	/***********************************************************************/
	/********************* Main Loop of the Simulation *********************/
	/***********************************************************************/

	while(tc <= tf){


		if(!i){
			// Compute forces using CELL list
			pe = forcecell(np,Lx,Ly,Lz,r,npc,nc,rcut2,sig2,kb,dim,eps,nspacBeads,sigma,celLenList);
			cout << "Initial potential energy per particle = " << pe/(double)np << endl;
			cout << "Steps    Time    P.E.   K.E.   Energy   Temperature" << endl;
		}

		/*********** Gaussian random deviates for the displacement ***********/
		// Random numbers Zeta(t) for each particle
		grand_disp(np,dim,zet,&seed);

		/*********** Gaussian random deviates for the displacement ***********/
		grand_disp(np,dim,tht,&seed1);

		// Move the particles using Velocity verlet with Langevin noise
		pe = move(np,r,zet,tht,dt,te,Lx,Ly,Lz,npc,nc,rcut2,
				sig2,&seed1,gama,kb,dim,eps,nspacBeads,sigma,celLenList);

		// Calculate kinetic energy of the system
		ke = kinEn(np,r);

		// Instantaneous temperature
		tinst = 2.0*ke/(3.0*((double)np - 1.0));

		// Per particle K.E.
		ke /= (double)np;

		// Total energy of the system
		etot = ke + pe;

		//if(i%100 == 0) cout <<" "<< i <<" "<< (double)i*dt <<" "<< pe <<" "<< ke <<" "<< etot <<" "<< tinst <<" "<< endl;

		// Storing trajectory and calculate correlations on FLY
		if(i > neq){ 

			// If enters into production runs, remove restart
			// INFO files stored during equilibration

			if((i - neq) == 1) removeRecreateTRJ();

			if(!(i%nsi)){

				// store trajectories
				traj(&tc,&i,np,r,0,fsize,dt);

				// Check if POLYMER BEADS interact with themselves through
				// their periodic images --> WHICH IS UNPHYSICAL
				pbcCheckPoly(nc,npc,r,Lx,Ly,Lz);
			}

			// calculate radial distribution functions
			if(!(i%grCal)) rdf(Lx,Ly,Lz,np,r,dim,nspacBeads,nStikers,nSpacers,npc);

			// Accumulating running variables
			pesum = pesum + pe;
			kesum = kesum + ke;
			etotsum = etotsum + etot;
			tmpsum = tmpsum + tinst;

			pesqsum = pesqsum + pe*pe;
			kesqsum = kesqsum + ke*ke;
			etotsqsum = etotsqsum + etot*etot;
			tmpsqsum = tmpsqsum + tinst*tinst;

			nsum_runvar += 1;

			if(!(i%nav)){

				peav = pesum/(double)nsum_runvar;
				keav = kesum/(double)nsum_runvar;
				etotav = etotsum/(double)nsum_runvar;
				tmpav = tmpsum/(double)nsum_runvar;

				pesqav = pesqsum/(double)nsum_runvar;
				kesqav = kesqsum/(double)nsum_runvar;
				etotsqav = etotsqsum/(double)nsum_runvar;
				tmpsqav = tmpsqsum/(double)nsum_runvar;

				flpe = sqrt(fabs(pesqav - peav*peav));
				flke = sqrt(fabs(kesqav - keav*keav));
				fletot = sqrt(fabs(etotsqav - etotav*etotav));
				fltmp = sqrt(fabs(tmpsqav - tmpav*tmpav));

				cout << "After " << i << " steps" << endl;
				cout << "Averaged over " << nsum_runvar << " steps" << endl;

				cout << "<PE> = " << peav << " +/-" << flpe << endl;
				cout << "<KE> = " << keav << " +/-" << flke << endl;
				cout << "<Etot> = " << etotav << " +/-" << fletot << endl;
				cout << "<T> = " << tmpav << " +/-" << fltmp << endl;

			}
		}else{
			// Save restart info and co-ordinates while equilibrating the system
			if(!(i%nsi)) saveResInfoEq(&tc,&i,np,fsize,r,dt);
		}

		tc += dt; // Update current time
		i++;  // Update step number

	} 

	// store trajectories for last step numbers those could not
	// be saved because of maximum no of records of a file
	traj(&tc,&i,np,r,1,fsize,dt);

	/************************************************************************/
	/****** DO NOT APPLY PBCs ON THE ACTUAL POSITIONS OF POLYMER BEADS ******/
	/*********** BECAUSE IT WILL BREAK COVALENT (HARMONIC BONDS) ************/
	/************************************************************************/
	// Print final configuration of the system
	out.open("crdOutLast.dat", ios::out);

	for(j = 0; j < np; j++) 
		out << r[j][0] <<" "<< r[j][1] <<" "<< r[j][2] <<" "<< r[j][3] <<" "<< r[j][4] <<" "<< r[j][5] << endl;

	out.close();

	delete [] r;
	delete [] zet;
	delete [] tht;
	delete [] eps;
	delete [] sigma;
	delete [] sig2;
	delete [] rcut2;

	return 0;
}
