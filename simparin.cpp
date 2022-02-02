void simparinit(double sigLJ, double sigWCA, double *sig2,
		double epsWCA, double epsLJ, double epsCI, double rcutLJ,
		double rcutWCA, double *rcut2, double *eps, double *sigma)
{

	double rcLJ,rcWCA,rcCI,sigCI;

	/*
	 * Cross interation LJ-WCA parameters are according to 
	 * the Lorentz-Bertholate (LB) rule
	 */

	// sigma_LJ-WCA particle diameters
	sigCI = 0.50*(sigLJ + sigWCA); 

	sigma[0] = sigLJ;
	sigma[1] = sigWCA;
	sigma[2] = sigCI;
	
	sig2[0] = sigLJ*sigLJ;
	sig2[1] = sigWCA*sigWCA;
	sig2[2] = sigCI*sigCI;

	rcLJ = rcutLJ*sigLJ;
	rcWCA = rcutWCA*sigWCA;

	// Mutual interaction between particles are repulsive,
	// which is a special case and must be modified according
	// to needs 
	rcCI = rcutWCA*sigCI;

	cout << "Warning: Mutual interaction between particles are repulsive (a special case), which must be modified as per needs" << endl;

	rcut2[0] = rcLJ*rcLJ;
	rcut2[1] = rcWCA*rcWCA;
	rcut2[2] = rcCI*rcCI;

	// Cross-interaction energies
	epsCI = sqrt(epsWCA*epsLJ);

	eps[0] = epsLJ;
	eps[1] = epsWCA;
	eps[2] = epsCI;

}
