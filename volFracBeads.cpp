double fracBeads(int nc, int npc, int nspacBeads, 
		double volume, double sigLJ, double sigWCA, 
		double *frac_WCA, int *nStickers, int *nSpacers)
{

	int i,j,nspac,nstik;
	double pi,frac_LJ;

	pi = 4.0*atan(1.0);

	nspac = nstik = 0;

	// Count sticker and spacer beads along a single chain

	for(i=0;i<npc;i++){
		if(i%nspacBeads == 0) nspac += 1;
		else nstik += 1;
	}

	// Multiply by number of chains in the system

	nspac = nspac*nc; // Total number of spacers in the system
	nstik = nstik*nc; // Total number of stickers in the system

	*nSpacers = nspac;
	*nStickers = nstik;

	*frac_WCA = nspac*sigWCA*sigWCA*sigWCA*pi/6.0/volume;
	frac_LJ = nstik*sigLJ*sigLJ*sigLJ*pi/6.0/volume;

	return frac_LJ;
}
