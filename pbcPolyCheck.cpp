// Function to calculate ened-to-end distance of 
// polymers and check if they cross PBCs
//
void pbcCheckPoly(int nc, int npc, double **r, 
		double Lx, double Ly, double Lz)
{
	int i,chain_end1,chain_end2;
	double dx,dy,dz;

#pragma omp parallel for firstprivate(nc,npc,Lx,Ly,Lz)\
	private(dx,dy,dz,chain_end1,chain_end2) schedule(dynamic)

	for(i = 0; i < nc; i++){

		// Indices of both ends of a chain
		chain_end1 = i*npc;
		chain_end2 = (i+1)*npc - 1;

		// Calculate end-to-end distance
		dx = r[chain_end1][0] - r[chain_end2][0];
		dy = r[chain_end1][1] - r[chain_end2][1];
		dz = r[chain_end1][2] - r[chain_end2][2];

		dx = sqrt(dx*dx); 
		dy = sqrt(dy*dy); 
		dz = sqrt(dz*dz); 

		if(dx > Lx/2.0 || dy > Ly/2.0 || dz > Lz/2.0){
			cout << "Polymer is too long for PBCs: Increase the box length(s)" << endl;
			cout << "chain index = " << i << " dx = " << dx <<" dy = " << dy <<" dz = " << dz << endl;
			exit(EXIT_FAILURE);
		} 
	}

}
