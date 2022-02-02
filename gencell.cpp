// Generate cells 
// Define cell neighbours
// generate a 3d map cell linear index from 3d integer index

void gencell(int ncx, int ncy, int ncz, int **cn)
{
	int i,j,k,l;

	l=0;

	for(i=0;i<ncx;i++){
		for(j=0;j<ncy;j++){
			for(k=0;k<ncz;k++){

				// Find the index of the 13 out of 26 neighbours

				// l - give the unique position of cell marked by
				// three indices (i,j,k)

				l = icell(i,j,k,ncx,ncy,ncz);

				// 13 neighbours of each central cell
				// throughout Simulation Box

				cn[l][0] = icell(i+1,j,k,ncx,ncy,ncz);
				cn[l][1] = icell(i+1,j+1,k,ncx,ncy,ncz);
				cn[l][2] = icell(i,j+1,k,ncx,ncy,ncz);
				cn[l][3] = icell(i-1,j+1,k,ncx,ncy,ncz);
				cn[l][4] = icell(i+1,j,k-1,ncx,ncy,ncz);
				cn[l][5] = icell(i+1,j+1,k-1,ncx,ncy,ncz);
				cn[l][6] = icell(i,j+1,k-1,ncx,ncy,ncz);
				cn[l][7] = icell(i-1,j+1,k-1,ncx,ncy,ncz);
				cn[l][8] = icell(i+1,j,k+1,ncx,ncy,ncz);
				cn[l][9] = icell(i+1,j+1,k+1,ncx,ncy,ncz);
				cn[l][10] = icell(i,j+1,k+1,ncx,ncy,ncz);
				cn[l][11] = icell(i-1,j+1,k+1,ncx,ncy,ncz);
				cn[l][12] = icell(i,j,k+1,ncx,ncy,ncz);
			}
		}
	}
}
