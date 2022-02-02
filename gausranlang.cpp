/* Give random displacements according to Maxwell-Boltzmann distributions */
void grand_disp(int np, int dim, double **r, long int *iseed)
{
	int i,j;

	for(i=0;i<np;i++)for(j=0;j<dim;j++) r[i][j] = gaus(iseed);

}
