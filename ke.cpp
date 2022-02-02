double kinEn(int np, double **r)
{
	int i,j;
	double ke;

	ke=0.0;

#pragma omp parallel for simd firstprivate(np) private(j)\
	reduction(+:ke) schedule(dynamic)

	for(i=0;i<np;i++)for(j=3;j<6;j++) ke += r[i][j]*r[i][j];
	
	ke *= 0.50;

	return ke;
}
