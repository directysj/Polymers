/*
 * Gaussian random deviates.
 *
 */


double gau(long int *idum)
{
    static int iset=0,i,j;
    static float gset;
    float fac,rsq,v1,v2;

    if (idum < 0) iset=0; 
    if (iset == 0){ 
        do {
            v1=2.0*ran3(idum)-1.0; 
            v2=2.0*ran3(idum)-1.0; 
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0); 

        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;

        iset=1; 
        return(v2*fac);

    } else { 
        iset=0; 
        return(gset); 
    }
}


void grand(int np, double  **c)
{
    static long idum;
    static int ctr=0;
    int x,is,i;

    if(ctr==0){
        idum=(long)time(NULL);
    }


    for(i=0;i<np;i++){
        c[i][3]=gau(&idum);
        c[i][4]=gau(&idum);
        c[i][5]=gau(&idum);
        //printf("%lf %lf %lf\n",c[i].x,c[i].y,c[i].z);
    }
}
