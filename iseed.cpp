long int seedrv()
{
	int i,is,x;
	long int idum;

	fstream fp;
	/* Generate seed of the rand */
	fp.open("seed.dat",ios::in);
	if(fp){
		fp >> is ;
		fp.close();
	}else{
		cout << "ERROR: seed.dat file does not exist" << endl;
		exit(EXIT_FAILURE);
	}

	for(i=0;i<is;i++) x = rand();

	fp.open("seed.dat",ios::out);
	fp << x/10000 << endl;
	fp.close();
	
	idum = x;

	return idum;
}
