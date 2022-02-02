// Program to store/generate restart file(s) while 
// equilibrating the system
void saveResInfoEq(double *tc, long *stp, int np, int fsize, 
		double **r, double dt)
{
	static int rsize,nr,k=0,m=0,n=0;
	static char fno[5],post[]=".bdt";
	static char rfname[20], rpre[]="trj/r";

	fstream filePointer;

	int i,j;

	// rsize - size of one record
	// fsize - file size in 256MB, 1024*1024*256
	// nr - number of records to be kept in the memory before writing
	// resi - restart index of the code, if restarted it will return
	// 1 otherwise 0

	if(k==0){	

		// size of one recored = size of 6n coordinates plus step number
		rsize = (int)sizeof(long) + (int)sizeof(double) + 6*np*(int)(sizeof(double));

		nr = fsize/rsize;
		cout << "trjrst.c: rsize = " << rsize << " nr= " << nr << endl;
		cout << "Save rastart INFO while equilibration" << endl;

		k++; // Increment the index

		// Check whether trj directory can be created
		i = system("if [ ! -d trj ]; then mkdir trj; fi");

		// If open restart flag file if it exists
		filePointer.open("trj/res.dat",ios::in);

		// If restart file does not exit 
		if(filePointer){

			filePointer >> i >> n >> m >> nr;
			filePointer.close();
			
			//Now reinstate the variables in case
			//of restart
			if(i == 1){
				//If restart is enabled and last file is 
				//complete no re-reading of saved trajectory is 
				//necessary
				//
				//Generate the file name where trajectory is stored
				sprintf(fno,"%d",n);
				strcpy(rfname,rpre);
				strcat(rfname,fno);
				strcat(rfname,post);
				cout << "trjrst.c: file to restore " << rfname << endl;

				//Open and restore current coordinates of simulation
				filePointer.open(rfname, ios::binary | ios::in);
				
				if(filePointer){
				filePointer.read((char*) stp,sizeof(long));
				filePointer.read((char*) tc,sizeof(double));
				for(i=0;i<np;i++) filePointer.read((char*) r[i],12*sizeof(double));
				filePointer.close();

				// Increment the step for starting from next step
				*stp++;
				*tc += dt;
			}else{
					cout << " trjrst.c:" << rfname << ":restart file not found exiting" << endl;
					exit(EXIT_FAILURE); //then exit
			}

				// start with a new storage file, increment its name and reset m
				n++;
				m = 0;
			}
			return;
		} else return;
	}

	// Increment to the next record
	m++;

	if(m == nr){

		// Restart file name
		sprintf(fno,"%d",n);
		strcpy(rfname,rpre);
		strcat(rfname,fno);
		strcat(rfname,post);
		cout << rfname << endl;

		// Record all coordinates in the file
		filePointer.open(rfname, ios::binary | ios::out);
		filePointer.write((char*) stp,sizeof(long));
		filePointer.write((char*) tc,sizeof(double));
		for(i=0;i<np;i++) filePointer.write((char*) r[i], 12*sizeof(double));
		filePointer.close();

		// Store restart flag file
		filePointer.open("trj/res.dat", ios::out);
		filePointer << 1 <<" "<< n <<" "<< m <<" "<< nr << endl;
		filePointer.close();

		//Increment the file name 
		n++;

		// Set the index of the storage memory to 0
		m = 0;

	}

	// Return the controlling variable back to main for restart 
	return;
}
