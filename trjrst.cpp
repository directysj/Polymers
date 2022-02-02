// Program to store trajectory for further analysis and generate restart file
void traj(double *tc, long *stp, int np, double **r,
		int term, int fsize, double dt)
{
	static int rsize,nr,k=0,m=0,n=0;
	static long *st;
	static double *ts;
	static double ***g;
	static char fname[20], fno[5],pre[]="trj/t",post[]=".bdt";
	static char rfname[20], rpre[]="trj/r";

	int i,j;
	fstream inp,chkf;

	// rsize - size of one record
	// fsize - file size in MB
	// nr - number of records to be kept in the memory before writing
	// resi - restart index of the code, if restarted it will return
	// 1 otherwise 0

	//printf("term = %d",term);

	if(!k){	

		//size of one recored = size of 6n coordinates plus step number
		rsize = (int)sizeof(long) + (int)sizeof(double) + 6*np*(int)(sizeof(double));

		nr = fsize/rsize;
		cout << "trjrst.c: rsize = " << rsize << " nr = " << nr << endl;

		//Create the directory trj

		//printf("trjrst.c:%d %d %d\n",nr,fsize,rsize);

		//Allocate space for storing the coordinates
		//
		g = new double**[nr];
		st = new long[nr];
		ts = new double [nr];
		for(i=0;i<nr;i++){
			g[i] = new double*[np];
			for(j=0;j<np;j++) g[i][j] = new double[6];
		}

		k++; // Increment the index

		// Check whether trj directory can be created
		i = system("if [ ! -d trj ]; then mkdir trj; fi");

		// If open restart flag file if it exists
		inp.open("trj/res.dat",ios::in);

		// If restart file does not exit 
		if(inp){
			inp >> i >> n >> m >> nr;
			inp.close();
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
				inp.open(rfname, ios::binary | ios::in);
				if(inp){
					inp.read((char*) stp,sizeof(long));
					inp.read((char*) tc,sizeof(double));
					for(i=0;i<np;i++) inp.read((char*) r[i],12*sizeof(double));
					inp.close();

					//Increment the step for starting from next step
					*stp++;
					*tc += dt;
				}else{
					cout << " trjrst.c: " << fname << ":restart file not found exiting" << endl;
					exit(EXIT_FAILURE);//then exit
				}

				if(m < nr){
					//Read in the trajectory buffer till it was
					//loaded earlier
					//
					sprintf(fno,"%d",n);
					strcpy(fname,pre);
					strcat(fname,fno);
					strcat(fname,post);

					cout << "trjrst.c: storing " << fname << endl;
					inp.open(fname, ios::binary | ios::in);
					for(i=0;i<m;i++){
						//Read the step number of the simulation
						inp.read((char*) &st[i],sizeof(long));
						//Read current time of the simulation
						inp.read((char*) &ts[i],sizeof(double));
						//Read the  phase space point coordinates + velocities from the file
						for(j=0;j<np;j++) inp.read((char*) g[i][j],6*sizeof(double));
					}

					inp.close();

				}else{

					//Else start with a new storage file
					//increment file name and reset m
					n++;
					m = 0;
				}
			}
			return;
		}else return;
	}

	//Store the phase space point to memory

	for(i=0;i<np;i++)for(j=0;j<6;j++) g[m][i][j] = (double)r[i][j];
	st[m] = *stp;
	ts[m] = *tc;

	//Increment to the next record
	m++;

	if(m == nr || term == 1 ){	

		//Generate the file name where trajectory to be stored
		sprintf(fno,"%d",n);
		strcpy(fname,pre);
		strcat(fname,fno);
		strcat(fname,post);

		cout << "trjrst.c: storing " << fname << endl;
		//Open the file to store file name
		//
		inp.open(fname,ios::binary | ios::out);
		chkf.open("out.dat",ios::app);
		for(i=0;i<m;i++){
			//Write the step number of the simulation
			inp.write((char*) &st[i],sizeof(long));

			//Write current time of the simulation
			inp.write((char*) &ts[i],sizeof(double));

			//Write the  phase space point coordinates + velocities to the file
			for(j=0;j<np;j++) inp.write((char*) g[i][j],6*sizeof(double));

			// Apppend step number for the checking purpose
			chkf << "step: " << st[i] << endl;

		}

		inp.close();
		chkf.close();

		//Section to  save restart file

		//Restart file name
		sprintf(fno,"%d",n);
		strcpy(rfname,rpre);
		strcat(rfname,fno);
		strcat(rfname,post);
		cout << rfname << endl;

		//Record all coordinates in the file
		inp.open(rfname, ios::binary | ios::out);
		inp.write((char*) stp,sizeof(long));
		inp.write((char*) tc,sizeof(double));
		for(i=0;i<np;i++) inp.write((char*) r[i],12*sizeof(double));
		inp.close();

		//Store restart flag file
		inp.open("trj/res.dat", ios::out);
		inp << 1 <<" "<< n <<" "<< m <<" "<< nr << endl;
		inp.close();

		//Increment the file name 
		n++;

		//Set the index of the storage memory to 0
		m = 0;

	}

	// Return the controlling variable back to main for restart 
	return;
}
