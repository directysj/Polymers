void removeRecreateTRJ()
{
	int i,j;
	// If enters into production runs, remove previously stored restart INFO files 
	cout << "Removing restart INFO files if succesfully enters into production runs" << endl;
	i = system("if [ -d trj ]; then rm -rf trj; fi");
	cout << "Removed restart INFO files" << endl;
	cout << "Recreate TRJ directory to store co-ordinates and restart INFO during production runs" << endl;
	j = system("if [ ! -d trj ]; then mkdir trj; fi");
}
