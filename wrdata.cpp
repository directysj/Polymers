#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

int main()
{

	int i,L = 5000,x,y;

	double rx,ry,rz;

	fstream out,inp;

	inp.open("conf.chains",ios::in);
	out.open("chains.out",ios::out);

	for(i=0;i<L;i++){
		inp >> x >> x >> x >> rx >> ry >> rz >> y >> y >> y;

		out << rx <<' '<< ry <<' '<< rz << endl;
	}

	inp.close();
	out.close();

	return 0;
}
