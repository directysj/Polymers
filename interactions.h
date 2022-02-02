#pragma once

/* Define class for potentials */
class Potentials
{

	double pLJ,pWCA,pBond;

	public:

	double peLJ(double eps, double dsi6)
	{
		pLJ = ((dsi6 - 1.00)*dsi6)*4.0*eps;
		return pLJ;
	}

	double peWCA(double eps, double dsi6)
	{
		pWCA = ((dsi6 - 1.00)*dsi6 + 0.250)*4.0*eps;
		return pWCA;
	}

	double peHarmonic(double springConstant, double bondDistance, double sigma)
	{
		pBond = springConstant*(bondDistance - sigma)*(bondDistance - sigma);
		return pBond;
	}

};

/* Define class for computing forces */
class Force
{
	double pair_force_nonBonding,pair_force_bond;

	public:

	double nonBonding(double eps, double ds2, double dsi6)
	{
		pair_force_nonBonding = 24.0*eps*(2.0*dsi6 - 1.00)*dsi6/ds2;
		return pair_force_nonBonding;
	}

	double harmonicBond(double springConstant, double sigmaByDistance)
	{
		pair_force_bond = springConstant*(sigmaByDistance - 1.0);
		return pair_force_bond;
	}
};

/* Define class for calculating distances */
class Distances
{

	double minImage,ds,ds2,dsi2,dsi6;

	public:

	double squareDistance(double dx, double dy, double dz)
	{
		ds2 = dx*dx + dy*dy + dz*dz;
		return ds2;
	}

	double inverseSquare(double sigma_square, double ds2)
	{
		dsi2 = sigma_square/ds2;
		return dsi2;
	}

	double inverseSquareCube(double dsi2)
	{
		dsi6 = dsi2*dsi2*dsi2;
		return dsi6;
	}

	double minimumImage(double boxLength, double ds)
	{
		minImage = ds - boxLength*rint(ds/boxLength);
		return minImage;
	}
};
