#include "DataHandler.h"
#include <cmath>

DataHandler::DataHandler(const char* filename)
{
	fstream file(filename);

	if (!file)
	{
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}

	file >> num_targets;

	file >> depot_x >> depot_y;

	radii = new double[num_targets];
	min_reward_bounds = new double[num_targets];
	max_risk_boounds = new double[num_targets];
	center_x = new double[num_targets];
	center_y = new double[num_targets];

	for (int i = 0; i < num_targets; i++)
	{
		file >> radii[i];
	}
	for (int i = 0; i < num_targets; i++)
	{
		file >> min_reward_bounds[i];
	}

	for (int i = 0; i < num_targets; i++)
	{
		file >> max_risk_boounds[i];
	}

	for (int i = 0; i < num_targets; i++)
	{
		file >> center_x[i] >> center_y[i];
	}
}

void DataHandler::getTangentLines(int C1, int C2)
{
	double a, b, c, d, r1, r2;
	a = center_x[C1 - 1];
	b = center_y[C1 - 1];
	c = center_x[C2 - 1];
	d = center_y[C2 - 1];
	r1 = radii[C1 - 1];
	r2 = radii[C2 - 1];

	if (r1 == r2) {
		
		
	
	
	}


	if (radii[C1 - 1] < radii[C1 - 1]) {
		a = center_x[C2 - 1];
		b = center_y[C2 - 1];
		c = center_x[C1 - 1];
		d = center_y[C1 - 1];
		r1 = radii[C2 - 1];
		r2 = radii[C1 - 1];
	}

	// outer tangent point on bigger circle
	double xp = (c * r1 - a * r2) / (r1 - r2);
	double yp = (d * r1 - b * r2) / (r1 - r2);

	double tmp1 = sqrt(pow(xp - a, 2) + pow(yp - b, 2) + pow(r1, 2));
	double xOTL_C1_1 = (pow(r1, 2) * (xp - a) + r1 * (yp - b) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + a;
	double xOTL_C1_2 = (pow(r1, 2) * (xp - a) - r1 * (yp - b) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + a;
	double yOTL_C1_1 = (pow(r1, 2) * (yp - b) - r1 * (xp - a) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + b;
	double yOTL_C1_2 = (pow(r1, 2) * (yp - b) + r1 * (xp - a) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + b;
	if ((b - yOTL_C1_1)*(yp - yOTL_C1_1) != (xOTL_C1_1 - a)*(xOTL_C1_1 - xp)) {
		swap(xOTL_C1_1, xOTL_C1_2);
	}
	Point outTanPoint_Ci_1 = { xOTL_C1_1, yOTL_C1_1 };
	Point outTanPoint_Ci_2 = { xOTL_C1_2, yOTL_C1_2 };


	// outer tangent point on circle j
	tmp1 = sqrt(pow(xp - c, 2) + pow(yp - d, 2) + pow(r1, 2));
	double xOTL3 = (pow(r1, 2) * (xp - c) + r1 * (yp - d) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + c;
	double xOTL4 = (pow(r1, 2) * (xp - c) - r1 * (yp - d) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + c;
	double yOTL3 = (pow(r1, 2) * (yp - d) - r1 * (xp - c) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + d;
	double yOTL4 = (pow(r1, 2) * (yp - d) + r1 * (xp - c) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + d;

	if ((d - yOTL3)*(yp - yOTL3) != (xOTL3 - c)*(xOTL3 - xp)) {
		swap(xOTL3, xOTL4);
	}
	Point outTanPoint_C1_1 = { xOTL1, yOTL1 };
	Point outTanPoint_C1_2 = { xOTL2, yOTL2 };

	//inner tangent points on circle i
	xp = (c * r0 + a * r1) / (r0 + r1);
	yp = (d * r0 + b * r1) / (r0 + r1);

	tmp1 = sqrt(pow(xp - a, 2) + pow(yp - b, 2) - pow(r0, 2));
	
	double xITL1 = (pow(r0, 2) * (xp - a) + r0 * (yp - b) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + a;
	double xITL2 = (pow(r0, 2) * (xp - a) - r0 * (yp - b) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + a;
	double yITL1 = (pow(r0, 2) * (yp - b) - r0 * (xp - a) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + b;
	double yITL2 = (pow(r0, 2) * (yp - b) + r0 * (xp - a) * tmp1) / (pow(xp - a, 2) + pow(yp - b, 2)) + b;
	if ((b - yITL1)*(yp - yITL1) != (xITL1 - a)*(xITL1 - xp)) { 
		swap(xITL1, xITL2);
	}
	// inner tangent points on circle j
	tmp1 = sqrt(pow(xp - c, 2) + pow(yp - d, 2) + pow(r1, 2));
	double xITL3 = (pow(r1, 2) * (xp - c) + r1 * (yp - d) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + c;
	double xITL4 = (pow(r1, 2) * (xp - c) - r1 * (yp - d) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + c;
	double yITL3 = (pow(r1, 2) * (yp - d) - r1 * (xp - c) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + d;
	double yITL4 = (pow(r1, 2) * (yp - d) + r1 * (xp - c) * tmp1) / (pow(xp - c, 2) + pow(yp - d, 2)) + d;

	if ((d - yITL3)*(yp - yITL3) != (xITL3 - c)*(xITL3 - xp)) {
		swap(xITL3, xITL4);
	}



}

DataHandler::~DataHandler()
{
	delete radii;
	delete min_reward_bounds;
	delete max_risk_boounds;
	delete center_x;
	delete center_y;
}

void DataHandler::printGraphInfo()
{
	cout << "1. depot:  " << depot_x << " " << depot_y << endl;
	cout << "2. # of cirlces:  " << num_targets << endl;
	cout << "3. the lengths of radius are:  ";
	for (int i = 0; i < num_targets; i++)
		cout << radii[i] << "  ";

	cout << "\n";
	cout << "4. lower bounds of rewards:  ";
	for (int i = 0; i < num_targets; i++)
		cout << min_reward_bounds[i] << "  ";

	cout << "\n";
	cout << "5. upper bounds of risks:  ";
	for (int i = 0; i < num_targets; i++)
		cout << max_risk_boounds[i] << "  ";

	cout << "\n";
	cout << "6. coordinates of targets' locations are:  \n";
	for (int i = 0; i < num_targets; i++)
		cout << center_x[i] << "\t" << center_y[i] << "\n";
	cout << "\n";
}