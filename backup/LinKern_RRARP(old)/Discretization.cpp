#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <math.h>
#include "Discretization.h"
#include "DataHandler.h"
#include "CostInfo.h"

using namespace std;
#define PI  3.1415926
#define INF numeric_limits<double>::infinity()


bool comparator(const pair<double, INDEX>& l, const pair<double, INDEX>& r) {
	return l.first < r.first;
}

Discretization::Discretization(int k, DataHandler *graph_arg, CostInfo* ci)
{

	this->CI = ci;
	para_ext_risk = 1;
	para_int_risk = 1;

	this->num_dstzn = k;
	this->num_targets = graph_arg->get_num_targets();

	depot_x = graph_arg->get_depot_x();
	depot_y = graph_arg->get_depot_y();

	this->center_x = graph_arg->get_center_x();
	this->center_y = graph_arg->get_center_y();

	subarc_angle = 2.0 * PI / k;
	num_nodes_sglIdx = num_targets * num_dstzn + 1;

	radii.push_back(0);
	for (int i = 0; i < num_targets; i++)
		radii.push_back(graph_arg->get_radii()[i]);
	radii.push_back(0);

	min_reward_bounds.push_back(0);
	for (int i = 0; i < num_targets; i++)
		min_reward_bounds.push_back(graph_arg->get_min_rewards()[i]);
	min_reward_bounds.push_back(0);

	max_risk_bounds.push_back(0);
	for (int i = 0; i < num_targets; i++)
		max_risk_bounds.push_back(graph_arg->get_max_risks()[i]);
	max_risk_bounds.push_back(0);
	
	//------------distance matrix with single-index--------------	
	nodes.resize(num_targets + 2); // coordinates of discretized points
	sglIdxMat.resize(num_nodes_sglIdx, vector<double>(num_nodes_sglIdx)); // single-index distance matrix
	intDist.resize(num_targets + 1, vector<double>(num_dstzn, -1)); // only interior distance vector 

	build();	
}

void Discretization::generate_nodes_crds()
{
	//1. first node: depot
	nodes[0].resize(1);
	nodes[0][0].x = depot_x;
	nodes[0][0].y = depot_y;

	//2. nodes on circles
	for (int pos = 1; pos <= num_targets; pos++) {
		nodes[pos].resize(num_dstzn);
		for (int i = 0; i < num_dstzn; i++)	{
			nodes[pos][i].x = radii[pos] * cos(i * subarc_angle) + center_x[pos - 1];
			nodes[pos][i].y = radii[pos] * sin(i * subarc_angle) + center_y[pos - 1];
		}
	}

	// 3. last node: depot
	nodes[num_targets + 1].resize(1);
	nodes[num_targets + 1][0].x = depot_x;
	nodes[num_targets + 1][0].y = depot_y;
}

void Discretization::cal_ext_dist()
{
	int s, t, i, j, idxmat_s, idxmat_t;
	double w = 0.0;
	//source <--> other nodes on circles
	for (t = 1; t <= num_targets; t++) {
		idxmat_t = (t - 1) * num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			w = sqrt(pow(nodes[t][i].x - nodes[0][0].x, 2) + pow(nodes[t][i].y - nodes[0][0].y, 2));
			sglIdxMat[0][idxmat_t + i] = w; // source is set to be connected with all entries 
			sglIdxMat[idxmat_t + i][0] = w;
		}
	}

	// s-th circle <--> t-th circle
	for (s = 1; s <= num_targets; s++){ // [target s]
		idxmat_s = (s - 1) * num_dstzn + 1; 
		for (i = 0; i < num_dstzn; i++)	{ // [target s, point i]
			for (t = 1; t <= num_targets; t++){ // [target t]
				if (s != t)	{
					idxmat_t = (t - 1) * num_dstzn + 1; // all target t's entries
					for (j = 0; j < num_dstzn; j++)	{
						w = sqrt(pow(nodes[s][i].x - nodes[t][j].x, 2) + pow(nodes[s][i].y - nodes[t][j].y, 2));
						sglIdxMat[idxmat_s + i][idxmat_t + j] = w; // target s' exit i points to target t's entry j;
						sglIdxMat[idxmat_t + j][idxmat_s + i] = w; // target t's exit j  points to target s' entry i;
					}
				}
			}
		}
	}

	/*
	for (i = 0; i < num_nodes_sglIdx; i++)
	{
		for (j = 0; j < num_nodes_sglIdx; j++)
		{
			cout << sglIdxMat[i][j] << '\t';
		}
		cout << "\n";
	}
	cout << "-----------------------------------------------" << endl;
	*/
}

void Discretization::cal_int_dist()
{
	int s, i, j, flag, idx;
	double w = 0.0;
	double risk_tmp, reward_tmp;
	vector<double> vec_risk(num_dstzn, -1);

	for (s = 1; s <= num_targets; s++) {
		// compute risk or reward for each interior path
		for (i = 0; i < num_dstzn; i++) {
			w = sqrt(pow(nodes[s][0].x - nodes[s][i].x, 2) + pow(nodes[s][0].y - nodes[s][i].y, 2));
			if (w > 0.000001) {
				risk_tmp = w;
				reward_tmp = w;
				if (risk_tmp < max_risk_bounds[s] && reward_tmp >= min_reward_bounds[s]) {
					vec_risk[i] = risk_tmp; // optimize later 					
				}
				else {
					vec_risk[i] = INF; // too risky or extremely less reward									  
				}
			}
			else {
				vec_risk[i] = INF;
			}
		}

		// interior distance
		intDist[s] = vec_risk;

		idx = (s - 1) * num_dstzn + 1;
		for (i = 0; i < num_dstzn; i++) {
			for (j = 0; j < num_dstzn; j++) {
				flag = (num_dstzn - i + j) % num_dstzn;
				sglIdxMat[idx + i][idx + j] = vec_risk[flag];
			}
		}
	}
	/*
	for (i = 0; i < num_nodes_sglIdx; i++)
	{
		for (j = 0; j < num_nodes_sglIdx; j++)
		{
			cout << sglIdxMat[i][j] << '\t';
		}
		cout << "\n";
	}
	cout << "-----------------------------------------------" << endl;
	*/
}

double Discretization::evaluate_risk(double radius, double w)
{
	/*
	double gama = 1.0;
	double R = 0.0;
	double a = w / 2;
	double d = sqrt(pow(radius, 2) - pow(a, 2));
	R = 2 * gama * (a*d / (pow(a, 2) + pow(d, 2)) + atan(a / d)) / pow(d, 3);
	*/

	double c1 = 1;
	double c2 = 1.0;

	double a = w / 2;
	double d = sqrt(pow(radius, 2) - pow(a, 2));
	double R = 0.0;
	double delta_t = 0.005;
	double t = 0;
	while (t <= a)
	{
		R = R + c1 * exp(-(t*t + d*d) / c2);
		t = t + delta_t;
	}
	return R;
}

double Discretization::evalute_reward(double radius, double w)
{
	double para_w = 1;
	double para_k = 1.0;

	double a = w / 2;
	double d = sqrt(pow(radius, 2) - pow(a, 2));

	double T = 0.0;
	double delta_t = 0.005;
	double t = 0.0;
	while (t <= a)
	{
		T = T + para_w * log2(1 + pow(para_k, 4) / (pow(d*d + t*t, 2) + 1));
		t = t + delta_t;
	}

	T = 2 * T;
	return T;

}

void Discretization::build() {
	generate_nodes_crds();
	cal_ext_dist();
	cal_int_dist();
}

Discretization::~Discretization(){

}

void Discretization::fast_sp2(vector<VERTEX>& spTour, double& optimal, vector<int> &seq)
{
	//SDS format (SD, short for shortest dist， shortest distances for all entries to source):
	//0: 0
	//1: 1-k
	//2: 1-k
	//3: ...
	//targets+1: 
	vector<vector<double>> SDS(num_targets + 2);

	int i, idxmat, idxmat_fr, idxmat_lat;

	vector<double> entry_dist(num_dstzn, INF);
	vector<double> exit_dist(num_dstzn, INF);
	vector<double> INF_dist(num_dstzn, INF);

	optimal = INF;
	double dist = 0.0;

	vector<vector<int>> backwardIndex(2 * num_targets + 2);
	backwardIndex[0].resize(1);
	backwardIndex[0][0] = -1;


	//1. depot: SDS is in default set to be 0.0.
	SDS[0].resize(1);
	SDS[0][0] = 0.0;

	//2: Entries on 1st circle
	int pos = 1;
	int id_circle = seq[pos];

	SDS[1].resize(num_dstzn);
	backwardIndex[1].resize(num_dstzn);
	backwardIndex[2].resize(num_dstzn);
	idxmat = (id_circle - 1) * num_dstzn + 1; // 1st entry at first target
	for (i = 0; i < num_dstzn; i++)	{
		entry_dist[i] = sglIdxMat[0][idxmat + i];
		SDS[1][i] = entry_dist[i];

		backwardIndex[1][i] = 0;
	}

	//2: Exits on 1st circle
	for (int i = 0; i < num_dstzn; i++)	{
		for (int j = 0; j < num_dstzn; j++)	{
			dist = entry_dist[i] + sglIdxMat[idxmat + i][idxmat + j];
			if (dist < exit_dist[j]) {
				exit_dist[j] = dist;
				backwardIndex[2][j] = i;
			}

		}
	}

	entry_dist = INF_dist;

	//3: nodes from seq[2] to seq[num_circles]
	int id_fr_circle = 0;
	int id_lat_circle = 0;
	for (pos = 2; pos <= num_targets; pos++) {// circle # seq[2], seq[3], seq[4], seq[5], seq[6]...
		SDS[pos].resize(num_dstzn);
		backwardIndex[2 * (pos - 1) + 1].resize(num_dstzn);
		backwardIndex[2 * pos].resize(num_dstzn);

		id_fr_circle = seq[pos - 1]; // the idx of front circle
		id_lat_circle = seq[pos]; // the idx of latter circle	

		idxmat_fr = (id_fr_circle - 1) * num_dstzn + 1;
		idxmat_lat = (id_lat_circle - 1) * num_dstzn + 1;

		//update entry
		for (int i = 0; i < num_dstzn; i++)	{
			for (int j = 0; j < num_dstzn; j++)	{
				dist = exit_dist[i] + sglIdxMat[idxmat_fr + i][idxmat_lat + j];
				if (dist < entry_dist[j]) {
					entry_dist[j] = dist;
					backwardIndex[2 * (pos - 1) + 1][j] = i;
				}
			}
		}

		// update SDS
		for (int i = 0; i < num_dstzn; i++)	{
			SDS[pos][i] = entry_dist[i];
		}

		exit_dist = INF_dist;
		// update exits of latter target
		for (int i = 0; i < num_dstzn; i++)	{
			for (int j = 0; j < num_dstzn; j++)	{
				dist = entry_dist[i] + sglIdxMat[idxmat_lat + i][idxmat_lat + j];
				if (dist < exit_dist[j]) {
					exit_dist[j] = dist;
					backwardIndex[2 * pos][j] = i;
				}

			}
		}

		entry_dist = INF_dist;
	}
	
	// last target
	id_circle = seq[num_targets];

	SDS[num_targets + 1].resize(1);
	backwardIndex[2 * num_targets + 1].resize(1);
	idxmat = (id_circle - 1) * num_dstzn + 1;
	for (int i = 0; i < num_dstzn; i++)	{
		dist = exit_dist[i] + sglIdxMat[idxmat + i][0];
		if (dist < optimal) {
			optimal = dist;
			backwardIndex[2 * num_targets + 1][0] = i;
		}

	}

	SDS[num_targets + 1][0] = optimal; // shortest distance

	// -------------------------------------------------------------------------------------------------------------
	// finding the shortest path by backtracking
	int lastIdx = backwardIndex[2 * num_targets + 1][0];
	vector<int> pathVec(2*num_targets);
	pos = 2 * num_targets - 1;
	pathVec[pos] = lastIdx;
	for (int i = 2 * num_targets; i > 1; i--) { // the total number of endpoints on boundaries is 2*n
		pos--;
		lastIdx = backwardIndex[i][lastIdx];
		pathVec[pos] = lastIdx;
	}

	// generate spTour for Lin-Kern searching
	VERTEX t;
	int tarID, subKey, matKey;
	for (int i = 0; i < pathVec.size(); i++) {
		tarID = seq[floor(i / 2) + 1];
		subKey = pathVec[i];
		matKey = (tarID - 1) * num_dstzn + subKey + 1;
		t = { tarID ,subKey, matKey };
		spTour[i + 2] = t;

	}

	spTour[0] = { 0, -1, 0 };
	spTour[1] = { 0, -1, 0 };
}

/*
void Discretization::generate_sortedCost()
{
double ** min_dist_nbds = CI->get_min_dist_nbds();

int s, t, i, j;

double w = 0.0;

//source <--> other nodes on circles
int idxmat_1, idxmat_2, idxmat_3, idxmat_4;

for (t = 1; t <= num_targets; t++) //t-th circle
{
idxmat_1 = (t - 1) * num_dstzn + 1;
idxmat_2 = (t - 1) * 2 * num_dstzn + 1;
for (i = 0; i < num_dstzn; i++) // (j+1)-th node in the t-th circle
{
w = adj_mat[0][idxmat_2 + i];
sortedCost[0][idxmat_1 + i] = { w, {t, i}}; // source is set to be connected with all entries
sortedCost[idxmat_1 + i][0] = { w, {0,-1}}; // all exits are set to be linked with source

matDist[0][idxmat_1 + i] = w;
matDist[idxmat_1 + i][0] = w;
}
}

// s-th circle <--> t-th circle
for (s = 1; s <= num_targets; s++)
{
idxmat_1 = (s - 1) * num_dstzn + 1;
idxmat_3 = (s - 1) * 2 * num_dstzn + num_dstzn + 1;
for (i = 0; i < num_dstzn; i++)
{
for (t = 1; t <= num_targets; t++)
{
if (s != t)
{
idxmat_2 = (t - 1) * num_dstzn + 1;
idxmat_4 = (t - 1) * 2 * num_dstzn + 1;
for (j = 0; j < num_dstzn; j++)
{
w = adj_mat[idxmat_3][idxmat_4]; // point {s,i} -> point {t,j}
w = para_ext_risk * (w - min_dist_nbds[s][t]);
sortedCost[idxmat_1 + i][idxmat_2 + j] = { w, {t, j}};
sortedCost[idxmat_2 + j][idxmat_1 + i] = { w, {s, i}};
matDist[idxmat_1 + i][idxmat_2 + j] = w;
matDist[idxmat_2 + j][idxmat_1 + i] = w;
}
}

}

}
}
int cnt = num_dstzn * num_targets + 1;

for (int i = 0; i < cnt; i++) {
sortedCost[i][i] = { INFINITY,{i, -1}};
}

// sorted
for (int i = 0; i < cnt; i++) {
sort(sortedCost[i].begin(), sortedCost[i].end(), comparator);
}

}
*/

