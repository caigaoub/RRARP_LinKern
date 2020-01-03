#include <iostream>
#include <math.h>
#include "costInfo.h"
#include <algorithm>

using namespace std;

bool comparator(const pair<double, int>& l, const pair<double, int>& r) {
	return l.first < r.first;
}



costInfo::costInfo(const char* filename) {

	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename
			<< "' for reading" << endl;
		throw (-1);
	}

	file >> n;// read the number of nodes in the graph

	x = new double[n]; // x-coordinate / created in heap
	y = new double[n]; // y-coordinates

	cost = new double*[n]; // create an array of type 'double*' in the heap

	for (int i = 0; i < n; i++) {
		cost[i] = new double[n];
		file >> x[i] >> y[i];
	}

	
	for (int i = 0; i < n; i++){
		for (int j = i + 1; j < n; j++) {
			cost[i][j] = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
			cost[j][i] = cost[i][j];
		}
	}
	
	for (int i = 0; i < n; i++) {
		cost[i][i] = INFINITY;
	}

	sortedCost.resize(n, vector < pair<double, int>>(n));
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			sortedCost[i][j] = make_pair(cost[i][j], j);
			sortedCost[j][i] = make_pair(cost[i][j], i);
		}
	}

	for (int i = 0; i < n; i++) {
		sortedCost[i][i] = make_pair(INFINITY, i);
	}
	// sorted
	for (int i = 0; i < n; i++) {
		sort(sortedCost[i].begin(), sortedCost[i].end(), comparator);
	}
	/*
	cout << "Cost: ----------------: " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << cost[i][j] << '\t';
		}
		cout << endl;
	}

	cout << "Sorted Cost: ----------------: " << endl;
	for (int i =0; i < n; i++)	{
		for (int j = 0; j < n; j++){
			cout << sortedCost[i][j].first << " (" << sortedCost[i][j].second<<") " << '\t';
		}
		cout << endl;
	}
	*/
}

double costInfo::getCost(int i, int j){
    if (i<j)
        return cost[i][j];
    else if(j<i)
        return cost[j][i];
    else
        return 0;
}

void costInfo::print(){
    for (int i=0; i<n; i++)
        for (int j=i+1; j<n; j++){
            cout<<i<<" "<<j<<" "<<cost[i][j]<<endl;
        }
}

costInfo::~costInfo(){
    
    for (int i=0; i<n; i++){
        delete(cost[i]);
    }
    delete x;
    delete y;
}


