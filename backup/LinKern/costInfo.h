
#ifndef _COSTINFO_H_
#define _COSTINFO_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

class costInfo{
    
private:
    int n;// number of nodes 
    double* x;
    double* y;
    double** cost;
	vector<vector<pair<double, int>>> sortedCost;
	
    
public:
    ~costInfo();
    costInfo(const char* filename);
	inline double** getCost() { return cost;};
	inline vector<vector<pair<double, int>>>& getSortedCost() { return sortedCost; };
    double getCost (int i, int j);
    void print();

    inline int getNodes(){return n;};
    inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

#endif