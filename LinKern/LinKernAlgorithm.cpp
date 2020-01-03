#pragma once
#include <ctime>
#include <iostream>
#include "LinKernAlgorithm.h"
/*
	Coded by Cai Gao, 4/26/2018. 
	***  Main Data Structure ***
	1. Tour: Current tour, vector<int>(N);
	2. Vertex: each vertex has ID and KEY
	3. 

*/


LinKern::LinKern(costInfo* ci) {

	this->N = ci->getNodes();
	this->D = ci->getCost();
	this->sortedCost = ci->getSortedCost();
	
	runTimes = 1;
	optimalLength = INFINITY;
}



void LinKern::solve() {

	kingsLanding(); // main function of Lin-Kern algorithm

}

/*  Main function of Lin-Kern Algorithm 
	a) Randomly create a feasible tour;
	b) Start Lin-Kern search for this tour.
	c) Once getting a positive gain, update the tour and go back to b)
	d) Stop util no positive gain exists and go back to a) for at most 'runTimes' times. 	
*/
void LinKern::kingsLanding() { 
//	vector<int> Tour = {4,8,2,9,3,0,1,7,6,5};

//	vector<int> Tour = {4,8,2,9,3,0,1,7,6,5 };
	
//	vector<int> Tour = {3,1,7,6,4,0,9,5,2,8 };

//	vector<int> Tour = {4,9,3,6,7,1,0,8,2,5};
//	vector<int> Tour = { 8,0,3,6,7,1,9,4,5,2 };
//	double curBestLength = LKSearch(Tour);
//	cout << "Algorithm stops!!" << endl;
	
	
	int iterations = 1;
	double curBestLength;
	vector<int> Tour = createRandomTour();
	
	while (true)
	{	
		curBestLength = LKSearch(Tour);

		if (curBestLength < optimalLength) {
			optimalLength = curBestLength;
			optimalTour = Tour;
		}
		else {
			iterations++;
			Tour = createRandomTour();
		}

		if (iterations > runTimes) {
			cout << "Lin-Kern Algorithm Stops!!" << endl;
			break;
		}
	}

	printSol(optimalTour);
	cout << "Tour length is " << optimalLength << endl;
	
}



/* initialize a random tour for LinKern algorithm */
vector<int> LinKern::createRandomTour() {
	vector<int> Tour;
	// create a sequence of increasing order
	for (int i = 0; i < N; i++) {
		Tour.push_back(i);
	}

	// randomly swap any two elements
	srand((unsigned)time(0));
	int count = 0;
	int randIdx1, randIdx2;
	while (true)
	{
		randIdx1 = (rand() % N);
		randIdx2 = (rand() % N);
		//	cout << randIdx1 << "  " << randIdx2 << endl;
		if (randIdx1 != randIdx2) {
			swap(Tour, randIdx1, randIdx2);
			count++;
		}

		if (count > N) {
			return Tour;
			break;
		}
	}
}

/*  swap i-th element with j-th element  */
void LinKern::swap(vector<int>& Tour, int i, int j) {
	if (i < 0 || i > N - 1 || j < 0 || j > N - 1) {
		cerr << " ERROR: index out of range!! " << endl;
	}
	int temp = Tour[i];
	Tour[i] = Tour[j];
	Tour[j] = temp;
}


/*  @@  Search current Tour and try to return a new tour with less length by Lin-Kern search 
	a) Choose some vertex 'randomly' as a start. 
	b) If some positive gain exists, go back to a)
	c) If no positive gain can be obtained from all vertices, stop searching for this Tour and go back to up-level. 
	*** Tour automatically got updated ***
	Return: total length of a new tour  @@*/
double LinKern::LKSearch(vector<int>& Tour) {
	
	double curMinLength;
	curMinLength = getTourLength(Tour);
	pair<double, double> ret_searchVtx; //return <minimal length, gain>
	
	while (true)  // stop until we search all possible cases and find no gain 
	{
		int vtx = 0;
		while (true) // start Lin-Kern search from vertex 'vtx' and break once finding some gain
		{
			ret_searchVtx = LKSearch(Tour, vtx);
			if (ret_searchVtx.first < curMinLength && ret_searchVtx.first != -1) {
				curMinLength = ret_searchVtx.first;
				break;
			}
			vtx++;

			if (vtx == N) {
				return curMinLength;
		
			}
		}
	}	
}


/*  @@ For the current tour, start Lin-Kern search from vertex 'v'. 
	a) search x1, return and go back to up-level if there exists some positive gain. Otherwise, go to b);
	b) search alternative x1, return  and go back to up-level if there exists some positive gain. Otherwise, go to c);
	c) no gain can be obtained from 'v', go back to up-level and choose another vertex.
	*** Notice that Tour get updated automatically if there exists a positive gain ***
    Return <current minimal length, gain>  @@ */

pair<double, double> LinKern::LKSearch(vector<int>& Tour, int key_t1) {
	pair<double, double> ret_suc;
	pair<double, double> ret_pred;
	pair<double, double> ret_none;
	while (true){
		 ret_suc = LKSearch(Tour, key_t1, trueKey(key_t1 + 1));
		 if (ret_suc.second > 0) { // choose x1
			 return ret_suc;
	
		 }
		 else {
			 ret_pred = LKSearch(Tour, key_t1, trueKey(key_t1 - 1));
			 if (ret_pred.second > 0) { // choose alternate x1
				 return ret_pred;
		
			 }
			 else { // no gain from current vertex v
				 ret_none.first = -1;
				 ret_none.second = -1;
				 return ret_none;
		
			 }
		 }
	}

}



/*  @@ In level 1, find available <t1, t2, t3, t4>
	a) if x2 gives a positive gain, return and go back to up-level; Otherwise go to b)
	b) if alternative x2 gives a positive gain, return and go back to up-lelve. Otherwise go to c)
	c) no gain from current t1 and t2;

*/
pair<double, double> LinKern::LKSearch(vector<int>& Tour, int key_t1, int key_t2) {
	
	vector<int> tempFlipTour = Tour;
	int k_flips = 0; // numbers of flips

	
	// get t1 and t2
	int id_t1 = Tour[key_t1];
	int id_t2 = Tour[key_t2];
	Vertex t1 = { id_t1, key_t1 };
	Vertex t2 = { id_t2, key_t2 };
	vector<pair<Vertex, Vertex>> X;
//	vector<pair<Vertex, Vertex>> Y;
	X.push_back(make_pair(t1, t2));


	double x1 = D[id_t1][id_t2];
	double x2, y1, g1;

	double G_star = 0.0;
	double Gi = 0.0;

	int id_t3, key_t3, id_t4, key_t4, id2_t4, key2_t4; // ID and Key
	int ngbor1, ngbor2;

	bool posGain_flag;
	for (int i = 0; i < N; i++) {
		// get t3
		id_t3 = sortedCost[id_t2][i].second;
		key_t3 = getKey(Tour, id_t3);
		Vertex t3 = { id_t3, key_t3 };
		
		y1 = D[id_t2][id_t3];
		ngbor1 = trueKey(key_t3 - 1);
		ngbor2 = trueKey(key_t3 + 1);
		// check whether g1 is positive. If not, break the for loop (sorted cost is used)
		g1 = x1 - y1;
		posGain_flag = isPositiveGain(Gi, g1);
		if (posGain_flag == false) { 
			break; 
		}
		// t3 must be unvisited and its two ngbors are unvisited as well. 
		if (!isInSet(X, key_t3) && !isInSet(X, ngbor1) && !isInSet(X, ngbor2)) { 
			// get the correct ngbor to close up the tour. The other t4 is stored as {id2_t4, key2_t4}
			id_t4 = getClosingUpNode(tempFlipTour, X, id_t3);
			// if t4 is adjacent to t1, then we can't close up the tour by connecting t1 and t4. 
			if (isNeighbor(Tour, key_t1, id_t4)) { 
				continue; 
			} 
			key_t4 = getKey(Tour, id_t4);			
			Vertex t4 = { id_t4, key_t4 };

			if (key_t4 == ngbor1) {
				key2_t4 = ngbor2;
				key2_t4 = Tour[key2_t4];
			}
			else {
				key2_t4 = ngbor1;
				key2_t4 = Tour[key2_t4];					
			}
			// We use tempFlipTour to find the closingup vertex. 
			flip(tempFlipTour, id_t1, id_t2, id_t3, id_t4); 
			
			// update Gain
			Gi += g1;

			// update G*. If G* gets updated at level 1, we set k_flips equal to 1. Otherwise, k_flips is still 0.
			bool imprvedG_Star = updateGStar(G_star, Gi, id_t1, id_t3, id_t4);
			if (imprvedG_Star == true) {
				k_flips = 1;
			}
			
			// push back <t3, t4> for level 2.
			X.push_back(make_pair(t3, t4));
			
			// level2 searching. (Note that these parameters may get updated)
			LKSearch_x2(Tour, tempFlipTour, X, k_flips, Gi, G_star);
			

			// if such pair of <t3,t4> returns positive gain, we FLIP and get new Tour. Otherwise, go back and try to find another <t3,t4>
			if (G_star > 0) { // 
				FLIP(Tour, X, k_flips);
				return make_pair(getTourLength(Tour), G_star);
			}
			else {
				X.pop_back(); //Current pair of <t3,t4> returns nothing. So we decide to choose another pair of <t3,t4>. 	
				Gi -= g1;     // current g1 is not chosen.
				tempFlipTour = Tour;  // no flip at level 1
			}

			// try alternative x2;
			/*
			vector<pair<Vertex, Vertex>> X2;
			X2.push_back(make_pair(t1, t2));
			Vertex alter_t4 = { id2_t4,key2_t4 };
			X2.push_back(make_pair(t3, alter_t4));
			*/

			
		}
	
	}
	// function need to return something even though G_star is 0 after checking all pairs of <t3,t4>

	return make_pair(getTourLength(Tour), G_star);
}

/*   @@ Level 2 search
	



*/
double LinKern::LKSearch_x2(vector<int>& Tour, vector<int>& tempFlipTour, vector<pair<Vertex, Vertex>>& X, int& k_flips, double& Gi, double& G_star)
{
	vector<int> copyTempFlipTour = tempFlipTour;
	vector<pair<Vertex, Vertex>> copyX = X;
	double y2, g2, y_star, g_star;

	
	int id_t1  = X[0].first.id;
	int key_t1 = X[0].first.key;
	int id_t4  = X[1].second.id;
	int key_t4 = X[1].second.key;
	
	double x2 = D[X.back().first.id][X.back().second.id];

	bool posGain_flag;
	
	int id_t5, key_t5, id_t6, key_t6;
	int ngbor1, ngbor2;


	for (int i = 0; i < N; i++) {
		//select t5
		id_t5 = sortedCost[id_t4][i].second; 
		// if t5 is a ngbor of t4, we can't add new y to connect t4 and t5.
		if (isNeighbor(Tour, key_t4, id_t5)){
			continue; 
		}
		key_t5 = getKey(Tour, id_t5);
		Vertex t5 = {id_t5, key_t5};

		y2 = D[id_t4][id_t5];
		ngbor1 = trueKey(key_t5 - 1);
		ngbor2 = trueKey(key_t5 + 1);
		// if G2 = g1+g2 >0, t5 is qualified to go into next step.
		g2 = x2 - y2;
		posGain_flag = isPositiveGain(Gi, g2);
		if (posGain_flag == false) { 
			break; 
		}
		// t5 must be unvisited
		if (!isInSet(X, key_t5) ) {
			// t6 is a closingup node. If t6 is a ngbor of t1, such pair <t5,t6> doesn't work.
			id_t6 = getClosingUpNode(tempFlipTour, X, id_t5);
			if (isNeighbor(Tour, key_t1, id_t6)) { 
				continue; 
			}	
			key_t6 = getKey(Tour, id_t6);
			if (!isInSet(X, key_t6)) {
							
				Vertex t6 = { id_t6, key_t6 };

				flip(tempFlipTour, id_t1, id_t4, id_t5, id_t6);

				// update gain
				Gi += g2;

				// stop condition
				if (Gi <= G_star) { 
					break; 
				}				

				// If G* got updated at current level, then update the best k as well. 
				bool improvedGain = updateGStar(G_star, Gi, id_t1, id_t5, id_t6);
				if (improvedGain == true) {
					k_flips = 2;
				}
				
				// push back <t5,t6> for continuing searching
				X.push_back(make_pair(t5, t6));

				// continue search t7, t8 and so on ...
				LKSearch_Continue(Tour, tempFlipTour, X, k_flips, Gi, G_star);

				// at level 2, such pair of <t5,t6> returns positive gain.
				if (G_star > 0) { 
					return G_star; // go back to up level. 
				}
				else {
					X = copyX; // choose another pair of <t5,t6>
					Gi -= g2;				
					tempFlipTour = copyTempFlipTour;
				}
			
			}			

		}
	}

	// If no "good" or "avaiable" <t5,t6>, go back to up level.
	return G_star;	
}

/*
	@@@ Level 3,4,5... Searching


*/
double LinKern::LKSearch_Continue(vector<int>& Tour, vector<int>& tempFlipTour, vector<pair<Vertex, Vertex>>& X, int& k_flips, double& Gi, double& G_star)
{
	bool posGain_flag, continueSearching;

	int key_t1 = X[0].first.key;
	int id_t1 = X[0].first.id;
	
	int key_last_t, id_last_t, key_next_t, id_next_t, id_new_t, key_new_t, ngbor1, ngbor2;
	double last_x, new_y, gi;

	int curLevel = 2;

	while (true)
	{
		continueSearching = false;

		curLevel++;
		key_last_t = X.back().second.key;
		id_last_t = X.back().second.id;
		last_x = D[X.back().first.id][X.back().second.id];

		for (int i = 0; i < N; i++) {
			id_new_t = sortedCost[id_last_t][i].second;
			if (isNeighbor(Tour, key_last_t, id_new_t)){
				continue; 
			}

			key_new_t = getKey(Tour, id_new_t);
			Vertex new_t = {id_new_t, key_new_t};

			new_y = D[id_last_t][id_new_t];
			ngbor1 = trueKey(id_new_t - 1);
			ngbor2 = trueKey(id_new_t + 1);

			gi = last_x - new_y;
			posGain_flag = isPositiveGain(Gi, gi);
			if (posGain_flag == false) { break; } // Once we are unable to find positive gain when reaching current vertex, no need to check other vertices in this level anymore.
			if (!isInSet(X, key_new_t) ) {
				
				id_next_t = getClosingUpNode(tempFlipTour, X, id_new_t);
				if (isNeighbor(Tour, key_t1, id_next_t)) {
					continue;
				}
				key_next_t = getKey(Tour, id_next_t);

				if (!isInSet(X, key_next_t)) {
					//update Gi
					Gi += gi;

					if (Gi <= G_star) { break; }

					continueSearching = true; // we continue search for next level.

					flip(tempFlipTour, id_t1, id_last_t, id_new_t, id_next_t);

					// update G*
					bool improvedGain = updateGStar(G_star, Gi, id_t1, id_new_t, id_next_t);
					if (improvedGain == true) {
						k_flips = curLevel;
					}

					Vertex next_t = { id_next_t, key_next_t };
					X.push_back(make_pair(new_t, next_t));

					break; // Once finding a pair of new t's, no need to find more available pairs at this level. 
				
				}
			}
		}

		//Once no further <t_2i, t_2i+1> satisfies the requirement of selection, or current best gain Gi <= G*, we stop such deep searching.
		if (continueSearching == false) {
			break;
		}
	}

	return G_star;

}
/*
	Like normal searching in level 2, all qualified y2's in this case can be checked.

*/

double LinKern::LKSearch_alternative_x2(vector<int>& Tour, vector<pair<Vertex, Vertex>>& X2, double& Gi, double& G_star)
{
	int id_t1 = X2[0].first.id;
	int key_t1 = X2[0].first.key;
	int id_t2 = X2[0].second.id;
	int id_t3 = X2[1].first.id;
	int id_t4 = X2[1].second.id;


	int id2_t4  = X2.back().second.id;
	int key2_t4 = X2.back().second.key;

	int id_t5, key_t5, id_t6, key_t6, ngbor1, ngbor2;
	double x2 = D[X2.back().first.id][X2.back().second.id];
	double y2, g2;

	bool posGain_flag;
	
	// get two segments
	pair<vector<int>, vector<int>> SEG = getTwoSegments(Tour, id_t1, id_t2, id_t3, id_t4);

	for (int i = 0; i < N; i++) {
		id_t5 = sortedCost[id2_t4][i].second;
		key_t5 = getKey(Tour, id_t5);

		Vertex t5 = { id_t5, key_t5 };
		y2 = D[id2_t4][id_t5];

		g2 = x2 - y2;
		posGain_flag = isPositiveGain(Gi, g2); // t5 need to result in positive Gi
		if (posGain_flag == false) {
			break;
		}

		if (!isInSet(X2, key_t5)) { // unvisited vertex
			ngbor1 = trueKey(key_t5 - 1);
			ngbor2 = trueKey(key_t5 + 1);

			if (isInSet(SEG.second, id_t5)) { // Case I of alternate x2
				// select t6
				if (!isInSet(X2, ngbor1)) {
					key_t6 = ngbor1; 
					if (!isInSet(X2,key_t6) && (D[id_t5][Tour[ngbor2]] > D[id_t5][Tour[ngbor1]])) {
						key_t6 = ngbor2;				
					}
				}
				else {
					if (!isInSet(X2, key_t6)) {
						key_t6 = ngbor2;
					}
					else{
						continue;
					}
				}
				id_t6 = Tour[key_t6];

				Vertex t6 = { id_t6, key_t6 };
				// update gain
				Gi += g2;

				
			//	update G*: updateGStar();
				
			//  get the new Tour: getNewTour();


			//  continue searching for more gain

			// if G* > 0, FLIP() and get gain. 
			}
			else { //Case II of alternative x2.
				// choose the closingup t6
				// choose t7 from segment t2~t3


				// update G*

				// get the new Tour: getNewTour2()

				// continue searching for more gain

				//if G* > 0, FLIP() and get gain
			
			}

		
		
		}



	
	}

	return -1;
}


pair<vector<int>, vector<int>> LinKern::getTwoSegments(vector<int>& Tour, int id_t1, int id_t2, int id_t3, int id_t4)
{
	int key_t1 = getKey(Tour, id_t1);
	int key_t2 = getKey(Tour, id_t2);
	int key_t3 = getKey(Tour, id_t3);
	int key_t4 = getKey(Tour, id_t4);
	
	vector<int> seg1, seg2; // two segments
	int temp;
	if (key_t2 == trueKey(key_t1 + 1)) { // t1 -> t2
		int k = 0;
		while (true)
		{
			temp = trueKey(key_t2 + k);
			seg2.push_back(Tour[temp]);
			if (temp == key_t3) {
				break;
			}
			k++;
		}
		
		k = 0;
		while (true)
		{
			temp = trueKey(key_t4 + k);
			seg1.push_back(Tour[temp]);
			if (temp == key_t1) {
				break;
			}
			k++;
		}
	}
	else { // t2 -> t1
	
		int k = 0;
		while (true)
		{
			temp = trueKey(key_t1 + k);
			seg1.push_back(Tour[temp]);
			if (temp == key_t4) {
				break;
			}
			k++;
		}

		k = 0;
		while (true)
		{
			temp = trueKey(key_t3 + k);
			seg2.push_back(Tour[temp]);
			if (temp == key_t2) {
				break;
			}
			k++;
		}
	
	}

	return make_pair(seg1, seg2);

}



bool LinKern::isInSet(vector<int>& seg, int id_t)
{
	for (int i = 0; i < seg.size(); i++) {
		if (id_t == seg[i]){
			return true;
		}
	}
	return false;
}


bool LinKern::isInSet(vector<pair<Vertex, Vertex>>& X, int key_ti)
{
	int lens = X.size();
	bool flag = false;
	for (int i = 0; i < lens; i++) {
		if (key_ti == X[i].first.key || key_ti == X[i].second.key) {
			flag = true;
			break;		
		}
	}
	return flag;
}

bool LinKern::isInTour(vector<int>& Tour, Edge e) {

	for (int i = 0; i < N - 1; i++) {
		Edge e1 = { Tour[i],Tour[i + 1]};
		if (isSameEdge(e1, e)) {
			return true;
		}	
	}
	Edge e1 = { Tour[N - 1],Tour[0] };
	if (isSameEdge(e1, e)) {
		return true;
	}
	return false;

}

bool LinKern::isSameEdge(Edge e1, Edge e2) {
	
	if (e1.u == e2.u && e1.v == e2.v) {
		return true;
	}
	if (e1.u == e2.v && e1.v == e2.u) {
		return true;
	}	
	return false;
}

bool LinKern::isNeighbor(vector<int>& Tour, int key_p, int id_t) {
	if (id_t == Tour[trueKey(key_p - 1)] || id_t == Tour[trueKey(key_p + 1)]) {
		return true;	
	}
	return false;
}


bool LinKern::isPositiveGain(double Gi, double g_local) {
	bool flag = false;
	if (Gi + g_local > 0) {
	//	Gi = Gi + g_local;
		flag = true;
	}
	return flag;
}

void LinKern::updateGain(double& Gi, double g_local) 
{
	Gi += g_local;
}


bool LinKern::updateGStar(double& G_star, double Gi,  int id_t1, int id_t3, int id_t4) {

	double g_star = D[id_t3][id_t4] - D[id_t1][id_t4];
	if (Gi + g_star > G_star) {
		G_star = Gi + g_star;
		return true;
	}
	return false;
}


int LinKern::getClosingUpNode(vector<int>& tempFlipTour, vector<pair<Vertex, Vertex>>& X, int id_t3) {
	int id_t1 = X[0].first.id;
	int id_t2 = X.back().second.id;

	int key_t1 = getKey(tempFlipTour, id_t1);
	int key_t2 = getKey(tempFlipTour, id_t2);
	int key_t3 = getKey(tempFlipTour, id_t3);

	int key_t4; //the new closing up node
	if (key_t2 == trueKey(key_t1 + 1)) { // if t2 the node after t1, then t4 is the node before t3;
		key_t4 = trueKey(key_t3 - 1); // then 
	}
	else { 
		key_t4 = trueKey(key_t3 + 1); // if t2 the node before t1, then t4 is the node after t3;
	}

	int id_t4 = tempFlipTour[key_t4];

//	flip(tempFlipTour, key_t1, key_t2, key_t3, key_t4);

	return id_t4;
}

void LinKern::FLIP(vector<int>& Tour, vector<pair<Vertex, Vertex>> X, int k)
{
	
	for (int i = 0; i < k; i++) {
		flip(Tour, X[0].first.id, X[i].second.id,X[i+1].first.id,X[i+1].second.id);
	}

}


void LinKern::flip(vector<int>& Tour, int id_a, int id_b, int id_c, int id_d) {
	vector<int> valVec;
	vector<int> keyVec;
	int key_a = getKey(Tour, id_a);
	int key_b = getKey(Tour, id_b);
	int key_c = getKey(Tour, id_c);
	int key_d = getKey(Tour, id_d);

	int temp;
	if (trueKey(key_a +1) == key_b) {
		int k = 0;
		while (true)
		{
			temp = trueKey(key_b + k);
			keyVec.push_back(temp);
			valVec.push_back(Tour[temp]);
			if (temp == key_d) {
				break;
			}
			k++;
		}
		for (int i = 0; i < keyVec.size(); i++) {//reverse
			Tour[keyVec[i]] = valVec[keyVec.size() - i - 1];
		}	
	}
	else {
		int k = 0;
		while (true)
		{
			temp = trueKey(key_a + k);
			keyVec.push_back(temp);
			valVec.push_back(Tour[temp]);
			if (temp == key_c) {
				break;
			}
			k++;
		}
		for (int i = 0; i < keyVec.size(); i++) {//reverse
			Tour[keyVec[i]] = valVec[keyVec.size() - i - 1];
		}
	
	
	}

}

int LinKern::getKey(vector<int>& Tour, int id) {
	int k = 0;
	while (true){
		if (id == Tour[k]){
			return k;
			break;
		}
		k++;
	}
}


void LinKern::printSol(vector<int>& Tour) {
	for (int i = 0; i < N; i++) {
		cout << Tour[i] << "   ";
	}
	cout << endl;
//	cout << "Optimal length is " << optimalLength << endl;
}

double LinKern::getTourLength(vector<int> T) {
	double lens = 0.0;
	for (int i = 0; i < N - 1; i++) {
		lens += D[T[i]][T[i + 1]];
	}
	lens += D[T[N - 1]][T[0]];

	return lens;
}

int LinKern::trueKey(int key_t)
{
	if (key_t > N - 1) 
		return key_t % N ;
	else if (key_t < 0)
		return N + key_t;
	else 
		return key_t;	
}



LinKern::~LinKern() {};