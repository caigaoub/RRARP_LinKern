#include "LinKern.h"
#include <iostream>
#include <ctime>
#include <algorithm>
#include <random>

bool comparator(pair<double, int>& l, pair<double, int>& r) {
	return l.first < r.first;
}

LinKern::LinKern(Discretization* disc, CostInfo* ci)
{
	this->disc = disc;
	this->N = disc->get_num_circles() + 1;
	this->num_dstzn = disc->get_num_dstzn();
	
	this->D = disc->getExtDist();
	this->intDist = disc->getIntDist();
	
	runTimes = 30;

	optimalLens = INFINITY;
}

void LinKern::solve()
{	
	kingsLanding();
}

void LinKern::kingsLanding()
{
	int iterations = 1;
	vector<VERTEX> spTour(2 * N);
	double curSPLens;
	double curTourOptLens;
	vector<int> tspSeq = createRandomTour();
	cout << "Create a New Tour: ";
//	for (int i = 0; i < tspSeq.size(); i++)
//		cout << tspSeq[i] << "  ";

	disc->fast_sp2(spTour, curSPLens, tspSeq);
	cout << curSPLens << endl;

	while (true)
	{
		LKSearch(spTour, curSPLens);
		cout << "Tour after LK-Searching:  " << curSPLens << endl;

		tspSeq = getTspTour(spTour);
		disc->fast_sp2(spTour, curTourOptLens, tspSeq);
		cout << "Tour Optimal Lens:  " << curTourOptLens << endl;
		cout << endl;

		if (curTourOptLens < optimalLens) {
			optimalLens = curTourOptLens;
			optimalSPTour = spTour;
		}
		
		iterations++;
		tspSeq = createRandomTour();
		disc->fast_sp2(spTour, curSPLens, tspSeq);
//		cout << "create a new tour: ";
//		for (int i = 0; i < tspSeq.size(); i++)
//			cout << tspSeq[i] << "  ";
		cout << "Create a New Tour:  " << curSPLens << endl;


		if (iterations > runTimes) {
			cout << "Lin-Kern Algorithm Stops!! and OPT is " << optimalLens << endl;
			break;
		}
	}

	
}


double LinKern::LKSearch(vector<VERTEX> & spTour, double& curSPLens) {

	double curMinLength;
	curMinLength = curSPLens;
	double bestGain; //return <minimal length, gain>

	while (true)  // stop until we search all possible cases and find no gain 
	{
		int vtx = 0;
		while (true) // start Lin-Kern search from vertex 'vtx' and break once finding some gain
		{
			bestGain = LKSearch(spTour, vtx);		
	//		bestGain = LKSearch(spTour, 0);
			if (bestGain > 0) {
				curMinLength = curSPLens - bestGain;
				curSPLens = curMinLength;
	//			double tm = getTourLens(spTour);

				break;
			}
			vtx++;

			if (vtx == N) { // if none of targets gives positive gain, we return to upper level
				curSPLens = curMinLength;
				return curMinLength;
			}
		}
	}
}

// <length after lk-opts, gain>

double LinKern::LKSearch(vector<VERTEX> & spTour, int key_t1)
{
	double ret_suc;
	double ret_pred;
	double ret_none;
	while (true) {
		int ngbor1 = trueKey(key_t1 + 1);
		ret_suc = LKSearch(spTour, key_t1, ngbor1);
	//	ret_suc = LKSearch(spTour, 0, 1);
		if (ret_suc > 0) { // choose x1
			return ret_suc;

		}
		else {
			int ngbor2 = trueKey(key_t1 - 1);
			ret_pred = LKSearch(spTour, key_t1, ngbor2);
			if (ret_pred > 0) { // choose alternate x1
				return ret_pred;

			}
			else { // no gain from current vertex v
				ret_none = -1;
				return ret_none;
			}
		}
	}

}


/*  
	@@ Given t1 and t2, we start LinKern searching for positive gain.
	return bestGain (or known as G* in original LinKern Algo.)

*/

double LinKern::LKSearch(vector<VERTEX>& spTour, int key_t1, int key_t2) {
	vector<vector<VERTEX>> newSPTourSet;

	vector<VERTEX>  tempFlipTour = spTour;	

	double G_star = 0.0;
	double Gi = 0.0;
	int k_flips;

	// ----------------------------------------------------------------------------------------------------------------------------
	// ===> Get x1: {point a, point b, length of x1}.
	LKEdgeX lkedge_x1 = get_x(spTour, key_t1, key_t2);
//	int spKey_p1 = lkedge_x1.spKeyA;
	VERTEX p1 = spTour[lkedge_x1.spKeyA];
//	int spKey_p2 = lkedge_x1.spKeyD;
	VERTEX p2 = spTour[lkedge_x1.spKeyD];
	double x1 = lkedge_x1.dist;

	// ----------------------------------------------------------------------------------------------------------------------------
	// ===> Compute the direction of t1->t2
	bool isCW; // 
	if (trueKey(key_t1 + 1) == key_t2)
		isCW = true;
	else
		isCW = false;
	
	// ===> select t3 
	vector<pair<int, int>> brokenSet;
	int id_t1 = p1.tarID;
	int id_t2 = p2.tarID;



	brokenSet.push_back(make_pair(id_t1, id_t2)); // This set is static based on spTour only 

	// ===> collect the candidates for t3
	int spKey_p3, key_t3, id_t3, ngborKey_1, ngborKey_2;
	double dist_p2p3;
	vector<pair<double, int>> approxDist;
	for (int key = 0; key < N; key++) { // key of target (or depot) of spTour
		key_t3 = key;
		ngborKey_1 = trueKey(key_t3 - 1); // consistent with original LinKern algo.
		ngborKey_2 = trueKey(key_t3 + 1);
		if (!isInSet(brokenSet, spTour[2 *key_t3].tarID) && !isInSet(brokenSet, spTour[2 * ngborKey_1].tarID) && !isInSet(brokenSet, spTour[2 * ngborKey_2].tarID)) {
			if (isCW == true) { // choose vertex whose spKey is odd, as p3
				spKey_p3 = 2 * key_t3 + 1; //odd
			}
			else {// choose vertex whose spKey is even, as p3
				spKey_p3 = 2 * key_t3; //even
			}

			// approximate length of y1 = |p2p3|
			dist_p2p3 = D[p2.matKey][spTour[spKey_p3].matKey];
			if (dist_p2p3 < x1) // add this vertex as p3 only if |p2p3| < |x1|
				approxDist.push_back(make_pair(dist_p2p3, spKey_p3));		
		}
	}
	sort(approxDist.begin(), approxDist.end(), comparator);

	// ===> Start choosing from the smallest |p2p3|
	for (int i = 0; i < approxDist.size(); i++) {
		Gi = 0.0;     // current g1 is not chosen.
		spKey_p3 = approxDist[i].second;
		key_t3 = floor(spKey_p3 / 2);
		id_t3 = spTour[spKey_p3].tarID;

		// ===> calculate y1
		VERTEX p3 = spTour[spKey_p3];
		LKEdgeY lkedge_y1 = get_y(spTour, p2, p3);
		double y1 = lkedge_y1.dist;
		

		// ===> check gain
		double g1 = x1 - lkedge_y1.dist;
		bool posGain_flag = isPositiveGain(Gi, g1); 
		if (posGain_flag == false) {
			continue; // break? Diff. from original Link-Kern Algo. 
		}
		
		// ===> update Gain (At this point, we know Gi is definitely positive)
		Gi += g1;

		// ===> find t4 and p4 
		int key_t4, spKey_p4;
		if (isCW) {
			key_t4 = trueKey(key_t3 - 1); // previous vertex
			spKey_p4 = 2 * key_t4; // even spKey
		}
		else {
			key_t4 = trueKey(key_t3 + 1); // succ vertex
			spKey_p4 = 2 * key_t4 + 1; // odd spKey
		}
		VERTEX p4 = spTour[spKey_p4];
		int id_t4 =p4.tarID;
		
		// ===> update G*: use key_t3 & key_t4 as parameters for computing last x and  use p1 and p4 as parameters for computing yStar
		LKEdgeX lkedge_x2 = get_x(spTour, key_t3, key_t4);
		double x2 = lkedge_x2.dist;
		LKEdgeY lkedge_yStar = get_y(spTour, p1, p4);
		double ystar_1 = lkedge_yStar.dist;

		bool isImproved = updateGStar(G_star, Gi, lkedge_x2.dist, lkedge_yStar.dist);
		if (isImproved == true) {
			k_flips = 1;
		}
		
		// ---------------------------------------------------------------------------------------------------------------
		// update four (or three) old middle points with new ones 
		updatePoints(tempFlipTour, lkedge_x1, lkedge_x2, lkedge_y1, lkedge_yStar);

		// flip
		flip(tempFlipTour, id_t1, id_t2, id_t3, id_t4);
		// ---------------------------------------------------------------------------------------------------------------


		// write this flipped tour into newSPTourSet
		newSPTourSet.push_back(tempFlipTour);
		
		//  ===> LKSearch_x2
		brokenSet.push_back(make_pair(id_t3, id_t4));
		pair<VERTEX, VERTEX> P = make_pair(p1, p4);
//		cout << "t1 = " << id_t1 << "  t2 = " << id_t2 << endl;
//		cout << "x1 = " << x1 << "y1 " << y1 << endl;;
//		cout << "t3 = " << id_t3 << "  t4 = " << id_t4 << endl;
//		cout << "x2 = " << lkedge_x2.dist << endl;
//		cout << "y*(2) = " << lkedge_yStar.dist << endl;
		LKSearch_x2(spTour, tempFlipTour, brokenSet, newSPTourSet, k_flips, Gi, G_star, P, lkedge_x2.dist);
		if (G_star > 0) {
			spTour = newSPTourSet[k_flips - 1];
			return G_star;
		}
		else {
			newSPTourSet.clear();
			brokenSet.pop_back(); //Current pair of <t3,t4> returns nothing. So choose another pair of <t3,t4>. 	
			
			tempFlipTour = spTour;  // no flip at level 1. 
		}
	}
	return G_star;
}
/* @@ LEVEL 2: search t5 and t6
	
*/

double LinKern::LKSearch_x2(vector<VERTEX>& spTour, vector<VERTEX>& tempFlipTour, vector<pair<int, int>>& brokenSet, vector<vector<VERTEX>>& newSPTourSet, int& k_flips, double& Gi, double& G_star, pair<VERTEX,VERTEX>& P, double x2)
{
	vector<vector<VERTEX>> copySPTourSet = newSPTourSet;
	vector<VERTEX> copyTempFlipTour = tempFlipTour;
	double copyGi = Gi;
	double g2;

	VERTEX p1 = P.first;
	VERTEX p4 = P.second;
	
//	int spKey_p1 = getSPKey(spTour, p1);
//	int spKey_p4 = getSPKey(spTour, p4);

	// ===> direction t1->t4 in tempFlipTour
	int key_t1_TFT = getKey(tempFlipTour, p1.tarID);
	int key_t4_TFT = getKey(tempFlipTour, p4.tarID);
	bool isCW;
	if (trueKey(key_t1_TFT + 1) == key_t4_TFT)
		isCW = true;
	else
		isCW = false;

	//---------------------------------------------------------------------------------------------------------------------------------
	// ===> select t5
	int id_t5, key_t5_TFT, spKey_p5_TFT, id_t6, key_t6_TFT, spKey_p6_TFT;
	double dist_p4p5;
	vector<pair<double, int>> approxDist;
	for (int key = 0; key < N; key++) { // key of tempFlipTour
		key_t5_TFT = key;
		id_t5 = tempFlipTour[2 * key_t5_TFT].tarID;
		if (!isInSet(brokenSet,id_t5)) { // this t5 should be unvisited
			if (isNeighbor(spTour, p4.tarID, id_t5)) { // t5 shouldn't be adjacent to t4 in spTour
				continue;
			}
			if (isCW == true) {
				spKey_p5_TFT = 2 * key_t5_TFT + 1;
			}
			else {
				spKey_p5_TFT = 2 * key_t5_TFT;
			}

			dist_p4p5 = D[p4.matKey][tempFlipTour[spKey_p5_TFT].matKey];
			approxDist.push_back(make_pair(dist_p4p5, spKey_p5_TFT));
		}	
	}
	sort(approxDist.begin(), approxDist.end(), comparator);


	//---------------------------------------------------------------------------------------------------------------------------------
	// ===> Start choosing from the smallest |p4p5|
	for (int i = 0; i < approxDist.size(); i++) {
		Gi = copyGi;     // current g2 is not chosen.

		spKey_p5_TFT = approxDist[i].second;
		key_t5_TFT = floor(spKey_p5_TFT / 2);
		id_t5 = tempFlipTour[spKey_p5_TFT].tarID;

		// compute y2
		VERTEX p5 = tempFlipTour[spKey_p5_TFT];
		LKEdgeY lkedge_y2 = get_y(spTour, p4, p5);
		double y2 = lkedge_y2.dist;

		// ===> check positivity of gain
		 g2 = x2 - lkedge_y2.dist;
		 bool posGain_flag = isPositiveGain(Gi, g2); 
		 if (posGain_flag == false) {
			 continue;
		 }

		 Gi += g2;

		 // find t6
		 if (isCW) { 
			 key_t6_TFT = trueKey(key_t5_TFT - 1); // previous vertex
			 spKey_p6_TFT = 2 * key_t6_TFT; // even spKey
			
		 }
		 else {
			 key_t6_TFT = trueKey(key_t5_TFT + 1); // succ vertex
			 spKey_p6_TFT = 2 * key_t6_TFT + 1; // odd spKey
		 }
		 VERTEX p6 = tempFlipTour[spKey_p6_TFT];
		 id_t6 = p6.tarID;
		
		 if (!isInSet(brokenSet, id_t6)) { // make sure t6 is unvisited
			 // ===> update G*
			 LKEdgeX lkedge_x3 = get_x(tempFlipTour, key_t5_TFT, key_t6_TFT);
			 LKEdgeY lkedge_yStar = get_y(spTour, p1, p6);

			 double x3 = lkedge_x3.dist;
			 double ystar_2 = lkedge_yStar.dist;
	//		 cout << "x3 = " << lkedge_x3.dist << "   " << " y*(2) = " << lkedge_yStar.dist << endl;

			 bool isImproved = updateGStar(G_star, Gi, lkedge_x3.dist, lkedge_yStar.dist);
			 if (isImproved == true) {
				 k_flips = 2;
			 }

			 // ===> stop condition
			 if (Gi <= G_star) {
				 // update four points 
				 LKEdgeX lkedge_x = get_x(tempFlipTour, key_t1_TFT, key_t4_TFT);
				 updatePoints(tempFlipTour, lkedge_x, lkedge_x3, lkedge_y2, lkedge_yStar);

				 // flip
				 flip(tempFlipTour, p1.tarID, p4.tarID, id_t5, id_t6);

				 // Push back this new sptour
				 newSPTourSet.push_back(tempFlipTour);

				 break;
			 }
		 
			 // update four points 
			 LKEdgeX lkedge_x =  get_x(tempFlipTour, key_t1_TFT,key_t4_TFT);
			 updatePoints(tempFlipTour, lkedge_x, lkedge_x3, lkedge_y2, lkedge_yStar);		 

			 // flip
			 flip(tempFlipTour, p1.tarID, p4.tarID, id_t5, id_t6);

			 // Push back this new sptour
			 newSPTourSet.push_back(tempFlipTour);

			 // <p1,p6>
			 brokenSet.push_back(make_pair(id_t5, id_t6));
			 pair<VERTEX, VERTEX> PP = make_pair(p1, p6);

			 // continue searching more gain
			 LKSearch_Continue(spTour, tempFlipTour, brokenSet, newSPTourSet, k_flips, Gi, G_star, PP, lkedge_x3.dist);

			 if (G_star > 0) {
		//		 spTour = newSPTourSet[k_flips - 1];
				 return G_star;
			 }
			 else {
				 newSPTourSet = copySPTourSet;
				 brokenSet.pop_back(); 
				 tempFlipTour = copyTempFlipTour;  
			 }

		 }
 	}
	return G_star;
}
/*  @@ LEVEl >=3


*/

double LinKern::LKSearch_Continue(vector<VERTEX>& spTour, vector<VERTEX>& tempFlipTour, vector<pair<int,int>>& brokenSet, vector<vector<VERTEX>>& newSPTourSet, int& k_flips, double& Gi, double& G_star, pair<VERTEX, VERTEX>& P, double last_x)
{
	bool posGain_flag, continueSearching;

	VERTEX p1 = P.first;
//	int spKey_p1 = getSPKey(spTour, p1);
	VERTEX last_p = P.second;
//	int spKey_last_p = getSPKey(spTour, last_p);

	double gi;

	int curLevel = 2;
	int key_t1_TFT, key_last_t_TFT, id_last_t, key_next_t_TFT, id_next_t, id_new_t, key_new_t_TFT, spKey_next_p_TFT, spKey_new_p_TFT;

	while (true)
	{
		continueSearching = false;
		curLevel++;

		// ===> get direction of t1->t2i in tempFlipTour
		key_t1_TFT  = getKey(tempFlipTour, p1.tarID);
		key_last_t_TFT = getKey(tempFlipTour, last_p.tarID);
		bool isCW;
		if (trueKey(key_t1_TFT + 1) == key_last_t_TFT)
			isCW = true;
		else
			isCW = false;
	
		// ===> get the candidates of next_p
		double dist_next_p_last_p;
		vector<pair<double, int>> approxDist;
		for (int key = 0; key < N; key++) { // key of tempFlipTour
			key_next_t_TFT = key;
			id_next_t = tempFlipTour[2 * key_next_t_TFT].tarID;
			if (!isInSet(brokenSet, id_next_t)) { // make sure next_t is unvisited 
				if (isNeighbor(spTour, last_p.tarID, id_next_t)) { 
					continue;
				}

				if (isCW == true) {
					spKey_next_p_TFT = 2 * key_next_t_TFT + 1;
				}
				else {
					spKey_next_p_TFT = 2 * key_next_t_TFT;
				}

				dist_next_p_last_p = D[last_p.matKey][tempFlipTour[spKey_next_p_TFT].matKey];
				approxDist.push_back(make_pair(dist_next_p_last_p, spKey_next_p_TFT));
			}
		}
		sort(approxDist.begin(), approxDist.end(), comparator);

		// ===> select next_p
		for (int i = 0; i < approxDist.size(); i++) {
			spKey_next_p_TFT = approxDist[i].second;
			key_next_t_TFT = floor(spKey_next_p_TFT / 2);
			id_next_t = tempFlipTour[spKey_next_p_TFT].tarID;
		
			// compute
			VERTEX next_p = tempFlipTour[spKey_next_p_TFT];
			LKEdgeY lkedge_new_y = get_y(spTour, last_p, next_p);
			
			gi = last_x - lkedge_new_y.dist;
			bool posGain_flag = isPositiveGain(Gi, gi); // Gi > 0
			if (posGain_flag == false) {
				break; 
			}
			
			// find new_t and new_p
			if (isCW) {
				key_new_t_TFT = trueKey(key_next_t_TFT - 1); // previous vertex
				spKey_new_p_TFT = 2 * key_new_t_TFT; // even spKey
			}
			else {
				key_new_t_TFT = trueKey(key_next_t_TFT + 1); // succ vertex
				spKey_new_p_TFT = 2 * key_new_t_TFT + 1; // odd spKey
			}

			VERTEX new_p = tempFlipTour[spKey_new_p_TFT];
			id_new_t = new_p.tarID;
			if (!isInSet(brokenSet, id_new_t)) { // the new_t must be unvisited
																
				// ===> update G*
				LKEdgeX lkedge_new_x = get_x(tempFlipTour, key_next_t_TFT, key_new_t_TFT);
				LKEdgeY lkedge_yStar = get_y(spTour, p1, new_p);

				Gi += gi;

				bool isImproved = updateGStar(G_star, Gi, lkedge_new_x.dist, lkedge_yStar.dist);
				if (isImproved == true) {
					k_flips = curLevel;
				}
			
				// stop condition
				if (Gi <= G_star) {

					// update four points 
					LKEdgeX lkedge_last_x = get_x(tempFlipTour, key_t1_TFT, key_last_t_TFT);
					updatePoints(tempFlipTour, lkedge_last_x, lkedge_new_x, lkedge_new_y, lkedge_yStar);

					// flip
					flip(tempFlipTour, p1.tarID, last_p.tarID, id_next_t, id_new_t);

					// Push back this new sptour
					newSPTourSet.push_back(tempFlipTour);

					break;
				}

				continueSearching = true; // we continue search for next level.
				
				VERTEX new_last_p = tempFlipTour[spKey_new_p_TFT]; // last_p should be changed before flipped. 

				// update four points 
				LKEdgeX lkedge_x = get_x(tempFlipTour, key_t1_TFT, key_last_t_TFT);
				updatePoints(tempFlipTour, lkedge_x, lkedge_new_x, lkedge_new_y, lkedge_yStar);
				// flip
				flip(tempFlipTour, p1.tarID, last_p.tarID, id_next_t, id_new_t);


				last_x = lkedge_new_x.dist;
				last_p = new_last_p;
				
				newSPTourSet.push_back(tempFlipTour);

				brokenSet.push_back(make_pair(id_next_t, id_new_t));

				break; // once finding a new pair of t's, no need to search others at  current level

			}


		}

		//Once no further <t_2i, t_2i+1> satisfies the requirement of selection, or current best gain Gi <= G*, we stop such deep searching.
		if (continueSearching == false) {
			break;
		}
	}

	return G_star;
}


bool LinKern::isPositiveGain(double Gi, double g_local) {
	bool flag = false;
	if (Gi + g_local > 0) {
		//	Gi = Gi + g_local;
		flag = true;
	}
	return flag;
}


bool LinKern::getDirection(int key_t1, int key_t2)
{
	bool clockwise;
	if (trueKey(key_t1 + 1) == key_t2)
		clockwise = true;
	else
		clockwise = false;
	return clockwise;
}


/* @@ return LKEdge x1 <length, four spKey's of this edge> 
   We assume the direction of spTour[0] -> spTour[1] is clockwise. If t1->t2 is cw, then the left-most in spTour[key_t1] is p1
   and right-most in spTour[key_t2] is p2.
*/
LKEdgeX LinKern::get_x(vector<VERTEX>& spTour, int key_t1, int key_t2)
{
	bool isCW; // if t1->t2 is clockwise, then let isCW be true; Otherwise, isCW = false;
	if (trueKey(key_t1 + 1) == key_t2) {
		isCW = true;
	}
	else {
		isCW = false;
	}

	int idx1, idx2, idx3, idx4;
	double dist;
	int spKeyA, spKeyB, spKeyC, spKeyD; 	

	if (isCW) { // isCW = true
		/*
		if (spTour[2*key_t1].tarID == 0) { //t1 is the depot, t2 is a circle
			idx1 = 0;
			idx2 = spTour[2 * key_t2].matKey;
			idx3 = spTour[2 * key_t2 + 1].matKey;
			dist = D[idx1][idx2] + D[idx2][idx3];
		}

		if (spTour[2 * key_t2].tarID == 0) { // t1 is a cirlce and t2 is a depot
			idx1 = spTour[2 * key_t1].matKey;
			idx2 = spTour[2 * key_t1 + 1].matKey;
			idx3 = 0;
			dist = D[idx1][idx2] + D[idx2][idx3];
		}

		if (spTour[2 * key_t2].tarID != 0 && spTour[2 * key_t2].tarID != 0) {
			idx1 = spTour[2 * key_t1].matKey;
			idx2 = spTour[2 * key_t1 + 1].matKey;
			idx3 = spTour[2 * key_t2].matKey;
			idx4 = spTour[2 * key_t2 + 1].matKey;
			dist = D[idx1][idx2] + D[idx2][idx3] + D[idx3][idx4];
		}
		*/
		idx1 = spTour[2 * key_t1].matKey;
		idx2 = spTour[2 * key_t1 + 1].matKey;
		idx3 = spTour[2 * key_t2].matKey;
		idx4 = spTour[2 * key_t2 + 1].matKey;


		dist = D[idx1][idx2] + D[idx2][idx3] + D[idx3][idx4];

		spKeyA = 2 * key_t1;
		spKeyB = (spKeyA + 1) % (2 * N);
		spKeyC = (spKeyB + 1) % (2 * N);
		spKeyD = (spKeyC + 1) % (2 * N);
	}
	else { // counter-clockwise
		/*
		if (spTour[2 * key_t1].tarID == 0) { //t1 is the depot
			idx1 = 0;
			idx2 = spTour[2 * key_t2 + 1].matKey;
			idx3 = spTour[2 * key_t2].matKey;
			dist = D[idx1][idx2] + D[idx2][idx3];
		}

		if (spTour[2 * key_t2].tarID == 0) { // t2 is the depot
			idx1 = spTour[2 * key_t1 + 1].matKey;
			idx2 = spTour[2 * key_t1].matKey;
			idx3 = 0;
			dist = D[idx1][idx2] + D[idx2][idx3];
		}

		if (spTour[2 * key_t2].tarID != 0 && spTour[2 * key_t2].tarID != 0) {
			idx1 = spTour[2 * key_t1 + 1].matKey;
			idx2 = spTour[2 * key_t1].matKey;
			idx3 = spTour[2 * key_t2 + 1].matKey;
			idx4 = spTour[2 * key_t2].matKey;
			dist = D[idx1][idx2] + D[idx2][idx3] + D[idx3][idx4];
		}
		*/
		idx1 = spTour[2 * key_t1 + 1].matKey;
		idx2 = spTour[2 * key_t1].matKey;
		idx3 = spTour[2 * key_t2 + 1].matKey;
		idx4 = spTour[2 * key_t2].matKey;
		dist = D[idx1][idx2] + D[idx2][idx3] + D[idx3][idx4];

		spKeyD = 2 * key_t2;
		spKeyC = (spKeyD + 1) % (2 * N);;
		spKeyB = (spKeyC + 1) % (2 * N);;
		spKeyA = (spKeyB + 1) % (2 * N);;
	}
	LKEdgeX E = { dist, spKeyA, spKeyB, spKeyC, spKeyD };
	return E;
}

/*  @@ To  

*/
LKEdgeY LinKern::get_y(vector<VERTEX>& spTour, VERTEX p2 , VERTEX p3)
{
	LKEdgeY  ret;
	if (p2.tarID == 0) { // t2 is a depot and t3 is a cirlce
		ret = get_y_type1(spTour, p2, p3);
	}

	if (p3.tarID == 0) {
		ret = get_y_type2(spTour, p2, p3);
	}

	if (p2.tarID != 0 && p3.tarID != 0) {
		ret = get_y_type3(spTour, p2, p3);
	}
	return ret;

}

/*
	y1: when t2 is the depot and t3 is a circle 
*/
LKEdgeY LinKern::get_y_type1(vector<VERTEX>& spTour, VERTEX p2, VERTEX p3)
{
	int subKey_pivot, matKey_pivot, subKey_bestPivot, matKey_bestPivot;

	double dist = INFINITY;
	double temp;

	// There are k points on circle t3, which means there are only at most k-1 different paths from p2 to p3.
	for (int i = 1; i < num_dstzn; i++) { // start counting from subKey_p3 + 1
		subKey_pivot = (p3.subKey + i) % num_dstzn;
		matKey_pivot = (p3.tarID - 1)* num_dstzn + subKey_pivot + 1;
		temp = D[0][matKey_pivot] + intDist[p3.tarID][i];
		if (temp < dist) {
			dist = temp;
			subKey_bestPivot = subKey_pivot;
			matKey_bestPivot = matKey_pivot;
		}	
	}

	VERTEX newPoint_t2 = { 0, -1, 0 };
	VERTEX newPoint_t3 = { p3.tarID, subKey_bestPivot, matKey_bestPivot };
	
	return{ dist, newPoint_t2, newPoint_t3 };

}
 
/* when t2 is a circle and t3 is the depot */
LKEdgeY LinKern::get_y_type2(vector<VERTEX>& spTour, VERTEX p2, VERTEX p3)
{
	int id_t2 = p2.tarID;
	int subKey_p2 = p2.subKey;
	int matKey_p2 = p2.matKey;

	int matKey_p3 = 0;

	int subKey_pivot, matKey_pivot, subKey_bestPivot, matKey_bestPivot;

	double dist = INFINITY;
	double temp;
	
	for (int i = 1; i < num_dstzn; i++) { // 
		subKey_pivot = (subKey_p2 + i) % num_dstzn;
		matKey_pivot = (id_t2 - 1)* num_dstzn + subKey_pivot + 1;
		temp = intDist[id_t2][i] + D[matKey_pivot][matKey_p3];
		if (temp < dist) {
			dist = temp;
			subKey_bestPivot = subKey_pivot;
			matKey_bestPivot = matKey_pivot;
		}
	}

	VERTEX newPoint_t2 = { id_t2, subKey_bestPivot, matKey_bestPivot };
	VERTEX newPoint_t3 = { 0, -1, 0 };
	
	return{ dist, newPoint_t2, newPoint_t3 };

}


/*
	when p2 and p3 are both located on circular boundaries
	Return a LKEdgeY. Here, LKEdgeY is a path of p2 -> newpoint_t2 -> newpoint_t3 -> p3
*/
LKEdgeY LinKern::get_y_type3(vector<VERTEX>& spTour, VERTEX p2, VERTEX p3)
{
	int matKey_p2 = p2.matKey;
	int matKey_p3 = p3.matKey;

	int id_t2 = p2.tarID;
	int id_t3 = p3.tarID;

	int subKey_pivot2, matKey_pivot2, matKey_bestPivot2;
	int subKey_pivot3, matKey_pivot3, matKey_bestPivot3;
	
	double temp;
	vector<int> backtrack(num_dstzn, -1);
	vector<double> vecSD(num_dstzn, INFINITY);

	for (int i = 0; i < num_dstzn; i++) { // t2
		subKey_pivot2 = (p2.subKey + i) % num_dstzn; 
		matKey_pivot2 = (id_t2 - 1)* num_dstzn + subKey_pivot2 + 1; // get the pivoting point on t2
		for (int j = 0; j < num_dstzn; j++) {
			subKey_pivot3 = (p3.subKey + j) % num_dstzn; 
			matKey_pivot3 = (id_t3 - 1)* num_dstzn + subKey_pivot3 + 1; 
			temp = intDist[id_t2][i] + D[matKey_pivot2][matKey_pivot3];
			if (temp < vecSD[j]) {
				vecSD[j] = temp;
				backtrack[j] = i;
			}
		}
	}

	int bkt;
	double sd = INFINITY;
	for (int k = 1; k < num_dstzn; k++) { // t3
		temp = vecSD[k] + intDist[id_t3][k];
		if (temp < sd) {
			sd = temp;
			subKey_pivot3 = (p3.subKey + k) % num_dstzn; // pivot's subkey
			matKey_pivot3 = (id_t3 - 1) * num_dstzn + subKey_pivot3 + 1; // compute this pivot 's matrix index
			bkt = k;
		}
	}

	subKey_pivot3 = (p3.subKey + bkt) % num_dstzn; // pivot's subkey
	matKey_pivot3 = (id_t3 - 1)* num_dstzn + subKey_pivot3 + 1; // compute this pivot 's matrix index
	
	subKey_pivot2 = (p2.subKey + backtrack[bkt]) % num_dstzn; // pivot's subkey
	matKey_pivot2 = (id_t2 - 1)* num_dstzn + subKey_pivot2 + 1; // compute this pivot 's matrix index

	VERTEX newPoint_t2 = { id_t2, subKey_pivot2, matKey_pivot2 };
	VERTEX newPoint_t3 = { id_t3, subKey_pivot3, matKey_pivot3 };

	return {sd, newPoint_t2, newPoint_t3 };

}

// search by key
bool LinKern::isInSet(vector<pair<int, int>>& brokenSet, int id_ti)
{
	int lens = brokenSet.size();
	bool flag = false;
	for (int i = 0; i < lens; i++) {
		if (id_ti == brokenSet[i].first || id_ti == brokenSet[i].second) {
			flag = true;
			break;
		}
	}
	return flag;
}

bool LinKern::updateGStar(double& G_star, double Gi, double x, double y_star) {

	double g_star = x - y_star;
	if (Gi + g_star > G_star) {
		G_star = Gi + g_star;
		return true;
	}
	return false;
}



void LinKern::updatePoints(vector<VERTEX>& tempSPTour, LKEdgeX x1, LKEdgeX x2, LKEdgeY y1, LKEdgeY yStar) 
{
	// get four points
	int q1 = x1.spKeyB;
	int q2 = x1.spKeyC;
	int q3 = x2.spKeyB;
	int q4 = x2.spKeyC;

	// update q2 and q3
	tempSPTour[q2] = y1.newPointA;
	tempSPTour[q3] = y1.newPointB;

	// update q1 and q4
	tempSPTour[q1] = yStar.newPointA;
	tempSPTour[q4] = yStar.newPointB;

}

int LinKern::getClosingUpNode(vector<VERTEX>& tempFlipTour, vector<pair<TARGET, TARGET>>& X, int id_t3) {
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

	int id_t4 = tempFlipTour[2 * key_t4].tarID;

	return id_t4;
}

void LinKern::flip(vector<VERTEX>& tempFlipTour, int id_a, int id_b, int id_c, int id_d) {
	vector<VERTEX> valVec;
	vector<int> keyVec;
	int key_a = getKey(tempFlipTour, id_a);
	int key_b = getKey(tempFlipTour, id_b);
	int key_c = getKey(tempFlipTour, id_c);
	int key_d = getKey(tempFlipTour, id_d);

	int tempKey;
	if (trueKey(key_a + 1) == key_b) { // cw
		int k = 0;
		while (true)
		{
			tempKey = trueKey(key_b + k);
			keyVec.push_back(2* tempKey);
			valVec.push_back(tempFlipTour[2* tempKey]);
			keyVec.push_back(2 * tempKey + 1);
			valVec.push_back(tempFlipTour[2 * tempKey + 1]);
			if (tempKey == key_d) {
				break;
			}
			k++;
		}
		for (int i = 0; i < keyVec.size(); i++) {//reverse the segment we find
			tempFlipTour[keyVec[i]] = valVec[keyVec.size() - i - 1];
		}
	}
	else { // ccw
		int k = 0;
		while (true)
		{
			tempKey = trueKey(key_a + k);
			keyVec.push_back(2 * tempKey);
			valVec.push_back(tempFlipTour[2 * tempKey]);
			keyVec.push_back(2 * tempKey + 1);
			valVec.push_back(tempFlipTour[2 * tempKey + 1]);
			if (tempKey == key_c) {
				break;
			}
			k++;
		}
		for (int i = 0; i < keyVec.size(); i++) {//reverse the other segment
			tempFlipTour[keyVec[i]] = valVec[keyVec.size() - i - 1];
		}
	}

}

int LinKern::getSPKey(vector<VERTEX>& spTour, VERTEX p) 
{
	for (int i = 0; i < spTour.size(); i++) {
		if (spTour[i].matKey == p.matKey)
			return i;
	}
}

int LinKern::getKey(vector<VERTEX>& spTour, int id) {
	int k = 0;
	while (true) {
		if (id == spTour[k].tarID) {
			return floor(k/2);
			break;
		}
		k++;
	}
}


int LinKern::trueKey(int key_t)
{
	if (key_t > N - 1)
		return key_t % N;
	else if (key_t < 0)
		return N + key_t;
	else
		return key_t;
}

int LinKern::matrixKey(int id, int i, bool entry_flag)
{
	if (entry_flag == true) {
		return 2 * (id - 1) * num_dstzn + i + 1;
	}
	else {
		return 2 * (id - 1) * num_dstzn + num_dstzn + i + 1;
	}
	
}





/*
@@ Initialize a random tour for LinKern algorithm
Note that: the first element is always 0, referring to depot
*/
vector<int> LinKern::createRandomTour() {
	vector<int> tspSeq;
	// create a sequence of increasing order
	for (int i = 0; i < N; i++) {
		tspSeq.push_back(i);
	}
	random_device rd;
	mt19937 gen(rd());
	srand(gen());
	int count = 0;
	int randIdx1, randIdx2;
	while (true)
	{
		// only swap a pair of positive index. (Keep zero at Tour[0])
		randIdx1 = (rand() % (N - 1)) + 1;
		randIdx2 = (rand() % (N - 1)) + 1;
		//	cout << randIdx1 << "  " << randIdx2 << endl;
		if (randIdx1 != randIdx2) {
			swap(tspSeq, randIdx1, randIdx2);
			count++;
		}

		if (count > N) {
			return tspSeq;
			break;
		}
	}
}

/*
@@ Swap i-th element with j-th element
*/
void LinKern::swap(vector<int>& Tour, int i, int j) {
	if (i < 0 || i > N - 1 || j < 0 || j > N - 1) {
		cerr << " ERROR: index out of range!! " << endl;
	}
	int temp = Tour[i];
	Tour[i] = Tour[j];
	Tour[j] = temp;
}

/*  @@ check whether id is a neighbor of Target p.

*/
bool LinKern::isNeighbor(vector<VERTEX>& spTour, int id_p, int id) {
	int key_p = getKey(spTour, id_p);
	if (id == spTour[2 * trueKey(key_p - 1)].tarID || id == spTour[2 * trueKey(key_p + 1)].tarID) {
		return true;
	}
	return false;
}

double LinKern::getTourLens(vector<VERTEX>& spTour) {
	double L = 0.0;
	for (int i = 0; i < 2 * N - 1; i++) {
		L += D[spTour[i].matKey][spTour[i + 1].matKey];
	}

	L += D[spTour[2*N-1].matKey][spTour[0].matKey];

	return L;
}

vector<int> LinKern::getTspTour(vector<VERTEX>& spTour) {

	double L = 0.0;
	vector<int> tspTour;
	int idx, id;
	for (int i = 0; i < 2 * N; i++) {
		if (spTour[i].tarID == 0) {
			idx = i;
			break;
		}
	}
	tspTour.push_back(0);
	while (true){
		idx += 2;
		id = spTour[(idx % (2 * N))].tarID;
		if (id == 0) {
			break;
		}

		tspTour.push_back(id);

	}
	return tspTour;
}

LinKern::~LinKern() {};