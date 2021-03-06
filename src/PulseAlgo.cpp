#include "PulseAlgo.h"
#define _USE_MATH_DEFINES
#include<bits/stdc++.h>

void Pulse::set_parameters(const int source, const int sink, const double demand_pct){
	this->_source = source;
	this->_sink = sink;
	this->_demand_pct = demand_pct;
	// this->_demand = this->_total_reward * this->_demand_pct;
}


void Pulse::read_PCCSP_instance(string filename){
	// this->_dataset = &dataset;
	// this->_taridx = tar;
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	int _nb_edges = -1;
	file >> _graphsize >> _nb_edges;

	_rewards_etd.resize(_graphsize, 0);
	string::size_type sz = 0;
	string tempstr = " ";
	for (int i = 0; i < _graphsize; i++){
		file >> tempstr;
		auto pos = tempstr.find_first_of(":");
		string temp = tempstr.substr(0, pos);
		int node = stoi(temp.substr(sz));
		pos = tempstr.find_last_of(":");
		temp = tempstr.substr(pos+1, tempstr.size());
		double rewval = stod(temp.substr(sz));
		// cout << node << ", " << rewval << endl;
		_rewards_etd[node] = rewval;
	}
	
	_total_reward = 0.0;
	for(int i=0; i< _graphsize; i++){
		_total_reward += _rewards_etd[i];
		// cout << _rewards_etd[i] << " ";
	}
	// exit(0);

	_graph.resize(_graphsize);
	for(int i=0; i<_graphsize; i++)
		_graph[i].resize(_graphsize, INF);

	for (int i = 0; i < _nb_edges; i++){
		file >> tempstr;
		auto pos = tempstr.find_first_of(":");
		string temp = tempstr.substr(0, pos);
		int s = stoi(temp.substr(sz));
		string temp_left = tempstr.substr(pos+1, tempstr.size());
		pos = temp_left.find_first_of(":");
		temp = temp_left.substr(0, pos);
		int t = stoi(temp.substr(sz));
		temp = temp_left.substr(pos+1, temp_left.size());
		double weight = stod(temp.substr(sz));
		_graph[s][t] = weight;
		_graph[t][s] = weight;
		// cout << s << ", " << t << ", " << weight << endl;
	}
	// exit(0);

	_domi_labels.resize(_graphsize);
	
}

void Pulse::read_PCCSP_instance2(string filename){
	// this->_dataset = &dataset;
	// this->_taridx = tar;
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	int _nb_edges = -1;
	file >> _graphsize >> _nb_edges;

	_rewards_etd.resize(_graphsize, 0);
	string::size_type sz = 0;
	string tempstr = " ";
	for (int i = 0; i < _graphsize; i++){
		file >> tempstr;
		auto pos = tempstr.find_first_of(":");
		string temp = tempstr.substr(0, pos);
		int node = stoi(temp.substr(sz));
		pos = tempstr.find_last_of(":");
		temp = tempstr.substr(pos+1, tempstr.size());
		double rewval = stod(temp.substr(sz));
		// cout << node << ", " << rewval << endl;
		_rewards_etd[node] = rewval;
	}
	
	_total_reward = 0.0;
	for(int i=0; i< _graphsize; i++){
		_total_reward += _rewards_etd[i];
		// cout << _rewards_etd[i] << " ";
	}
	// exit(0);

	_graph.resize(_graphsize);
	_rscgraph.resize(_graphsize);
	for(int i=0; i<_graphsize; i++){
		_graph[i].resize(_graphsize, INF);
		_rscgraph[i].resize(_graphsize, INF);
	}

	for (int i = 0; i < _nb_edges; i++){
		file >> tempstr;
		auto pos = tempstr.find_first_of(":");
		string temp = tempstr.substr(0, pos);
		int s = stoi(temp.substr(sz));
		string temp_left = tempstr.substr(pos+1, tempstr.size());
		pos = temp_left.find_first_of(":");
		temp = temp_left.substr(0, pos);
		int t = stoi(temp.substr(sz));
		temp = temp_left.substr(pos+1, temp_left.size());
		pos = temp.find_first_of(":");
		double weight = stod(temp.substr(0, pos).substr(sz));
		double resource = stod(temp.substr(pos+1, temp.size()).substr(sz));
		_graph[s][t] = weight;
		_graph[t][s] = weight;
		_rscgraph[s][t] = resource;
		_rscgraph[t][s] = resource;
		// cout << s << ", " << t << ", " << resource << endl;
	}
	// exit(0);

	_domi_labels.resize(_graphsize);
	
}
void Pulse::recursive_search(int curnode, vector<int> curpath, double pathreward, double pathrisk, ProgTime& timeout, bool log_on){
	if(_source == _sink)
		cerr << "Current code does not work when source is same as sink!! Regenerate the grpah "<< endl;
	if(true){
		timeout.end_prog();
		if(timeout._elapsed_secs >= 600.0){
			this->_status = false;
			return;
		}
	}
	_iterations++;
	// if(_iterations > 50)
	// 	exit(0);
	if (log_on){
		cout << '\n';
		cout << "@@" << _iterations << ". current path: [";
		for(unsigned i =0; i < curpath.size(); i++){
			cout << curpath[i] << ',';
		}
		cout << "], path reward: " << pathreward << ", path risk: " << pathrisk <<  ", consider moving to node " << curnode <<  endl;
	}
	if(curnode == _sink ){
		if(log_on){
			cout << pathreward << " ? " << _demand << " and "  << pathrisk + _graph[curpath.back()][_sink] << " ? " << _curbest_objval << endl;
		}
		if(pathreward >= _demand && pathrisk + _graph[curpath.back()][_sink] < _curbest_objval){
			_curbest_objval = pathrisk + _graph[curpath.back()][_sink];
			curpath.push_back(_sink);
			_curbest_path.clear();
			_curbest_path = curpath;
			if(log_on){
				cout << "!!!!!!!GET AN IMPROVEMENT: new obj: " << _curbest_objval << '\t' << " new path: [" << endl;
				for(unsigned i=0; i<_curbest_path.size(); i++){
					cout << _curbest_path[i] << ' ';
				}
				cout << "\n"; 
			}
		}else{
			if(log_on){
				cout << "... reach the sink but does not improve  " << endl;
			}
		}
		return;
	}


	if(curnode != _source && curnode != _sink){
		int lastnode = curpath.back();
		if(pathrisk + _lbrisk_sink[curnode] < _curbest_objval){ //primal bound check
			if(log_on)
				cout << "... still a valid partial path: check its budget " << endl;
			bool flag = is_intersected(_lbriskpaths_sink[curnode], curpath);
			if(log_on){
				cout << "curpath: ";
				for(unsigned j=0; j <curpath.size(); j++){
					cout << curpath[j] << ' ';
				}
				cout << "G(path~~>sink): ";
				for(unsigned j=0; j <_lbriskpaths_sink[curnode].size(); j++){
					cout << _lbriskpaths_sink[curnode][j] << ' ';
				}
				cout << "intersected? " << flag << '\n';
			}
			double add_reward = (flag)? 0.0:_rewards_lbrpaths_sink[curnode]-_rewards_etd[curnode];
			if(pathreward + _rewards_etd[curnode] + add_reward +0.000001 >= _demand){
				if(log_on){
					cout << "... ENOUGH BUDGET IS MET ... Pick extended path ";
				}
				pair<double,vector<int>> ret;
				if(flag == true){
					/*expand curpath to a full path by linking to shortest path to sink over subgraph */
					ret = quick_wrapup(curpath, curnode);
					if(log_on){
						cout << "on **Subgraph**: ";
						for(unsigned j=0; j <ret.second.size(); j++){
							cout << ret.second[j] << ' ';
						}
						cout << '\n';
					}
				}else{
					ret = make_pair(_lbrisk_sink[curnode], _lbriskpaths_sink[curnode]);
					if(log_on){
						cout << "on **orginal Graph**: ";
						for(unsigned j=0; j <ret.second.size(); j++){
							cout << ret.second[j] << ' ';
						}
						cout << '\n';
					}
				}

				if (pathrisk + _graph[lastnode][curnode] + ret.first < _curbest_objval){ // primal bound update 
					_curbest_objval = pathrisk + _graph[lastnode][curnode] + ret.first;
					_curbest_path.clear();
					_curbest_path.insert(_curbest_path.end(), curpath.begin(), curpath.end());
					_curbest_path.insert(_curbest_path.end(), ret.second.begin(), ret.second.end());
					if(log_on){
						cout << "!!!!!!!!find a better path. Total risk: " << _curbest_objval << " new path: ";
						for(unsigned i=0; i<_curbest_path.size(); i++){
							cout << _curbest_path[i] << ' ';
						}
						cout << "curbest_obj: " << _curbest_objval <<  "\n"; 
					}							
				}else{
					if(log_on)	
						cout << "... Propogate to sink. No better than cur OBJ." << endl;
				}
			}else{
				if(log_on){
					cout << "... NOT MEET BUDGET YET. NEED TO CONTINUE... " << endl;
				}
				vector<int> newpath(curpath);
				newpath.push_back(curnode);
	
				double newpathrisk = pathrisk + _graph[lastnode][curnode];
				double newpathreward = pathreward + _rewards_etd[curnode];
				// bool is_dominated = check_dominance(curnode, newpathrisk, newpathreward);
				if(true){
					if(log_on)	
						cout <<"... pass dominance check" << endl;
					// update_domilabels(curnode, newpathrisk, newpathreward);
					// update_domilabels(curnode, newpath, newpathrisk, newpathreward);

					// if(_iterations >9206 )
					// 	exit(0);
					for(int ngb=0; ngb <_graphsize; ngb++){
						if(_graph[curnode][ngb] < 1000.0){ //  && abs(_rewards_etd[ngb]) > 0.00001
							if(std::find(newpath.begin(), newpath.end(), ngb) == newpath.end()){			
								if(log_on){
									cout << "... " << ngb << " is not visited. *propogate* Input path: [";
									for(unsigned i=0; i<newpath.size(); i++)
										cout << newpath[i] << ' ';

									cout << ']' << endl;
								}								
								this->recursive_search(ngb, newpath, newpathreward, newpathrisk, timeout, log_on);
							}else{
								if(log_on)
									cout << "... " << ngb << " is already visited" << endl;
							}							
						}
					}
				}else{
					if(log_on)
						cout <<"... pruned by dominance" << endl;
				}
			}

		}else{
			if(log_on)
				cout << "... discarded by primal bound. risk of best full path: " << pathrisk + _lbrisk_sink[curnode] << " curbest:" << _curbest_objval << endl;
		}
	}else{ 
		/* current node is the source. Propogate to every neighbor */
		vector<int> curpath = {_source};
		for(int ngb=0; ngb<_graphsize; ngb++){
			if(_graph[curnode][ngb] < 10000.0 ){ // && abs(_rewards_etd[ngb]) > 0.00001
			// int ngb = 26;
				cout << "... Propogate from source " << _source << " to node " << ngb << endl;
				this->recursive_search(ngb, curpath, _rewards_etd[_source], 0.0, timeout, log_on);
			}
		}
	}
}




void Pulse::HeuRecSearch(int curnode, vector<int> curpath, double pathreward, double pathrisk, double pathrsc, ProgTime& timeout, bool log_on){
	if(_source == _sink)
		cerr << "Current code does not work when source is same as sink!! Regenerate the grpah "<< endl;
	if(true){
		timeout.end_prog();
		if(timeout._elapsed_secs >= 600.0){
			this->_status = false;
			return;
		}
	}
	_iterations++;
	// if (_iterations > 90800)
	// 	return;
	if (log_on){
		cout << '\n';
		cout << "@@" << _iterations << ". current path: [";
		for(unsigned i =0; i < curpath.size(); i++){
			cout << curpath[i] << ',';
		}
		cout << "], path reward: " << pathreward << ", path risk: " << pathrisk <<  ", consider moving to node " << curnode <<  endl;
	}
	if(curnode == _sink ){
		if(log_on){
			cout << pathreward << " ? " << _demand << " and "  << pathrisk + _graph[curpath.back()][_sink] << " ? " << _curbest_objval << endl;
		}
		if(pathreward +_rewards_etd[_sink] >= _demand && pathrisk + _graph[curpath.back()][_sink] < _curbest_objval && pathrsc + _rscgraph[curpath.back()][_sink] <= _budget){
			_curbest_objval = pathrisk + _graph[curpath.back()][_sink];
			curpath.push_back(_sink);
			_curbest_path.clear();
			_curbest_path = curpath;
			if(true){
				cout << "!!!!!!!GET AN IMPROVEMENT: new obj: " << _curbest_objval << '\t' << " new path: [" << endl;
				for(unsigned i=0; i<_curbest_path.size(); i++){
					cout << _curbest_path[i] << ' ';
				}
				cout << "\n"; 
			}
		}else{
			if(log_on){
				cout << "... reach the sink but does not improve  " << endl;
			}
		}
		return;
	}


	if(curnode != _source && curnode != _sink){
		int lastnode = curpath.back();
		// cout << pathrisk << ", " << _graph[lastnode][curnode] << ", " <<  _lbrisk_sink[curnode] << " ? " << _curbest_objval;
		// cout << " --- " << pathrsc  << ", " <<  _rscgraph[lastnode][curnode] << ", " << _lbrsc_sink[curnode] << " ? " << _budget << endl;
		if(pathrisk + _graph[lastnode][curnode] + _lbrisk_sink[curnode] < _curbest_objval && pathrsc + _rscgraph[lastnode][curnode] + _lbrsc_sink[curnode] <= _budget + 0.00001){ 
			if(log_on)
				cout << "... still a valid partial path: check its budget " << endl;

			// bool flag = true;
			// double add_reward = (flag)? 0.0:_rewards_lbrpaths_sink[curnode]-_rewards_etd[curnode];
			// if(pathreward + _rewards_etd[curnode] + 0.000001 >= _demand){
				// if(log_on){
				// 	cout << "... ENOUGH BUDGET IS MET ... Pick extended path ";
				// }
				// // pair<double,vector<int>> ret;
				// // if(flag){
				// 	/*expand curpath to a full path by linking to shortest path to sink over subgraph */
				// 	auto ret = quick_wrapup2(curpath, curnode);
				// 	if(log_on){
				// 		cout << "on **Subgraph**: ";
				// 		for(unsigned j=0; j <ret.second.size(); j++){
				// 			cout << ret.second[j] << ' ';
				// 		}
				// 		cout << '\n';
				// 	}

				// if (pathrisk + _graph[lastnode][curnode] + ret.first.first < _curbest_objval && pathrsc + _rscgraph[lastnode][curnode] + ret.first.second <= _budget){ // update primal bound  
				// 	_curbest_objval = pathrisk + _graph[lastnode][curnode] + ret.first.first;
				// 	_curbest_path.clear();
				// 	_curbest_path.insert(_curbest_path.end(), curpath.begin(), curpath.end());
				// 	_curbest_path.insert(_curbest_path.end(), ret.second.begin(), ret.second.end());
				// 	if(log_on){
				// 		cout << "!!!!!!!!find a better path. Total risk: " << _curbest_objval << " new path: ";
				// 		for(unsigned i=0; i<_curbest_path.size(); i++){
				// 			cout << _curbest_path[i] << ' ';
				// 		}
				// 		cout << "curbest_obj: " << _curbest_objval <<  "\n"; 
				// 	}							
				// }else{
				// 	if(log_on)	
				// 		cout << "... Propogate to sink. No better than cur OBJ." << endl;
				// }
			// }else{

				vector<int> newpath(curpath);
				newpath.push_back(curnode);
	
				double newpathrisk = pathrisk + _graph[lastnode][curnode];
				double newpathreward = pathreward + _rewards_etd[curnode];
				double	newpathrsc = pathrsc + _rscgraph[lastnode][curnode];
				// bool is_dominated = check_dominance(curnode, newpathrisk, newpathrsc);
				if(true){
					if(log_on)	
						cout <<"... pass dominance check" << endl;
					// update_domilabels(curnode, newpathrisk, newpathrsc);
					// update_domilabels(curnode, newpath, newpathrisk, newpathreward);

					for(int ngb=0; ngb <_graphsize; ngb++){
						if(_graph[curnode][ngb] < 1000.0){ //  && abs(_rewards_etd[ngb]) > 0.00001
							if(std::find(newpath.begin(), newpath.end(), ngb) == newpath.end()){			
								if(log_on){
									cout << "... " << ngb << " is not visited. *propogate* Input path: [";
									for(unsigned i=0; i<newpath.size(); i++)
										cout << newpath[i] << ' ';

									cout << ']' << endl;
								}
								// if(_lbrisk_sink[ngb] <= _lbrisk_sink[curnode] * 1.05){								
								this->HeuRecSearch(ngb, newpath, newpathreward, newpathrisk, newpathrsc, timeout, log_on);
								// }
							}else{
								if(log_on)
									cout << "... " << ngb << " is already visited" << endl;
							}							
						}
					}

				}else{
					if(log_on)
						cout <<"... pruned by dominance" << endl;
				}
			// }

		}else{
			if(log_on)
				cout << "... discarded by primal bound. risk of best full path: " << pathrisk + _lbrisk_sink[curnode] << " curbest:" << _curbest_objval << endl;
		}
	}else{ 
		/* current node is the source. Propogate to every neighbor */
		vector<int> curpath = {_source};
		for(int ngb=0; ngb<_graphsize; ngb++){
			if(_graph[curnode][ngb] < 10000.0 ){ // && abs(_rewards_etd[ngb]) > 0.00001
			 // ngb = 15;
				cout << "... Propogate from source " << _source << " to node " << ngb << endl;
				this->HeuRecSearch(ngb, curpath, _rewards_etd[_source], 0.0, 0.0, timeout, log_on);
			// break;
			}
		}
	}
}




void Pulse::update_domilabels(int node, double pathrisk, double pathrsc){
	RWP tmp;
	tmp._risk = pathrisk;
	tmp._resource = pathrsc;
	// tmp._path = vec_to_string(path);
	if(_domi_labels[node].size() == 0){
		_domi_labels[node].insert(tmp);
	}else{
		bool IN = true;
		// set<pair<double,double>>::iterator itr;
		for(auto itr=_domi_labels[node].begin(); itr !=_domi_labels[node].end();){
			if( (pathrisk <= itr->_risk && pathrsc < itr->_resource)|| (pathrisk < itr->_risk && pathrsc <= itr->_resource)){				
				itr = _domi_labels[node].erase(itr); // erase previous dominated path
			}else{
				itr++;
			}
			if(itr->_risk <= pathrisk && itr->_resource <= pathrsc){				
				IN = false;// if there is previous path that dominates or equals current path, not insert this path
			}
		}
		if(_domi_labels[node].size() < _maxlabsize && IN){
			_domi_labels[node].insert(tmp);
		}
	}
}


bool Pulse::check_dominance(int node, double pathrisk, double pathrsc){
	for(auto itr=_domi_labels[node].begin(); itr !=_domi_labels[node].end(); ){
		// cout << node << " cur: " << pathrisk << ", " << pathreward << " with: " << itr->_risk << ", " << itr->_reward << ", " << itr->_path << endl;
		// cout << node << " cur: " << pathrisk << ", " << pathreward << " with: " << itr->_risk << ", " << itr->_reward << endl;

		if((pathrisk >= itr->_risk && pathrsc >= itr->_resource)){		
			// cout << " -------unfinished---------" << endl;
			return true;
		}
		++itr;
	}
	// cout << " -------done---------" << endl;
	return false;

	// if((pathrisk >= _L1[node].first && pathreward <= _L1[node].second)||(pathrisk >= _L2[node].first && pathreward <= _L2[node].second)||(pathrisk >= _L3[node].first && pathreward <= _L3[node].second)){
	// 	return true;
	// }
	// return false;

}


bool Pulse::is_intersected(vector<int>& vec1, vector<int>& vec2){
	for(unsigned i=0; i<vec1.size(); i++){
		for(unsigned j=0; j<vec2.size(); j++){
			if(vec1[i] == vec2[j])
				return true;
		}
	}
	return false;
}


void Pulse::calc_leastrisk_sink(){
	// cout << "====>>>> Algorithm iniialization: least risky path to sink " << endl;
	vector<bool> visited(_graphsize, false);
	// int nb_voidpoints = 0;
	// for(int i=0; i<_graphsize; i++){
	// 	if(abs(_rewards_etd[i])<0.0000001){
	// 		visited[i] = true;
	// 		// nb_voidpoints++;
	// 	}
	// }
	// visited[_source] = false;
	// visited[_sink] = false;
	_lbrisk_sink.resize(_graphsize, INF);
	_lbrisk_sink[_sink] = 0.0;
	// cout << "sink: " << _sink;
	vector<int> prevnode(_graphsize, -2);

	for(int i=0; i< _graphsize; i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited[j]==false && _lbrisk_sink[j] < MIN){
				u = j;
				MIN = _lbrisk_sink[j];
			}
		}
		// cout << "u = " << u << endl;
		if(u==-1){
			cout << "the orginal graph is disconnected" << endl;
			exit(0);
		}
		visited[u] = true;
		for(int v=0; v<_graphsize; v++){
			if(visited[v] == false && _lbrisk_sink[u] + _graph[u][v] < _lbrisk_sink[v]){ 
				_lbrisk_sink[v] = _lbrisk_sink[u] + _graph[u][v];
				prevnode[v] = u;
			}
		}
	}
	// get the optimal path by backtracking
	_lbriskpaths_sink.resize(_graphsize);
	_rewards_lbrpaths_sink.resize(_graphsize, 0.0);
	for(int i=0; i<_graphsize; i++){
		int curnode = i;
		while(true){
			if(curnode == -2)
				break;
			_rewards_lbrpaths_sink[i] += _rewards_etd[curnode];
			_lbriskpaths_sink[i].push_back(curnode);
			curnode = prevnode[curnode];		
		}
	}
	// _rsc_lbrpath_source_sink = 0.0;
	// _rew_lbrpath_source_sink = 0.0;

	// for(int i = 0; i < (int)_lbriskpaths_sink[_source].size()-1; i++){
		// _rsc_lbrpath_source_sink += _rscgraph[_lbriskpaths_sink[_source][i]][_lbriskpaths_sink[_source][i+1]];
		// _rew_lbrpath_source_sink += _rscgraph[_lbriskpaths_sink[_source][i]][_lbriskpaths_sink[_source][i+1]];
	// }
	// _budget =  _rsc_lbrpath_source_sink * 2.0;

	if(false){
		for(int i=0; i<_graphsize; i++){
			cout << " ... node " << i << " to sink " << _sink << ": ";
			for(unsigned j=0; j<_lbriskpaths_sink[i].size();j++){
				cout << _lbriskpaths_sink[i][j] << ',';
			}
			cout << " ---- Total rewards: " << _rewards_lbrpaths_sink[i];
			cout << '\n';
		}
	}

	// _demand =  (_total_reward - _rewards_lbrpaths_sink[_source]) * _demand_pct + _rewards_lbrpaths_sink[_source];

	_demand =  _demand_pct * _rewards_lbrpaths_sink[_source];
	// _demand =  _demand_pct * _total_reward/_graphsize;

	// _demand =  (1.0 + _demand_pct) * _total_reward/_graphsize;

	// _demand =  (_total_reward) * _demand_pct;
	
	if (_demand <= _rewards_lbrpaths_sink[_source]){
		cout << "Reward on min-risk path already satisfies the demand!!!" << endl;
		// exit(0);
	}
}

void Pulse::calc_leastresource_sink(){
	vector<bool> visited(_graphsize, false);

	_lbrsc_sink.resize(_graphsize, INF);
	_lbrsc_sink[_sink] = 0.0;
	vector<int> prevnode(_graphsize, -2);

	for(int i=0; i< _graphsize; i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited[j]==false && _lbrsc_sink[j] < MIN){
				u = j;
				MIN = _lbrsc_sink[j];
			}
		}
		// cout << "u = " << u << endl;
		if(u==-1){
			cout << "the orginal graph is disconnected" << endl;
			exit(0);
		}
		visited[u] = true;
		for(int v=0; v<_graphsize; v++){
			if(visited[v] == false && _lbrsc_sink[u] + _rscgraph[u][v] < _lbrsc_sink[v]){ 
				_lbrsc_sink[v] = _lbrsc_sink[u] + _rscgraph[u][v];
				prevnode[v] = u;
			}
		}
	}
	// get the optimal path by backtracking
	_lbrscpaths_sink.resize(_graphsize);
	_rewards_lbrscpaths_sink.resize(_graphsize, 0.0);
	for(int i=0; i<_graphsize; i++){
		int curnode = i;
		while(true){
			if(curnode == -2)
				break;
			_rewards_lbrscpaths_sink[i] += _rewards_etd[curnode];
			_lbrscpaths_sink[i].push_back(curnode);
			curnode = prevnode[curnode];		
		}
	}
	// _budget =  _lbrsc_sink[_source] * 5;

	if(false){
		for(int i=0; i<_graphsize; i++){
			cout << " ... node " << i << " to sink " << _sink << ": ";
			for(unsigned j=0; j<_lbrscpaths_sink[i].size();j++){
				cout << _lbrscpaths_sink[i][j] << ',';
			}
			cout << " ---- consumed resource: " << _lbrsc_sink[i];
			cout << '\n';
		}
	}
	// _budget =  (_lbrscpaths_sink[_source] - _rewards_lbrscpaths_sink[_source]) * _demand_pct + _rewards_lbrscpaths_sink[_source];
	
	// if (_demand <= _rewards_lbrscpaths_sink[_source]){
		// cout << "Reward on min-risk path already satisfies the demand!!!" << endl;
		// exit(0);
	// }
}




pair<double, vector<int>> Pulse::quick_wrapup(const vector<int>& visited, int startnode){
	vector<bool> visited_etd(_graphsize, false);
	for(unsigned i=0; i<visited.size(); i++)
		visited_etd[visited[i]] = true;
	// for(unsigned i = 0; i < visited.size(); i++)
	// 	cout << visited[i] << ' ';
	// cout << '\n';
	vector<int> prevnode(_graphsize, -2);
	vector<double> dist(_graphsize, INF);
	dist[startnode] = 0.0;
	for(int i=0; i< _graphsize-(int)visited.size(); i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited_etd[j]==false && dist[j] < MIN){
				u = j;
				MIN = dist[j];
			}
		}
		if(u==-1 && visited_etd[_sink] == false){
			// cerr << "the orginal graph is disconnected" << endl;
			break;
		}
		visited_etd[u] = true;
		if(u == _sink){
			break;
		}
		for(int v=0; v<_graphsize; v++){
			if(visited_etd[v] == false && dist[u] + _graph[u][v] < dist[v]){
				dist[v] = dist[u] + _graph[u][v];
				prevnode[v] = u;
			}
		}


	}

	// cout << "prevnode: ";
	// for(int i=0; i<_graphsize; i++)
	// 	cout << prevnode[i] << ' ';
	// cout << '\n';
	// cout << "dist: ";
	// for(int i=0; i<_graphsize; i++)
	// 	cout << dist[i] << ' ';
	// cout << '\n';
	vector<int> optpath;
	if(visited_etd[_sink]){
		// get the optimal path by backtracking
		int curnode = _sink;
		while(curnode!= -2){
			optpath.insert(optpath.begin(),curnode);		
			curnode = prevnode[curnode];
		}
	}
	// cout << "path~~>sink: ";
	// for(unsigned i=0; i<optpath.size(); i++)
	// 	cout << optpath[i] << ' ';
	// cout << '\n';
	return make_pair(dist[_sink], optpath);
}



pair<pair<double,double>, vector<int>> Pulse::quick_wrapup2(const vector<int>& visited, int startnode){
	vector<bool> visited_etd(_graphsize, false);
	for(unsigned i=0; i<visited.size(); i++)
		visited_etd[visited[i]] = true;
	// for(unsigned i = 0; i < visited.size(); i++)
	// 	cout << visited[i] << ' ';
	// cout << '\n';
	vector<int> prevnode(_graphsize, -2);
	vector<double> dist(_graphsize, INF);
	dist[startnode] = 0.0;
	for(int i=0; i< _graphsize-(int)visited.size(); i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited_etd[j]==false && dist[j] < MIN){
				u = j;
				MIN = dist[j];
			}
		}
		if(u==-1 && visited_etd[_sink] == false){
			// cerr << "the orginal graph is disconnected" << endl;
			break;
		}
		visited_etd[u] = true;
		if(u == _sink){
			break;
		}
		for(int v=0; v<_graphsize; v++){
			if(visited_etd[v] == false && dist[u] + _graph[u][v] < dist[v]){
				dist[v] = dist[u] + _graph[u][v];
				prevnode[v] = u;
			}
		}


	}


	vector<int> optpath;
	if(visited_etd[_sink]){
		// get the optimal path by backtracking
		int curnode = _sink;
		while(curnode!= -2){
			optpath.insert(optpath.begin(),curnode);		
			curnode = prevnode[curnode];
		}
	}

	double future_rsc = 0.0;
	for(auto itr = optpath.begin(); itr < optpath.end()-1; itr++){
		future_rsc += _rscgraph[*itr][*(itr+1)];
	}

	return make_pair(make_pair(dist[_sink], future_rsc), optpath);
}



void Pulse::initialize_generalgraph(vector<vector<double>>& adjmatx, vector<double>& rewards, int size, double budget_pct){
	this->_graphsize = size;
	this->_source = 0;
	this->_sink = size -1;

	_rewards_etd = rewards;
	_total_reward = 0.0;
	for(int i=0; i< size; i++){
		// cout << _rewards_etd[i] << " ";
		_total_reward += _rewards_etd[i];
	}
	// cout << "total " << _total_reward << endl;
	_demand_pct = budget_pct;
	_demand = _total_reward * _demand_pct;
	_graph = adjmatx;

	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			if (_graph[i][j] == 0){
				_graph[i][j] = INF;
			}
			// cout << _graph[i][j] << " ";
		}
		// cout << '\n';
	}
	// _domi_labels.resize(_graphsize, make_pair(INF, 0));
	
}


void Pulse::initialize_triGraph(DataHandler& dataset, int tar, int lmt_obps, string graphfile){
	this->_dataset = &dataset;
	this->_taridx = tar;
	this->_graphsize = lmt_obps;

	_rewards_etd.resize(_graphsize, 0);
	_total_reward = 0.0;
	for(int i=0; i< _graphsize; i++){
		_rewards_etd[i] = _dataset->_all_obps[_taridx][i]._rewVal;
		_total_reward += _rewards_etd[i];
	}
	// _demand_pct = _dataset->_bdg_rewards_ratio[_taridx];
	// _demand = _total_reward * _demand_pct;

	fstream file(graphfile);
	if (!file) {
		cerr << "ERROR: could not open file '" << graphfile << "' for reading'" << endl;
		throw(-1);
	}

	_tmp_graph.resize(_graphsize);
	for(int i=0; i<_graphsize; i++)
		_tmp_graph[i].resize(_graphsize, false);

	string::size_type sz = 0;
	int nb_edges = -1;
	for (int i = 0; i < 6; i++){
		file >> nb_edges;
		vector<string> emptystr(nb_edges);
		for(int j = 0; j < nb_edges; j++){
			file >> emptystr[j];
			// cout << emptystr[j] << ' ';
		}
		// cout << "---------------------------" << endl;
		if(i == tar){
			for(int j = 0; j < nb_edges; j++){
				auto pos = emptystr[j].find_first_of(":");
				string temp = emptystr[j].substr(0, pos);
				int u = stoi(temp.substr(sz));
				string temp_left = emptystr[j].substr(pos+1, emptystr[j].size());
				int v = stoi(temp_left.substr(sz));
				// cout << " u, v = " << u << ", " << v << '\t';
				_tmp_graph[u][v] = true;
				_tmp_graph[v][u] = true;
			}
		}	
	}

	// exit(0);
	_graph.resize(_graphsize);
	for(int i=0; i < _graphsize; i++)
		_graph[i].resize(_graphsize, INF);
	Vertex p_i, p_j;
	for(int i=0; i < _graphsize; i++){
		for(int j=i+1; j < _graphsize; j++){
			if(_tmp_graph[i][j] == true){
				p_i._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i]._rad * cos(_dataset->_all_obps[_taridx][i]._angle);
				p_i._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i]._rad * sin(_dataset->_all_obps[_taridx][i]._angle);
				p_j._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][j]._rad * cos(_dataset->_all_obps[_taridx][j]._angle);
				p_j._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][j]._rad * sin(_dataset->_all_obps[_taridx][j]._angle);
				double dist = eucl_distance(p_i, p_j);
				_graph[i][j] = dist;
				_graph[j][i] = dist;	
			}	
		}
	}	
	
	_domi_labels.resize(_graphsize);
	
	_L1.resize(_graphsize,make_pair(INF, 0.0));
	_L2.resize(_graphsize,make_pair(INF, 0.0));
	_L3.resize(_graphsize,make_pair(INF, 0.0));

}


double Pulse::eucl_distance(const Vertex& p1, const Vertex& p2){
	return sqrt(pow(p1._x -p2._x,2) + pow(p1._y - p2._y,2));
}



double Pulse::calc_path_rewards(vector<int> & path){
	double total_reward = 0.0;
	for(unsigned i=0; i < path.size(); i++){
		total_reward += _rewards_etd[path[i]];
	}
	return total_reward;
}


void Pulse::print_opt_sol(){
	cout << " =====>>> total reward: " << _total_reward << ", Demand: " << _demand << ", Reward on min-risk path: " << _rewards_lbrpaths_sink[_source] << endl;
	// cout << " =====>>> Budget: " <<  _budget << endl;
	
	cout << " =====>>> Optimal path [";
	double verified_obj = 0;
	// double used_resouce = 0.0;
	for(unsigned i =0; i < _curbest_path.size(); i++){
		cout << _curbest_path[i] << ' ';
		if(i < _curbest_path.size()-1){
			verified_obj += _graph[_curbest_path[i]][_curbest_path[i+1]];
			// used_resouce += _rscgraph[_curbest_path[i]][_curbest_path[i+1]];
		}
	}
	cout << "]\n";
	cout << " =====>>> optimal obj value: " << _curbest_objval << ", collected reward: " << calc_path_rewards(_curbest_path) << endl; 
}


// void Pulse::write_opt_sol(){
// 	ofstream file;
// 	string cur_dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/ret/pulse_out/";
//     file.open(cur_dir + name_only + "_algo_" + to_string(algo_idx) + ".FisTest_" + to_string(fischetti_on));
// 	file << "---> instance_name: " << instance << '\n';
// 	file.close();
// }


string Pulse::vec_to_string(vector<int>& path){
	string str;
	for(auto itr=path.begin(); itr!=path.end(); itr++){
		str += to_string(*itr) + ",";
	}
	return str;
}