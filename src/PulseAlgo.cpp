#include "PulseAlgo.h"
#define _USE_MATH_DEFINES
#include<bits/stdc++.h>
void Pulse::initialize(DataHandler& dataset, int tar, int nb_dstzn, int lmt_obps){
	this->_dataset = &dataset;
	this->_taridx = tar;
	this->_graphsize = lmt_obps + 2;
	this->_source = 0;
	this->_sink = lmt_obps + 1;

	_rewards_etd.resize(_graphsize, 0);
	_total_reward = 0.0;
	for(int i=1; i<= lmt_obps; i++){
		_rewards_etd[i] = _dataset->_all_obps[_taridx][i-1]._rewVal;
		_total_reward += _rewards_etd[i];
		// cout << _rewards_etd[i] << " ";
	}

	// cout << '\n';
	// cout << _total_reward << endl;

	_budget_pct = _dataset->_bdg_rewards_ratio[_taridx];
	_budget_pct = 0.3;

	_budget = _total_reward * _budget_pct;
	// cout << "total reward: " << _total_reward << '\t' ;
	// cout << "_budget: " << _budget << endl;

	_graph.resize(_graphsize);
	for(int i=0; i < _graphsize; i++)
		_graph[i].resize(_graphsize, INF);
	Vertex p_i, p_j;
	for(int i=1; i <= lmt_obps; i++){
		for(int j=i+1; j <= lmt_obps; j++){
			p_i._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i-1]._rad * cos(_dataset->_all_obps[_taridx][i-1]._angle);
			p_i._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i-1]._rad * sin(_dataset->_all_obps[_taridx][i-1]._angle);
			p_j._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][j-1]._rad * cos(_dataset->_all_obps[_taridx][j-1]._angle);
			p_j._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][j-1]._rad * sin(_dataset->_all_obps[_taridx][j-1]._angle);
			_graph[i][j] = eucl_distance(p_i, p_j);
			_graph[j][i] = _graph[i][j];		
		}	
	}

	this->_nb_tnps = nb_dstzn;
	_tnps_tar.resize(_nb_tnps);
	double subarc_angle = 2.0 * M_PI /(double)_nb_tnps;
	for(int i =0; i < _nb_tnps; i++){
		_tnps_tar[i]._x = _dataset->_radii[_taridx] * cos(i * subarc_angle) + _dataset->_target_locs[_taridx]._x;
		_tnps_tar[i]._y = _dataset->_radii[_taridx] * sin(i * subarc_angle) + _dataset->_target_locs[_taridx]._y;
		// _tnps_tar[i].print();
	}

}


void Pulse::update_graph(Vertex& entry, Vertex& exit){
		// source and sink
	// cout << "entry and exit :" << endl;
	// entry.print();
	// exit.print();
	Vertex p;
	for(int i=1; i <= _graphsize-2; i++){
		p._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i-1]._rad * cos(_dataset->_all_obps[_taridx][i-1]._angle);
		p._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i-1]._rad * sin(_dataset->_all_obps[_taridx][i-1]._angle);
		// p.print();
		_graph[_source][i] = eucl_distance(entry, p);
		_graph[i][_source] = _graph[_source][i];
		_graph[_sink][i] = eucl_distance(exit, p);
		_graph[i][_sink] = _graph[_sink][i];
	}


	// for(int i=0; i < _graphsize ; i++){
	// 	for(int j=0; j < _graphsize; j++){
	// 		cout << _graph[i][j] << '\t';
	// 	}	
	// 	cout << '\n';
	// }
}

double Pulse::eucl_distance(const Vertex& p1, const Vertex& p2){
	return sqrt(pow(p1._x -p2._x,2) + pow(p1._y - p2._y,2));
}

void Pulse::recursive_search(int curnode, vector<int> curpath, double pathreward, double pathrisk, bool log_on){
	_iterations++;
	if (log_on){
		cout << '\n';
		cout << "@@" << _iterations << ". current path: [";
		for(unsigned i =0; i < curpath.size(); i++){
			cout << curpath[i] << ',';
		}
		cout << "], path reward: " << pathreward << ", path risk: " << pathrisk <<  ", consider moving to node " << curnode <<  endl;
	}
	if(curnode == _sink){
		if(pathreward >= _budget && pathrisk + _graph[curpath.back()][_sink] < _curbest_objval){
			_curbest_objval = pathrisk + _graph[curpath.back()][_sink];
			curpath.push_back(_sink);
			_curbest_path.clear();
			_curbest_path = curpath;
			// cout << "obj* (sink), " << _curbest_objval << " @iterations " << _iterations << endl;

			if(log_on){
				cout << "!!!!!!!GET AN IMPROVEMENT: new obj: " << _curbest_objval << '\t' << " new path: [" << endl;
				for(unsigned i=0; i<_curbest_path.size(); i++){
					cout << _curbest_path[i] << ' ';
				}
				cout << "\n"; 
			}
		}else{
			if(log_on){
				cout << "... reach the sink but no improvement " << endl;
			}
		}
		return;
	}
	if(curnode != _source){
		int lastnode = curpath.back();
		if(log_on)
			cout << "CYCLE CHECK:";
		if(std::find(curpath.begin(), curpath.end(), curnode) == curpath.end()){ 
			if(log_on)
				cout << "no cycle " << endl;
			if(pathrisk + _graph[lastnode][curnode] < _curbest_objval){
				if(log_on)
					cout << "... still a valid partial path: check its budget " << endl;
				if(pathreward + _rewards_etd[curnode] >= _budget){
					if(log_on)
						cout << "... ENOUGH BUDGET IS MET "<< endl;
					/*quick reach to sink */
					auto ret = quick_wrapup(curpath, curnode);
					// vector<int> path_sink({curnode,_sink});
					// pair<double,vector<int>> ret = make_pair(_graph[curnode][_sink], path_sink);
					if (pathrisk + _graph[lastnode][curnode] + ret.first < _curbest_objval){
						_curbest_objval = pathrisk + _graph[lastnode][curnode] + ret.first;
						_curbest_path.clear();
						_curbest_path.insert(_curbest_path.end(), curpath.begin(), curpath.end());
						_curbest_path.insert(_curbest_path.end(), ret.second.begin(), ret.second.end());
						if (ret.second.size() > 2)
							cout << "obj* (wrapup), " << _curbest_objval << " @iterations " << _iterations << endl;
						else
							cout << "obj* (sink), " << _curbest_objval << " @iterations " << _iterations << endl;

						if(log_on){
							cout << "!!!!!!!!find a better path. Total risk: " << _curbest_objval << " new path: ";
							for(unsigned i=0; i<_curbest_path.size(); i++){
								cout << _curbest_path[i] << ' ';
							}
							cout << "curbest_obj: " << _curbest_objval <<  "\n"; 
						}							
					}
				}else{
					if(log_on){
						cout << "... NOT ENOUGH BUDGET " << endl;
						cout << "... VOID DOMINANCE CHECK" << endl;
					}
					for(int ngb=1; ngb <_graphsize; ngb++){
						if(_graph[curnode][ngb] < 1000.0 ){
							if(std::find(curpath.begin(), curpath.end(), ngb) == curpath.end()){
								auto newpath = curpath;
								newpath.push_back(curnode);
								double newpathrisk = pathrisk + _graph[lastnode][curnode];
								double newpathreward = pathreward + _rewards_etd[curnode];
								if(log_on){
									cout << "... " << ngb << " is not visited. Move to it. newpath: [";
									for(unsigned i=0; i<newpath.size(); i++)
										cout << newpath[i] << ' ';
									cout << ']' << endl;
								}								
								this->recursive_search(ngb, newpath, newpathreward, newpathrisk, log_on);
							}else{
								if(log_on)
									cout << "... " << ngb << " is already visited" << endl;
							}							
						}
					}

				}

			}else{
				if(log_on)
					cout << "... deleted by primal bound " << endl;
			}

		}else{
			if(log_on)
				cout << "... cycle occurs " << endl;
		}


	}else{ 
		/* current node is the source. Propogate to every neighbor */
		vector<int> curpath = {0};
		for(int ngb=0; ngb<_graphsize; ngb++){
			if(_graph[curnode][ngb] < 10000.0){
				cout << "... Propogate from source to node " << ngb << endl;
				this->recursive_search(ngb, curpath, 0.0, 0.0, log_on);
			}
		}
	}
}


void Pulse::calc_leastrisk_sink(){
	vector<bool> visited(_graphsize, false);
	_lbrisk_sink.resize(_graphsize, INF);
	_lbrisk_sink[_sink] = 0.0;
	// cout << "sink: " << _sink;
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
			cerr << "the orginal graph is disconnected" << endl;
		}
		visited[u] = true;
		for(int v=0; v<_graphsize; v++){
			if(visited[v] == false && _lbrisk_sink[u] + _graph[u][v] < _lbrisk_sink[v]){
				_lbrisk_sink[v] = _lbrisk_sink[u] + _graph[u][v];
			}
		}
	}

	// cout << "least risk path to sink: " ;
	// for(int i=0; i<_graphsize; i++)
	// 	cout << _lbrisk_sink[i] << " ";
	// cout << endl;
}

pair<double, vector<int>> Pulse::quick_wrapup(const vector<int>& visited, int startnode){
	vector<bool> visited_etd(_graphsize, false);
	for(unsigned i=0; i<visited.size(); i++)
		visited_etd[visited[i]] = true;

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
		visited_etd[u] = true;
		if(u == _sink){
			break;
		}
		if(u==-1 && visited_etd[_sink] == false){
			// cerr << "the orginal graph is disconnected" << endl;
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

void Pulse::print_opt_sol(){
	cout << " =====>>> total reward: " << _total_reward << ", required budget: " << _budget << endl;
	cout << " =====>>> optimal obj value: " << _curbest_objval << "  optimal path: ";
	for(unsigned i =0; i < _curbest_path.size(); i++)
		cout << _curbest_path[i] << ' ';
	cout << '\n';
}

void Pulse::initialize(vector<vector<double>>& adjmatx, vector<double>& rewards, int size, double budget_pct){
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
	_budget_pct = budget_pct;
	_budget = _total_reward * _budget_pct;
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
}
