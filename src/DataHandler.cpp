#include "DataHandler.h"
#include "MacroName.h"
// #define is_server false


void DataHandler::read_instance(string filename) {
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	// Geographic Information of Targets
	auto pos = filename.find_last_of("/");
    _name = filename.substr(pos+1, filename.size());
    // size_t found = _name.find("c");
    double temp;
	file >> _nb_targets >> temp;
	file >> _depot1_loc._x >> _depot1_loc._y;
	file >> _depot2_loc._x >> _depot2_loc._y;
	_target_locs.resize(_nb_targets);
	_radii.resize(_nb_targets, -1.0);
	for (int i = 0; i < _nb_targets; i++) {
		file >> _target_locs[i]._x >> _target_locs[i]._y >> _radii[i];
	}
	file.close();
	// string dir = "";
	// if(is_server){
	// 	dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/dat/InnerPaths/";
	// }else{
	// 	dir = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerPaths/";
	// }
	pick_trigraphs(DIR);
}


void DataHandler::pick_trigraphs(string dir){
	vector<int> pickarr(_nb_targets, 0);
	for(int i=0; i < _nb_targets; i++)
		pickarr[i] = i+1;

	// for(auto itr=pickarr.begin(); itr != pickarr.end(); itr++){
	// 	cout << *itr << ',';
	// }
	// cout << '\n';


	_all_innermatr.resize(_nb_targets);
	string pathstr = " ";
	for(int t=0; t<_nb_targets; t++){
		_all_innermatr[t].resize(16);
		for(int i=0; i<16; i++){
			_all_innermatr[t][i].resize(16, make_pair(0.0, pathstr));
		}		
	}

	int index = 0;
	for(auto itr = pickarr.begin(); itr != pickarr.end(); itr++){
		string filename = dir + "1kCSP-"+to_string((*itr)) + ".dat";
		fstream file(filename);
		if (!file) {
			cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
			throw(-1);
		}
		int tmp = -1;
		int e1 = -1;
		int e2 = -1;
		double riskval = 0;
		for(int i = 0; i<120; i++){
			file >> tmp >> e1 >> e2 >> riskval >> pathstr; 
			// cout << e1 << " " << e2 << " " << riskval << ", " << pathstr << endl;
			_all_innermatr[index][e1][e2] = make_pair(riskval, pathstr);
			_all_innermatr[index][e2][e1] = make_pair(riskval, pathstr); 
		}
		index++;
	}
}

void DataHandler::build_riskgraph(){
	_points.resize(_nb_targets + 2);
	_points[0].resize(1); //departure point
	_points[0][0] = _depot1_loc;
	// _points[0][0].print();
	double _subarc_angle = 2.0 * M_PI / (double)_nb_dstzn;
	for (int pos = 0; pos < _nb_targets; pos++) {
		_points[pos+1].resize(_nb_dstzn);
		for (int i = 0; i < _nb_dstzn; i++) {
			_points[pos+1][i]._x = _radii[pos] * cos(i * _subarc_angle) + _target_locs[pos]._x;
			_points[pos+1][i]._y = _radii[pos] * sin(i * _subarc_angle) + _target_locs[pos]._y;
			// _points[pos+1][i].print();
		}
	}
	_points[_nb_targets + 1].resize(1); // landing depot
	_points[_nb_targets + 1][0] = _depot2_loc;
	// _points[_dataset->_nb_targets + 1][0].print();

	_graphsize = 2 * _nb_targets * _nb_dstzn + 2;
	_riskgraph.resize(_graphsize);
	for (int i = 0; i < _graphsize; i++) {
		_riskgraph[i].resize(_graphsize, make_pair(false,INF));
	}

	/* Generate the accumulated risk on outer paths*/
	int idxmat_1, idxmat_2, idxmat_3, idxmat_4;
	double val_risk;
	/* i) departure depot --> each entry turning point of a target */	
	for (int t = 1; t <= _nb_targets; t++) {	
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + 1; 
		for (int i = 0; i < _nb_dstzn; i++) {
			// val_risk = get_risk_outerTrajc(_points[t][i], _depot1_loc);
			val_risk = eucl_distance(_points[t][i], _depot1_loc);
			_riskgraph[0][idxmat_1 + i] = make_pair(true, val_risk);
			_riskgraph[idxmat_1 + i][0] = make_pair(true, val_risk);

		}
	}
	/* each pair of entry and exit between targets*/
	for (int s = 1; s <= _nb_targets; s++) { // ii) boundary s <-> boudary t
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // exit of target s
		idxmat_3 = (s - 1) * 2 * _nb_dstzn + 1; // entries  of target s
		for (int t = 1; t <= _nb_targets; t++) {
			if(s != t){
				idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1; // vertices on boundary t acted as entries
				idxmat_4 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // vertices on boundary t acted as exits
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						// val_risk = get_risk_outerTrajc(_points[s][i], _points[t][j]); 
						val_risk = eucl_distance(_points[s][i], _points[t][j]); 

						_riskgraph[idxmat_1 + i][idxmat_2 + j] = make_pair(true, val_risk);
						_riskgraph[idxmat_2 + j][idxmat_1 + i] = make_pair(true, val_risk);

						_riskgraph[idxmat_4 + j][idxmat_3 + i] = make_pair(true, val_risk);
						_riskgraph[idxmat_3 + i][idxmat_4 + j] = make_pair(true, val_risk);

						// cout<< s << " " << t << ": " << idxmat_1+i << "-->" << idxmat_2+j << " dist: "<< _G[idxmat_1 + i][idxmat_2 + j].second << endl;
						// cout << idxmat_1+i << "-->" << idxmat_2+j << " : " << _G[idxmat_1 + i][idxmat_2 + j].second << '\n';
					}
				}
			}
		}
	}
	/*all exits to arrival depot*/
	for (int s = 1; s <= _nb_targets; s++) { 	
		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			// val_risk = get_risk_outerTrajc(_points[s][i], _depot2_loc);
			val_risk = eucl_distance(_points[s][i], _depot2_loc);

			_riskgraph[idxmat_2 + i][_graphsize - 1] = make_pair(true, val_risk);
			_riskgraph[_graphsize - 1][idxmat_2 + i] = make_pair(true, val_risk);	
		}
	}


	int idx_row, idx_col;
	for (int s = 0; s < _nb_targets; s++) { 
		/* write admissible trajectories (meet minimum reward requirement) to matrix G */ 
		idx_row = s * 2 * _nb_dstzn + 1; 
		idx_col = s * 2 * _nb_dstzn + _nb_dstzn + 1;
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if(i != j){
					// int flag = (_nb_dstzn - i + j) % _nb_dstzn;
					_riskgraph[idx_row + i][idx_col + j] = make_pair(true, _all_innermatr[s][i][j].first);
					_riskgraph[idx_col + j][idx_row + i] = make_pair(true, _all_innermatr[s][i][j].first);
					
					// _G[idx_row + i][idx_col + j] = make_pair(false, INF);
				}
			}
		}
	}

	/*print for debugging: */
	if(false){
		for (int i = 0; i < _graphsize; i++) {
			for (int j = 0; j < _graphsize; j++) {
				cout << _riskgraph[i][j].second << ' ';
			}
			cout << '\n';
		}
	}

}

double DataHandler::eucl_distance(Vertex& v, Vertex& u) {
	return sqrt(pow(u._x - v._x, 2) + pow(u._y - v._y, 2));
}




double DataHandler::get_risk_outerTrajc(Vertex v, Vertex u) {
	// double total_risk = get_lineSeg_len(v, u);
	double total_risk = _par_outer * eucl_distance(v, u); // total risk over the line segment
	// cout << total_risk << "_*_";
	for (int i = 0; i < _nb_targets; i++) { // additional possible risk caused when passing through some region(s)
		tuple<bool, double, double> result = is_intersected(v, u, i);
		// cout << get<0>(result) << "**" << get<1>(result) << "**"<<get<2>(result) << endl;
		if (get<0>(result) == true){ // if current straight line intersects with region Ui, we add additional risk to this outer path
			// cout << "intersected target " << i+1 << endl;
			// total_risk = total_risk + get<1>(result) - get<2>(result);
			total_risk += get<1>(result); // add penalty
		}
	}
	// cout << total_risk << " ";
	// exit(0);
	return total_risk;
}

tuple<bool, double,double> DataHandler::is_intersected(Vertex v, Vertex u, int tar) {
	// current target center is denoted as o;
	Vertex o = {_target_locs[tar]._x, _target_locs[tar]._y };
	myVector v_2_o(v, o); // vector v->o
	myVector v_2_u(v, u); // vector v->u
	myVector o_2_u(o, u); // vector o->u
	double len_v_2_o = v_2_o.get_vecLen();
	double len_v_2_u = v_2_u.get_vecLen();
	double len_o_2_u = o_2_u.get_vecLen();	
	double angle_u = acos((len_o_2_u*len_o_2_u + len_v_2_u*len_v_2_u-len_v_2_o*len_v_2_o)/(2.0*len_o_2_u*len_v_2_u));
	double l2 = len_o_2_u * sin(angle_u); // vertical distance from center o to line segment vu
	double l1; // half length of the chord
	// double risk_vu_on_Ui; // risk of path vu when passing through target region Ui
	if (l2 + 0.0000001 < _radii[tar] && angle_u < (M_PI / 2.0) && len_o_2_u * cos(angle_u)+ 0.0000001 < len_v_2_u) {
		l1 = sqrt(_radii[tar] * _radii[tar] - l2 * l2);
		// risk_vu_on_Ui = inner_risk_function(l1, l2, tar);
		// cout << "vo:" << len_v_2_o << " vu:" << len_v_2_u << " ou:" << len_o_2_u << endl;
		// cout << "v:" << v._x << "," << v._y << endl;
		// cout << "u:" << u._x << "," << u._y << endl;
		// cout << "o:" << o._x << "," << o._y << endl;
		// cout << len_v_2_u << " " << l1 << " " << l2 <<  " r=" << _dataset->_radii[tar] << endl;
		// return make_tuple(true, risk_vu_on_Ui, outer_risk_function(2 * l1));

		return make_tuple(true, 2.0*l1, _par_outer* 2 * l1);

	}
	else {
		return make_tuple(false, NULL, NULL);
	}
}