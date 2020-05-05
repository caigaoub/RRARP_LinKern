#include "STEFormulation.h"


void STEFormulation::build_formul(GRBModel* model_MP_, DataHandler& dataset) {
	// _time = new ProgTime();
	// _time->start_prog();

	this->_dataset = &dataset;
	this->_G_dic = &(dataset._riskgraph_dic);
	this->_all_innerG = &(dataset._all_innermatr);	
	this->_nb_dstzn = dataset._nb_dstzn;
	this->_graphsize = dataset._graphsize;
	this->_nb_targets = _dataset->_nb_targets;

	this->_model = model_MP_;
	this->_size_var_x = _nb_targets + 2;
	this->_size_var_u = _graphsize;

	/* set up model parameters */
	_model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	_model->getEnv().set(GRB_IntParam_PreCrush, 1);
	_model->getEnv().set(GRB_IntParam_NumericFocus, 1);
	_model->getEnv().set(GRB_DoubleParam_TimeLimit, 7200);
	// _model->getEnv().set(GRB_DoubleParam_TimeLimit, GRB_INFINITY);

	// Step 0: create variables for ..
	_var_u = new GRBVar*[_size_var_u];
	for(int i = 0; i < _size_var_u; i++) {
		_var_u[i] = new GRBVar[_size_var_u];
	}
	for(int i = 0; i <_size_var_u; i++) {
		for(int j = 0; j < _size_var_u; j++) {
			// if((*_G_dic)[i][j].first == true)
				_var_u[i][j] = _model->addVar(0.0, 1.0, (*_G_dic)[i][j].second, GRB_CONTINUOUS, "u_" + itos(i) + "," + itos(j));
		}
	}
	_model->update();

	// Step 1: create variables for hamiltonian cycle problem 
	_var_x = new GRBVar*[_size_var_x];
	for(int i = 0; i < _size_var_x; i++) {
		_var_x[i] = new GRBVar[_size_var_x];
	}
	for(int i = 0; i <_size_var_x; i++) {
		for(int j = 0; j < _size_var_x; j++) {
			_var_x[i][j] = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(i) + "," + itos(j));
		}
	}
	_model->update();

	// Step 2: add constraints: degree constraints 
	GRBLinExpr expr1, expr2;
	for (int i = 0; i < _size_var_x; i++) {
		expr1 = 0;
		expr2 = 0;
		for (int j = 0; j < _size_var_x; j++) {
			expr1 += _var_x[i][j];
			expr2 += _var_x[j][i];
		}
		_model->addConstr(expr1 == 1, "deg1_row" + itos(i));
		_model->addConstr(expr2 == 1, "deg1_col" + itos(i));
		_model->update();
	}
	_model->update();

	_var_x[_size_var_x - 1][0].set(GRB_DoubleAttr_LB, 1);
	for (int i = 0; i < _size_var_x; i++) {
		_var_x[i][i].set(GRB_DoubleAttr_UB, 0);
	}

	// step 3: linking constraints
	int idxmat_1, idxmat_2, idxmat_3, idxmat_4;/* i) departure depot --> each entry turning point of a target */	
	for (int t = 1; t <= _nb_targets; t++) {	
		idxmat_1 = (t - 1) * 2 * _nb_dstzn + 1; 
		expr1 = 0;
		for (int i = 0; i < _nb_dstzn; i++) {
			 // _model->addConstr(_var_u[0][idxmat_1+i] <= _var_x[0][t]);
			 // _model->addConstr(_var_u[idxmat_1+i][0] <= _var_x[t][0]); 
			if((*_G_dic)[0][idxmat_1+i].first) 
				expr1 += _var_u[0][idxmat_1+i];
		}
		_model->addConstr(expr1 <= _var_x[0][t]);
	}
	_model->update();

	/* each pair of entry and exit between targets*/
	for (int s = 1; s <= _nb_targets; s++) { // ii) boundary s <-> boudary t
		idxmat_1 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // exit of target s
		idxmat_3 = (s - 1) * 2 * _nb_dstzn + 1; // entries  of target s
		for (int t = 1; t <= _nb_targets; t++) {
			if(s != t){
				idxmat_2 = (t - 1) * 2 * _nb_dstzn + 1; // vertices on boundary t acted as entries
				idxmat_4 = (t - 1) * 2 * _nb_dstzn + _nb_dstzn + 1; // vertices on boundary t acted as exits
				expr1 = 0;
				expr2 = 0;
				for (int i = 0; i < _nb_dstzn; i++) {
					for (int j = 0; j < _nb_dstzn; j++) {
						// _riskgraph[idxmat_1 + i][idxmat_2 + j] = make_pair(true, val_risk);
						// _riskgraph[idxmat_2 + j][idxmat_1 + i] = make_pair(true, val_risk);
						// _model->addConstr(_var_u[idxmat_1 + i][idxmat_2 + j] <= _var_x[s][t]);
						// _model->addConstr(_var_u[idxmat_2 + j][idxmat_1 + i] <= _var_x[t][s]);
						if((*_G_dic)[idxmat_1 + i][idxmat_2 + j].first) 
							expr1 += _var_u[idxmat_1 + i][idxmat_2 + j];
						// if((*G)[idxmat_2 + j][idxmat_1 + i].first) 
						// 	expr2 += _var_u[idxmat_2 + j][idxmat_1 + i];

						// _riskgraph[idxmat_4 + j][idxmat_3 + i] = make_pair(true, val_risk);
						// _riskgraph[idxmat_3 + i][idxmat_4 + j] = make_pair(true, val_risk);
						// _model->addConstr(_var_u[idxmat_4 + j][idxmat_3 + i] <= _var_x[t][s]);
						// _model->addConstr(_var_u[idxmat_3 + i][idxmat_4 + j] <= _var_x[s][t]);
						if((*_G_dic)[idxmat_4 + j][idxmat_3 + i].first) 
							expr2 += _var_u[idxmat_4 + j][idxmat_3 + i];
						// if((*G)[idxmat_3 + i][idxmat_4 + j].first) 
						// 	expr4 += _var_u[idxmat_3 + i][idxmat_4 + j];
					}
				}
				_model->addConstr(expr1 <= _var_x[s][t]);
				_model->addConstr(expr2 <= _var_x[t][s]);

			}
		}
	}
	_model->update();

	/*all exits to arrival depot*/
	for (int s = 1; s <= _nb_targets; s++) { 	
		idxmat_2 = (s - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		expr1 = 0;
		for (int i = 0; i < _nb_dstzn; i++) {
			// _riskgraph[idxmat_2 + i][_graphsize - 1] = make_pair(true, val_risk);
			// _riskgraph[_graphsize - 1][idxmat_2 + i] = make_pair(true, val_risk);	
			if((*_G_dic)[idxmat_2 + i][_graphsize - 1].first) 
				expr1 += _var_u[idxmat_2 + i][_graphsize - 1];
		}
		_model->addConstr(expr1 <= _var_x[s][_nb_targets+1]);
	}
	_model->update();


	// int idx_row, idx_col;
	// for (int s = 0; s < _nb_targets; s++) { 
	// 	/* write admissible trajectories (meet minimum reward requirement) to matrix G */ 
	// 	idx_row = s * 2 * _nb_dstzn + 1; 
	// 	idx_col = s * 2 * _nb_dstzn + _nb_dstzn + 1;
	// 	for (int i = 0; i < _nb_dstzn; i++) {
	// 		for (int j = 0; j < _nb_dstzn; j++) {
	// 			if(i != j){
	// 				// int flag = (_nb_dstzn - i + j) % _nb_dstzn;
	// 				_riskgraph[idx_row + i][idx_col + j] = make_pair(true, _all_innermatr[s][i][j].first);
	// 				_riskgraph[idx_col + j][idx_row + i] = make_pair(true, _all_innermatr[s][i][j].first);
					
	// 				// _G[idx_row + i][idx_col + j] = make_pair(false, INF);
	// 			}
	// 		}
	// 	}
	// }

	for(int i=1; i <_graphsize-1; i++){
		expr1 = 0;
		for(int j=0; j<_graphsize; j++){
			if((*_G_dic)[i][j].first)
				expr1 += _var_u[i][j];
			
			if((*_G_dic)[j][i].first)
				expr1 -= _var_u[j][i];
		}
		_model->addConstr(expr1 == 0);
		_model->update();
	}
	expr1 = 0;
	for(int j=0; j<_graphsize; j++){
		if((*_G_dic)[j][_graphsize-1].first)
			expr1 += _var_u[j][_graphsize-1];
	}
	_model->addConstr(expr1 == 1);
	_model->update();


	_model->write("model.lp");
}
/**
	solve master problem(TSP)
*/
void STEFormulation::solve_formul_wCB() {
	try {
		SubtourCuts * cb = new SubtourCuts(_var_x, _nb_targets);
		_model->setCallback(cb);
		_model->optimize();
		// _time->end_prog();
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL || _model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {	
			double obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			_total_nb_subtour_cuts = cb->get_nb_subtour_cuts();
			// double **sol = new double*[_size_var_x];
			// for (int i = 0; i < _size_var_x; i++) 
			// 	sol[i] = new double[_size_var_x];
			
			// for (int i = 0; i < _size_var_x; i++) {
			// 	for (int j = 0; j < _size_var_x; j++) {
			// 		sol[i][j] = _var_x[i][j].get(GRB_DoubleAttr_X);
			// 	}
			// }
			
			double **sol = new double*[_size_var_u];
			for (int i = 0; i < _size_var_u; i++) 
				sol[i] = new double[_size_var_u];
			
			for (int i = 0; i < _size_var_u; i++) {
				for (int j = 0; j < _size_var_u; j++) {
					sol[i][j] = _var_u[i][j].get(GRB_DoubleAttr_X);
				}
			}
			// double obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			cout << " ====>> objective of TSP solution: " << obj_val << endl;
			// cout << " ====>> optimal TSP Solution matrix: " << '\n';
			// for (int i = 0; i < _size_var_u; i++) {
			// 	for (int j = 0; j < _size_var_u; j++) {
			// 		// if(abs(sol[i][j]) < 1e-6){
			// 		// 		cout << 0 << "  ";
			// 		// }
			// 		if(abs(sol[i][j]-1) < 1e-6){
			// 			cout << i << " " << j <<  "1" << endl;
			// 		}
			// 		// cout << sol[i][j] << "   ";	
			// 	}
			// 	// cout << endl;
			// }
			// int len;
			// int *tour = new int[_size_var_x];
			// SubtourCuts::findsubtour(_size_var_x, sol, &len, tour);
			// _opt_seq.clear();
			// _opt_seq.resize(_size_var_x);
			// for (int i = 0; i < len; i++) {
			// 	_opt_seq.at(i) = tour[i];
			// 	// cout << tour[i] << ' ';
			// }
			// // cout << '\n';

			delete cb;
			// return make_pair(obj_val, v_val);
		}
		else {
			delete cb;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	// return make_pair(-INF, -INF);
}


void STEFormulation::get_optimal_sol(double ** sol) {
	// cout << " flag ***" << endl;
	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			sol[i][j] = _var_x[i][j].get(GRB_DoubleAttr_X);
		}
	}
	// cout << "end flag======" << endl;
}



STEFormulation::~STEFormulation(){
	for(int i = 0; i < _size_var_x; i++)
		delete[] _var_x[i];
	delete[] _var_x;

	for(int i = 0; i < _size_var_u; i++)
		delete[] _var_u[i];
	delete[] _var_u;
}
