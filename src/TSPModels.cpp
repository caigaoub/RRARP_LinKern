#include "TSPModels.h"

void TSPModel_STE::init_edge_weights(int size_node_arg, const vector<vector<double>> & costs_arg){
	if((int)costs_arg.size() != size_node_arg){
		cerr << "ERROR: size does not fit (TSPModels.cpp::line 4)" << endl;
	}
	this->_size_var_x = size_node_arg;
	this->_cost = new double*[_size_var_x];
	for(int i = 0; i < _size_var_x; i++)
		this->_cost[i] = new double[_size_var_x];

	/* symmetric cost */
	for(unsigned i = 0; i < costs_arg.size(); i++){
		for(unsigned j = i+1; j < costs_arg[i].size(); j++){
			_cost[i][j] = costs_arg[i][j];
			_cost[j][i] = costs_arg[i][j];
		}
	}
	_edge_cost_status = true;
}



void TSPModel_STE::create_formula(){
	if(_edge_cost_status == false){
		cerr << "ERROR: edge costs have been assigned yet!! " << endl;
	}
	// set up the model
	GRBEnv env = GRBEnv();
	_model = new GRBModel(env);
	_model->getEnv().set(GRB_IntParam_OutputFlag, 1);
	_model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	_model->getEnv().set(GRB_IntParam_PreCrush, 1);
	_model->getEnv().set(GRB_IntParam_MIPFocus, 1); // focus on feasibility
  	_model->getEnv().set(GRB_DoubleParam_TimeLimit, 7200);

  	_var_x = new GRBVar*[_size_var_x];
	for (int i = 0; i < _size_var_x; i++)
		_var_x[i] = new GRBVar[_size_var_x];

	for (int i = 0; i < _size_var_x; i++) {
		for (int j = 0; j < _size_var_x; j++) {
			_var_x[i][j] = _model->addVar(0.0, 1.0, _cost[i][j], GRB_BINARY, "x_" + itos(i) + "," + itos(j));
		}
	}
	_model->update();

	// Degree-1 constraints
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
	}
	_model->update();

	// Forbid edge from node back to itself
	for (int i = 0; i < _size_var_x; i++)	{
		_var_x[i][i].set(GRB_DoubleAttr_UB, 0);
	}
	_model->update();

	/********** Additional edge constraints **********/
	_var_x[_size_var_x - 1][0].set(GRB_DoubleAttr_LB, 1);
	_model->update();

}

void TSPModel_STE::solve()
{
	try	{
		SubtourElimCuts cb = SubtourElimCuts(_var_x, _size_var_x);
		_model->setCallback(&cb);
		_model->optimize();
		if (_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL){			
			double **sol = new double*[_size_var_x];
			for (int i = 0; i < _size_var_x; i++) 
				sol[i] = new double[_size_var_x];
			
			for (int i = 0; i < _size_var_x; i++) {
				for (int j = 0; j < _size_var_x; j++) {
					sol[i][j] = _var_x[i][j].get(GRB_DoubleAttr_X);
				}
			}
			double	obj_val = _model->get(GRB_DoubleAttr_ObjVal);
			_curbestobj = obj_val;
			cout << " ====>> objective of TSP solution: " << obj_val << endl;
			// cout << " ====>> optimal TSP Solution matrix: " << '\n';
			// for (int i = 0; i < _size_var_x; i++) {
			// 	for (int j = 0; j < _size_var_x; j++) {
			// 		if(abs(sol[i][j]) < 1e-6){
			// 				cout << 0 << "  ";
			// 		}
			// 		if(abs(sol[i][j]-1) < 1e-6){
			// 			cout << 1 << "  ";
			// 		}
			// 		// cout << sol[i][j] << "   ";	
			// 	}
			// 	cout << endl;
			// }
			int len;
			int *tour = new int[_size_var_x];
			SubtourElimCuts::find_subtour(_size_var_x, sol, &len, tour);
			// cout << "====>> current best visiting sequence: ";
			_opt_seq.clear();
			_opt_seq.resize(_size_var_x);
			for (int i = 0; i < len; i++) {
				_opt_seq.at(i) = tour[i];
				cout << tour[i] << ' ';
			}
			cout << '\n';

			for (int i = 0; i < _size_var_x; i++)
				delete[] sol[i];
			delete[] sol;
		}

	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
}


TSPModel_STE::~TSPModel_STE(){
	for (int i = 0; i < _size_var_x; i++)
		delete[] _var_x[i];
	delete[] _var_x;
	for (int i = 0; i < _size_var_x; i++)
		delete[] _cost[i];
	delete[] _cost;

	delete _model;
	
}
