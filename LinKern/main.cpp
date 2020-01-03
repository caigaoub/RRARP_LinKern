#include "gurobi_c++.h"
#include "costInfo.h"
#include "MTZFormulation.h"
#include "QAPFormulation.h"
#include "MCFormulation.h"
#include "SPFormulation.h"
#include "STEFormulation.h"
#include "LinKernAlgorithm.h"

int main()
{
	costInfo cost("TSP_instance_n_70_s_1.dat");
	GRBEnv* env = new GRBEnv();
	GRBModel model = GRBModel(*env);
	
	model.getEnv().set(GRB_DoubleParam_TimeLimit,1000);


	int choose_model_type = 4;


	if (choose_model_type == 1)
	{
		STEFormulation STEForm(&cost, &model);
		STEForm.solve(&model);
		STEForm.printSol(&model);
	}
	else if (choose_model_type == 2)
	{
		MCFormulation MCForm(&cost, &model);
		MCForm.solve(&model);
		MCForm.printSol(&model);
	}
	else if (choose_model_type == 3)
	{
		MTZFormulation MTZForm(&cost, &model);
		MTZForm.solve(&model);
		MTZForm.printSol(&model);

	}
	else if (choose_model_type == 4)
	{
		LinKern lk(&cost);

		lk.solve();
	}

	system("pause");
	return 0;

}





