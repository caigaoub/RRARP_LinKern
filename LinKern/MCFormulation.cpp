#include "gurobi_c++.h"
#include <vector>
#include "costInfo.h"
#include "MCFormulation.h"

MCFormulation::MCFormulation(costInfo* c, GRBModel* model){

    cost=c;
	n=cost->getNodes();

	// Build model
	model->set(GRB_IntAttr_ModelSense, 1);
    x = new GRBVar*[n];
    y = new GRBVar**[n];
    GRBLinExpr exp = 0;
    for (int i = 0; i < n; i++) {
        x[i]=model->addVars(n);
        y[i]=new GRBVar*[n];
        for (int j=0; j<n; j++) {
            y[i][j] = model->addVars(n);
        }
    }

    model->update();

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i!=j) {
                x[i][j].set(GRB_CharAttr_VType, GRB_BINARY);
                x[i][j].set(GRB_StringAttr_VarName, "x_"+cost->itos(i)+","+cost->itos(j));
                exp+=cost->getCost(i, j)*x[i][j];
                for (int s=1; s<n; s++)
                    y[i][j][s].set(GRB_StringAttr_VarName, "y_"+cost->itos(i)+","+cost->itos(j)+","+cost->itos(s));
            }
        }
    }

    model->setObjective(exp);

    for (int i = 0; i < n; i++){
        exp=0;
        for (int j = 0; j < n; j++){
            if (i!=j)
                exp+=x[i][j];
        }
        model->addConstr(exp == 1);
    }

    for (int i = 0; i < n; i++){
        exp=0;
        for (int j = 0; j < n; j++){
            if (i!=j)
                exp+=x[j][i];
        }
        model->addConstr(exp == 1);
    }

    for (int s = 1; s < n; s++){
        for (int i = 0; i < n; i++){
            exp=0;
            for (int j = 0; j < n; j++){
                if(i!=j){
                    exp+=y[i][j][s];
                    exp-=y[j][i][s];
                    model->addConstr(y[i][j][s]<=x[i][j]);
                }
            }
            if (i==0)
                model->addConstr(exp == 1);
            else
                if(s==i)
                    model->addConstr(exp == -1);
                else
                    model->addConstr(exp == 0);
        }
    }


    model->update();

    //model->write("model.lp");

}

void MCFormulation::solve(GRBModel *model){
    try{

        model->optimize();

    }catch(GRBException e){
        cout<<"Error number: "<<e.getErrorCode()<<endl;
        cout<<e.getMessage()<<endl;

    }catch (...){
        cout<<"Error during optimization"<<endl;
    }
}

void MCFormulation::printSol(GRBModel *model){
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            if (x[i][j].get(GRB_DoubleAttr_X)>0.5)
                cout<< i<<" "<<j<<endl;

    int i,k,t;
    i=0;
    k=-1;
    t=0;
    while(t<n){
      cout<<i<<endl;
      for (int j=0; j<n; j++){
        if(k!=j && x[i][j].get(GRB_DoubleAttr_X)>0.5){
            k=i;
            i=j;
          }
      }
      t++;
    }

}
