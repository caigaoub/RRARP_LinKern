#include "gurobi_c++.h"
#include <vector>
#include "costInfo.h"
#include "SPFormulation.h"

SPFormulation::SPFormulation(costInfo* c, GRBModel* model){

    cost=c;
	n=cost->getNodes();
    
	// Build model
	model->set(GRB_IntAttr_ModelSense, 1);
    x = new GRBVar**[n];
    GRBLinExpr exp = 0;
    for (int s = 0; s < n; s++) {
        x[s]=new GRBVar*[n];
        for (int i=0; i<n; i++) {
            x[s][i] = model->addVars(n);
        }
    }
    
    model->update();
    
    for (int s=0; s<n; s++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i!=j) {
                    x[s][i][j].set(GRB_CharAttr_VType, GRB_BINARY);
                    x[s][i][j].set(GRB_StringAttr_VarName, "x_"+cost->itos(s)+","+cost->itos(i)+","+cost->itos(j));
                    exp+=cost->getCost(i, j)*x[s][i][j];

            }            }
        }
    }
    
    model->setObjective(exp);
    
    
    exp=0;
    for (int j = 1; j < n; j++){
         exp+=x[0][0][j];
    }
    model->addConstr(exp == 1);
    
    for (int i = 1; i < n; i++){
        exp=0;
        for (int j = 1; j < n; j++){
            if(i!=j){
                exp+=x[1][i][j];
            }
        }
        exp-=x[0][0][i];
        model->addConstr(exp == 0);
    }
    
    for (int s = 1; s < n-2; s++){
        for (int i = 1; i < n; i++){
            exp=0;
            for (int j = 1; j < n; j++){
                if(i!=j){
                    exp+=x[s+1][i][j];
                    exp-=x[s][j][i];
                }
            }
            model->addConstr(exp == 0);
        }
    }
    
    for (int i = 1; i < n; i++){
        exp=0;
        for (int j = 1; j < n; j++){
            if(i!=j){
                exp-=x[n-2][j][i];
            }
        }
        exp+=x[n-1][i][0];
        model->addConstr(exp == 0);
    }
    
    exp=0;
    for (int j = 1; j < n; j++){
        exp+=x[n-1][j][0];
    }
    model->addConstr(exp == 1);
    
    for (int i = 1; i < n; i++){
        exp=0;
        exp+=x[n-1][i][0];
        for (int s = 1; s < n-1; s++){
            for (int j = 1; j < n; j++){
                if(i!=j){
                    exp+=x[s][i][j];
                }
            }
        }
        model->addConstr(exp <= 1);

    }
    
    model->update();
    
    //model->write("model.lp");
    
}

void SPFormulation::solve(GRBModel *model){
    try{
        
        model->optimize();
        
    }catch(GRBException e){
        cout<<"Error number: "<<e.getErrorCode()<<endl;
        cout<<e.getMessage()<<endl;
        
    }catch (...){
        cout<<"Error during optimization"<<endl;
    }
}

void SPFormulation::printSol(GRBModel *model){
    for (int s=0; s<n; s++)
        for (int i=0; i<n; i++)
            for (int j=0; j<n; j++)
                if (x[s][i][j].get(GRB_DoubleAttr_X)>0.5)
                    cout<<s<<" "<<i<<" "<<j<<endl;
        
}


