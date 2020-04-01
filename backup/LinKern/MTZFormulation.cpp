
#include "gurobi_c++.h"
#include <vector>
#include "costInfo.h"
#include "MTZFormulation.h"

MTZFormulation::MTZFormulation(costInfo* c, GRBModel* model){

    cost=c;
	n=cost->getNodes();
    
	// Build model
	model->set(GRB_IntAttr_ModelSense, 1);
    x = new GRBVar*[n];
    u = model->addVars(n);
    GRBLinExpr exp = 0;
    for (int i = 0; i < n; i++) {
        x[i]=model->addVars(n);
    }
    
    model->update();
    
    for (int i=0; i<n;i++){
        u[i].set(GRB_StringAttr_VarName, "u_"+cost->itos(i));
        for (int j=0; j<n; j++){
            if(i!=j){
                x[i][j].set(GRB_CharAttr_VType, GRB_BINARY);
                x[i][j].set(GRB_StringAttr_VarName, "x_"+cost->itos(i)+","+cost->itos(j));
                exp += cost->getCost(i, j)*x[i][j];
            }
        }
    }
    
    model->setObjective(exp);
    
    for (int i = 0; i < n; i++){
        exp=0;
        for (int j = 0; j < n; j++){
            if (i!=j){
                exp+=x[i][j];
                if (i>0&&j>0)
                    model->addConstr(u[i]-u[j]+n*x[i][j]<=n-1);
            }
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

    model->update();
    
    //model->write("model.lp");
    
}

void MTZFormulation::solve(GRBModel *model){
    try{
        
        model->optimize();
        
    }catch(GRBException e){
        cout<<"Error number: "<<e.getErrorCode()<<endl;
        cout<<e.getMessage()<<endl;
        
    }catch (...){
        cout<<"Error during optimization"<<endl;
    }
}

void MTZFormulation::printSol(GRBModel *model){
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            if (i!=j&&x[i][j].get(GRB_DoubleAttr_X)>0.5)
                cout<< i<<" "<<j<<endl;
        
}


