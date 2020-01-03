#include "gurobi_c++.h"
#include <vector>
#include "costInfo.h"
#include "QAPFormulation.h"

QAPFormulation::QAPFormulation(costInfo* c, GRBModel* model){

    cost=c;
	n=cost->getNodes();
    
	// Build model
	model->set(GRB_IntAttr_ModelSense, 1);
    x = new GRBVar*[n];
    u = new GRBVar**[n];
    GRBLinExpr exp = 0;
    for (int s = 0; s < n; s++) {
        x[s]=model->addVars(n);
        u[s]=new GRBVar*[n];
        for (int i=0; i<n; i++) {
            u[s][i] = model->addVars(n);
        }
    }
    
    model->update();
    
    for (int s=0; s<n; s++){
        for (int i=0; i<n; i++){
            x[s][i].set(GRB_CharAttr_VType, GRB_BINARY);
            x[s][i].set(GRB_StringAttr_VarName, "x_"+cost->itos(s)+","+cost->itos(i));
            for (int j=0; j<n; j++){
                u[s][i][j].set(GRB_StringAttr_VarName, "u_"+cost->itos(s)+","+cost->itos(i)+","+cost->itos(j));
                if (s<n-1){
                    exp+=cost->getCost(i, j)*u[s][i][j];
                    model->addConstr(u[s][i][j] >= x[s][i]+x[s+1][j]-1);
                }
            }
        }
    }

    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++){
            exp+=cost->getCost(i, j)*u[n-1][i][j];
            model->addConstr(u[n-1][i][j] >= x[n-1][i]+x[0][j]-1);
        }

    
    model->setObjective(exp);
    
    for (int s = 0; s < n; s++){
        exp=0;
        for (int i = 0; i < n; i++){
            exp+=x[s][i];
        }
        model->addConstr(exp == 1);
    }
    
    for (int i = 0; i < n; i++){
        exp=0;
        for (int s = 0; s < n; s++){
            exp+=x[s][i];
        }
        model->addConstr(exp == 1);
    }

    model->update();
    
    //model->write("model.lp");
    
}

void QAPFormulation::solve(GRBModel *model){
    try{
        
        model->optimize();
        
    }catch(GRBException e){
        cout<<"Error number: "<<e.getErrorCode()<<endl;
        cout<<e.getMessage()<<endl;
        
    }catch (...){
        cout<<"Error during optimization"<<endl;
    }
}

void QAPFormulation::printSol(GRBModel *model){
    for (int s=0; s<n; s++)
        for (int i=0; i<n; i++)
            if (x[s][i].get(GRB_DoubleAttr_X)>0.5)
                cout<< i<<" "<<s<<endl;
        
}


