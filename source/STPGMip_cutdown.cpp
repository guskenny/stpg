#include "STPGMip.h"

// callback class
class STPGCallback : public qol::Callback{
  // vectors for variables
  const std::vector<qol::Variable> &x;
  const std::vector<qol::Variable> &y;
  STPGModel *probModel;
  public:
  STPGCallback(STPGModel *_probModel,const std::vector<qol::Variable> &_x, const std::vector<qol::Variable> &_y):probModel(_probModel),x(_x),y(_y){};
  ~STPGCallback(){};
  Status callback(Progress where);
};


// solve method - creates the MIP model and solves it
void STPGMip::solve(){
  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);

	  mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    // initialise the MIP model
    initMIPModel(mip);

    // create callback object <--- DON'T THINK THIS IS RIGHT
    STPGCallback cb(probModel,x,y);

    // set callback function <--- DON'T THINK THIS IS RIGHT
    mip.setCallback(&cb);

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());
   
   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
}


// callback method - adds a constraint for each component with a cycle
qol::Callback::Status STPGCallback::callback(Progress where){
  
  std::vector<std::vector<int> > components;

  // get connected components with cycles, each component is a set of edges
  get_components(probModel,x,y,components);

  // for each component, add a cut to the model using addCut()
  for (int c = 0; c < components.size(); ++c){
    qol::Expression lhs;
    qol::Expression rhs;

    // vector to make sure we dont double up on node variables
    std::vector<int> added(probModel->n_nodes(),0);

    for (std::vector<int>::iterator e = components[c].begin() ; e != components[c].end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      lhs += x[*e];

      std::vector<int> endpoints{edge->getSrcID(),edge->getTgtID()}
      for (std::vector<int>::iterator i = endpoints.begin() ; i != endpoints.end(); ++i){
        if (!added[*i]){
          rhs += y[*i];
          added[*i] = 1;
        }
      }
    }
    rhs -= 1;
    addCut(lhs <= rhs);
    PE("cut added without exception")
  }
}

// initialises the MIP model variables and equality constraint
void STPGMip::initMIPModel(qol::MIPSolver &mip){
    x.clear();
    y.clear();

    // edge variables
    for (int e = 0; e < probModel->n_edges(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(e);
      x.push_back(mip.addVar(0,1,edge->getWt(),qol::Variable::CONTINUOUS));
      std::string var_name("x_" + e);
      mip.setVarName(x[e],var_name);
    }

    // node variables
    for (int n = 0; n < probModel->n_nodes(); ++n){
      int lb = 1 ? probModel->prob_graph.getNode(n)->isTerm() : 0;
      y.push_back(mip.addVar(lb,1,0,qol::Variable::BINARY));
      std::string var_name("y_" + n);
      mip.setVarName(y[n],var_name);
    }

    std::vector<double> x_start_vec(x.size(),0.0);
    std::vector<double> y_start_vec(y.size(),0.0);

    std::vector<int> sol_edges = best_sol.getSet();

    // set start values from best_sol
    for (std::vector<int>::iterator e = sol_edges.begin() ; e != sol_edges.end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      x_start_vec[*e] = 1;

      y_start_vec[edge->getSrcID()] = 1;
      y_start_vec[edge->getTgtID()] = 1;
    }

    // equality constraint
    qol::Expression lhs;
    qol::Expression rhs;
    for (int e = 0; e < probModel->n_edges(); ++e){
      lhs += x[e];
    }
    for (int n = 0; n < probModel->n_nodes(); ++n){
      rhs += y[n];
    }

    rhs -= 1;

    mip.addConstraint(lhs <= rhs).setName("init_lt");
    mip.addConstraint(lhs >= rhs).setName("init_gt");
}
