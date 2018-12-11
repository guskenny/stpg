#include "STPGMip.h"

STPGMip::STPGMip(const SettingsHandler sh, STPGModel *probModel, set_obj &best_sol) : probModel(probModel), best_sol(best_sol), sh(sh) {

  PE("\nMerge Solver initialised")
}

void STPGMip::solve(set_obj &sol){
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

    initMIPModel(mip);

    int count = 0;
    int *countPtr = &count;

    STPGCallback cb = STPGCallback(probModel,x,y,countPtr);

    mip.setCallback((qol::Callback *)&cb);
    // mip.setCallback(&cb);

    CPXENVptr env=(dynamic_cast<qol::CplexFormulation *>(mipPtr))->env;
    CPXsetintparam(env,CPX_PARAM_PREIND,0);
    //masterCplex.setParam(IloCplex::PreInd, IloFalse); 
    // Set the maximum number of threads to 1. 
    // This instruction is redundant: If MIP control callbacks are registered, 
    // then by default CPLEX uses 1 (one) thread only.
    // Note that the current example may not work properly if more than 1 threads 
    // are used, because the callback functions modify shared global data.
    // We refer the user to the documentation to see how to deal with multi-thread 
    // runs in presence of MIP control callbacks. 

    //masterCplex.setParam(IloCplex::Threads, 1); 
    CPXsetintparam(env,CPX_PARAM_THREADS,1);
  
    CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
    // Turn on traditional search for use with control callbacks
    //masterCplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
    CPXsetintparam(env,CPX_PARAM_MIPSEARCH,CPX_MIPSEARCH_TRADITIONAL);
  

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

qol::Callback::Status STPGCallback::callback(Progress where){
  // PF("\rCallback: " << (*countPtr)++)

  set_obj primal_x(probModel->n_edges());
  set_obj primal_y(probModel->n_nodes());

  for (int e = 0; e < primal_x.idx_size(); ++e){
    if (primalSolution[x[e]] > 0){
      primal_x.addElement(e);
    }
  }

  for (int n = 0; n < primal_y.idx_size(); ++n){
    if (primalSolution[y[n]] > 0){
      primal_y.addElement(n);
    }
  }

  std::vector<set_obj> components;

  // get connected components with cycles, each component is a set of edges
  get_components(probModel->prob_graph,primal_x,primal_y,components);

  // for each component, add a cut to the model using addCut()
  for (int c = 0; c < components.size(); ++c){
    std::vector<int> component = components[c].getSet();
    qol::Expression lhs;
    qol::Expression rhs;

    // vector to make sure we dont double up on node variables
    std::vector<int> added(probModel->n_nodes(),0);

    for (std::vector<int>::iterator e = component.begin() ; e != component.end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      lhs += x[*e];

      std::vector<int> endpoints{edge->getSrcID(),edge->getTgtID()};

      for (std::vector<int>::iterator i = endpoints.begin() ; i != endpoints.end(); ++i){
        if (!added[*i]){
          rhs += y[*i];
          added[*i] = 1;
        }
      }
    }
    rhs -= 1;
    addCut(lhs <= rhs);
  }
    return OK;
}

void STPGMip::initMIPModel(qol::MIPSolver &mip){
    PE("*** Initialising MIP model ***")

    PF("setting up edge variables.. ")

    x.clear();
    y.clear();

    for (int e = 0; e < probModel->n_edges(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(e);
      x.push_back(mip.addVar(0,1,edge->getWt(),qol::Variable::CONTINUOUS));
      std::string var_name("x_" + e);
      mip.setVarName(x[e],var_name);
    }

    PE("done!")
    PF("setting up vertex variables.. ")

    for (int n = 0; n < probModel->n_nodes(); ++n){
      int lb = 1 ? probModel->prob_graph.getNode(n)->isTerm() : 0;
      y.push_back(mip.addVar(lb,1,0,qol::Variable::BINARY));
      std::string var_name("y_" + n);
      mip.setVarName(y[n],var_name);
    }

    PE("done!")
    PF("setting start values.. ")

    std::vector<double> x_start_vec(x.size(),0.0);
    std::vector<double> y_start_vec(y.size(),0.0);

    std::vector<int> sol_edges = best_sol.getSet();

    for (std::vector<int>::iterator e = sol_edges.begin() ; e != sol_edges.end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      x_start_vec[*e] = 1;

      y_start_vec[edge->getSrcID()] = 1;
      y_start_vec[edge->getTgtID()] = 1;
    }

    // for (int e = 0; e < sol_edges.size(); ++e){
    //   PE(e)
    //   // Edge* edge = probModel->prob_graph.getEdge(sol_edges[e]);

    //   // x_start_vec[sol_edges[e]] = 1;

    //   // y_start_vec[edge->getSrcID()] = 1;
    //   // y_start_vec[edge->getTgtID()] = 1;
    // }

    PE("done!")
    PF("setting up initial constraints.. ")

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

    PE("done!")
}
