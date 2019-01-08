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
    param.setParamVal(qol::VERBOSITY,2);
    param.setParamVal(qol::RELGAP,0.001);

	  mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    initMIPModel(mip);

    int count = 0;
    int *countPtr = &count;

    STPGCallback cb = STPGCallback(probModel,x,y,countPtr,mipPtr);

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
  
    CHECK(mip.nConstr())  

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());
   
    sol.clear();

    CHECK(mip.getObjective())
    CHECK(mip.getObjectiveBound())
    CHECK(mip.nConstr())

    for (int e = 0; e < probModel->n_edges(); ++e){
      if (mip.getPrimal(x[e]) > 0){
        sol.addElement(e);
      }
    }    

    PE("*****************************************")
    PE("MIP solution " << ++(*countPtr) << " (final)")
    PE("*****************************************")

    int cost = 0;

    VF("x = { ")

    std::ostringstream latex_x;
    std::ostringstream latex_y;

    set_obj test_y(probModel->n_nodes());

    for (int e = 0; e < sol.size(); ++e){
      Edge * edge = probModel->prob_graph.getEdge(sol.get(e));
      
      cost += edge->getWt();

      VF("(" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") ")

      latex_x << std::to_string(edge->getSrcID()+1) << "/" << std::to_string(edge->getTgtID()+1)<< ",";

      test_y.addElement(edge->getTgtID());
      test_y.addElement(edge->getSrcID());
    }

    VE("}")
    VF("y = { ")

    for (int n = 0; n < test_y.idx_size(); ++n){
      if (test_y.is_element(n)){
        VF(n+1 << " ")
        latex_y << n+1 << ",";
      }
    }

    VE("}")
    PE("cost: " << cost << std::endl)

    VE("\nx = " <<latex_x.str())
    VE("y = " << latex_y.str()<<std::endl)

   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
}

qol::Callback::Status STPGCallback::callback(Progress where){
  // PE("Callback")

  set_obj primal_x(probModel->n_edges());
  set_obj primal_y(probModel->n_nodes());

  set_obj test_y(probModel->n_nodes());

  // PE("\n*****************************************")
  // PE("MIP solution " << ++(*countPtr))
  // PE("*****************************************")

  int cost = 0;

  // PF("x = { ")

  std::ostringstream latex_x;
  std::ostringstream latex_y;

  for (int e = 0; e < primal_x.idx_size(); ++e){
    if (primalSolution[x[e]] > 0){
      primal_x.addElement(e);

      Edge * edge = probModel->prob_graph.getEdge(e);

      cost += edge->getWt();

      // PF("(" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") ")

      latex_x << std::to_string(edge->getSrcID()+1) << "/" << std::to_string(edge->getTgtID()+1)<< ",";

      test_y.addElement(edge->getTgtID());
      test_y.addElement(edge->getSrcID());
    }
  }

  // PE("}")
  // PF("y = { ")

  for (int n = 0; n < primal_y.idx_size(); ++n){
    if (primalSolution[y[n]] > 0){
      primal_y.addElement(n);
      // PF(n+1 << " ")
      // latex_y << n+1 << ",";
    }
  }

  // PE("}")

  // PE("cost: " << cost << std::endl)
  // CHECK(test_y.size())
  // CHECK(primal_y.size())
  // CHECK(test_y.size()-primal_y.size())
  // test_y.subtract(primal_y);
  // CHECK(test_y.size())
  // PE("")

  std::vector<set_obj> components;

  // get connected components with cycles, each component is a set of edges
  get_components(probModel->prob_graph,primal_x,primal_y,components);

  // for each component, add a cut to the model using addCut()
  for (int c = 0; c < components.size(); ++c){
    std::vector<int> component = components[c].getSet();
    qol::Expression lhs;
    qol::Expression rhs;

    std::ostringstream lhs_os;
    std::ostringstream rhs_os;

    // vector to make sure we dont double up on node variables
    std::vector<int> added(probModel->n_nodes(),0);

    for (std::vector<int>::iterator e = component.begin() ; e != component.end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      lhs += x[*e];

      lhs_os << "x[" + std::to_string(edge->getSrcID()+1) + "," + std::to_string(edge->getTgtID()+1) + "]";

      if(e != component.end()-1)
        lhs_os << " + ";

      std::vector<int> endpoints{edge->getSrcID(),edge->getTgtID()};

      for (std::vector<int>::iterator i = endpoints.begin() ; i != endpoints.end(); ++i){
        if (!added[*i]){
          rhs += y[*i];
          rhs_os << "y[" + std::to_string((*i)+1) + "] + ";
          added[*i] = 1;
        }
      }
    }
    rhs -= 1;
    rhs_os << "(-1)";
    // PE("adding cut: " << lhs_os.str() << " <= " << rhs_os.str()<< std::endl)

    addCut(lhs <= rhs);

  }

  // CHECK(getBound())

  // PE("x = " <<latex_x.str())
  // PE("y = " << latex_y.str()<<std::endl)
    return OK;
}

void STPGMip::initMIPModel(qol::MIPSolver &mip){
    PE("*** Initialising MIP model ***")

    PF("setting up edge variables.. ")

    x.clear();
    y.clear();

    VE("Edges: ")

    for (int e = 0; e < probModel->n_edges(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(e);

      VE(edge->getID() << " (" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") = " << edge->getWt())

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

    std::vector<int> sol_edges = best_sol.getSet();

    std::ostringstream x_start;
    std::ostringstream y_start;

    x_start << "x = { ";

    y_start << "y = { ";

    for (std::vector<int>::iterator e = sol_edges.begin() ; e != sol_edges.end() && sh.WARM_START; ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);


      y_start << edge->getSrcID()+1 << ", " << edge->getTgtID()+1 << ", ";

      x_start << "(" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") ";     
    }

    x_start << "}";
    y_start << "}";

    VE("\n" << x_start.str())
    VE("\n" << y_start.str())

    PE("done!")
    PF("setting up initial constraints.. ")

    qol::Expression lhs;
    qol::Expression rhs;
    for (int e = 0; e < probModel->n_edges(); ++e){
      lhs += x[e];

      // edge constraints
      Edge * edge = probModel->prob_graph.getEdge(e);
      qol::Expression e_lhs = 2*x[e];
      qol::Expression e_rhs = y[edge->getTgtID()] + y[edge->getSrcID()];
      std::string var_name("e_" + e);
      mip.addConstraint(e_lhs <= e_rhs).setName(var_name);      
    }
    for (int n = 0; n < probModel->n_nodes(); ++n){
      rhs += y[n];

      // // node constraints
      // Node * node = probModel->prob_graph.getNode(n);
      // qol::Expression n_lhs = y[n];
      // qol::Expression n_rhs;

      // std::vector<Edge*> adj_edges; 
      // node->getEdges(adj_edges);
      // for (std::vector<Edge*>::iterator e = adj_edges.begin(); e != adj_edges.end(); ++e){
      //   n_rhs += x[(*e)->getID()];
      // }
      // std::string var_name("n_" + n);
      // mip.addConstraint(n_lhs <= n_rhs).setName(var_name); 
    }

    rhs -= 1;

    mip.addConstraint(lhs <= rhs).setName("init_lt");
    mip.addConstraint(lhs >= rhs).setName("init_gt");

    PE("done!")
}
