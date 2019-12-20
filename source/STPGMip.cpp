#include "STPGMip.h"

STPGMip::STPGMip(const SettingsHandler sh, STPGModel *probModel, set_obj &best_sol) : probModel(probModel), best_sol(best_sol), sh(sh) {

  PE("\nMerge Solver initialised")
}

void STPGMip::solve(set_obj &sol){

  // PE("Testing terminals:")
  // std::vector<int> terms = probModel->getTerms();
  // PF("getTerms terminals: ")
  // for (int i = 0; i < terms.size(); ++i){
  //   PF(terms[i] << " ")
  // }
  // PE(" ")
  // PF("isTerm terminals: ")
  // for (int i = 0; i < probModel->n_nodes(); ++i){
  //   if (probModel->prob_graph.getNode(i)->isTerm()){
  //     PF(i << " ")
  //   }
  // }
  // PE(" ")

  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);
    param.setParamVal(qol::RELGAP,0.001);
    if (sh.MIP_TIME){
      param.setParamVal(qol::TIMELIMIT,sh.MIP_TIME);
    }

    PE("initialising MIP..")
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

  VE("\n*****************************************")
  VE("MIP solution " << ++(*countPtr))
  VE("*****************************************")

  int cost = 0;
  double frac_cost = 0;

  VF("x = { ")

  std::ostringstream latex_x;
  std::ostringstream latex_y;
  std::ostringstream latex_cut_x;
  std::ostringstream latex_cut_y;

  for (int e = 0; e < primal_x.idx_size(); ++e){
    if (primalSolution[x[e]] > 0){
      primal_x.addElement(e);

      Edge * edge = probModel->prob_graph.getEdge(e);

      cost += edge->getWt();
      frac_cost += primalSolution[x[e]]*edge->getWt();

      VF(primalSolution[x[e]]<<"(" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") ")

      latex_x << std::to_string(edge->getSrcID()+1) << "/" << std::to_string(edge->getTgtID()+1)<< ",";

      test_y.addElement(edge->getTgtID());
      test_y.addElement(edge->getSrcID());
    }
  }

  VE("}")
  VF("y = { ")

  for (int n = 0; n < primal_y.idx_size(); ++n){
    if (primalSolution[y[n]] > 0){
      primal_y.addElement(n);
      VF(n+1 << " ")
      latex_y << n+1 << ",";
    }
  }

  VE("}")

  VE("full cost: " << cost << std::endl)
  VE("solution cost: " << frac_cost << std::endl)
  // CHECK(test_y.size())
  // CHECK(primal_y.size())
  // CHECK(test_y.size()-primal_y.size())
  // test_y.subtract(primal_y);
  // CHECK(test_y.size())
  // VE("")

  VE("x = " <<latex_x.str())
  VE("y = " << latex_y.str()<<std::endl)

  std::vector<set_obj> comp_edges;
  std::vector<set_obj> comp_nodes;
  std::vector<int> comp_type; // if component contains a terminal

  // get connected components with cycles, each component is a set of edges
  // CHECK(x.size())
  // PF("getting components.. ")
  get_components(probModel->prob_graph,primal_x,primal_y,comp_edges,comp_nodes,comp_type);
  // PE("done!")

  // CHECK(comp_edges.size())
  // CHECK(comp_nodes.size())
  // CHECK(is_term.size())

  // for each component, add a cut to the model using addCut()
  for (int c = 0; c < comp_edges.size(); ++c){
    std::vector<int> comp_edge = comp_edges[c].getSet();
    std::vector<int> comp_node = comp_nodes[c].getSet();
    
    qol::Expression lhs;
    qol::Expression rhs;

    std::ostringstream lhs_os;
    std::ostringstream rhs_os;

    // vector to make sure we dont double up on node variables
    std::vector<int> added(probModel->n_nodes(),0);

    int coeff = 1;

    // CHECK(comp_type[c])
    // switch (comp_type[c]) {
    //     case NOTERM_CYCLE: comp_node.pop_back(); break;
    //     case TERM_CYCLE:  rhs += -1;rhs_os << "(-1)"; break;
    //     case SINGLE_EDGE: rhs += -1;rhs_os << "(-1)"; break;
    //     default: {rhs += -1;rhs_os << "(-1) +";}
    // }

    if (comp_type[c] == NOTERM){
      comp_node.pop_back();
    }
    else{
      rhs += -1;
      rhs_os << "(-1)";
    }


    // add edges to cut
    for (std::vector<int>::iterator e = comp_edge.begin() ; e != comp_edge.end(); ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      lhs += coeff*x[*e];

      lhs_os << "x[" + std::to_string(edge->getSrcID()+1) + "," + std::to_string(edge->getTgtID()+1) + "]";
      latex_cut_x << std::to_string(edge->getSrcID()+1) << "/" << std::to_string(edge->getTgtID()+1)<< ",";

      if(e != comp_edge.end()-1)
        lhs_os << " + ";

    }


    // bool term = false;

    // for (int i = 0; i < comp_node.size(); ++i){
    //   if (probModel->prob_graph.getNode(comp_node[i])->isTerm())
    //     term = true;
    // }

    // if no terminals in component, remove one arbitrary node (back node)
     // add nodes to cut
    for (std::vector<int>::iterator n = comp_node.begin() ; n != comp_node.end(); ++n){
      rhs += y[*n];
      rhs_os << " + y[" + std::to_string((*n)+1) + "]";
      latex_cut_y << *n+1 << ",";
    }

    //   std::vector<int> endpoints{edge->getSrcID(),edge->getTgtID()};

    //   for (std::vector<int>::iterator i = endpoints.begin() ; i != endpoints.end(); ++i){
    //     if (!added[*i]){
    //       rhs += y[*i];
    //       rhs_os << "y[" + std::to_string((*i)+1) + "] + ";
    //       added[*i] = 1;
    //     }
    //   }
    // }

    // rhs -= 1;
    // rhs_os << "(-1)";
    VE("adding cut: " << lhs_os.str() << " <= " << rhs_os.str()<< std::endl)


    addCut(lhs <= rhs);

  }

  // CHECK(getBound())
  VE("\ncut_x = " <<latex_cut_x.str())
  VE("cut_y = " << latex_cut_y.str()<<std::endl)


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

    double init_sum = 0;

    std::ostringstream x_start;
    std::ostringstream y_start;

    x_start << "x = {";

    y_start << "y = { ";

    for (std::vector<int>::iterator e = sol_edges.begin() ; e != sol_edges.end() && sh.WARM_START; ++e){
      Edge* edge = probModel->prob_graph.getEdge(*e);

      mip.setPrimalStart(x[*e],1);
      mip.setPrimalStart(y[edge->getSrcID()],1);
      mip.setPrimalStart(y[edge->getTgtID()],1);

      y_start << edge->getSrcID()+1 << ", " << edge->getTgtID()+1 << ", ";

      x_start << " (" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << "): " << edge->getWt();

      init_sum += edge->getWt();
    }

    x_start << "}";
    y_start << "}";


    PE("\n" << x_start.str())
    PE("\n" << y_start.str())

    PE("done!")
    CHECK(init_sum)
    PF("setting up initial constraints.. ")

    qol::Expression lhs;
    qol::Expression rhs;
    for (int e = 0; e < probModel->n_edges(); ++e){
      lhs += x[e];

      // edge constraints
      Edge * edge = probModel->prob_graph.getEdge(e);
      qol::Expression e_lhs = x[e];
      qol::Expression n1_rhs = y[edge->getTgtID()]; 
      qol::Expression n2_rhs = y[edge->getSrcID()];

      std::string var_name1("e1_" + e);
      std::string var_name2("e2_" + e);
      mip.addConstraint(e_lhs <= n1_rhs).setName(var_name1);
      mip.addConstraint(e_lhs <= n2_rhs).setName(var_name2);      
    }
    for (int n = 0; n < probModel->n_nodes(); ++n){
      rhs += y[n];
    }

    rhs += -1;

    mip.addConstraint(lhs == rhs).setName("init_main");
    // mip.addConstraint(lhs >= rhs).setName("init_gt");

    PE("done!")
}
