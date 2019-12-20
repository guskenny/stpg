#include "STPGMergeMip.h"

STPGMergeMip::STPGMergeMip(const SettingsHandler sh, STPGModel *probModel, merge_sol &groups, set_obj &best_sol) : probModel(probModel), groups(groups), best_sol(best_sol), sh(sh) {

  PE("\nMerge Solver initialised")
}

void STPGMergeMip::solve(set_obj &sol){

  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);
    // param.setParamVal(qol::RELGAP,0.001);
    if (sh.MIP_TIME){
      param.setParamVal(qol::TIMELIMIT,sh.MIP_TIME);
    }

	  mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    initMIPModel(mip);

    int count = 0;
    int *countPtr = &count;

    STPGMergeCallback cb = STPGMergeCallback(probModel,z,countPtr,mipPtr,groups);

    mip.setCallback((qol::Callback *)&cb);
    // mip.setCallback(&cb);

    CPXENVptr env=(dynamic_cast<qol::CplexFormulation *>(mipPtr))->env;
    
    CPXsetintparam(env,CPX_PARAM_PREIND,sh.PRESOLVE);

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

    // if (sh.STOP_OPTIMAL){
    //   CPXsetdblparam(env,CPX_PARAM_CUTUP,probModel->getOptimal());
    // }

    if (sh.HEURISTIC_SEARCH){
      CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR,1.0);
      CPXsetintparam(env,CPX_PARAM_HEURFREQ,1);
      // CPXsetintparam(env,CPX_PARAM_POLISHAFTEREPAGAP,0.1);
      CPXsetintparam(env,CPX_PARAM_RINSHEUR,1);
      CPXsetintparam(env,CPX_PARAM_POLISHAFTERINTSOL,1);
    }    
  
    CHECK(mip.nConstr()) 
    CHECK(mip.nVar())  

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());
   
    sol.clear();

    // extract solution from z variables
    for (int g = 0; g < groups.size(); ++g){
      if (mip.getPrimal(z[g]) > 0){
        for (int v = 0; v < groups[g].size(); ++v){
          if (groups[g][v] < probModel->n_edges()){
            sol.addElement(groups[g][v]);
          }
        }
      }
    }

   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
}

qol::Callback::Status STPGMergeCallback::callback(Progress where){

  // PE("CALLBACK!")

  set_obj primal_x(probModel->n_edges());
  set_obj primal_y(probModel->n_nodes());


  // extract primal x (edge) and primal y (node) information from z variables
  for (int g = 0; g < groups.size(); ++g){
    if (primalSolution[z[g]] > 0){
      for (int v = 0; v < groups[g].size(); ++v){
        if (groups[g][v] < probModel->n_edges()){
          primal_x.addElement(groups[g][v]);
        }
        else{
          primal_y.addElement(groups[g][v] - probModel->n_edges()); 
        }
      }
    }
  }

  std::vector<set_obj> comp_edges;
  std::vector<set_obj> comp_nodes;
  std::vector<int> comp_type; // if component contains a terminal (not used for merge formulation)

  // uses the same get_components function as previously to return a vector of edge sets and
  // node sets for each component
  get_components(probModel->prob_graph,primal_x,primal_y,comp_edges,comp_nodes,comp_type);

  // for each component, add a cut to the model using addCut()
  for (int c = 0; c < comp_edges.size(); ++c){
    // PE("component "<<c << ", size: " << comp_edges[c].size())
    // vector to compute all coefficients at once
    std::vector<int> coeffs(groups.size(), 0);


    // iterate over all edges in component and subtract edges from coefficients
    for (int e = 0; e < comp_edges[c].size(); ++e){
      // PE("c: "<< c << ", e: " << e << ", comp_edges[c][e]: " << comp_edges[c][e])
      // PE("groups.find_group("<<comp_edges[c][e]<<") = " << groups.find_group(comp_edges[c][e]))
      coeffs[groups.find_group(comp_edges[c][e])]--;
    }

    // iterate over all nodes in component and add nodes to coefficients (skip first node)
    for (int n = 1; n < comp_nodes[c].size(); ++n){
      coeffs[groups.find_node_group(comp_nodes[c][n])]++;
    }

    // iterate over groups and if coefficient isnt zero then add coeff[g] * z_g to cut
    qol::Expression lhs;
    for (int g = 0; g < groups.size(); ++g){
      if (coeffs[g] != 0){
        lhs += coeffs[g] * z[g];
      }
    }

    addCut(lhs >= 0);
  }
    return OK;
}

void STPGMergeMip::initMIPModel(qol::MIPSolver &mip){
    PE("*** Initialising MIP model ***")

    z.clear();

    PF("setting up binary variables..")
    // set up binary variables for groups, iterates over each group and computes the
    // cost of each edge within that group
    for (int g = 0; g < groups.size(); ++g){
      int cost = 0;
      int count = 0;
      // PF("\nSetting vars for group " << g << ":")
      // iterate over all variables in the group
      for (int v = 0; v < groups[g].size(); ++v){
        int var_idx = groups[g][v];
        // if variable index indicates edge, add to total cost for group
        if (var_idx < probModel->n_edges()){
          cost += probModel->prob_graph.getEdge(var_idx)->getWt();
        }
        // PF(" " << var_idx)
        count++;
      }

      // create binary variable for group
      z.push_back(mip.addVar(0,1,cost,qol::Variable::BINARY));
      std::string var_name("z_" + g);
      mip.setVarName(z[g],var_name);
      // PF(" - group: " << z.size()-1 << ", cost: "<< cost << ", count: " << count)
    }

    // set warm start values
    if (sh.WARM_START){
      PF("done!\nsetting warm start..")
      for (int i = 0; i < best_sol.size(); ++i){
        mip.setPrimalStart(z[groups.find_group(best_sol[i])], 1);
      }
    }

    PF("done!\nsetting groups with terminals to 1..")
    // set all groups containing terminals to 1
    std::vector<int> terms = probModel->getTerms();
    for (int t = 0; t < terms.size(); ++t){
      // add index for node to number of edges to get position in group vector
      mip.setVarLB(z[groups.find_node_group(terms[t])], 1);
    }

    PF("done!\nsetting tree constraints..")

    // set initial "tree" constraint
    // iterates over all groups and computes a coefficient of |N| - |E| for each group
    // the sum of these coefficients times the binary variables must equal 1
    // i.e., there must be exactly one less edge than node in the solution

    qol::Expression lhs;

    for (int g = 0; g < groups.size(); ++g){
      int coeff = 0;
      for (int v = 0; v < groups[g].size(); ++v){
        if (groups[g][v] < probModel->n_edges()){
          coeff--;
        }
        else{
          coeff++;
        }
      }
      lhs += coeff*z[g];
    }
    mip.addConstraint(lhs == 1).setName("init_tree");

    PF("done!\nsetting group constraints..")

    // set initial group constraints
    // z_p <= z_q if p contains edge (i,j) and q contains node i and p != q
    // iterate over groups (z_p) and find the endpoints of all edges in that group
    // then find the related groups for all of those endpoints (z_q) and add constraint
    // z_p <= z_q

    for (int p = 0; p < groups.size(); ++p){
      qol::Expression z_p = z[p];

      // use set object data structure to avoid duplicates
      set_obj end_groups(groups.size());

      for (int v = 0; v < groups[p].size(); ++v){
        if (groups[p][v] >= probModel->n_edges()){
          continue;
        }
        // add the groups of each endpoint to the list of groups to make
        // constraints about
        Edge * edge = probModel->prob_graph.getEdge(groups[p][v]);
        end_groups.addElement(groups.find_node_group(edge->getTgtID()));
        end_groups.addElement(groups.find_node_group(edge->getSrcID()));
      }

      // iterate over all endpoint groups and if not the same as z_p, add constraint
      // z_p <= z_q
      for (int q_idx = 0; q_idx < end_groups.size(); ++q_idx){
        if (end_groups[q_idx] == p){
          continue;
        }
        qol::Expression z_q = z[end_groups[q_idx]];

        mip.addConstraint(z_p <= z_q);//.setName(std::string("z_")+std::string(p)+"_"+std::string(end_groups[q_idx]));
      }
    }
    PE("done!\n")
}
