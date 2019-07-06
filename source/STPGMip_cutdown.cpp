#include "STPGMip.h"

/* 
This is a cutdown version of what i am using, to make it a bit easier to follow
I have also put some "helper functions" such as get_components() from a utility
file at the bottom..

the problem instances are converted into a custom undirected graph class, and
then put into a class STPGModel which holds all the instance information..

solutions are the custom class set_obj (set object) which is a data structure that
allows for fast adding, removal and retrieval of objects in a set of integers

SettingsHandler is a class that reads a settings.ini file at the start and sets
all the parameters so it doesnt have to be recompiled every time
*/

// callback class for STPG
class STPGCallback : public qol::Callback{
    STPGModel *probModel;
    const std::vector<qol::Variable> &x;
    const std::vector<qol::Variable> &y;
    int *countPtr;
    qol::MIPSolver *mipPtr;
    public:
    STPGCallback(STPGModel *_probModel,const std::vector<qol::Variable> &_x, const std::vector<qol::Variable> &_y,int *_countPtr,qol::MIPSolver *_mipPtr):probModel(_probModel),x(_x),y(_y),countPtr(_countPtr),mipPtr(_mipPtr){};
    ~STPGCallback(){};
    Status callback(Progress where);
};

STPGMip::STPGMip(const SettingsHandler sh, STPGModel *probModel, set_obj &best_sol) : probModel(probModel), best_sol(best_sol), sh(sh) {}

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

    // initialise mip model
    initMIPModel(mip);

    int count = 0;
    int *countPtr = &count;

    // set up callbacks
    STPGCallback cb = STPGCallback(probModel,x,y,countPtr,mipPtr);

    // set callback function in solver
    mip.setCallback((qol::Callback *)&cb);

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
   
    sol.clear();

    // extract solution from MIP model
    for (int e = 0; e < probModel->n_edges(); ++e){
      if (mip.getPrimal(x[e]) > 0){
        sol.addElement(e);
      }
    }    
   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
}

qol::Callback::Status STPGCallback::callback(Progress where){

  // sets for nodes and edges in primal solution
  set_obj primal_x(probModel->n_edges());
  set_obj primal_y(probModel->n_nodes());

  set_obj test_y(probModel->n_nodes());

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

    // create cut
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

    // start values
    std::vector<int> sol_edges = best_sol.getSet();

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
    }

    rhs -= 1;

    mip.addConstraint(lhs <= rhs).setName("init_lt");
    mip.addConstraint(lhs >= rhs).setName("init_gt");
}


// here are the "helper functions"


// function to extract all connected components in a solution and include "bridge" edges
// within components
void get_components(const U_Graph &graph, const set_obj &x,const set_obj &y, std::vector<set_obj> &components){

  // clear components vector
  components.clear();

  // set_obj for visited vertices
  set_obj visited(graph.getNumNodes());

  // find all nodes induced by the edge set
  for (int i = 0; i < x.size(); ++i){
    Edge *e = graph.getEdge(x.get(i));
    int src = e->getSrcID();
    int tgt = e->getTgtID();

    set_obj temp_comp(graph.getNumEdges());
    if (!y.is_element(src) || !y.is_element(tgt)){
      temp_comp.addElement(x.get(i));
      components.push_back(temp_comp);
    }
    visited.addElement(e->getSrcID());
    visited.addElement(e->getTgtID());
  }

  // find components
  while (!visited.empty()){
    std::vector<int> stack;
    // check if vertex is in solution
    stack.push_back(visited.get(0));
    
    bool cycle = false;

    // make set_obj for temporary vertices and edges
    set_obj temp_x(x.idx_size());
    set_obj temp_y(y.idx_size());

    while(!stack.empty()){
      // get current node
      int curr_idx = stack.back();
      stack.pop_back();

      // mark current as visited
      visited.removeElement(curr_idx);

      // add current to temp_y
      temp_y.addElement(curr_idx);

      Node *node = graph.getNode(curr_idx);
      std::vector<Edge*> adj_edges;
      node->getEdges(adj_edges);
      
      // iterate over adjacent edges and if edge is in x, add it to component
      for (std::vector<Edge*>::iterator e = adj_edges.begin(); e != adj_edges.end(); ++e){
        // skip edge if not in x, else add it to component
        if (!x.is_element((*e)->getID())){
          continue;
        }
        else{
          temp_x.addElement((*e)->getID());
        }

        // find endpoints of the edge
        std::vector<int> endpoints;   
        endpoints.push_back((*e)->getSrcID());
        endpoints.push_back((*e)->getTgtID());

        // if both endpoints already in temp_y, flag that cycle exists, otherwise add endpoints
        if (temp_y.is_element(endpoints[0]) && temp_y.is_element(endpoints[1])){
          cycle = true;
        }
        else{
          for (int i = 0; i < endpoints.size(); ++i){
            if (!temp_y.is_element(endpoints[i])){
              stack.push_back(endpoints[i]);  
            }
          }
        }
      } // end iterator
    } // end stack while

    // prune the subgraph to remove all but the cycle
    prune_subgraph_noterm(graph,temp_x);

    // reset temp_y to store nodes of pruned subgraph
    temp_y.clear();

    // add all nodes from new subgraph
    for (int e=0; e < temp_x.size(); ++e){
      Edge * edge = graph.getEdge(temp_x.get(e));
      temp_y.addElement(edge->getSrcID());
      temp_y.addElement(edge->getTgtID());
    }

    // find all extra edges that bridge component
    for (int i=0; i < temp_y.size()-1; ++i){
      for (int j=i+1; j < temp_y.size(); ++j){
        Edge * e = graph.getEdgePair(temp_y.get(i),temp_y.get(j));
        if (e){
          temp_x.addElement(e->getID());
        }
      }
    }   

    if (temp_x.size() > 0){
      components.push_back(temp_x);
    }

  } // end visit while
}

// removes "dangling" edges from subgraph, regardless of terminal status (for creating
// cuts based on components)
void prune_subgraph_noterm(const U_Graph &graph, set_obj &subgraph){

  int subgraph_size = subgraph.size() + 1;
  
  std::vector<int> degree(graph.getNumNodes(), 0);

  // populate degree vector
  for (int e = 0; e < subgraph.size(); ++e){
    int src = graph.getEdge(subgraph.get(e))->getSrcID();
    int tgt = graph.getEdge(subgraph.get(e))->getTgtID();
    degree[src]++;
    degree[tgt]++;
  }

  // loop while there are still edges to remove
  while (subgraph_size - subgraph.size() != 0){
    subgraph_size = subgraph.size();

    // delete edges with endpoint with only one degree
    for (int e = 0; e < subgraph.size(); ++e){
      int src = graph.getEdge(subgraph.get(e))->getSrcID();
      int tgt = graph.getEdge(subgraph.get(e))->getTgtID();
      if (degree[src] < 2 || degree[tgt] < 2){
        subgraph.removeElement(subgraph.get(e));
        degree[src]--;
        degree[tgt]--;
        break;
      }
    }
  }
}