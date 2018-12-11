#include <STPGSolver.h>

STPGSolver::STPGSolver(int argc, const char **argv){

	sh = SettingsHandler();
	probModel = new STPGModel(argv[argc-1],sh);

}

void STPGSolver::solve(){
	
	set_obj init_sol = set_obj(probModel->n_edges());

	if (sh.INIT_SOL_TYPE==0){
		getInitSol(init_sol);
	}
	else{
		getRandSol(init_sol);
	}

	int init_val = get_sol_value(probModel->prob_graph, init_sol);

	// print_tree(probModel->prob_graph, sol, std::string(probModel->getFPath()+probModel->getName()+".tree"));
	// PAUSE

	PE("verifying initial solution...")
	verify(probModel->prob_graph, probModel->getTerms(), init_sol);


	STPGMip MIPsolver(sh,probModel,init_sol);

	set_obj temp_sol = set_obj(probModel->n_edges());


	MIPsolver.solve(temp_sol);

	LocalSearch ls(sh,probModel);
	
	qol::CpuTimer solve_timer;

	std::vector<int> obj_vals;
	std::vector<int> sol_vals;

	int best_val = init_val;

	std::vector<set_obj> sols;

	for (int swap = 0; swap < sh.NUM_MOVES; ++swap){
		PE("\niteration: " << swap << " of " << sh.NUM_MOVES)
		set_obj sol = init_sol;

		// ls.doSwaps(sol);
		int curr_val = ls.singleSwaps(sol,obj_vals);

		sols.push_back(sol);

		best_val = std::min(curr_val, best_val);

		sol_vals.push_back(curr_val);
	}

	PE("\nbefore swapping: " << init_val) 
	PE("after swapping: " << best_val)

	SolutionMerger sm(sh, probModel);

	std::vector<std::vector<int> > groups;
	std::vector<int> group_map;

	sm.merge(sols, groups, group_map);

	std::ofstream _out;

	if (sh.RECORD_DATA){
		_out.open("tracking_data/" + probModel->getInstGrp() + probModel->getName() + ".out");
		_out << "OBJECTIVE TRACKING - " << probModel->getName() << std::endl;
	}
	
	_out << "! NAME " << probModel->getName() << std::endl;
	_out << "! CPU_TIME " << solve_timer.elapsedSeconds() << std::endl;
	_out << "! SWAPS " << sh.DEAD_MOVES << std::endl;
	_out << "! JUMP_FREQ " << sh.NUM_JUMPS << std::endl;
	_out << "! POP_SIZE " << sh.NUM_MOVES;

	for (int i = 0; i < sol_vals.size(); ++i){
		_out << sol_vals[i] << std::endl;		
	}

	_out.close();

}

// int STPGSolver::getSolValue(const set_obj &sol){
// 	std::vector<int> edges = sol.getSet();

// 	int total_wt = 0;

// 	for (int e = 0; e < edges.size(); ++e){
// 		total_wt += probModel->prob_graph.getEdge(edges[e])->getWt();
// 	}

// 	return total_wt;
// }


// int STPGSolver::getSolValue(const set_obj &sol, const U_Graph &graph){
// 	std::vector<int> edges = sol.getSet();

// 	int total_wt = 0;
	
// 	for (int e = 0; e < edges.size(); ++e){
// 		total_wt += graph.getEdge(edges[e])->getWt();
// 		PE2("edge (" << graph.getEdge(edges[e])->getSrcID() << ", " << graph.getEdge(edges[e])->getTgtID() << ")")
// 	}

// 	return total_wt;
// }

void STPGSolver::getRandSol(set_obj &sol){

	PE("Getting random solution")
	// initialise ants on terminal nodes
	std::vector<int> ants = probModel->getTerms();

	std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 rng = std::mt19937(seed);


	PF("ants: ")
	for (int i = 0 ; i < ants.size(); ++i){
		PF(ants[i] << " ")
	}
	PE("")

	std::shuffle(std::begin(ants), std::end(ants), rng);

	PF("ants: ")
	for (int i = 0 ; i < ants.size(); ++i){
		PF(ants[i] << " ")
	}
	PE("")

	std::vector<int> visited(probModel->prob_graph.getNumNodes(),-1);

	bool ants_alive = true;
	while (ants_alive){
		ants_alive = false;
		for (int ant = 0; ant < ants.size(); ++ant){
			if (ants[ant] == -1){
				continue;
			}
			ants_alive = true;
			Node *ant_node = probModel->prob_graph.getNode(ants[ant]);
			std::vector<Node*> connected;
			ant_node->getConnected(connected);
			// std::shuffle(std::begin(connected), std::end(connected), rng);
			for (int n = 0; connected.size(); ++n){
				if (visited[connected[n]->getID()] == ant){
					continue;
				}
				if (visited[connected[n]->getID()] == -1){
					visited[connected[n]->getID()] = ant;
					Edge *edge = probModel->prob_graph.getEdgePair(ants[ant],connected[n]->getID());
					ants[ant] = connected[n]->getID();
					sol.addElement(edge->getID());
					break;
				}
				if (visited[connected[n]->getID()] > -1){
					Edge *edge = probModel->prob_graph.getEdgePair(ants[ant],connected[n]->getID());
					ants[ant] = -1;
					sol.addElement(edge->getID());
					break;
				}
			}
		}
	}
	PE("fixing solution")
	// remove any cycles
	get_mst(probModel->prob_graph, sol);

	// prune tree of any dangling vertices
	prune_subgraph(probModel->prob_graph, sol);	
}

void STPGSolver::getInitSol(set_obj &sol){

	// get terminal nodes
	std::vector<int> terms = probModel->getTerms();

	// generate distange graph for terminals
	U_Graph dist(terms.size());

	std::vector<int> dist_terms;
	for (int t = 0; t < terms.size(); ++t){
		for (int s = 0; s < t; ++s){
			dist.addEdge(s,t,probModel->fw_dist[terms[s]][terms[t]]);

		}
		PF2(terms[t] << " ")
		dist_terms.push_back(t);
	}
	PE2("")

	// get MST of distance graph
	set_obj dist_sol(dist.getNumEdges());
	dist_sol.setAll();

	get_mst(dist, dist_sol);

	std::vector<int> dist_edges = dist_sol.getSet();

	sol.clear();

	int total_path=0;

	// add all shortest paths to solution
	for (int e = 0; e < dist_edges.size(); ++e){
		std::vector<int> path_edges;
		int u = terms[dist.getEdge(dist_edges[e])->getSrcID()];
		int v = terms[dist.getEdge(dist_edges[e])->getTgtID()];
		int path_length = probModel->get_path(path_edges, u, v);
		sol.addElements(path_edges);
		PE2("path length: " << path_length)
		total_path += path_length;
	}

	// remove any cycles
	get_mst(probModel->prob_graph, sol);

	// prune tree of any dangling vertices
	prune_subgraph(probModel->prob_graph, sol);	
}