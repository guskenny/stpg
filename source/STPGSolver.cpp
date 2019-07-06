#include <STPGSolver.h>

STPGSolver::STPGSolver(int argc, const char **argv){

	sh = SettingsHandler();
	probModel = new STPGModel(argv[argc-1],sh);

}

void STPGSolver::solve(){
	
	set_obj init_sol = set_obj(probModel->n_edges());
	int init_val = 0;

	bool sol_valid = false;

	while (!sol_valid){

		if (sh.INIT_SOL_TYPE==0){
			getInitSol(init_sol);
		}
		else{
			getRandSol(init_sol);
		}

		init_val = get_sol_value(probModel->prob_graph, init_sol);

		// print_tree(probModel->prob_graph, sol, std::string(probModel->getFPath()+probModel->getName()+".tree"));
		// PAUSE

		PE("verifying initial solution...")
		sol_valid = verify(probModel->prob_graph, probModel->getTerms(), init_sol);

		if (!sol_valid){
			PE("solution not valid! trying again!")
			continue;
		}

		CHECK(get_sol_value(probModel->prob_graph, init_sol))

	}
	// STPGMip MIPsolver(sh,probModel,init_sol);

	// set_obj temp_sol = set_obj(probModel->n_edges());

	// MIPsolver.solve(init_sol);

	// PE("verifying MIP solution...")
	// verify(probModel->prob_graph, probModel->getTerms(), init_sol);
	
	// CHECK(get_sol_value(probModel->prob_graph, init_sol))

	
	qol::CpuTimer solve_timer;

	std::vector<int> obj_vals;
	std::vector<int> sol_vals;

	int best_val = init_val;

	std::vector<set_obj> full_sols;

	set_obj best_sol = init_sol;

	set_obj best_full_sol;

	for (int iter = 0; iter < sh.NUM_ITER; ++iter){
		LocalSearch ls(sh,probModel);
		CHECK(iter)
		full_sols.clear();
		init_val = get_sol_value(probModel->prob_graph, best_sol);
		for (int swap = 0; swap < sh.NUM_MOVES; ++swap){
			PF("\riteration: " << swap << " of " << sh.NUM_MOVES)
			set_obj sol = best_sol;
			// ls.doSwaps(sol);
			int curr_val = ls.singleSwaps(sol,obj_vals);
			
			// PE("verifying swap solution...")
			// verify(probModel->prob_graph, probModel->getTerms(), sol);

			// create full solution with nodes and edges for merging
			set_obj full_sol(probModel->n_edges() + probModel->n_nodes());
			for (int i = 0; i < sol.size(); ++i){
				full_sol.addElement(sol[i]);
				Edge *edge = probModel->prob_graph.getEdge(sol[i]);
				// add end nodes for each edge, index for nodes is |E| + node index
				full_sol.addElement(probModel->n_edges() + edge->getSrcID());
				full_sol.addElement(probModel->n_edges() + edge->getTgtID());
			}

			full_sols.push_back(full_sol);

			if (best_val > curr_val){
				best_val = curr_val;
				best_sol = sol;
				best_full_sol = full_sol;
			}

			sol_vals.push_back(curr_val);
		}

		PE("\nbefore swapping: " << init_val) 
		PE("after swapping: " << best_val << std::endl)


		std::ostringstream b;

		b << "Best sol: " << std::endl << "x = { ";
		for (int e = 0; e < best_sol.size(); ++e){
			Edge *edge = probModel->prob_graph.getEdge(best_sol.get(e));
			b << "(" << edge->getSrcID()+1 << "," << edge->getTgtID()+1 << ") "; 
		}

		b << "}" << std::endl;

		VE(b.str())

		SolutionMerger sm = SolutionMerger(sh);

		merge_sol merged_sol(probModel->n_edges(), probModel->n_nodes());

		sm.merge(full_sols, merged_sol);

		set_obj post_merge_sol(probModel->n_edges());

		STPGMergeMip MIPsolver(sh,probModel,merged_sol,best_full_sol);
		MIPsolver.solve(post_merge_sol);


		best_full_sol.clear();

		for (int i = 0; i < post_merge_sol.size(); ++i){
				best_full_sol.addElement(post_merge_sol[i]);
				Edge *edge = probModel->prob_graph.getEdge(post_merge_sol[i]);
				// add end nodes for each edge, index for nodes is |E| + node index
				best_full_sol.addElement(probModel->n_edges() + edge->getSrcID());
				best_full_sol.addElement(probModel->n_edges() + edge->getTgtID());
		}

		PE("\nbefore merge: " << get_sol_value(probModel->prob_graph, best_sol)) 
		PE("after merge: " << get_sol_value(probModel->prob_graph, post_merge_sol) << std::endl)

		// if (get_sol_value(probModel->prob_graph, post_merge_sol) < get_sol_value(probModel->prob_graph, best_sol)){
		// 	std::cin.get();
		// }
		
		best_sol = post_merge_sol;
		best_val = get_sol_value(probModel->prob_graph, best_sol);
	}
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
	std::vector<int> terms = probModel->getTerms();
	std::vector<std::vector<int>> ants(terms.size(), std::vector<int>());
	for (int ant = 0; ant < terms.size(); ++ant){
		ants[ant].push_back(terms[ant]);
	}
	std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 rng = std::mt19937(seed);


	// PF("ants: ")
	// for (int i = 0 ; i < ants.size(); ++i){
	// 	PF(ants[i] << " ")
	// }
	// PE("")

	// std::shuffle(std::begin(ants), std::end(ants), rng);

	// PF("ants: ")
	// for (int i = 0 ; i < ants.size(); ++i){
	// 	PF(ants[i] << " ")
	// }
	// PE("")

	std::vector<int> nodes_visited(probModel->prob_graph.getNumNodes(),-1);

	bool ants_alive = true;
	while (ants_alive){
		ants_alive = false;
		std::vector<int> connected;
		for (int ant = 0; ant < ants.size(); ++ant){
			if (ants[ant].empty()){
				continue;
			}
			ants_alive = true;
			int ant_node_id = ants[ant].back();
			Node *ant_node = probModel->prob_graph.getNode(ant_node_id);
			connected.clear();
			ant_node->getConnectedIDs(connected);
			std::shuffle(std::begin(connected), std::end(connected), rng);
			bool moved = false;
			for (int n = 0; n < connected.size(); ++n){
				if (nodes_visited[connected[n]] == ant){
					continue;
				}
				else if (nodes_visited[connected[n]] < 0){
					Edge *edge = probModel->prob_graph.getEdgePair(ant_node_id,connected[n]);
					nodes_visited[connected[n]] = ant;
					ants[ant].push_back(connected[n]);
					sol.addElement(edge->getID());
					moved = true;
					break;
				}
				else if (nodes_visited[connected[n]] > -1){
					Edge *edge = probModel->prob_graph.getEdgePair(ant_node_id,connected[n]);
					ants[ant].clear();
					sol.addElement(edge->getID());
					moved = true;
					break;
				}
			}
			if (!moved){
				ants[ant].pop_back();
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

	// generate distance graph for terminals
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