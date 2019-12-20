#include <LocalSearch.h>

int LocalSearch::singleSwaps(set_obj &sol, std::vector<int> &obj_vals){
	int dead_count = 0;
	int curr_obj = 0;
	int best_obj = get_sol_value(probModel->prob_graph, sol);

	obj_vals.push_back(best_obj);

	set_obj best_sol = sol;

	// PE(" ")

	// initial jump away from the start
	for (int i = 0; i < sh.NUM_INIT_JUMPS; ++i){
		// get keypath
		std::vector<int> kp_edges;
		std::vector<int> kp_ends;
		get_keypath(sol, kp_edges, kp_ends);
		// PE("keypath: (" << kp_ends[0] <<  ", " << kp_ends[1] << ")")
		// remove keypath from solution
		sol.removeElements(kp_edges);
		// PF("blerp.. ")
		repairRandom(sol,kp_ends);
		// PF(" done.. ")
	}

	while (dead_count < sh.DEAD_MOVES){
		// PF("\rdead count: " << dead_count << " of " << sh.DEAD_MOVES)
		// get keypath
		std::vector<int> kp_edges;
		std::vector<int> kp_ends;
		get_keypath(sol, kp_edges, kp_ends);
		// remove keypath from solution
		sol.removeElements(kp_edges);
		repairBest(sol,kp_ends);
		
		int new_obj = get_sol_value(probModel->prob_graph, sol);

		if (new_obj < best_obj){
			best_obj = new_obj;
			best_sol = sol;
		}
		
		if (curr_obj != new_obj){
			dead_count = 0;
			curr_obj = new_obj;
		}
		else{
			dead_count++;
		}
		obj_vals.push_back(new_obj);
	}

	// PE("\nbest_obj: " << best_obj)

	return best_obj;
}

int LocalSearch::doSwaps(set_obj &sol, std::vector<int> &obj_vals){

	int dead_count = 0;
	int curr_obj = 0;
	int best_obj = get_sol_value(probModel->prob_graph, sol);

	// std::ofstream _out;

	// if (sh.RECORD_DATA){
	// 	_out.open("tracking_data/" + probModel->getInstGrp() + probModel->getName() + ".out");
	// 	_out << "OBJECTIVE TRACKING - " << probModel->getName() << std::endl;
	// }

	// qol::CpuTimer solve_timer;

	// std::vector<int> obj_vals;

	obj_vals.push_back(best_obj);

	set_obj best_sol = sol;

	// PE(" ")
	int new_obj = 0;

	for (int swap = 0; swap < sh.NUM_SWAPS; ++swap){

		// initial jump away from the start
		for (int i = 0; i < sh.NUM_INIT_JUMPS; ++i){
			// get keypath
			std::vector<int> kp_edges;
			std::vector<int> kp_ends;
			get_keypath(sol, kp_edges, kp_ends);
			// PE("keypath: (" << kp_ends[0] <<  ", " << kp_ends[1] << ")")
			// remove keypath from solution
			sol.removeElements(kp_edges);
			// PF("blerp.. ")
			repairRandom(sol,kp_ends);
			// PF(" done.. ")
		}

		// PF("\rmove: " << swap << " of " << sh.NUM_SWAPS)

		while (dead_count < sh.DEAD_MOVES){

			// // jump and repair solution
			// if (sh.JUMP_FREQ > 0 && dead_count > (sh.NUM_SWAPS/sh.JUMP_FREQ)){
			// 	sol = best_sol;
			// 	for (int i = 0; i < sh.NUM_JUMPS; ++i){
			// 		// get keypath
			// 		std::vector<int> kp_edges;
			// 		std::vector<int> kp_ends;
			// 		get_keypath(sol, kp_edges, kp_ends);
			// 		// remove keypath from solution
			// 		sol.removeElements(kp_edges);
			// 		repairRandom(sol,kp_ends);
			// 		// PF(" done.. ")
			// 	}
			// 	// PE("jump")
			// }
			// else{
				// get keypath
			std::vector<int> kp_edges;
			std::vector<int> kp_ends;
			get_keypath(sol, kp_edges, kp_ends);
			// remove keypath from solution
			sol.removeElements(kp_edges);
			repairBest(sol,kp_ends);
			// }		

			new_obj = get_sol_value(probModel->prob_graph, sol);

			if (new_obj < best_obj){
				best_obj = new_obj;
				best_sol = sol;
			}

			if (curr_obj != new_obj){
				dead_count = 0;
				curr_obj = new_obj;
			}
			else{
				dead_count++;
			}
		// }


		// // repair solution
		// if (sh.JUMP_FREQ > 0 && (move % (sh.NUM_MOVES/sh.JUMP_FREQ) == 0)){
		// 	repairRandom(sol,kp_ends);
		// }
		// else{
		// 	repairBest(sol,kp_ends);
		// }

		// VE("\nverifying repaired solution...")
		// verify(probModel->prob_graph, probModel->getTerms(), sol);

		// repair is not working properly.. check node sets


		}
		obj_vals.push_back(new_obj);
	}

	// _out << "! NAME " << probModel->getName() << std::endl;
	// _out << "! CPU_TIME " << solve_timer.elapsedSeconds() << std::endl;
	// _out << "! SWAPS " << sh.NUM_SWAPS << std::endl;
	// _out << "! JUMP_FREQ " << sh.JUMP_FREQ << std::endl;

	// for (int i = 0; i < obj_vals.size(); ++i){
	// 	_out << obj_vals[i] << std::endl;		
	// }

	// _out.close();

	// PE("\nbest_obj: " << best_obj)

	// PE("solve timer: " << solve_timer.elapsedSeconds())

	return best_obj;

}

void LocalSearch::getNodeSets(const set_obj &sol, const std::vector<int> &kp_ends, std::vector<std::vector<int> >&node_sets){

	// each of the keypath ends will correspond to the disjoint trees
	for (int i = 0; i < kp_ends.size(); ++i){
		int n_root = kp_ends[i];
		std::vector<int> n_stack;
		n_stack.push_back(n_root);

		std::vector<int> visited(probModel->prob_graph.getNumNodes(),0);

		while (!n_stack.empty()){
			int node_idx = n_stack.back();
			n_stack.pop_back();

			node_sets[i].push_back(node_idx);

			visited[node_idx] = 1;

			Node *node = probModel->prob_graph.getNode(node_idx);

			std::vector<Edge*> edges;
			const int n_edges = node->getEdges(edges);

			for (int e = 0; e<n_edges; ++e){
				if (!sol.is_element(edges[e]->getID())){
					continue;
				}
				std::vector<int> src_tgt;
				src_tgt.push_back(edges[e]->getSrcID());
				src_tgt.push_back(edges[e]->getTgtID());
				for (int st = 0; st < src_tgt.size(); ++st){
					if (visited[src_tgt[st]]){
						continue;
					}
					n_stack.push_back(src_tgt[st]);
				}
			}
		}
	}
	// for (int i = 0; i < node_sets.size(); ++i){
	// 	PF("set " << i << ": ")
	// 	for (int j = 0; j < node_sets[i].size(); ++j){
	// 		PF(node_sets[i][j] << " ")
	// 	}
	// 	PE(" ")
	// }
}

void LocalSearch::repairBest(set_obj &sol, const std::vector<int> &kp_ends){

	std::vector<std::vector<int> > node_sets(2,std::vector<int>());

	getNodeSets(sol, kp_ends, node_sets);

	// find shortest path connecting two node sets
	int best_u = 0;
	int best_v = 0;
	int best_path = INF;

	for (int u = 0; u < node_sets[0].size(); ++u){
		for (int v =0; v < node_sets[1].size(); ++v){
			if (probModel->fw_dist[node_sets[0][u]][node_sets[1][v]] < best_path && node_sets[0][u] != node_sets[1][v]){
				best_u = node_sets[0][u];
				best_v = node_sets[1][v];
				best_path =probModel->fw_dist[node_sets[0][u]][node_sets[1][v]]; 
			}
		}
	}

	VE("shortest path found: (" << best_u << ", " << best_v << "), value: " << best_path)

	// get edge set for shortest path
	std::vector<int> uv_edges;
	probModel->get_path(uv_edges, best_u, best_v);

	// add edge set to solution
	sol.addElements(uv_edges);

	// remove any cycles
	get_mst(probModel->prob_graph, sol);

	// prune tree of any dangling vertices
	prune_subgraph(probModel->prob_graph, sol);	
}

void LocalSearch::repairRandom(set_obj &sol, const std::vector<int> &kp_ends){

	std::vector<std::vector<int> > node_sets(2,std::vector<int>());

	getNodeSets(sol, kp_ends, node_sets);

	std::vector<int> rand_uv;	

	for (int i = 0; i < node_sets.size(); ++i){
		// find random elements of the two node sets
		std::uniform_int_distribution<int> uni(0,node_sets[i].size()-1);
	    rand_uv.push_back(node_sets[i][uni(rng)]);

	}

	VE("random path found: (" << rand_uv[0] << ", " << rand_uv[1] << "), value: " << probModel->fw_dist[rand_uv[0]][rand_uv[1]])


	// get edge set for shortest path
	std::vector<int> uv_edges;
	// PF(" get path.. " <<rand_uv[0]<< ","<< rand_uv[1]<<" ")
	probModel->get_path(uv_edges, rand_uv[0], rand_uv[1]);

	// PF(" test.. ")
	// add edge set to solution
	sol.addElements(uv_edges);


	// remove any cycles
	get_mst(probModel->prob_graph, sol);

	// prune tree of any dangling vertices
	prune_subgraph(probModel->prob_graph, sol);	
}

void LocalSearch::get_keypath(const set_obj &sol, std::vector<int> &path_edges, std::vector<int> &kp_ends){
	path_edges.clear();
	kp_ends.clear();

	int start_edge = sol.getRandomElement(rng);

	std::vector<int> visited(probModel->prob_graph.getNumEdges(),0);

	std::vector<int> edge_stack;
	edge_stack.push_back(start_edge);

	while (!edge_stack.empty()){
		int edge_idx = edge_stack.back();
		edge_stack.pop_back();
		Edge *edge = probModel->prob_graph.getEdge(edge_idx);
		path_edges.push_back(edge_idx);
		visited[edge_idx] = 1;

		std::vector<int> src_tgt;
		src_tgt.push_back(edge->getSrcID());
		src_tgt.push_back(edge->getTgtID());

		for (int i = 0; i < src_tgt.size(); ++i){
			Node *node = probModel->prob_graph.getNode(src_tgt[i]);
			if (!node->isKey()){
				std::vector<Edge*> edges;
				const int n_edges = node->getEdges(edges);
				for (int e = 0; e < n_edges; ++e){
					if (sol.is_element(edges[e]->getID()) && (visited[edges[e]->getID()] == 0)){
						edge_stack.push_back(edges[e]->getID());
						break;
					}
				}
			}
			else{
				kp_ends.push_back(node->getID());
			}
		}
	}
}

LocalSearch::~LocalSearch(){
}