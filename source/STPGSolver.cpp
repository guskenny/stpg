#include <STPGSolver.h>

STPGSolver::STPGSolver(int argc, const char **argv){

	sh = SettingsHandler();
	probModel = new STPGModel(argv[argc-1],sh);

}

void STPGSolver::solve(){

	time_t now = time(0);
	tm *ltm = localtime(&now);

	int size_count = 100;
	std::vector<int> pop_sizes={1500};
	// while(size_count < 5000){
	// 	pop_sizes.push_back(size_count);
	// 	size_count += 100;
	// }
	// pop_sizes.push_back(1500);

	int num_groups = 0;

	for(int _run = 0; _run < sh.FULL_RUNS*pop_sizes.size(); ++_run){

	PE("\n++++++++++++++++++++++++++++++++")
	PE("     STARTING RUN "<<_run+1)
	PE("++++++++++++++++++++++++++++++++\n")

	bool end_early = false;

	std::vector<int> inits;
	set_obj init_sol = set_obj(probModel->n_edges());
	int init_val = 0;

	// for(int j = 0; j < 20; ++j){
	bool sol_valid = false;
	// set_obj _init_sol = set_obj(probModel->n_edges());
	// init_val=0;	
	while (!sol_valid){

		switch (sh.INIT_SOL_TYPE) {
        case 0: {PE("HEURISTIC!")getInitSol(init_sol);break;}
        case 1: {PE("RANDOM!")getRandSol(init_sol);break;}
        default: {
			PE("getting worst feasible solution")
			for (int i = 0; i < probModel->n_edges(); ++i){
				init_sol.addElement(i);
			}
			get_mst(probModel->prob_graph, init_sol);

			// prune tree of any dangling vertices
			// prune_subgraph(probModel->prob_graph, init_sol);	
			}
    	}

		// if (sh.INIT_SOL_TYPE==0){
		// 	getInitSol(init_sol);
		// }
		// else if(sh.INIT_SOL_TYPE==1){
		// 	getRandSol(init_sol);
		// }
		// else{
		// 	PE("getting worst feasible solution")
		// 	for (int i = 0; i < probModel->n_edges(); ++i){
		// 		init_sol.addElement(i);
		// 	}
		// 	get_mst(probModel->prob_graph, init_sol);

		// 	// prune tree of any dangling vertices
		// 	prune_subgraph(probModel->prob_graph, init_sol);	
		// }

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
	// inits.push_back(init_val);
	// }

	// std::ofstream _inits;
	// _inits.open("inits_log.csv",std::ios_base::app);
	// _inits << probModel->getName() ;
	// for (int i=0; i < inits.size(); ++i){
	// 	_inits << "," << inits[i];
	// }
	// _inits << std::endl;
	// _inits.close();
	// return;

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

	set_obj best_full_sol(probModel->n_edges() + probModel->n_nodes());
	for (int i = 0; i < init_sol.size(); ++i){
		best_full_sol.addElement(init_sol[i]);
		Edge *edge = probModel->prob_graph.getEdge(init_sol[i]);
		// add end nodes for each edge, index for nodes is |E| + node index
		best_full_sol.addElement(probModel->n_edges() + edge->getSrcID());
		best_full_sol.addElement(probModel->n_edges() + edge->getTgtID());
	}

	int iter;

	LocalSearch ls(sh,probModel);

	for (iter = 0; iter < sh.NUM_ITER; ++iter){
		CHECK(iter)
		full_sols.clear();
		init_val = get_sol_value(probModel->prob_graph, best_sol);
		// for (int move = 0; move < sh.NUM_MOVES; ++move){
		for (int move = 0; move < pop_sizes[_run%pop_sizes.size()]; ++move){
			PF("\riteration: " << move << " of " << sh.NUM_MOVES)
			set_obj sol = best_sol;
			
			// int curr_val = ls.singleSwaps(sol,obj_vals);
			int curr_val = ls.doSwaps(sol,obj_vals);
			
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
			
			if (best_val <= probModel->getOptimal() && best_val > 0){
				// best_val = probModel->getOptimal();
				// std::cout << "*** FOUND OPTIMAL SOLUTION - ENDING SEARCH EARLY ***" << std::endl;
				end_early = true;
				// break;
			}
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

		if (sh.RAND_SPLIT){
			for (int i = 0; i < sh.NUM_SPLITS; ++i){
				set_obj bitstring(probModel->n_edges() + probModel->n_nodes());

				for (int j = 0; j < bitstring.size();++j){
					if (rand()%2) bitstring.addElement(j);
				}

				full_sols.push_back(bitstring);
			}
		}
		SolutionMerger sm = SolutionMerger(sh);

		merge_sol merged_sols(probModel->n_edges(), probModel->n_nodes());

		if (sh.SOLVE_MIP){
			sm.merge(full_sols, merged_sols);
		}

		num_groups = merged_sols.size();
		CHECK(num_groups)

		// std::ofstream _plot_grp;
		// _plot_grp.open(probModel->getName()+"_pop_plot.dat",std::ios_base::app);
		// _plot_grp << pop_sizes[_run%pop_sizes.size()] << "," << num_groups << std::endl;
		// _plot_grp.close();
		// continue;

		if (!sh.RAND_SPLIT){
			if (sh.SPLIT_FACTOR > 0){
				splitFactor(merged_sols);
			}
			else{
				split(merged_sols);
			}
		}

		set_obj post_merge_sol(probModel->n_edges());

		PE("optimal: " << probModel->getOptimal())
		PE("instance " << probModel->getName() << std::endl)

		if (sh.SOLVE_MIP){
			STPGMergeMip MIPsolver(sh,probModel,merged_sols,best_full_sol);
			MIPsolver.solve(post_merge_sol);
		}
		else{
			post_merge_sol = best_sol;
		}

		best_full_sol.clear();

		for (int i = 0; i < post_merge_sol.size(); ++i){
				best_full_sol.addElement(post_merge_sol[i]);
				Edge *edge = probModel->prob_graph.getEdge(post_merge_sol[i]);
				// add end nodes for each edge, index for nodes is |E| + node index
				best_full_sol.addElement(probModel->n_edges() + edge->getSrcID());
				best_full_sol.addElement(probModel->n_edges() + edge->getTgtID());
		}

		PE("\nbefore merge: " << get_sol_value(probModel->prob_graph, best_sol)) 
		PE("after merge: " << get_sol_value(probModel->prob_graph, post_merge_sol))
		PE("optimal: " << probModel->getOptimal())
		PE("instance " << probModel->getName() << std::endl)


		// if (get_sol_value(probModel->prob_graph, post_merge_sol) < get_sol_value(probModel->prob_graph, best_sol)){
		// 	std::cin.get();
		// }
		
		best_sol = post_merge_sol;
		best_val = get_sol_value(probModel->prob_graph, best_sol);

		if (best_val <= probModel->getOptimal() && best_val > 0){
			// best_val = probModel->getOptimal();
			std::cout << "*** FOUND OPTIMAL SOLUTION - ENDING SEARCH EARLY ***" << std::endl;
			end_early = true;
			break;
		}
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

	std::cout << "Instance " << probModel->getName() << std::endl;
	std::cout << "Final objective: " << best_val << std::endl;
	std::cout << "Total CPU time taken: " << solve_timer.elapsedSeconds() << std::endl;

	// time_t now = time(0);
	// tm *ltm = localtime(&now);

	std::string year = std::to_string(1900+ltm->tm_year);
	std::string month = std::to_string(1+ltm->tm_mon);
	std::string day = std::to_string(ltm->tm_mday);
	std::string hour = std::to_string(ltm->tm_hour);
	std::string minute = std::to_string(ltm->tm_min);
	std::string second = std::to_string(ltm->tm_sec);

	if (ltm->tm_year < 10)
		year = "0" + std::to_string(1+ltm->tm_year);

	if (ltm->tm_mon < 10)
		month = "0" + std::to_string(1+ltm->tm_mon);

	if (ltm->tm_mday < 10)
		day = "0" + std::to_string(ltm->tm_mday);

	if (ltm->tm_sec < 10)
		second = "0" + std::to_string(ltm->tm_sec);

	if (ltm->tm_min < 10)
		minute = "0" + std::to_string(ltm->tm_min);

	if (ltm->tm_hour < 10)
		hour = "0" + std::to_string(ltm->tm_hour);
	
	std::ofstream _log;
	_log.open("results_log.csv",std::ios_base::app);
	_log << "," <<year << "-" << month << "-" << day << ",";
	_log << hour << ":" << minute << ":" << second << ",";
	_log << probModel->getName() << "," << best_val << "," << probModel->getOptimal() << ",";
	if (best_val <= probModel->getOptimal()){
		_log << "*" << ",";
	}
	else
	{
		_log << ",";
	}
	_log << solve_timer.elapsedWallTime() << "," << solve_timer.elapsedSeconds() << ",";
	// _log << iter << "," << sh.NUM_ITER << "," << sh.NUM_MOVES << ",";
	_log << iter << "," << sh.NUM_ITER << "," << pop_sizes[_run%pop_sizes.size()] << ",";
	_log << sh.NUM_SPLITS << "," << sh.SPLIT_FACTOR << "," << sh.PREPROCESS << "," << sh.HEURISTIC_SEARCH << "," << sh.MIP_TIME << "," << sh.NUM_SWAPS << "," << sh.NUM_INIT_JUMPS << "," << sh.NUM_JUMPS << "," << sh.JUMP_FREQ << "," << init_val << "," << probModel->n_vars() << "," << num_groups << std::endl;
	_log.close();



	}//end full run
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

void STPGSolver::DFSsplit(const merge_sol &merged_sols, const int group, std::vector<int> &connected_vars, std::mt19937 &rng){

	connected_vars.clear();

	int start_var = merged_sols.groups[group].getRandomElement(rng);

	std::vector<int> visited(probModel->n_edges() + probModel->n_nodes(),0);

	// CHECK(group)

	std::vector<int> var_stack;
	var_stack.push_back(start_var);
	connected_vars.push_back(start_var);

	visited[start_var]=1;

	while (!var_stack.empty()){
		int var_idx = var_stack.back();
		var_stack.pop_back();

		// if v is edge find nodes
      	if (var_idx < probModel->n_edges()){
      		Edge *edge = probModel->prob_graph.getEdge(var_idx);
      		std::vector<int> src_tgt;
			src_tgt.push_back(edge->getSrcID());
			src_tgt.push_back(edge->getTgtID());
			for (int i = 0; i < src_tgt.size(); ++i){
				Node *node = probModel->prob_graph.getNode(src_tgt[i]);
				// if node not in current group then continue
				if ((visited[src_tgt[i]+probModel->n_edges()] == 0) && merged_sols.groups[group].is_element(src_tgt[i]+probModel->n_edges())){
					visited[src_tgt[i]+probModel->n_edges()] = 1;
						// group == merged_sols.group_map[src_tgt[i]+probModel->n_edges()]){
					// PE("adding node: " << src_tgt[i]+probModel->n_edges())
					var_stack.push_back(src_tgt[i]+probModel->n_edges());
					connected_vars.push_back(src_tgt[i]+probModel->n_edges());
				}
				else{
					visited[src_tgt[i]+probModel->n_edges()] = 1;
				}
			}
      	}
      	// if v is node find edges
      	else{
      		Node *node = probModel->prob_graph.getNode(var_idx - probModel->n_edges());
      		std::vector<Edge*> edges;
      		node->getEdges(edges);
			for (int edge_ptr_idx = 0; edge_ptr_idx < edges.size(); ++edge_ptr_idx){
				int edge_idx = edges[edge_ptr_idx]->getID();
				if ((visited[edge_idx] == 0) && merged_sols.groups[group].is_element(edge_idx)){
					visited[edge_idx] = 1;
					// group == merged_sols.group_map[edge_idx]){
					// PE("adding edge: " << edge_idx)
					var_stack.push_back(edge_idx);
					connected_vars.push_back(edge_idx);
				}
				else{
					visited[edge_idx] = 1;
				}
	      	}
		}
	}

	// std::vector<int> test_seen(probModel->n_edges() + probModel->n_nodes(),0);

	// for (int i = 0; i < connected_vars.size(); ++i){
	// 	if (test_seen[connected_vars[i]]){
	// 		PE("duplicate! " << connected_vars[i])
	// 	}
	// 	else{
	// 		// PE("not duplicate") 
	// 		test_seen[connected_vars[i]]=1;
	// 	}
	// }

}

void STPGSolver::splitFactor(merge_sol &merged_sols){
	// set fixed group number (no point in splitting already split groups)
	// for number of splits (default = all groups){
	// 		choose random group
	// 		choose random variable (if size of group = 2 then just split it in half)
	// 		get all variables connected to variable that are in the same group
	//		if size of connected variables = size of group, then remove last 
	//			variable from array of connected variables
	// 		remove connected variables from original group
	// 		add to new group
	// }
	// check groups for 0 element groups

	std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 rng = std::mt19937(seed);

    int num_vars = probModel->n_edges() + probModel->n_nodes();

    set_obj temp_group(num_vars);
      	
  	std::vector<int> connected_vars;

    int fixed_group_size = merged_sols.size();

    PE(fixed_group_size << " groups before split")
    PF("splitting groups by factor " << sh.SPLIT_FACTOR << " or adding " << sh.NUM_SPLITS << " groups..")

    int tries = 0;

    set_obj splittable(num_vars);
    for (int i = 0; i < merged_sols.size(); ++i){
    	if (merged_sols.groups.size() > 1){
    		splittable.addElement(i);
    	}
    }

	while ((merged_sols.size() < double(num_vars)/sh.SPLIT_FACTOR) || (((merged_sols.size()-fixed_group_size) < sh.NUM_SPLITS) && (merged_sols.size() < num_vars))){
	// while ((merged_sols.size() < double(num_vars)/sh.SPLIT_FACTOR) || ((merged_sols.size() < sh.NUM_SPLITS) && (merged_sols.size() < num_vars))){
		// choose random group

		int g = splittable.getRandomElement(rng);

   //    	// if group is size 0 or 1 continue
   //    	if (merged_sols.groups[g].size() == 0){
   // //    		tries++;
			// // if (tries < 100000){
	  // //     		break;
   // //    		}
   //    		// PE("ZERO GROUP!!")
   //    		continue;
   //    	}
   //    	if (merged_sols.groups[g].size() == 1){
   //    		// decrement s to try again
   // //    		tries++;
			// // if (tries < 100000){
	  // //     		break;
   // //    		}
   //    		continue;
   //    	}

      	temp_group.clear();
      	connected_vars.clear();

		tries = 0;
  		
  		// recursively do DFS until all connected vars in group are found
      	DFSsplit(merged_sols,g,connected_vars,rng);

      	// PE("before - orig: " << merged_sols.groups[g].size())
      	// CHECK(connected_vars.size())
      	if (connected_vars.size() == 0){
      		CHECK(connected_vars.size())
      		CHECK(merged_sols.groups[g].size())
      	}

      	// if group is size 2 split one variable and continue
      	if (connected_vars.size() == 2){
      		merged_sols.group_map[connected_vars[0]] = merged_sols.groups.size();
			temp_group.addElement(connected_vars[0]);
			merged_sols.groups[g].removeElement(connected_vars[0]);
			merged_sols.groups.push_back(temp_group);
			splittable.removeElement(g);
      		continue;
      	}
      	// CHECK(connected_vars.size())
      	// if all variables from group were added, remove last variable in set so that a split is maintained
      	if (connected_vars.size() == merged_sols.groups[g].size()){
      		connected_vars.pop_back();
      	}

       	// add group to back of list and remove from current group - update group map
      	for (std::vector<int>::iterator var = connected_vars.begin(); var != connected_vars.end(); ++var){
      		temp_group.addElement(*var);
      		merged_sols.group_map[*var] = merged_sols.groups.size();
      		merged_sols.groups[g].removeElement(*var);
      	}
      	// CHECK(temp_group.size())

      	if(temp_group.size() > 1){
      		splittable.addElement(merged_sols.groups.size());
      	}
      	if(merged_sols.groups[g].size() < 2){
      		splittable.removeElement(g);
      	}

      	merged_sols.groups.push_back(temp_group);

      	// CHECK(merged_sols.groups[merged_sols.groups.size()-1].size())
      	// PE("after - orig: " << merged_sols.groups[g].size() << ", new: " << merged_sols.groups[merged_sols.groups.size()-1].size() << ", group: " << merged_sols.groups.size()-1)
		
	}

	PE(" done!")
	PE(merged_sols.groups.size() << " groups after split")
	int big_group_count = 0;
	int zero_group_count = 0;
	int one_group_count = 0;
	for (int i = 0; i < merged_sols.groups.size(); ++i){
		if (merged_sols.groups[i].size() > 1){
			big_group_count++;
			continue;
		}
		if (merged_sols.groups[i].size() == 0){
			zero_group_count++;
			continue;
		}
		if (merged_sols.groups[i].size() == 1){
			one_group_count++;
			continue;
		}
	}
	CHECK(zero_group_count)
	CHECK(one_group_count)
	CHECK(big_group_count)
}

void STPGSolver::split(merge_sol &merged_sols){
	// set fixed group number (no point in splitting already split groups)
	// for number of splits (default = all groups){
	// 		choose random group
	// 		choose random variable (if size of group = 2 then just split it in half)
	// 		get all variables connected to variable that are in the same group
	//		if size of connected variables = size of group, then remove last 
	//			variable from array of connected variables
	// 		remove connected variables from original group
	// 		add to new group
	// }
	// check groups for 0 element groups

	std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 rng = std::mt19937(seed);

    set_obj temp_group(probModel->n_edges() + probModel->n_nodes());
      	
  	std::vector<int> connected_vars;

    int fixed_group_size = merged_sols.size();

    PE(fixed_group_size << " groups before split")
    PF("splitting groups " << sh.NUM_SPLITS << " times..")

    int tries = 0;

	for (int s = 0; s < sh.NUM_SPLITS; ++s){
		// choose random group
		std::uniform_int_distribution<int> uni(0,merged_sols.size()-1);
      	int g = uni(rng);
      	// if group is size 0 or 1 continue
      	if (merged_sols.groups[g].size() == 0){
      		tries++;
			if (tries < 10000){
	      		s--;
      		}
      		// PE("ZERO GROUP!!")
      		continue;
      	}
      	if (merged_sols.groups[g].size() == 1){
      		// decrement s to try again
      		tries++;
			if (tries < 10000){
	      		s--;
      		}
      		continue;
      	}

      	temp_group.clear();
      	connected_vars.clear();

		tries = 0;
  		
  		// recursively do DFS until all connected vars in group are found
      	DFSsplit(merged_sols,g,connected_vars,rng);

      	// PE("before - orig: " << merged_sols.groups[g].size())
      	// CHECK(connected_vars.size())
      	if (connected_vars.size() == 0){
      		CHECK(connected_vars.size())
      		CHECK(merged_sols.groups[g].size())
      	}

      	// if group is size 2 split one variable and continue
      	if (connected_vars.size() == 2){
      		merged_sols.group_map[connected_vars[0]] = merged_sols.groups.size();
			temp_group.addElement(connected_vars[0]);
			merged_sols.groups[g].removeElement(connected_vars[0]);
			merged_sols.groups.push_back(temp_group);
      		continue;
      	}
      	// CHECK(connected_vars.size())
      	// if all variables from group were added, remove last variable in set so that a split is maintained
      	if (connected_vars.size() == merged_sols.groups[g].size()){
      		connected_vars.pop_back();
      	}

       	// add group to back of list and remove from current group - update group map
      	for (std::vector<int>::iterator var = connected_vars.begin(); var != connected_vars.end(); ++var){
      		temp_group.addElement(*var);
      		merged_sols.group_map[*var] = merged_sols.groups.size();
      		merged_sols.groups[g].removeElement(*var);
      	}
      	// CHECK(temp_group.size())
      	merged_sols.groups.push_back(temp_group);
      	// CHECK(merged_sols.groups[merged_sols.groups.size()-1].size())
      	// PE("after - orig: " << merged_sols.groups[g].size() << ", new: " << merged_sols.groups[merged_sols.groups.size()-1].size() << ", group: " << merged_sols.groups.size()-1)
		
	}

	PE(" done!")
	PE(merged_sols.groups.size() << " groups after split")
	int big_group_count = 0;
	int zero_group_count = 0;
	int one_group_count = 0;
	for (int i = 0; i < merged_sols.groups.size(); ++i){
		if (merged_sols.groups[i].size() > 1){
			big_group_count++;
			continue;
		}
		if (merged_sols.groups[i].size() == 0){
			zero_group_count++;
			continue;
		}
		if (merged_sols.groups[i].size() == 1){
			one_group_count++;
			continue;
		}
	}
	CHECK(zero_group_count)
	CHECK(one_group_count)
	CHECK(big_group_count)
}

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

	PF("getting heuristic solution.. ")
	// get terminal nodes
	std::vector<int> terms = probModel->getTerms();

	CHECK(terms.size())

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
	PE("done! ")
}