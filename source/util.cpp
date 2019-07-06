// utility functions
#include <util.h>

void floyd_warshall(const U_Graph &graph, std::vector<std::vector<int> > &dist, std::vector<std::vector<int> > &next){
	

	const int n_nodes = graph.getNumNodes();
	const int n_edges = graph.getNumEdges();


	dist = std::vector<std::vector<int> > (n_nodes, std::vector<int>(n_nodes,INF));
	next = std::vector<std::vector<int> > (n_nodes, std::vector<int>(n_nodes,-1));

	for (int i = 0; i < n_nodes; ++i){
		dist[i][i] = 0;
	}

	for (int e = 0; e < n_edges; ++e){
		const Edge *edge = graph.getEdge(e);
		int src = edge->getSrc()->getID();
		int tgt = edge->getTgt()->getID();
		if (src >= n_nodes || tgt >= n_nodes){
			PE("error found in ("<<src<<", "<<tgt<<")")
		}
	}


	for (int e = 0; e < n_edges; ++e){
		const Edge *edge = graph.getEdge(e);
		int src = edge->getSrc()->getID();
		int tgt = edge->getTgt()->getID();
		// PF("setting weight for ("<<src<<", "<<tgt<<").. ")
		dist[src][tgt] = edge->getWt();
		dist[tgt][src] = edge->getWt();
		// PE("done!")
		next[src][tgt] = tgt;
		next[tgt][src] = src;
	}

	for (int k = 0; k < n_nodes; ++k){
		for (int i = 0; i < n_nodes; ++i){
			for (int j = 0; j < n_nodes; ++j){
				if (dist[i][j] > dist[i][k] + dist[k][j]){
					dist[i][j] = dist[i][k] + dist[k][j];
					next[i][j] = next[i][k];
				}
			}
		}
	}

	// PE("")
	// for (int i = 0; i < dist.size(); ++i){
	// 	for (int j = 0; j < dist.size(); ++j){
	// 		if (dist[i][j] < 10){
	// 			PF(dist[i][j] << "   ")
	// 			continue;
	// 		}
	// 		if (dist[i][j] < 100){
	// 			PF(dist[i][j] << "  ")
	// 			continue;
	// 		}
	// 		if (dist[i][j] < 1000){
	// 			PF(dist[i][j] << " ")
	// 		}
	// 	}
	// 	PE("")
	// }
}

void get_mst(const U_Graph &graph, set_obj &sol){

	std::vector<int> edge_set = sol.getSet();

	std::vector<std::pair <int, int> > sorted_edges;

	for (int i = 0; i < edge_set.size(); ++i){
		sorted_edges.push_back(std::make_pair(graph.getEdge(edge_set[i])->getWt(),edge_set[i]));
	}

	// sort edge set in descending weight order so that can pop from back
	std::sort(sorted_edges.rbegin(), sorted_edges.rend());

	// vector to keep track of covered vertices (no cycles)
	std::vector<int> covered(graph.getNumNodes(),-1);

	std::vector<set_obj> subsets;

	sol.clear();

	while(!sorted_edges.empty()){
		int edge_idx = sorted_edges.back().second;
		sorted_edges.pop_back();
		Edge *edge = graph.getEdge(edge_idx);
		std::vector<int> src_tgt;
		src_tgt.push_back(edge->getSrcID());
		src_tgt.push_back(edge->getTgtID());

		std::vector<int> src_tgt_set;

		for (int i = 0; i < src_tgt.size(); ++i){
			bool found = false;
			// search subsets to find node
			for (int s = 0; s<subsets.size();++s){
				if (subsets[s].is_element(src_tgt[i])){
					src_tgt_set.push_back(s);
					found = true;
					break;
				}
			}
			// if not found in current subsets add new subset
			if (!found){
				subsets.push_back(set_obj(graph.getNumNodes()));
				subsets[subsets.size()-1].addElement(src_tgt[i]);
				src_tgt_set.push_back(subsets.size()-1);
			}
		}

		// PF("sets: " << src_tgt_set[0] <<  ", " << src_tgt_set[1])

		// test if cycle created
		if (src_tgt_set[0] != src_tgt_set[1]){
			// PE(" -> adding edge!")
			sol.addElement(edge_idx);
			// union two subsets
			subsets[src_tgt_set[0]].makeUnion(subsets[src_tgt_set[1]]);
			// remove 2nd subset
			subsets[src_tgt_set[1]] = subsets.back();
			subsets.pop_back();
		}
		else{
			// PE(" -> skipping edge!")
		}
	}
}

void get_components(const U_Graph &graph, const set_obj &x,const set_obj &y, std::vector<set_obj> &comp_edges, std::vector<set_obj> &comp_nodes, std::vector<int> &comp_type){


	// clear components vector
	comp_edges.clear();
	comp_nodes.clear();
	comp_type.clear();

	// set_obj for visited vertices
	set_obj visited(graph.getNumNodes());

	// visited.fill();

	// set_obj y(graph.getNumNodes());
	// for (int i = 0; i < x.size(); ++i){
	// 	Edge *e = graph.getEdge(x.get(i));
	// 	y.addElement(e->getSrcID());
	// 	y.addElement(e->getSrcID());
	// }

	// find all nodes induced by the edge set
	for (int i = 0; i < x.size(); ++i){
		Edge *e = graph.getEdge(x.get(i));
		int src = e->getSrcID();
		int tgt = e->getTgtID();

		set_obj temp_comp_e(graph.getNumEdges());
		set_obj temp_comp_n(graph.getNumNodes());

		// // if edge doesnt have both endpoints, add as separate component
		if (!y.is_element(src) || !y.is_element(tgt)){
			temp_comp_e.addElement(x.get(i));
			temp_comp_n.addElement(src);
			temp_comp_n.addElement(tgt);
			comp_nodes.push_back(temp_comp_n);
			comp_edges.push_back(temp_comp_e);
			comp_type.push_back(graph.getNode(src)->isTerm() || graph.getNode(tgt)->isTerm());
		}

		visited.addElement(src);
		visited.addElement(tgt);
	}

		while (!visited.empty()){
		std::vector<int> stack;
		// check if vertex is in solution
		// if (y.is_element(visited.get(0))){
		stack.push_back(visited.get(0));
		// }
		// else{
		// 	visited.removeElement(visited.get(0));
		// 	continue;
		// }
		
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

		// if no cycles skip component
		// if (!cycle){
		// 	continue;
		// }


		// prune the subgraph to remove all but the cycle
		prune_subgraph_noterm(graph,temp_x);

		for (int e=0; e < temp_x.size(); ++e){
			if (temp_x.get(e) < 0){
				CHECK(temp_x.get(e))
				std::cin.get();
			}
		}

		// reset temp_y to store nodes of pruned subgraph
		temp_y.clear();

// TODO: check if any src or tgt on any edges are -1.. otherwise find what the -1 is!

		bool term = false;
		// add all nodes from new subgraph
		for (int e=0; e < temp_x.size(); ++e){
	// PF("setting up.. " << e << " from size " << temp_x.size() << "/" << temp_x.set_data.size() << ": " << temp_x.get(e) << ".. ")
			Edge * edge = graph.getEdge(temp_x.get(e));
	// PE("done!")
			temp_y.addElement(edge->getSrcID());
			temp_y.addElement(edge->getTgtID());
			// check if there is a terminal in the component
			if (edge->getSrc()->isTerm() || edge->getTgt()->isTerm()){
				term = true;
			}
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

		// CHECK(temp_x.size())
		// CHECK(temp_y.size())
		if (temp_x.size() > 0){
			comp_nodes.push_back(temp_y);
			comp_edges.push_back(temp_x);
			comp_type.push_back(term);
		}

	} // end visit while
}

// removes "dangling" edges from subgraph, regardless of terminal status
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

	std::vector<int> edges_to_remove;
	std::vector<int> edges_to_skip(subgraph.size(),0);

	// loop while there are still edges to remove
	while (subgraph_size - subgraph.size() != 0){
		subgraph_size = subgraph.size();

		// delete edges with endpoint with only one degree
		for (int e = 0; e < subgraph.size(); ++e){
			if (edges_to_skip[e]){
				continue;
			}
			int src = graph.getEdge(subgraph.get(e))->getSrcID();
			int tgt = graph.getEdge(subgraph.get(e))->getTgtID();
			if (degree[src] < 2 || degree[tgt] < 2){
				int edge = subgraph.get(e);
				edges_to_remove.push_back(edge);
				edges_to_skip[e]++;
				// subgraph.removeElement(subgraph.get(e));
				degree[src]--;
				degree[tgt]--;
				break;
			}
		}
	}
	// remove edges
	subgraph.removeElements(edges_to_remove);
}

// removes "dangling" edges that arent terminals
void prune_subgraph(const U_Graph &graph, set_obj &subgraph){

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
	while (subgraph_size - subgraph.size() != 0 && !subgraph.empty()){
		subgraph_size = subgraph.size();

		// delete edges with endpoint with only one degree that is not a terminal
		for (int e = 0; e < subgraph.size(); ++e){
			int src = graph.getEdge(subgraph.get(e))->getSrcID();
			int tgt = graph.getEdge(subgraph.get(e))->getTgtID();
			if ((degree[src] < 2 && !graph.getNode(src)->isTerm()) || (degree[tgt] < 2 && !graph.getNode(tgt)->isTerm())){
				subgraph.removeElement(subgraph.get(e));
				degree[src]--;
				degree[tgt]--;
				break;
			}
		}
	}
}

bool verify(const U_Graph &graph, const std::vector<int> &terms, const set_obj &sol){

	if (sol.empty()){
		PE("SOLUTION EMPTY!")
		return false;
	}

	int path_cost =0;
	std::vector<int> v_visited(graph.getNumNodes(),0);
	std::vector<int> e_visited(graph.getNumEdges(),0);
	std::vector<int> edge_stack;
	edge_stack.push_back(sol.get(0));

	while (!edge_stack.empty()){
		int edge_idx = edge_stack.back();
		edge_stack.pop_back();
		e_visited[edge_idx] = 1;
		Edge *edge = graph.getEdge(edge_idx);
		path_cost += edge->getWt();
		std::vector<int> endpoints;		
		endpoints.push_back(edge->getSrcID());
		endpoints.push_back(edge->getTgtID());

		for (int i = 0; i < endpoints.size(); ++i){
			Node *node = graph.getNode(endpoints[i]);
			if (v_visited[endpoints[i]]){
				continue;
			}

			std::vector<Edge*> adj_edges; 
			node->getEdges(adj_edges);
			for (int e = 0; e < adj_edges.size();++e){
				int adj_idx = adj_edges[e]->getID(); 
				if (!sol.is_element(adj_idx) || e_visited[adj_idx]){
					continue;
				}
				edge_stack.push_back(adj_idx);
			}
			v_visited[endpoints[i]] = 1;
		}
	}

	VE("Total solution cost: " << path_cost)

	std::vector<int> missed_terms; 
	for (int t = 0; t< terms.size(); ++t){
		if (!v_visited[terms[t]]){
			missed_terms.push_back(terms[t]);
		}
	}

	if (missed_terms.empty()){
		VE("\nAll terminals visited!")
	}
	else{
		PE("Error! Terminals missed: ")
		for (int t=0;t<missed_terms.size();++t){
			PF(missed_terms[t] << " ")
		}
		PE("")
		return false;
	}

	VE("\n*****************************************")
    VE("tested solution:")
    VE("*****************************************")

    std::ostringstream latex_x;
	std::ostringstream latex_y;

    int cost = 0;

    VF("x = { ")

    set_obj test_y(graph.getNumNodes());

    for (int e = 0; e < sol.size(); ++e){
      Edge * edge = graph.getEdge(sol.get(e));
      
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
    VE("cost: " << cost << std::endl)
    VE("\nx = " <<latex_x.str())
    VE("y = " << latex_y.str()<<std::endl)
	return true;
}


int get_sol_value(const U_Graph &graph, const set_obj &sol){
	std::vector<int> edges = sol.getSet();

	int total_wt = 0;

	for (int e = 0; e < edges.size(); ++e){
		total_wt += graph.getEdge(edges[e])->getWt();
	}

	return total_wt;
}

void print_tree(const U_Graph &graph, const set_obj &sol, const std::string &fname){
	std::ofstream ofile;
  	ofile.open (fname);
  	for (int i = 0; i < sol.size(); ++i){
  		ofile << graph.getEdge(sol.get(i))->getSrcID() << " " << graph.getEdge(sol.get(i))->getTgtID() << std::endl;
  	}
  	ofile.close();
}