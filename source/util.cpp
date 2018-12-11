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
		dist[src][tgt] = edge->getWt();
		dist[tgt][src] = edge->getWt();
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

void get_components(const U_Graph &graph, const set_obj &x, const set_obj &y, std::vector<set_obj> &components){

	// clear components vector
	components.clear();

	// filled set_obj for visited vertices
	set_obj visited.fill(y.idx_size());

	while (!visited.empty()){
		set_obj temp_x(x.idx_size());
		set_obj temp_y(y.idx_size());
	}

}

// removes "dangling" edges
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
		VF("Error! Terminals missed: ")
		for (int t=0;t<missed_terms.size();++t){
			VF(missed_terms[t] << " ")
		}
		VE("")
		return false;
	}

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