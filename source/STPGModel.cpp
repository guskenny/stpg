#include "STPGModel.h"

STPGModel::STPGModel(const char *argv, const SettingsHandler &_sh) : sh(_sh){

	PF("building graph..")

	std::string f_name(argv);

	std::size_t botDirPos = f_name.find_last_of("/");
	_fpath = f_name.substr(0, botDirPos+1);
	std::string file_name = f_name.substr(botDirPos+1,f_name.length());
	std::size_t dotPos = file_name.find_last_of(".");
	_name = file_name.substr(0,dotPos);

	std::size_t topDirPos = _fpath.find_first_of("/");
	_inst_grp = _fpath.substr(topDirPos+1, _fpath.length());

	std::ifstream in(f_name);

	if (!in){
		PE("\nFile: " << f_name <<" not found!")
		exit(0);
	}

	std::string in_line;

	int _n_terms=0;

	while ((std::getline(in, in_line))) {
		std::vector<std::string> split_str;
		boost::split(split_str, in_line, [](char c){return c == ' ';});



		if (split_str[LABEL] == "Nodes"){
			// std::cout << "resizing graph.." << std::flush;
			orig_graph.resizeNodes(stoi(split_str[N1]));
			_n_nodes = stoi(split_str[N1]);
			_n_full_nodes = stoi(split_str[N1]);
			// std::cout << " done!" << std::endl;
        	continue;
		}
		else if (split_str[LABEL] == "Edges"){
			// std::cout << "adding number of edges.." << std::flush;
			_n_edges = stoi(split_str[N1]);
			_n_full_edges = stoi(split_str[N1]);
			// std::cout << " done!" << std::endl;
        	continue;
		}
		else if (split_str[LABEL] == "Optimal"){
			_optimal = stoi(split_str[N1]);
        	continue;
		}
		else if (split_str[LABEL] == "Terminals"){
			// std::cout << "adding number of terminals.." << std::flush;
			_n_terms = stoi(split_str[N1]);
			// std::cout << " done!" << std::endl;
        	continue;
		}
		else if (split_str[LABEL] == "E"){
			// std::cout << "adding edge (" << stoi(split_str[N1])-1 <<","<< stoi(split_str[N2])-1 << "), wt: " << split_str[WT] << ".." <<std::flush;
			orig_graph.addEdge(stoi(split_str[N1])-1, stoi(split_str[N2])-1, stoi(split_str[WT]));
			// std::cout << " done!" << std::endl;
        	continue;
		}
		else if (split_str[LABEL] == "T"){
			// std::cout << "adding terminal.." << std::flush;
			terms.push_back(stoi(split_str[N1])-1);
			full_terms.push_back(stoi(split_str[N1])-1);
			orig_graph.getNode(stoi(split_str[N1])-1)->setTerm(true);
			orig_graph.getNode(stoi(split_str[N1])-1)->setKey(true);
			// std::cout << " done!" << std::endl;
        	continue;
		}
	}

	PE(" done!")
	PE("\n***** GRAPH STATS *****")
	PE("Nodes: " << _n_nodes*(_n_nodes == orig_graph.getNumNodes()) << ", Edges: " << _n_edges*(_n_edges == orig_graph.getNumEdges()) << ", Terminals: " << n_terms()*(_n_terms == n_terms()))

	qol::CpuTimer run_timer;

	if (sh.PREPROCESS){
		PE("\n***** PRE-PROCESSING GRAPH *****")
		pre_process();
	}
	else{
		prob_graph = orig_graph;
	}
	PF("running floyd warshall algorithm... ")
	floyd_warshall(prob_graph,fw_dist,fw_next);
	PE("done!")
	findKeys();

	PE("timer: " << run_timer.elapsedSeconds())
	// reconstruct();

}

void STPGModel::pre_process(){
	std::vector<std::vector<int> > d_dash;
	std::vector<std::vector<int> > n_dash;

	std::vector<EdgeID> del_edges;

	// save copy of original graph
	// U_Graph temp_graph = orig_graph;
	
	bool keep_removing = true;

	int nodes_to_remove = 0;

	// while we can still remove edges from graph
	while (keep_removing){
		keep_removing = false;

		floyd_warshall(orig_graph,d_dash,n_dash);
		
		int e_count = orig_graph.getNumEdges();


		// std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;

		// remove edges that have a shorter path otherwise
		PF1("removing edges that have shorter path otherwise... ")
		PE1(" ")
		for (int e = 0; e < e_count; ++e){
			const Edge *edge = orig_graph.getEdge(e);
			int src = edge->getSrc()->getID();
			int tgt = edge->getTgt()->getID();
			if (d_dash[src][tgt] < edge->getWt()){
				del_edges.push_back(e);
				keep_removing = true;
				// pp_moves.push_back(Move(src,tgt,edge->getWt()));
				PE1("removing edge (" << src << ", " << tgt << ", " << edge->getWt()<< ")")
			}
		}

		orig_graph.deleteEdges(del_edges);
		PE1("done!")

		// std::cout << "removed " << del_edges.size() << " edges" << std::endl;

		// // std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
		// std::cout << "Nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()*(prob_graph.getNumEdges()==(e_count-del_edges.size()))<< std::endl;

		del_edges.clear();
		// del_nodes.clear();

		// e_count = prob_graph.getNumEdges();

		PF1("removing nodes of degree 1... ")
		PE1(" ")
		// removing vertices of degree 1 or 2
		for (int n = 0; n < orig_graph.getNumNodes(); ++n){
			const Node *node = orig_graph.getNode(n);
			if (node->isTerm()){
				continue;
			}
			if (node->getDegree() == 1){
				PE1("node " << node->getID()+1 << " has degree 1")
				std::vector<Edge*> edges;
				node->getEdges(edges);
				del_edges.push_back(edges[0]->getID());
				keep_removing = true;
				nodes_to_remove++;
				PE1("removing edge (" << edges[0]->getSrc()->getID()+1 << ", " << edges[0]->getTgt()->getID()+1 << ", " << edges[0]->getWt() << ")")
				// pp_moves.push_back(Move(edges[0]->getSrc()->getID(),edges[0]->getTgt()->getID(), edges[0]->getWt()));
				// del_nodes.push_back(n);
			}
		}

		// std::cout << "-> removing " <<del_edges.size() << " edges... " << std::flush; 
		orig_graph.deleteEdges(del_edges);
		// std::cout << "-> removing " <<del_nodes.size() << " nodes... " << std::flush; 
		// prob_graph.deleteNodes(del_nodes);
		// for (int i = 0; i < del_nodes.size(); ++i){
		// 	PE1("removing node " << del_nodes[i])
		// 	pp_moves.push_back(Move(del_nodes[i],-1, -1));
		// }

		PE1("done!")

		// std::cout << "removed " << del_edges.size() << " edges" << std::endl;
		// std::cout << "removed " << del_nodes.size() << " nodes" << std::endl;

		// // std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
		// std::cout << "Nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()*(prob_graph.getNumEdges()==(e_count-del_edges.size()))<< std::endl;

		del_edges.clear();
		// del_nodes.clear();

		e_count = orig_graph.getNumEdges();

		////////////////////////////

		PF1("removing nodes of degree 2... ")
		PE1(" ")
		// removing vertices of degree 1 or 2
		for (int n = 0; n < orig_graph.getNumNodes(); ++n){
			const Node *node = orig_graph.getNode(n);
			if (node->isTerm()){
				continue;
			}
			if (node->getDegree() == 2){
				std::vector<Edge*> edges;
				node->getEdges(edges);
				bool bail = false;
				// if any of the nodes are terms, then bail!
				if (edges[0]->getSrc()->isTerm() || edges[0]->getTgt()->isTerm() || edges[1]->getSrc()->isTerm() ||edges[1]->getTgt()->isTerm()){
					bail = true;
				}
				if (bail){
					continue;
				}

				int src_0 = edges[0]->getSrc()->getID();
				int tgt_0 = edges[0]->getTgt()->getID(); 
				int wt_0 = edges[0]->getWt();
				int src_1 = edges[1]->getSrc()->getID();
				int tgt_1 = edges[1]->getTgt()->getID();
				int wt_1 = edges[1]->getWt();

				int n1 = (src_0 != n) ? src_0 : tgt_0;
				int n2 = (src_1 != n) ? src_1 : tgt_1;

				int new_src = std::min(n1,n2);
				int new_tgt = std::max(n1,n2);
				

				edges[0]->setSrc(orig_graph.getNode(new_src));
				edges[0]->setTgt(orig_graph.getNode(new_tgt));
				edges[0]->setWt(wt_0 + wt_1);

				// delete edges and nodes
				PE1("removing edge (" << src_0 << ", " << tgt_0 << ", " << wt_0 << ")")
				PE1("removing edge (" << src_1 << ", " << tgt_1 << ", " << wt_1 << ")")
				// pp_moves.push_back(Move(src_0,tgt_0,wt_0));
				// pp_moves.push_back(Move(src_1,tgt_1,wt_1));

				del_edges.push_back(edges[0]->getID());
				del_edges.push_back(edges[1]->getID());
				orig_graph.deleteEdges(del_edges);
				del_edges.clear();
				orig_graph.addEdge(new_src,new_tgt,wt_0 + wt_1);
				nodes_to_remove++;
				keep_removing = true;


				// prob_graph.addEdge(edges[0]->getID(),new_src,new_tgt,wt_0+wt_1);				
				// pp_moves.push_back(Move(new_src,new_tgt,-1));
				PE1("adding edge (" << new_src << ", " << new_tgt << ", " << wt_0+wt_1 << ")")
				// del_nodes.push_back(n);
				// n=0;
			}
		}

		// orig_graph.deleteEdges(del_edges);	

		/////////////////////////////


		prob_graph.deleteEdges(del_edges);

		// prob_graph.deleteNodes(del_nodes);		
		// for (int i = 0; i < del_nodes.size(); ++i){
		// 	pp_moves.push_back(Move(del_nodes[i],-1, -1));
		// }

		// std::cout << "done!" << std::endl;

		// std::cout << "removed " << del_edges.size() << " edges" << std::endl;
		// std::cout << "added " << del_edges.size()/2 << " edges" << std::endl;
		// std::cout << "removed " << del_nodes.size() << " nodes" << std::endl;

		// std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
		
		PE1("\nOld: nodes: " << _n_nodes << ", edges: " << _n_edges)
		PE1("New: nodes: " << orig_graph.getNumNodes() << ", Edges: " << orig_graph.getNumEdges()<< std::endl)


		del_edges.clear();
		// del_nodes.clear();

	}

	CHECK(nodes_to_remove)

	PE1("should end up with: "<<_n_full_nodes - nodes_to_remove)

	// set up node maps for pre-processed graphs
	node_pp_map = std::vector<NodeID>(_n_full_nodes,-1);
	pp_node_map.clear();

	terms.clear();

	// get all nodes that have more than 1 degree or are terminals and make maps
	for (int n = 0; n < _n_full_nodes; ++n){
		Node* node = orig_graph.getNode(n);
		if (node->getDegree() > 0){
			node_pp_map[n] = pp_node_map.size();
			pp_node_map.push_back(n);
		}
	}

	PE("testing maps:")
	for (int n = 0; n < _n_full_nodes; ++n){
		if (node_pp_map[n] == -1){
			continue;
		}
		if (n != pp_node_map[node_pp_map[n]]){
			PE("ERROR IN MAP!")
		}
	}

	// CHECK(temp_graph.getNumNodes())
	// CHECK(temp_graph.getNumEdges())
	CHECK(orig_graph.getNumNodes())
	CHECK(orig_graph.getNumEdges())

	// prob_graph = U_Graph(pp_node_map.size());
	prob_graph.resizeNodes(pp_node_map.size());

	CHECK(prob_graph.getNumNodes())
	CHECK(prob_graph.getNumEdges())
	int edge_count = 0;

	CHECK(pp_node_map.size())

	// reconstruct pre-processed graph
	for (int n = 0; n < pp_node_map.size(); ++n){
		// PE("node: " << n << " - " << pp_node_map[n])
		Node* orig_node = orig_graph.getNode(pp_node_map[n]);
		NodeID orig_id = orig_node->getID();
		// set terminal status
		if (orig_node->isTerm()){
			prob_graph.setTerm(n,true);
			// terms.push_back(n);
		}
		for (int e = 0; e < orig_node->getOutDegree(); ++e){
			int tgt = orig_node->getOutEdge(e)->getTgtID();
			int wt = orig_node->getOutEdge(e)->getWt();
			// PF1("adding edge: "<< new_graph.getNumEdges() << " = (" << orig_id+1 <<","<< tgt+1 <<","<<wt<<") = (" << n <<","<< node_pp_map[tgt] <<","<<wt<<").. ")
			prob_graph.addEdge(node_pp_map[orig_id],node_pp_map[tgt],wt);
			// PE1("done!")
		}
	}


	// prob_graph = new_graph;

	terms.clear();
	for (int n = 0; n < prob_graph.getNumNodes(); ++n){
		if (prob_graph.getNode(n)->isTerm()){
			terms.push_back(n);
		}
	}

	PE("testing terms:")
	for (int t = 0; t < full_terms.size(); ++t){
		bool found = false;
		for (int n = 0; n < terms.size(); ++n){
			if (full_terms[t] == pp_node_map[terms[n]]){
				found = true;
				break;
			}
		}
		if (!found){
			PE("term " << full_terms[t] << " not matched!")
		}
	}

	// PE("testing edges:")
	// for (int e = 0; e < prob_graph.getNumEdges(); ++e){
	// 	bool found = false;
	// 	Edge* temp_edge = prob_graph.getEdge(e);

	// 	if (!found){
	// 		PE("term " << terms[e] << " not matched!")
	// 	}
	// }


	// for (int i = prob_graph.getNumNodes()-1; i>-1; --i){
	// 	if (prob_graph.getNode(i)->getDegree() ==0){
	// 		del_nodes.push_back(i);
	// 		pp_moves.push_back(Move(i,-2, -2));
	// 	}
	// }
	// prob_graph.deleteNodes(del_nodes);		

	PE("\nOld: nodes: " << orig_graph.getNumNodes() << ", edges: " << orig_graph.getNumEdges())
	PE("New: nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()<< std::endl)

	// fw_dist.clear();
	// fw_next.clear();

	_n_nodes = prob_graph.getNumNodes();
	_n_edges = prob_graph.getNumEdges();

	_n_full_nodes = _n_nodes;
	_n_full_edges = _n_edges;
	_n_terms = terms.size();
	full_terms = terms;
}

// void STPGModel::pre_process(){
// 	std::vector<std::vector<int> > d_dash;
// 	std::vector<std::vector<int> > n_dash;

// 	// int e_count=prob_graph.getNumEdges()+1;
// 	std::vector<NodeID> del_nodes;
// 	std::vector<EdgeID> del_edges;

// 	// save copy of original graph
// 	orig_graph = prob_graph;
// 	_n_full_nodes = orig_graph.getNumNodes();
// 	_n_full_edges = orig_graph.getNumEdges();

// 	bool keep_removing = false;

// 	int nodes_to_remove = 0;

// 	// while we can still remove edges from graph
// 	while (keep_removing){
// 		keep_removing = false;

// 		floyd_warshall(prob_graph,d_dash,n_dash);
		
// 		int e_count = prob_graph.getNumEdges();


// 		// std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;

// 		// remove edges that have a shorter path otherwise
// 		PF1("removing edges that have shorter path otherwise... ")
// 		PE1(" ")
// 		for (int e = 0; e < e_count; ++e){
// 			const Edge *edge = prob_graph.getEdge(e);
// 			int src = edge->getSrc()->getID();
// 			int tgt = edge->getTgt()->getID();
// 			if (d_dash[src][tgt] < edge->getWt()){
// 				del_edges.push_back(e);
// 				keep_removing = true;
// 				// pp_moves.push_back(Move(src,tgt,edge->getWt()));
// 				PE1("removing edge (" << src << ", " << tgt << ", " << edge->getWt()<< ")")
// 			}
// 		}

// 		prob_graph.deleteEdges(del_edges);
// 		PE1("done!")

// 		std::cout << "removed " << del_edges.size() << " edges" << std::endl;

// 		// // std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
// 		// std::cout << "Nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()*(prob_graph.getNumEdges()==(e_count-del_edges.size()))<< std::endl;

// 		del_edges.clear();
// 		// del_nodes.clear();

// 		// e_count = prob_graph.getNumEdges();

// 		PF1("removing nodes of degree 1... ")
// 		PE1(" ")
// 		// removing vertices of degree 1 or 2
// 		for (int n = 0; n < prob_graph.getNumNodes(); ++n){
// 			const Node *node = prob_graph.getNode(n);
// 			if (node->isTerm()){
// 				continue;
// 			}
// 			if (node->getDegree() == 1){
// 				PE1("node " << node->getID()+1 << " has degree 1")
// 				std::vector<Edge*> edges;
// 				node->getEdges(edges);
// 				del_edges.push_back(edges[0]->getID());
// 				keep_removing = true;
// 				nodes_to_remove++;
// 				PE1("removing edge (" << edges[0]->getSrc()->getID()+1 << ", " << edges[0]->getTgt()->getID()+1 << ", " << edges[0]->getWt() << ")")
// 				// pp_moves.push_back(Move(edges[0]->getSrc()->getID(),edges[0]->getTgt()->getID(), edges[0]->getWt()));
// 				// del_nodes.push_back(n);
// 			}
// 		}

// 		// std::cout << "-> removing " <<del_edges.size() << " edges... " << std::flush; 
// 		prob_graph.deleteEdges(del_edges);
// 		// std::cout << "-> removing " <<del_nodes.size() << " nodes... " << std::flush; 
// 		// prob_graph.deleteNodes(del_nodes);
// 		// for (int i = 0; i < del_nodes.size(); ++i){
// 		// 	PE1("removing node " << del_nodes[i])
// 		// 	pp_moves.push_back(Move(del_nodes[i],-1, -1));
// 		// }

// 		PE1("done!")

// 		std::cout << "removed " << del_edges.size() << " edges" << std::endl;
// 		// std::cout << "removed " << del_nodes.size() << " nodes" << std::endl;

// 		// // std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
// 		// std::cout << "Nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()*(prob_graph.getNumEdges()==(e_count-del_edges.size()))<< std::endl;

// 		del_edges.clear();
// 		// del_nodes.clear();

// 		// e_count = prob_graph.getNumEdges();

// 		// PF1("removing nodes of degree 2... ")
// 		// PE1(" ")
// 		// // removing vertices of degree 1 or 2
// 		// for (int n = 0; n < prob_graph.getNumNodes(); ++n){
// 		// 	const Node *node = prob_graph.getNode(n);
// 		// 	if (node->isTerm()){
// 		// 		continue;
// 		// 	}
// 		// 	if (node->getDegree() == 2){
// 		// 		std::vector<Edge*> edges;
// 		// 		node->getEdges(edges);
// 		// 		int src_0 = edges[0]->getSrc()->getID();
// 		// 		int tgt_0 = edges[0]->getTgt()->getID(); 
// 		// 		int wt_0 = edges[0]->getWt();
// 		// 		int src_1 = edges[1]->getSrc()->getID();
// 		// 		int tgt_1 = edges[1]->getTgt()->getID();
// 		// 		int wt_1 = edges[1]->getWt();

// 		// 		int n1 = (src_0 != n) ? src_0 : tgt_0;
// 		// 		int n2 = (src_1 != n) ? src_1 : tgt_1;

// 		// 		int new_src = std::min(n1,n2);
// 		// 		int new_tgt = std::max(n1,n2);
				
// 		// 		edges[0]->setSrc(prob_graph.getNode(new_src));
// 		// 		edges[0]->setTgt(prob_graph.getNode(new_tgt));
// 		// 		edges[0]->setWt(wt_0 + wt_1);

// 		// 		// delete edges and nodes
// 		// 		PE1("removing edge (" << src_0 << ", " << tgt_0 << ", " << wt_0 << ")")
// 		// 		PE1("removing edge (" << src_1 << ", " << tgt_1 << ", " << wt_1 << ")")
// 		// 		// pp_moves.push_back(Move(src_0,tgt_0,wt_0));
// 		// 		// pp_moves.push_back(Move(src_1,tgt_1,wt_1));

// 		// 		// del_edges.push_back(edges[0]->getID());
// 		// 		del_edges.push_back(edges[1]->getID());


// 		// 		// prob_graph.addEdge(edges[0]->getID(),new_src,new_tgt,wt_0+wt_1);				
// 		// 		// pp_moves.push_back(Move(new_src,new_tgt,-1));
// 		// 		PE1("adding edge (" << new_src << ", " << new_tgt << ", " << wt_0+wt_1 << ")")
// 		// 		// del_nodes.push_back(n);
// 		// 		n=0;
// 		// 	}
// 		// }

// 		// prob_graph.deleteEdges(del_edges);	

// 		// prob_graph.deleteEdges(del_edges);

// 		// prob_graph.deleteNodes(del_nodes);		
// 		// for (int i = 0; i < del_nodes.size(); ++i){
// 		// 	pp_moves.push_back(Move(del_nodes[i],-1, -1));
// 		// }

// 		// std::cout << "done!" << std::endl;

// 		std::cout << "removed " << del_edges.size() << " edges" << std::endl;
// 		// std::cout << "added " << del_edges.size()/2 << " edges" << std::endl;
// 		// std::cout << "removed " << del_nodes.size() << " nodes" << std::endl;

// 		// std::cout << "node 19 degree: " << prob_graph.getNode(19)->getDegree()<< std::endl;
		
// 		PE1("\nOld: nodes: " << _n_nodes << ", edges: " << _n_edges)
// 		PE1("New: nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()<< std::endl)


// 		del_edges.clear();
// 		// del_nodes.clear();

// 	}

// 	CHECK(nodes_to_remove)

// 	PE1("should end up with: "<<_n_full_nodes - nodes_to_remove)

// 	// set up node maps for pre-processed graphs
// 	node_pp_map = std::vector<NodeID>(_n_full_nodes,-1);
// 	pp_node_map.clear();

// 	prob_graph.clear();
// 	terms.clear();

// 	// get all nodes that have more than 1 degree or are terminals and make maps
// 	for (int n = 0; n < _n_full_nodes; ++n){
// 		Node* node = orig_graph.getNode(n);
// 		if (node->getDegree() > 0 || node->isTerm()){
// 			node_pp_map[n] = pp_node_map.size();
// 			pp_node_map.push_back(n);
// 		}
// 	}

// 	CHECK(prob_graph.getNumNodes())
// 	CHECK(prob_graph.getNumEdges())
// 	CHECK(orig_graph.getNumNodes())
// 	CHECK(orig_graph.getNumEdges())

// 	// prob_graph = U_Graph(pp_node_map.size());
// 	U_Graph new_graph(pp_node_map.size());

// 	CHECK(new_graph.getNumNodes())
// 	CHECK(new_graph.getNumEdges())
// 	int edge_count = 0;

// 	CHECK(pp_node_map.size())

// 	// reconstruct pre-processed graph
// 	for (int n = 0; n < pp_node_map.size(); ++n){
// 		// PE1("node: " << n)
// 		Node* orig_node = orig_graph.getNode(pp_node_map[n]);
// 		NodeID orig_id = orig_node->getID();
// 		// set terminal status
// 		if (orig_node->isTerm()){
// 			new_graph.setTerm(n,true);
// 			// terms.push_back(n);
// 		}
// 		for (int e = 0; e < orig_node->getOutDegree(); ++e){
// 			int tgt = orig_node->getOutEdge(e)->getTgtID();
// 			int wt = orig_node->getOutEdge(e)->getWt();
// 			// PF1("adding edge: "<< new_graph.getNumEdges() << " = (" << orig_id+1 <<","<< tgt+1 <<","<<wt<<") = (" << n <<","<< node_pp_map[tgt] <<","<<wt<<").. ")
// 			new_graph.addEdge(n,node_pp_map[tgt],wt);
// 			// PE1("done!")
// 		}
// 	}


// 	prob_graph = orig_graph;

// 	terms.clear();
// 	for (int n = 0; n < prob_graph.getNumNodes(); ++n){
// 		if (prob_graph.getNode(n)->isTerm()){
// 			terms.push_back(n);
// 		}
// 	}

// 	PE("testing terms:")
// 	for (int t = 0; t < terms.size(); ++t){
// 		bool found = false;
// 		for (int n = 0; n < full_terms.size(); ++n){
// 			if (terms[t] == node_pp_map[full_terms[n]]){
// 				found = true;
// 				break;
// 			}
// 		}
// 		if (!found){
// 			PE("term " << terms[t] << " not matched!")
// 		}
// 	}

// 	PE("testing edges:")
// 	for (int e = 0; e < prob_graph.getNumEdges(); ++e){
// 		bool found = false;
// 		Edge* temp_edge = prob_graph.getEdge(e);

// 		if (!found){
// 			PE("term " << terms[e] << " not matched!")
// 		}
// 	}


// 	// for (int i = prob_graph.getNumNodes()-1; i>-1; --i){
// 	// 	if (prob_graph.getNode(i)->getDegree() ==0){
// 	// 		del_nodes.push_back(i);
// 	// 		pp_moves.push_back(Move(i,-2, -2));
// 	// 	}
// 	// }
// 	// prob_graph.deleteNodes(del_nodes);		

// 	PE("\nOld: nodes: " << orig_graph.getNumNodes() << ", edges: " << orig_graph.getNumEdges())
// 	PE("New: nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges()<< std::endl)

// 	// fw_dist.clear();
// 	// fw_next.clear();

// 	std::vector<std::vector<int> > pp_dist;
// 	std::vector<std::vector<int> > pp_next;

// 	PF("running floyd warshall algorithm... ")
// 	floyd_warshall(prob_graph,fw_dist,fw_next);
// 	PE("done!")

// 	fw_dist = pp_dist;
// 	fw_next = pp_next;

// 	findKeys();


// 	_n_nodes = prob_graph.getNumNodes();
// 	_n_edges = prob_graph.getNumEdges();

// }

void STPGModel::reconstruct(){
	std::vector<EdgeID> del_edges;

	PE("\n***** RECONSTRUCTING GRAPH *****")
	PF("Rebuilding nodes and edges... ")
	PE1(" ")
	while(pp_moves.size() > 0){
		Move move = pp_moves.back();
		PF1("Move(" << move._src << ", " << move._tgt << ", " << move._wt << ")")

		pp_moves.pop_back();

		switch (move._wt) {
          case -1:      // remove edge
          	PF1(" removing edge.. ")
          	del_edges.push_back(prob_graph.getEdgePair(move._src,move._tgt)->getID());
          	PE1("done!")
            continue;
          case -2:       // replace node
          	PF1(" inserting node.. ")
            prob_graph.insertNode(move._src);
            PE1("done!")
            continue;
          default:      // add edge
          	PF1(" adding edge.. ")
          	prob_graph.addEdge(move._src, move._tgt,move._wt);
          	PE1("done!")
            continue;
        }
	}

	PE("done!")
	prob_graph.deleteEdges(del_edges);
	PE("Nodes: " << prob_graph.getNumNodes() << ", Edges: " << prob_graph.getNumEdges())

	PF("running floyd warshall algorithm... ")
	floyd_warshall(prob_graph,fw_dist,fw_next);
	PE("done!")
}

// use next matrix from floyd warshall algorithm to find shortest path
int STPGModel::get_path(std::vector<int> &path_edges, int u, int v){
	int path_length=0;
	path_edges.clear();

	PE2("finding path " << u <<  " --> " << v)

	Edge *curr = prob_graph.getEdgePair(u,fw_next[u][v]);
	path_edges.push_back(curr->getID());
	path_length += curr->getWt();
	while (true){
		u = fw_next[u][v];
		if (u == v){
			break;
		}
		PF2("test get edge pair (" << u << ", " << fw_next[u][v] << ").. ")
		curr = prob_graph.getEdgePair(u,fw_next[u][v]);
		PE2("done!")
		path_edges.push_back(curr->getID());
		path_length += curr->getWt();
	}
	return path_length;
}

void STPGModel::findKeys(){
	for (int n = 0; n < prob_graph.getNumNodes(); ++n){
		if (prob_graph.getNode(n)->getDegree() > 2 || prob_graph.getNode(n)->isTerm()){
			prob_graph.getNode(n)->setKey(true);
		}
	}
}

STPGModel::~STPGModel(){
}