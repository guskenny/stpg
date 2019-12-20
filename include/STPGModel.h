#ifndef __STPGModel_H__
#define __STPGModel_H__
#include <debug.h>
#include "u_graph.h"
#include "SettingsHandler.h"
#include <string>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <util.h>
#include <set_obj.h>
#include "QOL/CpuTimer.h"


const int LABEL = 0;
const int N1 = 1;
const int N2 = 2;
const int WT = 3;

struct Move
{
  int  _src, _tgt, _wt;
  Move(int src, int tgt, int wt) : _src(src),_tgt(tgt),_wt(wt){}
};


class STPGModel{
	private: 
		int _n_full_nodes;
		int _n_nodes;
		int _n_full_edges;
		int _n_edges;
		int _n_terms;
		int _optimal;
		std::string _name;
		std::string _fpath;
		std::string _inst_grp;
		std::vector<int> terms;
		std::vector<int> full_terms;
		std::vector<Move> pp_moves;
		const SettingsHandler &sh;
		// index is original graph node, value is pp graph node 
		std::vector<NodeID> node_pp_map;
		// index is pp graph node, value is original graph node 
		std::vector<NodeID> pp_node_map;

	public:
		std::vector<std::vector<int> > fw_dist;
		std::vector<std::vector<int> > fw_next;
		U_Graph orig_graph;
		U_Graph prob_graph;
		STPGModel();
		STPGModel(const char *argv, const SettingsHandler &_sh);
		void pre_process();
		void reconstruct();
		std::vector<int> getTerms(){return terms;}
		std::string getName(){return _name;}
		int getOptimal(){return _optimal;}
		std::string getFPath(){return _fpath;}
		std::string getInstGrp(){return _inst_grp;}
		int get_path(std::vector<int> &path,int u, int v);
		void findKeys();
		~STPGModel();

		const int n_nodes(){return _n_nodes;};
		const int n_edges(){return _n_edges;};
		const int n_vars(){return _n_nodes + _n_edges;};
		const int n_terms(){return terms.size();};
		
};

#endif