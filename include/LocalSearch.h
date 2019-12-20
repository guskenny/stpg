#ifndef __LocalSearch_H__
#define __LocalSearch_H__
#include <debug.h>
#include "u_graph.h"
#include <STPGModel.h>
#include "SettingsHandler.h"
#include <string>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <util.h>
#include <set_obj.h>
#include <random>
#include "QOL/CpuTimer.h"

class LocalSearch{
	private: 
		STPGModel *probModel;
		const SettingsHandler &sh;
		std::mt19937 rng;

	public:
		LocalSearch(const SettingsHandler sh, STPGModel *model) : sh(sh), probModel(model) {
	      std::cout << "LocalSearch initialised" << std::endl;
	      std::random_device r;
	      std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
	      rng = std::mt19937(seed);
    	};

		LocalSearch(const char *argv, const SettingsHandler &_sh);
		void get_keypath(const set_obj &sol, std::vector<int> &path_edges, std::vector<int> &kp_ends);
		void getNodeSets(const set_obj &sol, const std::vector<int> &kp_ends, std::vector<std::vector<int> >&node_sets);
		int doSwaps(set_obj &sol, std::vector<int> &obj_vals);
		int singleSwaps(set_obj &sol, std::vector<int> &obj_vals);
		void repairBest(set_obj &sol, const std::vector<int> &kp_ends);
		void repairRandom(set_obj &sol, const std::vector<int> &kp_ends);
		~LocalSearch();

};

#endif