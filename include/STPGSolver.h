#ifndef __STPGSolver_H__
#define __STPGSolver_H__

#include <debug.h>
#include <STPGModel.h>
#include <STPGMip.h>
#include <STPGMergeMip.h>
#include <SettingsHandler.h>
#include <util.h>
#include <random>
#include <algorithm>
#include <LocalSearch.h>
#include <SolutionMerger.h>
#include <sstream>
#include <ctime>

class STPGSolver{
	private:
		STPGModel *probModel;
		SettingsHandler sh;

	public:
		STPGSolver(int argc, const char **argv);
		~STPGSolver(){ free(probModel);}
		void solve();
		void DFSsplit(const merge_sol &merged_sols, const int group, std::vector<int> &connected_vars, std::mt19937 &rng);
		void split(merge_sol &merged_sols);
		void splitFactor(merge_sol &merged_sols);
		void getInitSol(set_obj &sol);
		void getRandSol(set_obj &sol);
		// int getSolValue(const set_obj &sol);
		// int getSolValue(const set_obj &sol, const U_Graph &graph);


};

#endif