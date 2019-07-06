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

class STPGSolver{
	private:
		STPGModel *probModel;
		SettingsHandler sh;

	public:
		STPGSolver(int argc, const char **argv);
		~STPGSolver(){ free(probModel);}
		void solve();
		void getInitSol(set_obj &sol);
		void getRandSol(set_obj &sol);
		// int getSolValue(const set_obj &sol);
		// int getSolValue(const set_obj &sol, const U_Graph &graph);


};

#endif