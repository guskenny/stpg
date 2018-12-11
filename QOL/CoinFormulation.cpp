#ifdef _USE_COIN_
#include "CoinFormulation.h"
#include <boost/thread.hpp>

namespace qol {

	void  QolCoinMessageHandler::checkSeverity() {
		boost::this_thread::interruption_point(); // allow thread interruption
		if (currentMessage_.severity_=='S') {
			throw qol::Exception("Severe error in coin solver");
			//CoinMessageHandler does abort() - a really bad idea!
		}
	}
	int QolCoinMessageHandler::print() {
		int result = 0;
		if(qolLevel > 0)
		   result = CoinMessageHandler::print(); // do normal printing
		boost::this_thread::interruption_point(); // insert an interruption point
		return result;
	}


	CoinFormulation::CoinFormulation(double timelimit, qol::SolutionType type) : MIPSolver(type), model(new QolColFormulation)
	{
		//OsiClpSolverInterface solver = OsiClpSolverInterface();
		solver = OsiClpSolverInterface();
		solver.getModelPtr()->setMaximumSeconds(timelimit);
		//solver2 = &solver;
		cbcmodel.reset(new CbcModel(solver));
		cbcmodel->setMaximumSeconds(timelimit);
		cbcmodel->solver()->passInMessageHandler(&msgHandler);
		cbcmodel->passInMessageHandler(&msgHandler);
	};

	Status CoinFormulation::solveExact()
	{
		if (model != 0) { load();  }
		try{
			int numcols = cbcmodel->getNumCols();
			//cbcmodel->getEventHandler()->setAction(CbcEventHandler::CbcEvent::
			// Set up cuts
			
			CglProbing generator1;
			generator1.setUsingObjective(true);
			generator1.setMaxPass(1);
			generator1.setMaxPassRoot(3);
			// Number of unsatisfied variables to look at
			generator1.setMaxProbe(10);
			generator1.setMaxProbeRoot(100);
			// How far to follow the consequences
			generator1.setMaxLook(5);
			generator1.setMaxLookRoot(50);
			// Only look at rows with fewer than this number of elements
			//generator1.setMaxElements(200);
			generator1.setRowCuts(3);
			
			CglGomory generator2;
			// try larger limit
			//generator2.setLimit(30);
			
			CglKnapsackCover generator3;

			CglRedSplit generator4;
			// try larger limit
			//generator4.setLimit(200);
			
			CglClique generator5;
			generator5.setStarCliqueReport(false);
			generator5.setRowCliqueReport(false);
			
			CglMixedIntegerRounding2 mixedGen;
			CglFlowCover flowGen;
			
			// Add in generators
			// Experiment with -1 and -99 etc
			if (numcols < 7000){
				// try larger limit
				generator2.setLimit(300);
				generator4.setLimit(200);
				cbcmodel->addCutGenerator(&generator1, 99, "Probing");
				//cbcmodel->addCutGenerator(&generator2, -1, "Gomory");
				cbcmodel->addCutGenerator(&generator3, 99, "Knapsack");
				//cbcmodel->addCutGenerator(&generator4, -1, "RedSplit");
				cbcmodel->addCutGenerator(&generator5, 99, "Clique");
				cbcmodel->addCutGenerator(&flowGen, 99, "FlowCover");
				cbcmodel->addCutGenerator(&mixedGen, 99, "MixedIntegerRounding");
				cbcmodel->addCutGenerator(&generator2, 99, "Gomory");
				cbcmodel->addCutGenerator(&generator4, 99, "RedSplit");
			}
			else{
				generator2.setLimit(30);
				cbcmodel->addCutGenerator(&generator2, -99, "Gomory");
			}
			//
			// Say we want timings
			
			int numberGenerators = cbcmodel->numberCutGenerators();
			int iGenerator;
			for (iGenerator = 0; iGenerator<numberGenerators; iGenerator++) {
				CbcCutGenerator * generator = cbcmodel->cutGenerator(iGenerator);
				generator->setTiming(true);
			}
			
			// Set up heuristics
			CbcHeuristicDiveFractional heur1(*cbcmodel);
			heur1.setWhen(3);
			heur1.setHeuristicName("Dive");
			cbcmodel->addHeuristic(&heur1);

			if (numcols > 7000){
				cbcmodel->setNumberStrong(0);
				CbcHeuristicFPump heur2(*cbcmodel);
				heur2.setWhen(1);
				heur2.setHeuristicName("FPUMP");
				heur2.setMaximumTime(cbcmodel->getDblParam(CbcModel::CbcMaximumSeconds)*0.9);
				cbcmodel->addHeuristic(&heur2);
			}
			// Allow rounding heuristic
			
			CbcRounding heuristic1(*cbcmodel);
			heuristic1.setHeuristicName("Round");
			cbcmodel->addHeuristic(&heuristic1);
			
			// And local search when new solution found

			CbcHeuristicLocal heuristic2(*cbcmodel);
			heuristic2.setHeuristicName("Local");
			cbcmodel->addHeuristic(&heuristic2);
			
			//Preprocess
			/*
			if (numcols < 10000){
				CbcStrategyDefault strategy(true, 1, 1);
				// Set up pre-processing to find sos if wanted
				strategy.setupPreProcessing(2);
				cbcmodel->setStrategy(strategy);
			}
			*/

			// Set up branching priorities (basically say branch on the assignment varaibles first)

			int iColumn;
			int numberColumns = cbcmodel->getNumCols();
			/* We are going to add a single follow-on object but we
			want to give low priority to existing integers
			As the default priority is 1000 we don't actually need to give
			integer priorities but it is here to show how.
			*/
			// Normal integer priorities
			int * priority = new int[numberColumns];
			int numberIntegers = 0;
			for (iColumn = 0; iColumn<numberColumns; iColumn++) {
				if (cbcmodel->isBinary(iColumn)) {
					priority[numberIntegers++] = 1; // high priority
				}
				else if (cbcmodel->isInteger(iColumn)){
					priority[numberIntegers++] = 100; // low priority
				}
			}
			/* Second parameter is set to true for objects,
			and false for integers. This indicates integers */
			cbcmodel->passInPriorities(priority, false);
			delete[] priority;


			// do the solve and get result status
			cbcmodel->initialSolve();

			initialsolvevalue = cbcmodel->getBestPossibleObjValue();

			if (cbcmodel->solver()->isAbandoned()){ 
				return FAILED; 
			}
			/// Is optimality proven?
			else if (cbcmodel->solver()->isProvenOptimal()){
				//std::cout << "Optimal" << std::endl;
				// optimal! good to keep going
			}
			// Is primal infeasiblity proven?
			else if (cbcmodel->solver()->isProvenPrimalInfeasible()) {
				return INFEASIBLE; // if infeasible now (relaxed), gonna always be
			}
			// Is it unbounded?
			else if (cbcmodel->solver()->isProvenDualInfeasible()) {
				return UNBOUNDED;
			}
			else{
				// Did you end the solve early?
				return FAILED; // probably a timeout
			}

			cbcmodel->branchAndBound();
		}catch(boost::thread_interrupted){
			return ABORTED;
		}catch(...){
			if(cbcmodel->bestSolution() == 0)
				return FAILED;
			else
				return HEURISTIC; // got some solution at least
		}
			/// Is primal infeasiblity proven?
			if (cbcmodel->isProvenInfeasible()) {
				if (cbcmodel->solver()->isAbandoned()){
					return FAILED; // numerical difficulties
				}
				if (cbcmodel->solver()->isProvenDualInfeasible()){
					return FAILED;
				}
				if (cbcmodel->solver()->isProvenPrimalInfeasible()){
					return INFEASIBLE; // actually infeasible?
				}
				if (cbcmodel->bestSolution() == 0){
					return FAILED;
				}
				return HEURISTIC;
			}

			// Did you end the solve early?
			if (cbcmodel->isSecondsLimitReached()){
				if (cbcmodel->bestSolution() == 0){
					return FAILED;
				}
				return HEURISTIC;
			}

			/// Are there a numerical difficulties?
			if (cbcmodel->isAbandoned()) return FAILED;
			
			// Is it unbounded?
			if (cbcmodel->isProvenDualInfeasible()) return UNBOUNDED;

			/// Is optimality proven?
			if (cbcmodel->isProvenOptimal()){
				if (cbcmodel->solver()->isProvenOptimal()) // can think its optimal when its not...
					return OPTIMAL;
				else
					return HEURISTIC;
			}
			return (cbcmodel->bestSolution() == 0) ? FAILED : HEURISTIC;

		} // end solveExact()

	void CoinFormulation::load(QolColFormulation *mip, std::string name){
		if (mip == 0){
			if (model == 0) return;	// nothing to load
			mip = model;
		}
		const size_t nCol = mip->collb.size();
		if (nCol == 0){
			std::cout << "Error in " << __FILE__ << ":" << __LINE__ << ": Trivial problem, no columns provided." << std::endl;
		}
		size_t nz = 0; // number of non-zeros
		
		// load matrix in one hit
		for (size_t i = 0; i < mip->matrix.size(); ++i) nz += mip->matrix[i].size();
		std::vector<int> index;
		std::vector<CoinBigIndex> start;
		std::vector<double> value;
		index.reserve(nz);
		value.reserve(nz);
		start.reserve(nCol + 1);
		start.push_back(0);
		for (size_t i = 0; i < nCol; ++i){
			std::sort(mip->matrix[i].begin(), mip->matrix[i].end());
			for (size_t j = 0; j<mip->matrix[i].size(); ++j){
				value.push_back(mip->matrix[i][j].first);
				index.push_back(mip->matrix[i][j].second);
			}
			start.push_back((int)value.size());
		}
		int numrows = (int)mip->rhs.size();
		std::vector<double> range(mip->rhs.size(), 0.0); // no ranged constraints

		cbcmodel->solver()->loadProblem((int)nCol, (int)mip->rhs.size(),
			&start[0], &index[0], &value[0], &mip->collb[0], &mip->colub[0],
			&mip->cost[0], &mip->sense[0], &mip->rhs[0], &range[0]);

		// set other stuff (like names and variable type)
		cbcmodel->solver()->setIntParam(OsiNameDiscipline, 2); // full names for everything
		//solver->setColNames(colname,0,(int)colname.size(),0); // set all column names
		for (size_t i = 0; i<mip->colname.size(); ++i)
			if (mip->colname[i] != "") cbcmodel->solver()->setColName((int)i, mip->colname[i]);
		for (size_t i = 0; i<mip->rowname.size(); ++i)
			if (mip->rowname[i] != "") cbcmodel->solver()->setRowName((int)i, mip->rowname[i]);

		std::vector<int> integerCols;
		for (size_t i = 0; i<mip->coltype.size(); ++i)
		if (mip->coltype[i] == 'I' || mip->coltype[i] == 'B') integerCols.push_back((int)i);
		for (int i = (int)integerCols.size() - 1; i >= 0; --i) {
			cbcmodel->solver()->setInteger(integerCols[i]);
		}

		
		if (model){					// indicate that model is loaded
			delete model;
			model = 0;
		}

	} // end load()

	/* // Fine control over message handling and verbosity - not needed for now
	class CoinMessageHandlerQOL : public CoinMessageHandler {
	public: 
		virtual int print() { 
			std::cout << "XXX " << this->currentMessage().message() << std::endl; 
			return 0;
		}
	};
	*/
	void CoinFormulation::setParameters(qol::Parameters params) {
		double val;
		val = params.getParamValue(VERBOSITY);
		if (val != InvalidParam){
			msgHandler.qolLevel=int(val);
			// don't stop Coin/Cbc from stopping to print as this is the only
			// way we can insert interruption points :-(
			//cbcmodel->setLogLevel(int(val)); 
			//cbcmodel->passInMessageHandler(new CoinMessageHandlerQOL);
			//cbcmodel->setPrintingMode(int(val));
			//if (val < 1){
			//	cbcmodel->setPrintFrequency(99999); // avoid printing
			//	cbcmodel->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
			//}
		}
		/*
		val = params.getParamValue(THREADS);
		if (val != InvalidParam){
			lastErrorCode = CPXsetintparam(env, CPX_PARAM_THREADS, int(val));
			qolAssert(lastErrorCode == 0, getErrorMessage(lastErrorCode));
		}
		*/
		val = params.getParamValue(PRINTFREQ);
		if (val != InvalidParam){
			cbcmodel->setPrintFrequency(std::max(1,int(val)));
		}
		val = params.getParamValue(ABSGAP);
		if (val != InvalidParam){
			bool OK=cbcmodel->setDblParam(CbcModel::CbcAllowableGap, val);
			qolAssert(OK, "Coin CBC: Failed to set absolute gap");
		}
		val = params.getParamValue(RELGAP);
		if (val != InvalidParam){
			bool OK = cbcmodel->setDblParam(CbcModel::CbcAllowableFractionGap, val);
			qolAssert(OK, "Coin CBC: Failed to set relative gap");
		}
		val = params.getParamValue(EPINT);
		if (val != InvalidParam){
			bool OK = cbcmodel->setDblParam(CbcModel::CbcIntegerTolerance, val);
			qolAssert(OK, "Coin CBC: Failed to set integer tolerance");
		}
	
		val = params.getParamValue(ITLIMIT);
		if (val != InvalidParam){
			bool OK = cbcmodel->solver()->setIntParam(OsiMaxNumIteration,int(val));
			qolAssert(OK, "Coin CBC: Failed to set iteration limit");
		} 

		val = params.getParamValue(TIMELIMIT);
		if (val != InvalidParam){
			bool OK = cbcmodel->setDblParam(CbcModel::CbcMaximumSeconds, val);
			qolAssert(OK, "Coin CBC: Failed to set time limit");
		}
		// Note: CBC is single threaded so CPUTIME & WALLTIME should be about the same
		val = params.getParamValue(CPUTIME);
		if (val != InvalidParam){
			bool OK = cbcmodel->setDblParam(CbcModel::CbcMaximumSeconds, val);
			qolAssert(OK, "Coin CBC: Failed to set time limit");
		}
		val = params.getParamValue(WALLTIME);
		if (val != InvalidParam){
			bool OK = cbcmodel->setDblParam(CbcModel::CbcMaximumSeconds, val);
			qolAssert(OK, "Coin CBC: Failed to set time limit");
		}
		/* CBC has HeuristicFractionGap but no direct control over heuristic frequency
		val = params.getParamValue(HEURISTIC_FREQ);
		if (val != InvalidParam){
			lastErrorCode = CPXsetlongparam(env, CPX_PARAM_HEURFREQ, long(val));
			qolAssert(lastErrorCode == 0, getErrorMessage(lastErrorCode));
		}
		*/
	} // end setParam



}
#endif
