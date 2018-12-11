/***************************************************************
* Constraint Programming / Beam Search / Ant Colony Optimisation Hybrid
* Author: Andreas.Ernst@monash.edu  December 2015
****************************************************************/

#ifndef __QOL_CPFORMULATION_H__
#define __QOL_CPFORMULATION_H__

#include "QolColFormulation.h"
#include "QolMIPSolver.h"
#include <iostream>
#include <limits>
#include <set>
#include  <boost/ptr_container/ptr_vector.hpp>

namespace qol {

  // for future restructuring:
  // Should set up a CProgHeurFormulation class that can form a plug & play framework
  // with components that replace:
  // CPstate - with lagrangian relaxation or other more complex constraint propagation
  // SearchObject - try various solve/optimse methods
  // Perhaps some additional cut generation as optional extra

  // the intention of this class is to do the constraint propagation
  // the CPstate also captures the current state when doing backtracking
  // or other heuristic changes to the state that we may want to unwind
  // by copying an old state - should probably go somewhere else eventually
	class CPstate { // current state of variables
	public:
		enum Consistency { Unknown = -1, Infeasible = 0, OK = 1, Fixed=2 }; // Unknown prior to propagate
	protected:
		std::vector<double> collb, colub;  // lowest/highest possible value 
		std::vector<double> rowlb, rowub;  // lowest/highest row value
		std::vector<double> maxChg; // max change in LHS for any row r due to a single variable moving
		std::vector<int> freeCnt; // number of non-fixed variables for each constraint (zero if non-active)
		std::vector<int> constrCnt; // number of active constraints for a variable (zero if fixed)
		std::vector<Index> rows; // rows that are still worth checking
    	const RowColFormulation *mip;
		bool changed;//std::set<qol::Variable> changed; // for faster propagation
		double cutoff; // heuristic (upper) bound for objective
		enum BoundType {LOWER=1,UPPER=2,BOTH=3};
		Consistency updateRowBounds(Index v,BoundType which=BOTH); // update contribution to max/min row values for variable
		void reverseRowBounds(Index v,BoundType which=BOTH); // reverse contribution
		void deactivateRow(Index r,size_t pos=0); // update counters - this row no longer active, optional position in rows[]
		bool checkRowInactive(Index r); // recalculate maxChg[r],freeCnt[r] and return true if constraint can't be violated
	public:
		CPstate() : mip(0),cutoff(qol::inf),state(Unknown) {}
		// CP state does NOT copy the MIP - it must not be deleted by the caller
		// until after the CPstate destructor has been called (or at 
		// least after the last call to propagate() )
		CPstate(const RowColFormulation &MIP) { init(&MIP); }
		void init(const RowColFormulation *mip);
		// basic variable query/modification
		virtual bool isFixed(qol::Variable v) const { return collb[v] == colub[v]; }
		void checkRowBnds() const { // for debugging consistency checks
			std::vector<double> lb(mip->nConstr(),0),ub(mip->nConstr(),0);
			std::vector<int> free(mip->nConstr(),0);
			for(Index v=0;v<mip->matrix.size();++v)
				for (Index c = 0; c < mip->matrix[v].size(); ++c) {
					const double coeff = mip->matrix[v][c].first;
					const Index r = mip->matrix[v][c].second;
					if(!isFixed(v)) ++free[r];
					if (coeff > 0) {
						lb[r] += coeff*collb[v];
						ub[r] += coeff*colub[v];
					} else {
						lb[r] += coeff*colub[v];
						ub[r] += coeff*collb[v];
					}
				}
			for(Index r=0;r<mip->nConstr();++r){
				if( fabs(rowlb[r]-lb[r]) > eps){
					std::cout<<"ERRRO "<<r<<": " << mip->getConstrName(r) << "lb="<<rowlb[r] <<"!="<<lb[r]<<std::endl;
				}
				if( fabs(rowub[r]-ub[r]) > eps){
					std::cout<<"ERRRO "<<r<<": " << mip->getConstrName(r) << "ub="<<rowub[r] <<"!="<<ub[r]<<std::endl;
				}
				if(freeCnt[r] != 0 && free[r] != freeCnt[r]){
					std::cout<<"ERRRO "<<r<<": " << mip->getConstrName(r) << "free="<<freeCnt[r]<<"!="<<free[r]<<std::endl;
				}
			}
		}
		virtual Consistency fixVar(qol::Variable v, double val) {
			BoundType bnds= BOTH;
			if(fabs(collb[v]-val) <= eps){
				bnds=UPPER;
			}
			if(fabs(colub[v]-val) <= eps){
				if(bnds==UPPER) return OK; // nothing is changing
				bnds=LOWER;
			}
			//std::cout << "Fixing "<<mip->getVarName(v) << "==" <<val<<std::endl;
			reverseRowBounds(v,bnds);
			collb[v] = colub[v] = val;
			const Consistency status=updateRowBounds(v,bnds);
			//checkRowBnds(); // debug check 
			constrCnt[v]=0; // should be the case anyway but just in case
			changed = true;//changed.insert(v);
			return status;
		}
		virtual Consistency setVarUB(qol::Variable v, double ub) {
			if(mip->coltype[v] != qol::Variable::CONTINUOUS)
				ub = floor(ub + eps); // make sure we have an integer
			if (ub < colub[int(v)]) { // only if tighter
				if( ub < collb[int(v)]-eps ){ return (state = Infeasible);  }
				//std::cout << "Fixing "<<mip->getVarName(v) << "<=" <<ub<<std::endl;
				reverseRowBounds(v,UPPER);
				colub[int(v)] = std::max(collb[int(v)],ub); // max is just to prevent rounding error
				changed = true;//changed.insert(v);
				return updateRowBounds(v,UPPER);
				//checkRowBnds(); // debug check 
			}
			return OK;
		}
		const std::vector<double> getVarUB() const { return colub;  }
		virtual double getVarUB(qol::Variable v) const {
			return colub[int(v)];
		}
		// update lower bound on a variable. Returns OK or Infeasible
		virtual Consistency setVarLB(qol::Variable v, double lb) {
			if (mip->coltype[v] != qol::Variable::CONTINUOUS)
				lb = ceil(lb - eps); // make sure we have an integer
			if (lb > collb[(int)v]) { // only if tighter
				if( lb > colub[int(v)]+eps ){ return (state = Infeasible); }
				//std::cout << "Fixing "<<mip->getVarName(v) << ">=" <<lb<<std::endl;
				reverseRowBounds(v,LOWER);
				collb[int(v)] = std::min(lb,colub[int(v)]); // min to prevent rounding error 
				changed = true;//changed.insert(v); 
				return updateRowBounds(v,LOWER);
				//checkRowBnds(); // debug check V
			}
			return OK;
		}
		const std::vector<double> getVarLB() const { return collb; }
		virtual double getVarLB(qol::Variable v) const {
			return collb[int(v)];
		}
		// Warning: setting cutoff means any solution with higher value will 
		// be considered "Infeasible"
		virtual void setCutoff(double heuristicBound); // need to set cutoff & do some propagation
		virtual double getCutoff() const { return cutoff; }
		Consistency propagate(); // returns true if is feasible 
		Consistency getState() const { return state; }
		/// getObjectiveBound() returns a (weak) lower/upper bound for
		/// a minimisation/maximisation problem
		double getObjectiveBound() const;

		// check whether solution satisfies constraint and optionally
		// return the product of mip->matrix * solution
		Consistency checkSolution(const std::vector<double> &solution,
			std::vector<double> *LHS = 0) const;
		// check that current variable bounds match consistency state
		Consistency verify(int verbose=0) const;

	protected:
		Consistency state;
	};

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Variable Group: want to detect places where we set (at most) 1 of several binaries to 1
	// This is used to get a better ACO pheremone representation 
	///////////////////////////////////////////////////////////////////////////////////////////////
	class VariableGroup {
	protected:
		Index idx;
	public:
		void setIndex(Index i) { idx = i; }
		Index getIndex() const { return idx; }
		virtual Index nValues() const { return 0; } // number of possible values
		virtual Index nFeasible(const CPstate &s) const { return 0; } // number of possible values
		// returns number of feasible values (like nFeasible) together with list of 
		// feasible values (for use in setVariables)
		virtual Index getFeasible(const CPstate &s,std::vector<Index> &feas) const {return 0;}
		virtual void setVariables(CPstate &state, Index val) const {}
		virtual void setFeasible(CPstate &state, Index val) const {}
		virtual double getObjectiveBound(CPstate &state) const {	return 0;}
		virtual void print(const CPstate &state) const {
			std::cout << idx << ":" << std::endl;
		}
	};
	class SOSequalOne : public VariableGroup { // group of binaries with exactly one variable = 1
	protected:
		std::vector<Index> vars;
	public:
		SOSequalOne() {}
		virtual void addVar(Index v) { vars.push_back(v); }
		virtual Index nValues() const { return vars.size(); }
		virtual Index nFeasible(const CPstate &state) const {
			Index cnt=0;
			for(size_t i = 0; i < vars.size(); ++i) { if( !state.isFixed(vars[i]) ) ++cnt; }
			return cnt;
		}
		virtual Index getFeasible(const CPstate &state,std::vector<Index> &feas) const {
			feas.clear(); feas.reserve(nValues()); 
			for(size_t i = 0; i < vars.size(); ++i) { 
				if( !state.isFixed(vars[i]))  feas.push_back(i); 
			}
			return (Index)feas.size();
		}
		virtual void setVariables(CPstate &state, Index val) const
		{
			for (size_t i = 0; i < vars.size() && state.getState()==CPstate::OK; ++i)
				if (i == val) state.fixVar(vars[i],1);
				else state.fixVar(vars[i],0);
		}
		virtual void setFeasible(CPstate &state, Index val) const
		{  Index cnt=0;
			for (size_t i = 0; i < vars.size() && state.getState()==CPstate::OK; ++i){
				if(state.isFixed(vars[i])) continue;
				if (val == cnt) state.fixVar(vars[i],1);
				else state.fixVar(vars[i],0);
				++cnt; 
			}
		}
		virtual void print(const CPstate &state) const {
			std::cout << idx << ":";
			for(size_t i=0;i<vars.size();++i)
				std::cout << " " << vars[i] << (state.isFixed(vars[i]) ? '!' : '~');
			std::cout << " = " << nFeasible(state) << std::endl;
		}
	};
	class SOSatmostOne: public SOSequalOne { // group of binaries with exactly one variable <= 1
	public:
		SOSatmostOne() {}
		//void addVar(Index v) { vars.push_back(v); }
		virtual Index nValues() const { return vars.size()+1; }
		virtual Index nFeasible(const CPstate &state) const {
			Index cnt=0;
			for(size_t i = 0; i < vars.size(); ++i) { if( !state.isFixed(vars[i]) ) ++cnt; }
			return (cnt > 0 ? (cnt+1) : 0);
		}
		virtual Index getFeasible(const CPstate &state,std::vector<Index> &feas) const {
			feas.clear(); feas.reserve(nValues()); 
			for(size_t i = 0; i < vars.size(); ++i) { 
				if( !state.isFixed(vars[i]))  feas.push_back(i); 
			}
			if(! feas.empty() ) // InvalidIndex: all binaries == 0
				feas.push_back(qol::InvalidIndex); 
			return (Index)feas.size();
		}		
		//virtual void setVariables(CPstate &state, Index val) const
		//{ for (size_t i = 0; i < vars.size(); ++i)
		//	if (i == val) state.fixVar(1);	else state.fixVar(0);
		//}s
	};
	class SingleVarGroup : public VariableGroup {
	protected:
		Index var;
		double minVal, maxVal;
		int steps;
	public:
		SingleVarGroup(Index v=InvalidIndex,
			double minV = 0, double maxV = 1, int stp = 1)
			: var(v), minVal(minV), maxVal(maxV), steps(stp) {}
		virtual Index nValues() const { return steps + 1; }
		virtual Index nFeasible(const CPstate &state) const {
			const Index minStep = (Index)ceil( (std::max(minVal,state.getVarLB(var)-eps)-minVal)/steps);
			const Index maxStep = (Index)floor((state.getVarUB(var)+eps-minVal)/steps)+1;
			return maxStep - minStep;
		}
		virtual Index getFeasible(const CPstate &state,std::vector<Index> &feas) const {
			const Index minStep = (Index)ceil( (std::max(minVal,state.getVarLB(var)-eps)-minVal)/steps);
			const Index maxStep = (Index)floor((state.getVarUB(var)+eps-minVal)/steps)+1;
			feas.resize(maxStep-minStep);
			for(Index i=0;i<feas.size();++i) feas[i] = i+minStep;
			return (Index)feas.size();
		}		
		virtual void setFeasible(CPstate &state,Index val){
		    const Index minStep = (Index)ceil( (std::max(minVal,state.getVarLB(var)-eps)-minVal)/steps);
			state.fixVar(var,minVal + (val+minStep)*(maxVal - minVal) / steps);	
		}
		virtual void setVariables(CPstate &state, Index val)
		{   
			state.fixVar(var,minVal + val*(maxVal - minVal) / steps);	
		}
		virtual void print(const CPstate &state) const {
			std::cout << idx << ":" << " " << var << (state.isFixed(var) ? '!' : '~') << std::endl;
		}
	};




  /// CPFormulation is like a column formulation but provides methods for 
  /// constraint propagation 
  class CPBeamACO  : public FormulationMIPSolver<RowColFormulation> {
  protected:
	  double bestObj;
	  std::vector<double> bestSolution;
	  std::vector<double> bestLHS; // = mip.matrix * bestSolution
	  void updateBest(CPstate *s = 0) {
		  if (!s) s = state;
		  if (s->getState() == CPstate::Fixed && s->getObjectiveBound() < bestObj) {
			  if(s->checkSolution(s->getVarLB(), &bestLHS) != CPstate::OK){
				  std::cerr<< "ERROR: infeasible solution with value " << s->getObjectiveBound() <<std::endl;
				  s->verify(1);
			  }else{
				bestSolution = s->getVarLB(); // == getVarUB()
				bestObj = s->getObjectiveBound();
			  }
		  }
	  }
	  Parameters params;
	// define variables groups by running through all constraints
	void initVarGroups(); // called by init()
	void initPheromone(); // init PheromoneMatrix, called by init();
  public: 

	CPstate *state; // current state (after a solve)
	boost::ptr_vector<VariableGroup> varGroup; 
	typedef std::vector<DblVec> PheromoneMatrix;
	PheromoneMatrix pheromone;

	void init(); // initialise state, rowMatrix etc
	CPBeamACO() : bestObj(qol::inf), state(0) {}
	~CPBeamACO() { if(state) delete state;}

	
    // methods needed to implement a MIPsolver
	virtual Status solveRelaxed() { 
		init(); // initalise and propagate
		switch (state->getState()) {
		case CPstate::Unknown: return qol::FAILED;
		case CPstate::Infeasible: return qol::INFEASIBLE;
		case CPstate::OK: 
			return (state->getObjectiveBound() <= -qol::inf) ?
				qol::UNBOUNDED : qol::HEURISTIC;
		case CPstate::Fixed: updateBest();
			return qol::HEURISTIC;
		}
		return qol::HEURISTIC; // cannot reach here
	}
	virtual Status solveBeam(); // alternative solveExact
	virtual Status solveExact();
	virtual double getObjective() const { return bestObj;  }
	virtual double getObjectiveBound() const {
		qolAssert(state != 0, "Cannot get object bound for before solving CPBeamACO");
		return state->getObjectiveBound();
	}

	virtual double getPrimal(const Variable &v) const
	{
		if(bestSolution.size()==0) throw qol::Exception("No primal solution available");
		if (bestSolution.size() == 0) return 0; // no solution available
		qolAssert(v < bestSolution.size(), "no solution for variable in CPBeamACO");
		return bestSolution[v];
	}
	virtual double getPrimal(const Constraint &c) const
	{  
		if(bestLHS.size()==0) throw qol::Exception("No primal solution available");
		if (bestLHS.empty()) return 0; // no solution available
		qolAssert((Index)c < bestLHS.size(), "no solution for constraint in CPBeamACO");
		return bestLHS[c];
	}

	/// get the reduced cost of the variable (only after relaxed solve)
	virtual double getDual(const Variable &v) const {
		qolAssert(false, "getDual(Variable) but no dual available in CPBeamACO");
		throw qol::Exception("No dual solution available");
		return 0;
	}

	/// get the dual value (only after relaxed solve)
	virtual double getDual(const Constraint &c) const {
		qolAssert(false, "getDual(constraint) but no dual available in CPBeamACO");
		throw qol::Exception("No dual solution available");
		return 0;
	}
	virtual void setPriority(const Variable &v, int pri){}; // 0 = default, higher means fix/branch earlier
	// addSOS adds a special ordered set of type 1 (ie at most 1 variable = 1) to the model
	// the optional priority provides a priority for branching on this set
	virtual void addSOS(const std::vector<Variable> &vars, int pri = 0){}; // declare the variables to be an SOS set

	virtual void setParameters(qol::Parameters params_) { params = params_; } // no parameters to set
  protected:
	  std::vector<int> priority; // priority for each variable;
	  std::vector<std::vector<Index> > SOS; // special ordered sets

  }; // end QolColFormulation
}
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
