#ifndef __QOL_COINFORMULATION__
#define __QOL_COINFORMULATION__

#ifdef _USE_COIN_

#include "QolMIPSolver.h"
#include "QolColFormulation.h"
#include <coin/OsiClpSolverInterface.hpp>
#include <coin/CbcModel.hpp>
#include <coin/CoinFinite.hpp>
// Heuristics
#include <coin/CbcHeuristicDiveFractional.hpp>
#include <coin/CbcHeuristicFPump.hpp>
#include <coin/CbcHeuristicLocal.hpp>
// Cuts
#include <coin/CbcStrategy.hpp>
#include <coin/CbcCutGenerator.hpp>
#include <coin/CglGomory.hpp>
#include <coin/CglProbing.hpp>
#include <coin/CglKnapsackCover.hpp>
#include <coin/CglRedSplit.hpp>
#include <coin/CglClique.hpp>
#include <coin/CglFlowCover.hpp>
#include <coin/CglMixedIntegerRounding2.hpp>
// Preprocessing
#include <coin/CglPreProcess.hpp>



namespace qol {

	class QolCoinMessageHandler : public CoinMessageHandler {
		public:
		 int qolLevel; // only print messages below this level (0=off)
		 QolCoinMessageHandler () : qolLevel(0) {}
		 virtual void checkSeverity(); // NEVER just abort (as the default implementation does)
		 virtual int print();
   };

	class CoinFormulation : public MIPSolver {
	public:
		// Interface to the COIN solver
		//std::auto_ptr<OsiSolverInterface> solver;
		//OsiClpSolverInterface* solver2;
		OsiClpSolverInterface solver;
		QolCoinMessageHandler msgHandler; // Use our own message handler to allow interruptions
		double initialsolvevalue;
		std::auto_ptr<CbcModel> cbcmodel;
		// QolColFormulation to hold the model until it is loaded, means I only need to write the
		// load function right now, plus a few others
		QolColFormulation *model;



		CoinFormulation(double timelimit, qol::SolutionType type = qol::PRIMAL);
		~CoinFormulation() { if(model) delete model; model=0;}

		virtual void setObjCoeff(qol::Variable v, double obj){
			if (model) { model->setObjCoeff(v, obj); return; }
			cbcmodel->solver()->setObjCoeff(v, obj);
		}

		virtual double getObjCoeff(qol::Variable v) const {
			if (model) return model->getObjCoeff(v);
			return cbcmodel->solver()->getObjCoefficients()[(int)v];
		}

		virtual void setVarUB(qol::Variable v, double ub) {
			if (model) { model->setVarUB(v, ub); return; }
			cbcmodel->solver()->setColUpper(v, ub);
		}

		virtual double getVarUB(qol::Variable v) const {
			if (model) return model->getVarUB(v);
			return cbcmodel->solver()->getColUpper()[(int)v];
		}

		virtual void setVarLB(qol::Variable v, double lb) {
			if (model){ model->setVarLB(v, lb); return; }
			cbcmodel->solver()->setColLower(v, lb);
		}

		virtual double getVarLB(qol::Variable v) const {
			if (model) return model->getVarLB(v);
			return cbcmodel->solver()->getColLower()[(int)v];
		}

		virtual void setVarType(qol::Variable v, char type) {
			qolAssert(type == qol::Variable::CONTINUOUS ||
				type == qol::Variable::BINARY ||
				type == qol::Variable::INTEGER,
				"Invalid variable type " << type << " for variable " << v);
			if (model) {
				model->setVarType(v, type); 
				return;
			}
			int va = v;
			int * var = &va;
			if (type == qol::Variable::CONTINUOUS)
				cbcmodel->solver()->setContinuous(var, 1);
			else if ( type == qol::Variable::INTEGER)
				cbcmodel->solver()->setInteger(var, 1);
			else{ // Binary
				cbcmodel->solver()->setInteger(var, 1);
				setVarLB(v, 0.0);
				setVarUB(v, 1.0);
			}
		}

		virtual char getVarType(qol::Variable var) const {
			if (model) return model->getVarType(var);
			switch (cbcmodel->solver()->getColType(true)[(int)var]){
				case 0: return qol::Variable::CONTINUOUS; break;
				case 1: return qol::Variable::BINARY; break;
				case 2: return qol::Variable::INTEGER; break;
				default: return qol::Variable::CONTINUOUS;
			}
		}

		virtual void setVarName(qol::Variable var, const std::string &name) {
			if (model){ model->setVarName(var, name); return; }
			cbcmodel->solver()->setColName(var, name);
		}

		virtual const std::string getVarName(qol::Variable var) const {
			if (model) return model->getVarName(var);
			return cbcmodel->solver()->getColName(var);
		}

		// constraint modification
		virtual void setSense(qol::Constraint cnstr, char s) {
			qolAssert(s == qol::Constraint::LE ||
				s == qol::Constraint::EQ ||
				s == qol::Constraint::GE, "Invalid constraint sense " << s);
			if (model){ model->setSense(cnstr, s); return; }
			//TODO: this for cbcmodel, a little fiddly
		}

		virtual char getSense(const qol::Constraint cnstr) const {
			if (model) return model->getSense(cnstr);
			return cbcmodel->solver()->getRowSense()[(int)cnstr];
		}

		virtual void setRHS(qol::Constraint cnstr, double rhs) {
			if (model){ model->setRHS(cnstr, rhs); return; }
			//TODO: this for cbcmodel, a little fiddly
		}

		virtual double getRHS(qol::Constraint cnstr) const {
			if (model) return model->getRHS(cnstr);
			return 0.0;
			//TODO: this for cbcmodel, a little fiddly
		}

		virtual void setConstrName(qol::Constraint cnstr, const std::string &name) {
			if (model){ model->setConstrName(cnstr, name); return; }
			cbcmodel->solver()->setRowName(cnstr, name);
		}

		virtual const std::string getConstrName(qol::Constraint cnstr) const {
			if (model) return model->getConstrName(cnstr);
			return cbcmodel->solver()->getRowName(cnstr);
		}

		virtual Index nVar() const {
			if (model) return model->nVar();
			return cbcmodel->getNumCols();
		}

		virtual Index nConstr() const {
			if (model) return model->nConstr();
			return cbcmodel->getNumRows();
		}

		void setNumVar(Index n){
			if (model) {
				model->setNumVar(n);
				return;
			}
			//TODO: Don't know if you can /w cbc / osi, at least the function isn't just sitting there
		}

		virtual ConstraintMIP addConstraint(const Row &row){
			if (model) return model->addConstraint(row);
			//TODO
			//cbcmodel->solver()->addRow(row.lhs, row.op, row.rhs, 0);
			//return getConstr(this->nConstr() - 1);
		}

		virtual bool relaxedIsExact(){ return false; }

		virtual void writeLP(const char *filename){
			if (model) {
				model->writeLP(filename); return;
			}
			cbcmodel->solver()->writeLp(filename);
		}

		virtual Status solveRelaxed() { return FAILED; }

		virtual Status solveExact(); 

		virtual double getObjectiveBound() const {
			if (model != 0) { return 0.0; }
			if (cbcmodel->solver()->isProvenOptimal())
				return cbcmodel->getBestPossibleObjValue();
			else
				return initialsolvevalue;
		}

		virtual double getObjective() const {
			if (model == 0) return cbcmodel->getObjValue();
			return 0.0;
		}

		virtual double getPrimal(const Variable &v) const {
			qolAssert(isLoaded(),
				"problem not loaded in CoinFormulation::getPrimal()");
			if(cbcmodel->bestSolution()==0) throw qol::Exception("No primal solution available");
			return cbcmodel->bestSolution()[(int)v];
			//return cbcmodel->getCbcColSolution()[(int)v];
			//return cbcmodel->getColSolution()[(int)v];
		}

		virtual double getDual(const Variable &v) const { 
			if (model != 0) return 0.0;
			if(cbcmodel->getReducedCost()==0) throw qol::Exception("No dual solution available");
			return cbcmodel->getReducedCost()[(int)v];
		}

		virtual double getPrimal(const Constraint &c) const { 
			if (model != 0) return 0.0; 
			return 0.0;
		}

		virtual double getDual(const Constraint &c) const { 
			if (model != 0)return 0.0; 
			if(cbcmodel->getRowPrice()==0) throw qol::Exception("No dual solution available");
			return cbcmodel->getRowPrice()[(int)c];
		}

		virtual void extractSolution() {}

		// New stuff so it compiles

		/// load data into coin (will load from this->model by default)
		virtual void load(qol::QolColFormulation *mip = 0,
			std::string name = "CoinFormulation");

		virtual void load(const qol::QolColFormulation *mip,
			std::string name = "CoinFormulation"){};

		/// is all data loaded into COIN?
		virtual bool isLoaded() const { return model == 0; }

		virtual std::string getErrorMessage(int error = -1) const{ return 0; };

		virtual void getRay(DblVec &ray) const
		{
			ray.resize(nVar());
			cbcmodel->solver()->getDualRays(nConstr()+nVar(), true);
			//int status = CPXgetray(env, LP, &ray[0]);
			//if (status) qolAssert(status == 0, getErrorMessage(status));
		}
		virtual double getRay(const Variable &v) const
		{
			DblVec ray; getRay(ray); return ray[ray.size()-1 - (nVar()-1-(int)v)];
		} // crude implementation

		virtual void setParameters(qol::Parameters p);

		virtual double getNodeBound() const
		{
			double obj;
			//int lastErrorCode = CPXgetcallbackinfo(env, callbackData, callbackWherefrom,
			//	CPX_CALLBACK_INFO_BEST_REMAINING, &obj);
			//if (lastErrorCode) qolAssert(lastErrorCode == 0, getErrorMessage(lastErrorCode));
			return obj;
		}
		virtual double getNodeBest() const
		{
			double obj;
			//int lastErrorCode = CPXgetcallbackinfo(env, callbackData, callbackWherefrom,
			//	CPX_CALLBACK_INFO_BEST_INTEGER, &obj);
			//if (lastErrorCode) qolAssert(lastErrorCode == 0, getErrorMessage(lastErrorCode));
			return obj;
		}

		virtual void addCut(const Row &row, bool isGlobal = true){};

		void convertSenseToBound(const char sense, const double right,
			const double range,
			double& lower, double& upper) const
		{
			double inf = COIN_DBL_MAX;
			switch (sense) {
			case 'E':
				lower = upper = right;
				break;
			case 'L':
				lower = -inf;
				upper = right;
				break;
			case 'G':
				lower = right;
				upper = inf;
				break;
			case 'R':
				lower = right - range;
				upper = right;
				break;
			case 'N':
				lower = -inf;
				upper = inf;
				break;
			}
		}
	};
}
#endif   //_USE_COIN_

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
 */
