#ifndef __QOLSOLVER_FORMULATION_H__
#define __QOLSOLVER_FORMULATION_H__

#include "QolMIP.h"
#include "QolSolutionPrimal.h"
#include "QolSolutionFull.h"
#include "QolParams.h"
#include "QolColFormulation.h"

namespace qol {
  class Callback;
  
  class MIPSolver: public MIP {
  protected:
	  /* what is the point of this? - each solver stores current solution in its own way */
    SolutionBase * sln;
    qol::SolutionType _type;
	Callback *callback;
  public:
	  MIPSolver():sln(0),_type(qol::PRIMAL),callback(0) {}
	  
	  MIPSolver(SolutionType type):MIP(), _type(type), callback(0){
		if (type==qol::FULL) sln = new SolutionFull();
		else sln = new SolutionPrimal();
	  }
	  virtual ~MIPSolver() {} // { if (sln) delete (sln);}
    virtual void writeLP (const char * filename) {};
	/// work out what solver is really being used. Default is QOLFORM
	/// (ie no solver)
	virtual FormulationType getSolverType() const { return QOLFORM; } 
	
	virtual void load(const qol::QolColFormulation *mip=0,
					  std::string name="MIP"
					  ) = 0;	// ensure all data is loaded into the solver
	virtual bool isLoaded() const {return false; } // is all data loaded?
    virtual Status solveRelaxed() = 0;
    virtual Status solveExact() = 0;
	/// optimize() matches naming convention of cplex/gurobi
	virtual Status optimize() { return solveExact(); }


	/// initial primal solution - no need to set this for every variable normally
 	virtual void setPrimalStart(const Variable &v,double val)
		  { std::cerr << "ERROR: warmstart with initial solution not implemented for solver "
					  << getSolverType() << std::endl; }
	/// val must be of length number of variables - gives complete solution
	virtual void setPrimalStart(const std::vector<double> &val)
		  { std::cerr << "ERROR: warmstart with initial solution not implemented for solver "
					  << getSolverType() << std::endl; }

	/// error messages: translate solver specific error codes into something
	/// meaningful
	/// The actual error code can be provided or will return the message
	/// of the last error by default
	/// Note: the ability to get the last error code is somewhat solver
	///       specific (eg with CPLEX we need to store this explicitly so
	///       it can't be kept for constant methods)
	///       on the other hand Gurobi only really gives a string message
	///       for the last error.
	virtual std::string getErrorMessage(int error=-1) const =0;
    /// inspect solution
    /// get value of current solution (this may be a lower bound if only the
    /// relaxed problem has been solved)
    virtual double getObjective() const = 0;
    /// get the variable value.
    virtual double getPrimal(const Variable &v) const = 0;
    /// get the reduced cost of the variable (only after relaxed solve)
    virtual double getDual(const Variable &v) const = 0;
    /// get the constraint activity level (row * solution_vec)
    virtual double getPrimal(const Constraint &c) const = 0;
    /// get the dual value (only after relaxed solve)
    virtual double getDual(const Constraint &c) const = 0;

	virtual double getObjectiveBound() const = 0;
	/// get an extremal ray in the unbounded case
	virtual void getRay(DblVec &ray) const = 0;
	/// get an extramal ray single value
	/// Note: for CPLEX can only do this so better to get whole ray
	virtual double getRay(const Variable &v) const
	{ DblVec ray; getRay(ray); return ray[v]; } // crude implementation
	  

    virtual void extractSolution ();

    virtual SolutionBase* getSolution () const {return sln;}

    virtual void setParameters (qol::Parameters p) = 0;
	
	/// the following allows directly setting a Gurobi/CPLEX/... parameter
	/// requires that you know what kind of formulation this is
// 	virtual void setNativeIntParam(int param,int value)=0;
// 	virtual void setNativeDblParam(int param,double value)=0;
// 	virtual void setNativeStrParam(int param,std::string value)=0;
// 	virtual int getNativeIntParam(int param)=0;
// 	virtual double getNativeDblParam(int param)=0;
// 	virtual std::string getNativeStrParam(int param)=0;

	/// set a callback function (for cuts/lazy-constraints)
	/// Could also be used for infomation callback but would be quite
	/// inefficient
	/// Installs the callback for the actual solver used
	/// Note: the callback implementation from the actual solver
	///       needs to set the callback->x solution and then
	///       call callback->callback(where) with the right 'where'
	virtual void setCallback(Callback *cb) {callback= cb;}
	virtual Callback *getCallback() { return callback;}
	virtual const Callback *getCallback() const { return callback;}
	
	/// methods for callback to inspect solution
  protected:
	friend class Callback;
	/// get current relaxed solution value at a branch & bound node
	virtual double getNodeBound() const = 0;
	/// get current integer solution value at a branch & bound node
	virtual double getNodeBest() const = 0;
	/// add a global cut (or local cut if isGlobal=false)
	virtual void addCut(const Row &row,bool isGlobal=true)=0;
	
  };


  /// Callback class to enable monitoring and cut generation during branch &
  /// bound.
  class Callback {
  protected:
	/// mip is a pointer to the solver implementation.
	/// This can in principle be used to access some of the native solver
	/// calls (by downcasting to the correct solver type and calling the
	/// the native methods with the env/LP pointers contained here) but that
	/// should only be done as a last resort and if you know what you are
	/// doing.
	MIPSolver *mip;				
  public:
	Callback() : mip(0) {}
	
	/// this method only required by setCallback();
	void setSolver(MIPSolver *solver){mip=solver;}
	
	/// Current fractional solution
	DblVec primalSolution;

	void addCut(const Row &row,bool isGlobal=true) {
	  assert(mip!=0);
	  mip->addCut(row,isGlobal);
	}
	
	enum Status {OK,ABORT};
	/// where are we in the solve process
	/// PRESOLVE (not sure if this ever happens
	/// SIMPLEX - solving root node relaxation
	/// ROOT - root node cut loop
	/// MIP - normal branch & bound node
	/// MIP_SOL - integer feasible branch & bound node (leaf node)
	/// MIP_UNBD - unbounded branch & bound node (leaf node)
	enum Progress {PRESOLVE,SIMPLEX,ROOT,MIP,MIP_SOL,MIP_UNBD};
	/// callback function must be implemented by the user and return either OK
	/// to continue or ABORT to stop the solving process
	/// if cuts is non-empty they will be added.
	virtual Status callback(Progress where) = 0;		
	/// current relaxed bound
	double getBound() const { assert(mip!=0); return mip->getNodeBound(); }
	/// current best feasible solution
	double getBest() const { assert(mip!=0); return mip->getNodeBest(); }
	
  };

  //###
  //### TODO: Use template version to have RowFormulation or RowColFormulation with both descriptions
  //###

  // abstract base class that delegates all of the work of storing the formulation into a ColFormulation object
  // However no actual solving stuff is implemented
  template<class FORMULATION> class FormulationMIPSolver : public MIPSolver {
  public:
	  /// All data is stored in the model
	  FORMULATION model;
	  std::string name;

	  FormulationMIPSolver() {}

	  /// create new formulation, optionally with existing environment
	  FormulationMIPSolver(const qol::QolColFormulation &mip
		  //,std::string name = "gurobi", qol::SolutionType type = qol::PRIMAL
		  ) : model(mip) {}

	  virtual FormulationType getSolverType() const { return qol::QOLFORM; } // no solver

	  /// load data into this object (will load from this->mip by default)
	  virtual void load(const qol::QolColFormulation *mip = 0,std::string _name="") {
		  if (mip) { model.load(mip); }
		  name = _name;
	  }

	  /// is all data loaded into internal memory?
	  virtual bool isLoaded() const { return true; } // no separate load operation

	  virtual Status solveRelaxed() = 0;

	  virtual Status solveExact() = 0;

	  /// error messages: translate solver specific error codes into something
	  /// meaningful
	  /// The actual error code can be provided or will return the message
	  /// of the last error by default
	  /// Note: the ability to get the last error code is somewhat solver
	  ///       specific (eg with CPLEX we need to store this explicitly so
	  ///       it can't be kept for constant methods)
	  virtual std::string getErrorMessage(int error = -1) const {
		  return "UNKNOWN ERROR";
	  }

	  virtual Index nVar() const { return model.nVar(); }

	  virtual Index nConstr() const { return model.nConstr(); }

	  virtual void setNumVar(Index n) { model.setNumVar(n); }
	  virtual void setObjCoeff(qol::Variable v, double obj) { model.setObjCoeff(v, obj); }
	  virtual double getObjCoeff(qol::Variable v) const { return model.getObjCoeff(v); }
	  virtual void setVarUB(qol::Variable v, double ub) { model.setVarUB(v, ub); }
	  virtual double getVarUB(qol::Variable v) const { return model.getVarUB(v); }
	  virtual void setVarLB(qol::Variable v, double lb) { model.setVarLB(v, lb); }
	  virtual double getVarLB(qol::Variable v) const { return model.getVarLB(v); }
	  virtual void setVarType(qol::Variable v, char type) { model.setVarType(v, type); }
	  virtual char getVarType(qol::Variable v) const { return model.getVarType(v); }
	  virtual void setVarName(qol::Variable v, const std::string &name) { model.setVarName(v, name); }
	  virtual const std::string getVarName(qol::Variable v) const { return model.getVarName(v); }
	  // constraint modification
	  virtual ConstraintMIP addConstraint(const Row &row) { return model.addConstraint(row); }
	  virtual void setSense(qol::Constraint cnstr, char s) { model.setSense(cnstr, s); }
	  virtual char getSense(const qol::Constraint cnstr) const { return model.getSense(cnstr); }
	  virtual void setRHS(qol::Constraint cnstr, double r) { model.setRHS(cnstr, r); }
	  virtual double getRHS(qol::Constraint cnstr) const { return model.getRHS(cnstr); }
	  virtual void setConstrName(qol::Constraint cnstr, const std::string &name) { model.setConstrName(cnstr, name); }
	  virtual std::string getConstrName(qol::Constraint cnstr) const {
		  return model.getConstrName(cnstr);}
	  virtual bool relaxedIsExact() { return false; }

	  /// write formulation out to file as plain .lp file
	  virtual void writeLP(const char *filename) { model.writeLP(filename); }

	  // inspect solution
	  /// get value of current solution (this may be a lower bound if only the
	  /// relaxed problem has been solved)
	  virtual double getObjective() const = 0;
	  /// return lower bound in case branch & bound hasn't finished
	  virtual double getObjectiveBound() const = 0;

	  ///<summary>
	  ///  If the LP is unbounded there must be a direction in which the
	  ///  objective keeps improving infinitely. This call gives one such
	  ///  direction
	  ///</summary>
	  ///<return>Extremal ray, a vector of length nVar()</return>
	  ///WARNING: this only works when GRB_INT_PARAM_INFUNBDINFO (InfUnbdInfo)
	  ///         is set to 1 before solving (it's 0 by default)
	  virtual void getRay(DblVec &ray) const { qolAssert(false, "getRay() not implemented"); }
	  /// get an extramal ray single value
	  ///WARNING: this only works when GRB_INT_PARAM_INFUNBDINFO (InfUnbdInfo)
	  ///         is set to 1 before solving (it's 0 by default)
	  virtual double getRay(const Variable &v) const {
		  qolAssert(false, "getRay() not implemented");
		  return 0;
	  }

	  /// get the variable value. 
	  virtual double getPrimal(const Variable &v) const {
		  qolAssert(false, "getPrimal(Variable) for class without a solver");
		  return 0;
	  }

	  /// get the reduced cost of the variable (only after relaxed solve)
	  virtual double getDual(const Variable &v) const {
		  qolAssert(false, "getDual(Variable) for class without a solver");
		  return 0;
	  }

	  /// get the constraint activity level (row * solution_vec)
	  virtual double getPrimal(const Constraint &c) const {
		  qolAssert(false, "getPrimal(Constraint) for class without a solver");
		  return 0.0;
	  }

	  /// get the dual value (only after relaxed solve)
	  virtual double getDual(const Constraint &c) const {
		  qolAssert(false, "getDual(constraint) for class without a solver");
		  return 0;
	  }

	  //virtual void extractSolution (); // defined in QolMIPSolver.h

	  virtual void setParameters(qol::Parameters params) {} // no parameters to set
	  //----------------------- callback related ------------------------
	  /// Install callback - there can only ever be one callback
	  /// (note in future we could allow multiple callbacks to be called
	  /// in sequence)
	  virtual void setCallback(Callback *cb) {
		  qolAssert(false, "setCallback() for solver without callback");
	  }
	protected:
		/// get current relaxed solution value at a branch & bound node
		virtual double getNodeBound() const { qolAssert(false, "getNodeBound() not implemented"); return 0; };
		/// get current integer solution value at a branch & bound node
		virtual double getNodeBest() const { qolAssert(false, "getNodeBest() not implemented"); return 0; }
	  /// add a global cut (or local cut if isGlobal=false)
	  virtual void addCut(const Row &row, bool isGlobal = true)
	  {
		  qolAssert(false, "addCut() for solver without cutting mechanism");
	  }
  }; // end class ColMIPSolver

  
} // end namespace qol

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
