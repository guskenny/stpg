/// Solve a MIP with Gurobi
/// Note: Gurobi's C++ interface appears to be just a wrapper around the C
/// interface which is not always that efficient so this implementation uses
/// the raw C interface directly (same as for the CplexFormulation
/// implementation)
#ifndef __QOL_GUROBIFORMULATION__
#define __QOL_GUROBIFORMULATION__
#ifndef _USE_GUROBI_
  class GurobiFormulation;
#else
extern "C" {
#include <gurobi_c.h>
};
#include <deque>
#include <iostream>
#include <limits>
#include "QolMIPSolver.h"
#include "QolColFormulation.h"


namespace qol {


  class GurobiFormulation : public MIPSolver {
  public:
	/// small helper class to manage a single license for use with multiple
	/// models or GurobiFormulations
	class GurobiLicenseManager {
	  int cnt;					// how many models are registered
	  GRBenv *env;
	public:
	  GurobiLicenseManager() : cnt(0),env(0) {}
	  ~GurobiLicenseManager() { if(env) GRBfreeenv(env); env=0;}
	  GRBenv *grabGurobiEnv()
	  { if(env == 0){
		  int error = GRBloadenv(&env,"gurobi.log");
		  qolAssert(error==0,"No gurobi licence (" << error << ")");
		}
		++cnt; return env; }
	  void releaseGurobiEnv()
	  { if(--cnt == 0){ GRBfreeenv(env); env=0; } }
	}; 
	static GurobiLicenseManager gurobiLicenseManager;
	/// Warning: the env pointer changes once a model is loaded
	/// Initially it points to a global environment, when isLoaded()
	/// env points to the LP's local copy of env!
    GRBenv *env; 
    GRBmodel *LP;
	/// loading directly into GUROBI is slow, so load into mip until we need
	/// to load into gurobi. This pointer is NULL once we have loaded the mip
	/// into gurobi's internal memory
	QolColFormulation *model;
  private:
    /// don't allow copying of a formulation
    GurobiFormulation(const GurobiFormulation &c, qol::SolutionType type=qol::PRIMAL)
	  : MIPSolver (type), env(gurobiLicenseManager.grabGurobiEnv()), model (0)
	{}
    std::vector<char> varType; // column or variable type
  public:
    /// create new formulation, optionally with existing environment
    GurobiFormulation(qol::SolutionType type=qol::PRIMAL)
	  : MIPSolver(type), model(new QolColFormulation)
	{ env = gurobiLicenseManager.grabGurobiEnv();}
	GurobiFormulation(const qol::QolColFormulation &mip,
					  std::string name="gurobi",
					  qol::SolutionType type=qol::PRIMAL);
	virtual ~GurobiFormulation() {	  
	  if(LP) GRBfreemodel(LP); // free's LP->env I think
	  gurobiLicenseManager.releaseGurobiEnv(); // free's global env if required
	  if(model) delete model;
	}
	
	virtual FormulationType getSolverType() const { return GUROBI; }	  
	
	/// load data into gurobi (will load from this->mip by default)
	virtual void load(const qol::QolColFormulation *mip=0,
					  std::string name="gurobi");	
	
	/// is all data loaded into CPLEX's internal memory?
	virtual bool isLoaded() const { return model == 0; } 

	  virtual bool isLP() const {
		  std::vector<char> type(nVar());
		  GRBgetcharattrarray(LP,GRB_CHAR_ATTR_VTYPE,0,nVar(),type.data());
		  for(int v=0;v<nVar();++v)
			  if(type[v] != qol::Variable::CONTINUOUS) return false;
		  return true;
	  }
		  
	  
    virtual Status solveRelaxed() {
		if( ! isLoaded() ) load();
		if(isLP() ){
			if(GRBoptimize(LP)) return FAILED;
		}else {
			// no easy way to do this. Either change variable types or
			// try just solving the root node.
			double oldNodeLim ;
			if(GRBgetdblparam(env,GRB_DBL_PAR_NODELIMIT,&oldNodeLim))
				qolAssert(false,getErrorMessage());
			if(GRBsetdblparam(env,GRB_DBL_PAR_NODELIMIT,1))
				qolAssert(false,getErrorMessage());
			if(GRBupdatemodel(LP))
				qolAssert(false,getErrorMessage());
			// Solve the model with slacks
			/*int error =*/
			if(GRBoptimize(LP)) return FAILED;
			GRBsetdblparam(env,GRB_DBL_PAR_NODELIMIT,oldNodeLim);
		}
		int optimstatus;
		if(GRBgetintattr(LP,GRB_INT_ATTR_STATUS,&optimstatus))
			qolAssert(false,getErrorMessage());

      switch (optimstatus) {
	  case GRB_INFEASIBLE:
	  case GRB_INF_OR_UNBD:
		return INFEASIBLE;
	  case GRB_UNBOUNDED:
		return UNBOUNDED;
	  case GRB_OPTIMAL:
		return OPTIMAL;
	  default:
		return FAILED;
      }
    }

    virtual Status solveExact() {
	  if( ! isLoaded()) load();
	  int error=GRBupdatemodel(LP);
	  if(error) qolAssert(error==0,getErrorMessage());	  
      error = GRBoptimize(LP);
	  if(error) return FAILED;
      int optimstatus;
	  error =GRBgetintattr(LP,GRB_INT_ATTR_STATUS,&optimstatus);
	  qolAssert(error==0,getErrorMessage());	  
      switch (optimstatus) {
        case GRB_INFEASIBLE:
        case GRB_UNBOUNDED:
        case GRB_INF_OR_UNBD:
          return INFEASIBLE;
        case GRB_OPTIMAL:
          return OPTIMAL;
        default:
          return FAILED;
      }
    }

	/// initial primal solution - no need to set this for every variable normally
 	virtual void setPrimalStart(const Variable &v,double val)
		  {  if(! isLoaded() ) load();
			  GRBsetdblattrelement(LP,GRB_DBL_ATTR_START,int(v),val);
		  }
	/// val must be of length number of variables - gives complete solution
	virtual void setPrimalStart(const std::vector<double> &val)
		  { for(Index i=0;i<val.size();++i) setPrimalStart(Variable(i),val[i]);}

	  

	/// error messages: translate solver specific error codes into something
	/// meaningful
	/// The actual error code can be provided or will return the message
	/// of the last error by default
	/// Note: the ability to get the last error code is somewhat solver
	///       specific (eg with CPLEX we need to store this explicitly so
	///       it can't be kept for constant methods)
	virtual std::string getErrorMessage(int error=-1) const;
	  
    virtual Index nVar() const
    { if(model) return model->nVar();
	  int n;
	  if(GRBgetintattr(LP,GRB_INT_ATTR_NUMVARS,&n))
		qolAssert(false,getErrorMessage());		
	  return n;
	}

    virtual Index nConstr() const
    { if(model) return model->nConstr();
	  int n;
	  if(GRBgetintattr(LP,GRB_INT_ATTR_NUMCONSTRS,&n))
		qolAssert(false,getErrorMessage());
	  return n;
	}

    void setNumVar(Index n){
	  if(model) { model->setNumVar(n); return;}
	  int change = n-nVar();
	  int err=0;
	  if(change > 0){
		std::vector<double> ub(n,1.0);
		std::vector<char> vtype(n,GRB_BINARY);
		err=GRBaddvars(LP,change,0, // numnz=0
				   0,0,0, // vbeg,vind,vval = 0
				   0,0,&ub[0],&vtype[0],0); // obj,lb,ub,vtype,varnames
		qolAssert(err==0,getErrorMessage());
		err=GRBupdatemodel(LP);
	  }else if(change < 0) {
		std::vector<int> varstodel(-change);;
		for(int i=0;i<-change;++i) { varstodel[i] = n+i;}
		err=GRBdelvars(LP,-change,&varstodel[0]);
		qolAssert(err==0,getErrorMessage());		
		err=GRBupdatemodel(LP);
	  }
	  if(err) qolAssert(false,getErrorMessage());
    }

    virtual void setObjCoeff(qol::Variable v,double obj) {
	  if(model){ model->setObjCoeff(v,obj); return; }
	  if(GRBsetdblattrelement(LP,GRB_DBL_ATTR_OBJ,v,obj))
		qolAssert(false,getErrorMessage());		
	}

    virtual double getObjCoeff(qol::Variable v) const {
	  if(model) return model->getObjCoeff(v);
	  double obj;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_OBJ,v,&obj))
		qolAssert(false,getErrorMessage());
	  return obj;}

    virtual void setVarUB(qol::Variable v,double ub) {
	  if(model) return model->setVarUB(v,ub);
	  if(GRBsetdblattrelement(LP,GRB_DBL_ATTR_UB,v,ub))
		qolAssert(false,getErrorMessage());
	}

    virtual double getVarUB(qol::Variable v) const {
	  if(model) return model->getVarUB(v);
	  double ub;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_UB,v,&ub))
		qolAssert(false,getErrorMessage());		
      return ub;
	}

    virtual void setVarLB(qol::Variable v,double lb) {
	  if(model) return model->setVarLB(v,lb);
	  if(GRBsetdblattrelement(LP,GRB_DBL_ATTR_LB,v,lb))
		qolAssert(false,getErrorMessage());		
	}
	
	
    virtual double getVarLB(qol::Variable v) const {
	  if(model) return model->getVarLB(v);
	  double lb;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_LB,v,&lb))
		qolAssert(false,getErrorMessage());		
	  return lb;
	}

    virtual void setVarType(qol::Variable v,char type) {	
	  if(model) return model->setVarType(v,type);
      qolAssert(type==qol::Variable::CONTINUOUS ||
				type==qol::Variable::BINARY ||
				type==qol::Variable::INTEGER,
				"Invalid variable type " << type << " for variable " << v);
	  if(GRBsetcharattrelement(LP,GRB_CHAR_ATTR_VTYPE,v,type))
		qolAssert(false,getErrorMessage());		
      //varType[int(v)] = type;
    }
	  
	virtual char getVarType(qol::Variable v) const {
	  if(model) return model->getVarType(v);
	  char type;
	  if(GRBgetcharattrelement(LP,GRB_CHAR_ATTR_VTYPE,v,&type))
		qolAssert(false,getErrorMessage());		
	  return type;
    }

    virtual void setVarName (qol::Variable v, const std::string &name)
	{ if(model) return model->setVarName(v,name);
	  if(GRBsetstrattrelement(LP,GRB_STR_ATTR_VARNAME,v,name.c_str()))
		qolAssert(false,getErrorMessage());		
	}

    virtual const std::string getVarName (qol::Variable v) const
	{ if(model) return model->getVarName(v);
	  char *name;
	  if(GRBgetstrattrelement(LP,GRB_STR_ATTR_VARNAME,v,&name))
		qolAssert(false,getErrorMessage());		
	  return std::string(name);
	}

    // constraint modification
    virtual ConstraintMIP addConstraint(const Row &row) ;

    virtual void setSense(qol::Constraint cnstr,char s) {
	  if(model) return model->setSense(cnstr,s);
      qolAssert(s==qol::Constraint::LE ||
				s==qol::Constraint::EQ ||
				s==qol::Constraint::GE,"Invalid constraint sense " << s);
	  if(GRBsetcharattrelement(LP,GRB_CHAR_ATTR_SENSE,cnstr,s))
		qolAssert(false,getErrorMessage());		
    }

    virtual char getSense(const qol::Constraint cnstr) const {
	  if(model) return model->getSense(cnstr);
	  char s;
	  if(GRBgetcharattrelement(LP,GRB_CHAR_ATTR_SENSE,cnstr,&s))
		qolAssert(false,getErrorMessage());		
	  return s; 
    }

    virtual void setRHS(qol::Constraint cnstr,double r) {
	  if(model) return model->setRHS(cnstr,r);
	  if(GRBsetdblattrelement(LP,GRB_DBL_ATTR_RHS,cnstr,r))
		qolAssert(false,getErrorMessage());		
    }

    virtual double getRHS(qol::Constraint cnstr) const {
	  if(model) return model->getRHS(cnstr);
	  double r;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_RHS,cnstr,&r))
		qolAssert(false,getErrorMessage());
	  return r; 
    }

    virtual void setConstrName (qol::Constraint cnstr, const std::string &name)
	{ if(model) return model->setConstrName(cnstr,name);
	  if(GRBsetstrattrelement(LP,GRB_STR_ATTR_CONSTRNAME,cnstr,name.c_str()))
		qolAssert(false,getErrorMessage());
	}

    virtual std::string getConstrName (qol::Constraint cnstr) const {
	  if(model) return model->getConstrName(cnstr);
	  char *name;
	  if(GRBgetstrattrelement(LP,GRB_STR_ATTR_CONSTRNAME,cnstr,&name))
		qolAssert(false,getErrorMessage());
	  return std::string(name);
	}

    virtual bool relaxedIsExact (){
      return false;
    }

    /// write formulation out to file
    /// cplex will write a range of file types which are infered from the filename
    /// this includes .lp.gz or .lp.bz2 (compressed LP files), .mps/.lp files
    /// and .rmp/.rlp (MPS or LP format with all names changed to generic names)
    virtual void writeLP(const char *filename) {
		if(model) model->writeLP(filename);
		else{
			if(GRBupdatemodel(LP)) qolAssert(false,getErrorMessage());
			if(GRBwrite(LP,filename)) qolAssert(false,getErrorMessage());
		}
    }

    // inspect solution
    /// get value of current solution (this may be a lower bound if only the
    /// relaxed problem has been solved)
    virtual double getObjective() const {
	  double obj;
	  if(GRBgetdblattr(LP,GRB_DBL_ATTR_OBJVAL,&obj))
		throw Exception(getErrorMessage());
	  return obj;
    }

    /// return lower bound in case branch & bound hasn't finished
    virtual double getObjectiveBound() const {
	  double obj;
	  if(GRBgetdblattr(LP,GRB_DBL_ATTR_OBJBOUND,&obj))
		throw Exception(getErrorMessage());
	  return obj;
    }
	  
	///<summary>
	///  If the LP is unbounded there must be a direction in which the
	///  objective keeps improving infinitely. This call gives one such
	///  direction
	///</summary>
	///<return>Extremal ray, a vector of length nVar()</return>
	///WARNING: this only works when GRB_INT_PARAM_INFUNBDINFO (InfUnbdInfo)
	///         is set to 1 before solving (it's 0 by default)
	virtual void getRay(DblVec &ray) const
	{ ray.resize(nVar());
	  if(GRBgetdblattrarray(LP,GRB_DBL_ATTR_UNBDRAY,0,nVar(),&ray[0]))
		qolAssert(false,getErrorMessage());
	}
	/// get an extramal ray single value
	///WARNING: this only works when GRB_INT_PARAM_INFUNBDINFO (InfUnbdInfo)
	///         is set to 1 before solving (it's 0 by default)
	virtual double getRay(const Variable &v) const
	{ double ray;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_UNBDRAY,v,&ray))
		qolAssert(false,getErrorMessage());
	  return ray;
	} // crude implementation

	  
    /// get the variable value. 
    virtual double getPrimal(const Variable &v) const {
	  double x;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_X,v,&x))
		throw Exception(getErrorMessage());		
	  return x;
    }

    /// get the reduced cost of the variable (only after relaxed solve)
    virtual double getDual(const Variable &v) const {
	  double rc;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_RC,v,&rc))
		throw Exception(getErrorMessage());		
	  return rc;
    }

    /// get the constraint activity level (row * solution_vec)
    virtual double getPrimal(const Constraint &c) const {
	  double slack;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_SLACK,c,&slack))
		throw Exception(getErrorMessage());
	  // IS THIS RIGHT??? or is it RHS+slack
	  // not sure how slack is defined in gurobi (RHS-LHS, LHS-RHS or even
	  // worse |RHS-LHS| depending on whether we have <= or >= constraint)??
	  return getRHS(c) - slack; 
    }

    /// get the dual value (only after relaxed solve)
    virtual double getDual(const Constraint &c) const {
	  double pi;
	  if(GRBgetdblattrelement(LP,GRB_DBL_ATTR_PI,c,&pi))
		throw Exception(getErrorMessage());		
	  return pi;
    }

    //virtual void extractSolution (); // defined in QolMIPSolver.h

    virtual void setParameters (qol::Parameters params);
	//----------------------- callback related ------------------------
	/// Install callback - there can only ever be one callback
	/// (note in future we could allow multiple callbacks to be called
	/// in sequence)
	virtual void setCallback(Callback *cb);
	/// add a global cut (or local cut if isGlobal=false)
	virtual void addCut(const Row &row,bool isGlobal=true);
	
	// the following function should really be private but needs to be public
	// based on the way the CPLEX callable library interface is organised
	void setGurobiCallbackInfo(void *data,int wherefrom)
	{ callbackData = data; callbackWherefrom = wherefrom; }
  protected:
	//current values as passed to callback function by cplex
	void *callbackData;			// as provided by CPLEX in callback
	int callbackWherefrom; 		// what CPLEX is currently doing
	virtual double getNodeBound() const
	{ double obj=0;
	  int err=0;
	  switch(callbackWherefrom){
	  case GRB_CB_MIP:
		err=GRBcbget(callbackData,callbackWherefrom,
				 GRB_CB_MIP_OBJBND,&obj);				break;
	  case GRB_CB_MIPSOL:
		err=GRBcbget(callbackData,callbackWherefrom,
				 GRB_CB_MIPSOL_OBJBND,&obj);				break;
	  case GRB_CB_MIPNODE:	
		err=GRBcbget(callbackData,callbackWherefrom,
				 GRB_CB_MIPNODE_OBJBND,&obj);				break;
	  default:
		return -qol::inf;		// no solution available
	  }
	  if(err){
		qolAssert(false,getErrorMessage());
		return -qol::inf;
	  }
	  return obj;
	}
	virtual double getNodeBest() const
	{ double obj;
	  int err=0;
	  switch(callbackWherefrom){
	  case GRB_CB_MIP:
		err=GRBcbget(callbackData,callbackWherefrom,
					 GRB_CB_MIP_OBJBST,&obj);				break;
	  case GRB_CB_MIPSOL:
		err=GRBcbget(callbackData,callbackWherefrom,
					 GRB_CB_MIPSOL_OBJBST,&obj);			break;
	  case GRB_CB_MIPNODE:	
		err=GRBcbget(callbackData,callbackWherefrom,
					 GRB_CB_MIPNODE_OBJBST,&obj);			break;
	  default:
		return qol::inf;		// no solution available
	  }
	  if(err){
		qolAssert(false,getErrorMessage());
		return qol::inf;
	  }
	  return obj;
	}	
  };
}
#endif
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
