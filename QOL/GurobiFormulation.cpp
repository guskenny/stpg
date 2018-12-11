#include "GurobiFormulation.h"
#ifdef _USE_GUROBI_
namespace qol {

	// static/global object for holding a gurobi license/environment:
	GurobiFormulation::GurobiLicenseManager GurobiFormulation::gurobiLicenseManager;
	
	std::string GurobiFormulation::getErrorMessage(int error) const
	{
		if(error == -1){
			if(env) return GRBgeterrormsg(env);
		}
		// No generic method for looking up any error message
		// This is a basic hacked version
		switch(error){
			case 0: return "";
			case GRB_ERROR_OUT_OF_MEMORY: return "OUT_OF_MEMORY";
			case GRB_ERROR_NULL_ARGUMENT: return "NULL_ARGUMENT";
			case GRB_ERROR_INVALID_ARGUMENT: return "INVALID_ARGUMENT";
			case GRB_ERROR_UNKNOWN_ATTRIBUTE: return "UNKNOWN_ATTRIBUTE";
			case GRB_ERROR_DATA_NOT_AVAILABLE: return "DATA_NOT_AVAILABLE";
			case GRB_ERROR_INDEX_OUT_OF_RANGE: return "INDEX_OUT_OF_RANGE";
			case GRB_ERROR_UNKNOWN_PARAMETER: return "UNKNOWN_PARAMETER";
			case GRB_ERROR_VALUE_OUT_OF_RANGE: return "VALUE_OUT_OF_RANGE";
			case GRB_ERROR_NO_LICENSE: return "NO_LICENSE";
			case GRB_ERROR_SIZE_LIMIT_EXCEEDED: return "SIZE_LIMIT_EXCEEDED";
			case GRB_ERROR_CALLBACK: return "CALLBACK";
			case GRB_ERROR_FILE_READ: return "FILE_READ";
			case GRB_ERROR_FILE_WRITE: return "FILE_WRITE";
			case GRB_ERROR_NUMERIC: return "NUMERIC";
			case GRB_ERROR_IIS_NOT_INFEASIBLE: return "IIS_NOT_INFEASIBLE";
			case GRB_ERROR_NOT_FOR_MIP: return "NOT_FOR_MIP";
			case GRB_ERROR_OPTIMIZATION_IN_PROGRESS: return "OPTIMIZATION_IN_PROGRESS";
			case GRB_ERROR_DUPLICATES: return "DUPLICATES";
			case GRB_ERROR_NODEFILE: return "NODEFILE";
			case GRB_ERROR_Q_NOT_PSD: return "Q_NOT_PSD";
			case GRB_ERROR_QCP_EQUALITY_CONSTRAINT: return "QCP_EQUALITY_CONSTRAINT";
			case GRB_ERROR_NETWORK: return "NETWORK";
			case GRB_ERROR_JOB_REJECTED: return "JOB_REJECTED";
			case GRB_ERROR_NOT_SUPPORTED: return "NOT_SUPPORTED";
			default:{
				std::stringstream tmp;
				tmp << "Unknown Gurobi error " << error;
				return tmp.str();
			}
		}
	}
	
  GurobiFormulation::GurobiFormulation(const qol::QolColFormulation &mip,
											std::string name,
											qol::SolutionType type)
	: MIPSolver(type), env(0),LP(0),model(0)
  {
	  env = gurobiLicenseManager.grabGurobiEnv();
	  qolAssert(env!=0,"GRBloadenv failed: " << getErrorMessage());
	  load(&mip);
  }
  
  void GurobiFormulation::load(const qol::QolColFormulation *mip,
							   std::string name)
  {
	if(mip==0){
	  if(model == 0) return;	// nothing to load
	  mip = model;
	}
    Index nz=0;
    const Index matrixSize = mip->matrix.size();
    const Index numVars = mip->nVar();
    const Index numConstr = mip->nConstr();

    for(Index i=0;i<matrixSize;++i) nz+=mip->matrix[i].size();

    std::vector<int> matbeg(numVars+1),matcnt(numVars),matind(nz);
    DblVec matval(nz);
    nz = 0;
    for(Index i=0;i<matrixSize;++i){
      matbeg[i]=int(nz);
      const Index matrixI = mip->matrix[i].size();
      for(Index j=0;j<matrixI;++j){
        if( mip->matrix[i][j].first == 0.0) continue;
        matval[nz] = mip->matrix[i][j].first;
        matind[nz] = int(mip->matrix[i][j].second);
        ++nz;
      }
      matcnt[i] = int(nz-matbeg[i]);
    }
    matbeg[mip->nVar()] = int(nz); // needed by some solvers but perhaps not by CPLEX
    std::vector<char *> colname(numVars,0),rowname(numConstr,0);
    for(Index i=0;i<mip->colname.size();++i)
      colname[i] = (char*) mip->colname[i].c_str();
    for(Index i=0;i<mip->rowname.size();++i)
	  rowname[i] = (char*) mip->rowname[i].c_str();
	// gurobi's declaration of the load function is a bit stuffed in that
	// it fails to declar all of the pointers to be const (even though
	// everything bar the LP is input/read-only for this function)
	// The cast's below try to shut up some compiler warnings
	int error =
		GRBloadmodel(env,&LP,name.c_str(),mip->nVar(),mip->nConstr(),
				 GRB_MINIMIZE,	// always minimising (should add senese to MIP)
				 0.0,(double *)&mip->cost[0], // objective constant & coefficients
				 (char *)&mip->sense[0],(double *)&mip->rhs[0],
				 &matbeg[0],&matcnt[0],
				 &matind[0],&matval[0],
				 (double *)&mip->collb[0],(double *)&mip->colub[0],
				 (char *)&mip->coltype[0],
				 &colname[0],&rowname[0]);				 
	qolAssert(error==0,getErrorMessage());
	if(model){					// indicate that model is loaded
	  delete model;
	  model = 0;
	}
	// LP has it's own copy of env once loaded
	env = GRBgetenv(LP);
  }
	
  ConstraintMIP GurobiFormulation::addConstraint(const Row &row)
  {
      if(model){model->addConstraint(row); return getConstr(model->nConstr()-1); }
      std::vector<int> cind;
      DblVec cval;
      /*
	cind.reserve(row.lhs.size());
	cval.reserve(row.lhs.size());
	for(Index i=0;i<row.lhs.size();++i){
	  if (fabs(row.lhs[i].first)<=eps) continue;
	  cind.push_back(row.lhs[i].first);
	  cval.push_back(row.lhs[i].second);
	  }*/
      row.setCoeffVec(cind,cval);
      //std::cout << "setCoeffVec got " << cind.size() << " values\n";
      //for(size_t i=0;i<cind.size();++i)
      //  std::cout <<"\t"<<cind[i] << ": " << cval[i] << std::endl;
      int err=GRBaddconstr(LP,cind.size(),&cind[0],&cval[0],row.op,row.rhs,0);
      //std::cout << "Error code = " << err << std::endl;
      qolAssert(err==0,getErrorMessage());
      err=GRBupdatemodel(LP); // is this a good idea? perhaps only do occasionally
      qolAssert(err==0,getErrorMessage());
      
      return getConstr(this->nConstr()-1);
  }
  

  void GurobiFormulation::setParameters (qol::Parameters params) {
    double val;
	int err=0;

    val = params.getParamValue (VERBOSITY);
    if (val!=InvalidParam){
	  err = GRBsetintparam (env, GRB_INT_PAR_OUTPUTFLAG, std::min(1,int(val)));
	  qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (THREADS);
    if (val!=InvalidParam){
      err = GRBsetintparam(env, GRB_INT_PAR_THREADS, int(val));
	  qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (PRINTFREQ);
    if (val!=InvalidParam){
      err = GRBsetintparam(env,GRB_INT_PAR_DISPLAYINTERVAL,int(val));
	  qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (ABSGAP);
    if (val!=InvalidParam){
	  err = GRBsetdblparam(env,GRB_DBL_PAR_MIPGAPABS,val);
	  qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (EPS);
    if (val!=InvalidParam){
	  err = GRBsetdblparam(env,GRB_DBL_PAR_OPTIMALITYTOL,val);
	  qolAssert(err==0,getErrorMessage());
	}

    val = params.getParamValue (ITLIMIT);
    if (val!=InvalidParam){
	  err = GRBsetdblparam(env,GRB_DBL_PAR_ITERATIONLIMIT,val);
	  qolAssert(err==0,getErrorMessage());
	}	  

    val = params.getParamValue (TIMELIMIT);
    if (val!=InvalidParam){
	  err = GRBsetdblparam(env,GRB_DBL_PAR_TIMELIMIT,val);
	  qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (CPUTIME);
    if (val!=InvalidParam){
		// actually can't limit CPU time - at least limit walltime
		err = GRBsetdblparam(env,GRB_DBL_PAR_TIMELIMIT,val);
		qolAssert(err==0,getErrorMessage());
	}
    val = params.getParamValue (WALLTIME);
    if (val!=InvalidParam){
	if(params.getParamValue(THREADS) != InvalidParam)
	    val *= params.getParamValue(THREADS);
	else
	    val *= 32; // guess how many threads we might have
	err = GRBsetdblparam(env,GRB_DBL_PAR_TIMELIMIT,val);
	qolAssert(err==0,getErrorMessage());
    }
	if(err); // just so we don't get warnings about unused variables
  }

  //========================================================================
  // CALLBACKS
  //========================================================================
# if ! defined(_MSC_VER) || (_MSC_VER < 12)
  extern "C" {					// give this function C linkage
	  int GurobiGenericCallback(GRBmodel *LP, void *cbdata, int where,
		  void *userdata);
  }
  int GurobiGenericCallback(GRBmodel *LP, void *cbdata, int wherefrom, void *gurobiFormulation)
# else
  int __stdcall
	  GurobiGenericCallback(GRBmodel *LP, void *cbdata, int wherefrom, void *gurobiFormulation)
# endif
  {
	try{
	  GurobiFormulation *mip=(GurobiFormulation *)gurobiFormulation;
	  mip->setGurobiCallbackInfo(cbdata,wherefrom);
	  Callback::Progress progress;
	  Callback *callback = mip->getCallback();
	  if(callback == 0) return 2; // should never happen
	  if(callback->primalSolution.size() == 0) return 2; // should never happen
	  int err=0;
	  switch (wherefrom) {
	  case GRB_CB_POLLING: return 0; // ignore polling
	  case GRB_CB_PRESOLVE: return 0; // ignore
	  case GRB_CB_SIMPLEX: return 0;	// root node solve
	  case GRB_CB_BARRIER:	   return 0;	// root node solve
	  case GRB_CB_MESSAGE:	   return 0;	// what?
	  case GRB_CB_MIP:  // what does MIP mean as opposed to MIPNODE or MIPSOL?
		return 0; //probably best to ingore as we can't add cuts/lazyconstraints
		//progress = Callback::ROOT; // not sure if this is right??
	  case GRB_CB_MIPNODE: {
		int cbStatus;
		err=GRBcbget(cbdata,wherefrom,GRB_CB_MIPNODE_STATUS,&cbStatus);
		qolAssert(err==0,mip->getErrorMessage());		
		bool haveSoln=false;
		switch(cbStatus){
		case GRB_OPTIMAL:
		case GRB_SUBOPTIMAL:
		  haveSoln=true;
		case GRB_INFEASIBLE:
		case GRB_INF_OR_UNBD:	// no dual ray. set DualReductions=0 to work
								// out which it is
		  progress = Callback::MIP;
		case GRB_UNBOUNDED:
		  progress = Callback::MIP_UNBD;
		default:
		  progress = Callback::MIP; // hit some limit, interupted or whatever
		}
		if(haveSoln)
		  err=GRBcbget(cbdata,wherefrom,GRB_CB_MIPNODE_REL,
					   &callback->primalSolution[0]);
		break;}
	  case GRB_CB_MIPSOL:
		progress = Callback::MIP_SOL;		
		err=GRBcbget(cbdata,wherefrom,GRB_CB_MIPSOL_SOL,
					 &callback->primalSolution[0]);
		break;
	  default:
		qolAssert(false,"Unhandled 'where' = " << wherefrom
				  << " in gurobi callback" );
		return GRB_ERROR_CALLBACK;
	  }
	  if(err) qolAssert(err==0,mip->getErrorMessage());		
	  if(callback->callback(progress) == Callback::ABORT)
		  return GRB_ERROR_CALLBACK;
	}catch(qol::Exception error){
		throw error;
	}catch(...){ // is this a good idea? maybe should just pass exception on
	  return GRB_ERROR_CALLBACK;					// failed;
	}
	return 0;					// OK
  }


  void GurobiFormulation::setCallback(Callback *cb)
  {
	callback = cb;
	callback->setSolver(this);
	callback->primalSolution.resize(nVar(),0.0);	// make sure we have memory for solution
	if(!isLoaded())
		load();							// can't set callback without LP
	int error = GRBsetcallbackfunc(LP,&GurobiGenericCallback,this);
	qolAssert(error==0,getErrorMessage(error));
  }

  void GurobiFormulation::addCut(const Row &row,bool isGlobal)
  {
	std::vector<int> ind;
	DblVec val;
	ind.reserve(row.lhs.size());
	val.reserve(row.lhs.size());
	for(size_t i=0;i<row.lhs.size();++i)
	  if(row.lhs[i].first < -eps || row.lhs[i].first > eps){
		val.push_back(row.lhs[i].first);
		ind.push_back(row.lhs[i].second);
	  }
	if(val.empty()) return;		// can't add
	qolAssert(callbackWherefrom==GRB_CB_MIPNODE ||
			  callbackWherefrom==GRB_CB_MIPSOL,
			  "Gurobi cannot add cut when where is " << callbackWherefrom
			  << " only GRB_CB_MIPNODE or GRB_CB_MIPSOL allowed");
	qolAssert(isGlobal,"Gurobi cannot accept local cuts");
	int err=0;
	// not really clear what the difference is between cut & lazy
	// (cut perhaps always just add a cut?)
	if(callbackWherefrom==GRB_CB_MIPSOL){
		err=GRBcblazy(callbackData,(int)val.size(),&ind[0],&val[0],row.op,row.rhs);
	}else// if(callbackWherefrom == GRB_CB_MIPNODE)
		err=GRBcbcut(callbackData,(int)val.size(),&ind[0],&val[0],row.op,row.rhs);
	if(err)
		qolAssert(false,getErrorMessage());
  }
	

} //namespace
#endif
