// Solve a MIP with CPLEX.
// this module uses the basic CPLEX C API because
// (a) there are less linker issues with different variants of C++ compilers
// (b) we already have most of the Concert functionality provided by QolFormulation
#include "CplexFormulation.h"
#include "QolUtil.h"
#ifdef _USE_CPLEX_

namespace qol {
	
	
  // create a formulation by copying all of the data in the Formulation
  // into CPLEX's internal solver memory

  CplexFormulation::CplexFormulation (const qol::QolColFormulation &mip, std::string name,
                        CPXENVptr e, qol::SolutionType type) : MIPSolver(type)
  {
    int status;
    if(e)
      env = e;
    else
      env=CPXopenCPLEX(&status);
	LP=0;
    //if(env) LP=CPXcreateprob(env,&status,name.c_str());
    //if(!LP) return;
	model=0;
	load(&mip,name);
  }
  void CplexFormulation::load(const qol::QolColFormulation *mip,
							  std::string name)
  {
	if(mip==0){
	  if(model == 0) return;	// nothing to load
	  mip = model;
	}
	int status;
	if(!env)  // should never happen really but just in case
	  env=CPXopenCPLEX(&status);
	LP=CPXcreateprob(env,&status,name.c_str());
	
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
    CPXcopylpwnames(env,LP,mip->nVar(),mip->nConstr(),
					CPX_MIN, //always minimising to keep things simple for now
					&mip->cost[0],&mip->rhs[0],&mip->sense[0],&matbeg[0],&matcnt[0],
					&matind[0],&matval[0],&mip->collb[0],&mip->colub[0],
					0, // no ranged constraints allowed (in current implementation)
					&colname[0],&rowname[0]);
	//CPXcopyctype(env,LP,&mip->coltype[0]);
	varType = mip->coltype;		// do we really need this??
	CPXcopyctype(env,LP,&varType[0]);
	if(model){					// indicate that model is loaded
	  delete model;
	  model = 0;
	}
  } // end constructor converting Formulation to CplexFormulation

	

  void CplexFormulation::setNumVar(Index n)
  {   //< resize to given number of variables
	  if(model) {
		  model->setNumVar(n);
		  return;
	  }
    if(n == nVar()) return;
    varType.resize(n, Variable::BINARY);
    if(n < nVar()){
      CPXdelcols(env,LP,int(n),int(nVar())-1);
      return;
    }
    Index nCol = n-nVar(); // number of new columns
    DblVec obj(nCol,0.0),lb(nCol,0.0),ub(nCol,1.0),cmatval(1,0.0);
    std::vector<int> cmatbeg(nCol,0), cmatind(1,0), newCols(nCol);
    CPXaddcols(env,LP,nCol,0, // no coefficients
      &obj[0],&cmatbeg[0],&cmatind[0],&cmatval[0],&lb[0],&ub[0],
      0); // default names created if required
    //std::vector<char> ctype(nCol,qol::Variable::BINARY);
    //for(int i=0;i<nCols;++i) newCols[i] = i+nVar();
    //CPXchgctype(env,LP,nCol,&newCols[0],&ctype[0]);
    //matrix.resize(n);
  }

  ConstraintMIP CplexFormulation::addConstraint(const Row &row) {
      if(model){ model->addConstraint(row); return getConstr(model->nConstr()-1);}
      std::vector<int> rmatbeg(2,0),rmatcnt(1),rmatind;//(row.lhs.size());
      DblVec rmatval;//(row.lhs.size());
      int nz=row.setCoeffVec(rmatind,rmatval); // originally = 0
    /*
    for(Index i=0;i<row.lhs.size();++i){
      if (fabs(row.lhs[i].first)<=eps) continue;
      //matrix[row.lhs[i].second].push_back (CoeffIdxPair (row.lhs[i].first, rhs.size()));
      //                if(row.lhs[i].first == 0.0) continue;
      rmatval[nz] = row.lhs[i].first;
      rmatind[nz++] = int(row.lhs[i].second);
      }
    sense.push_back (row.op);
    rhs.push_back (row.rhs);
    */
    rmatcnt[0] = rmatbeg[1] = nz;
    CPXaddrows(env,LP,0, // no new columns
      1, // one new row
      nz,&row.rhs,&row.op,&rmatbeg[0],&rmatind[0],&rmatval[0],
	       0,0); // no row or column names
    return getConstr (this->nConstr()-1);
  }

    std::string CplexFormulation::getConstrName (qol::Constraint cnstr) const {
	if(model) return model->getConstrName(cnstr);
	int size;
	int lastErrorCode=CPXgetrowname(env, LP, 0, 0, 0, &size, cnstr, cnstr);
	if (size == 0) return std::string("");
	//if(lastErrorCode) throw Exception(getErrorMessage(lastErrorCode));

	size *= -1;
	std::vector<char> buf(size+1);
	char *cname[2];
	int tmp;
	lastErrorCode=CPXgetrowname(env, LP, cname, &buf.front(), size,
				    &tmp, cnstr, cnstr);
	if(lastErrorCode) throw Exception(getErrorMessage(lastErrorCode));
	return std::string(cname[0]);
    }

    void CplexFormulation::setPrimalStart(const Variable &v,double val)
    {   int vIdx=v;
	if((int)primalStartIdx.size() > vIdx && primalStartIdx[vIdx] == vIdx)
	    primalStartVal[vIdx]=val;
	else{
	    primalStartIdx.push_back(vIdx);
	    primalStartVal.push_back(val);
	}
    }

    void CplexFormulation::setPrimalStart(const std::vector<double> &val)
    {
	qolAssert(nVar()==(int)val.size(),"Incorrect size for primal start vector");
	primalStartIdx.resize(nVar());
	for(int i=0;i<nVar();++i) primalStartIdx[i]=i;
	primalStartVal = val;	    
    }


    
  std::string CplexFormulation::getErrorMessage(int error) const
  {
	char buffer[CPXMESSAGEBUFSIZE];
	const char *errstr;
	if(error == -1) error = lastErrorCode;
	if(error == 0) return "";	// no error
	errstr = CPXgeterrorstring(env,error,buffer);
	if(errstr == NULL)
	  sprintf(buffer,"Unknown CPLEX error %d",error);
	return std::string(buffer);
  }

//   void CplexFormulation::extractSolution () {
//     if (sln) delete(sln);
//     if (_type==qol::FULL)
//       sln = new qol::SolutionFull();
//     else
//       sln = new qol::SolutionPrimal();
//     sln->setObjective (getObjective ());
//     const Index numVars = nVar();
//     for (Index v=0; v<numVars; ++v) {
//       if (_type==qol::PRIMAL) {
//         ((SolutionPrimal*)sln)->push_back (getPrimal (qol::Variable(v)));
//       }
//       else {
//         qol::ExtendedVariable variable (qol::Variable(v), getVarName(qol::Variable(v)),
//           getPrimal(qol::Variable(v)), getDual(qol::Variable(v)));
//         ((SolutionFull*)sln)->push_back (variable);
//       }
//     }
//   }

  void CplexFormulation::setParameters (qol::Parameters params) {
    double val;
    val = params.getParamValue (VERBOSITY);
    if (val!=InvalidParam){
	  lastErrorCode = CPXsetintparam (env, CPX_PARAM_SCRIND, std::min(1,int(val)));
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
    val = params.getParamValue (THREADS);
    if (val!=InvalidParam){
      lastErrorCode = CPXsetintparam (env, CPX_PARAM_THREADS, int(val));
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
    val = params.getParamValue (PRINTFREQ);
    if (val!=InvalidParam){
      lastErrorCode = CPXsetlongparam (env, CPX_PARAM_MIPINTERVAL, long(val));
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
	
    val = params.getParamValue (ABSGAP);
    if (val!=InvalidParam){
      lastErrorCode=CPXsetdblparam (env, CPX_PARAM_EPAGAP, val);
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));	  
	}
    val = params.getParamValue (RELGAP);
    if (val!=InvalidParam){
      lastErrorCode=CPXsetdblparam (env, CPX_PARAM_EPGAP, val);
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));	  
	}
    val = params.getParamValue (EPINT);
    if (val!=InvalidParam){
      lastErrorCode=CPXsetdblparam (env, CPX_PARAM_EPINT, val);
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}	  

    val = params.getParamValue (ITLIMIT);
    if (val!=InvalidParam){
      lastErrorCode = CPXsetintparam (env, CPX_PARAM_STRONGITLIM
									  /*CPX_PARAM_ITLIMIT*/, int(val));
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
	
    val = params.getParamValue (TIMELIMIT);
    if (val!=InvalidParam){
		lastErrorCode = CPXsetdblparam (env, CPX_PARAM_TILIM, val);	
		qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
    val = params.getParamValue (CPUTIME);
	if (val != InvalidParam ){
		lastErrorCode = CPXsetintparam(env,CPX_PARAM_CLOCKTYPE,1);
		lastErrorCode = CPXsetdblparam(env,CPX_PARAM_TILIM, val);
		qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}		
    val = params.getParamValue (WALLTIME);
	if (val != InvalidParam ){
		lastErrorCode = CPXsetintparam(env,CPX_PARAM_CLOCKTYPE,2);
		lastErrorCode = CPXsetdblparam(env,CPX_PARAM_TILIM, val);
		qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}		
		
    if (val!=InvalidParam){
		lastErrorCode = CPXsetdblparam (env, CPX_PARAM_TILIM, val);	
		qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
    val = params.getParamValue (HEURISTIC_FREQ);
    if (val!=InvalidParam){
		lastErrorCode = CPXsetlongparam (env, CPX_PARAM_HEURFREQ, long(val));	
		qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	}
  }

  //========================================================================
  // CALLBACKS
  //========================================================================
  //extern "C" {					// give this function C linkage
  //	int CplexGenericCallback(CALLBACK_CUT_ARGS);
  //}
  int CPXPUBLIC CplexGenericCallback(CPXCENVptr xenv, void *cbdata,
							 int wherefrom, void *cplexFormulation,
							 int *useraction_p)
  {
	try{
	  CplexFormulation *mip=(CplexFormulation *)cplexFormulation;
	  mip->setCplexCallbackInfo(cbdata,wherefrom);
	  Callback::Progress progress;
	  //std::cout << "Cplex cutcallback - " << wherefrom << std::endl;
	  switch (wherefrom) {
      case CPX_CALLBACK_MIP_CUT_FEAS:
		progress = Callback::MIP_SOL; 
		break;
	  case CPX_CALLBACK_MIP_CUT_UNBD:
		progress = Callback::MIP_UNBD;
		break;
      case CPX_CALLBACK_MIP_CUT_LAST: 
		progress = Callback::MIP; // end of node's cut loop
		break;
      case CPX_CALLBACK_MIP_CUT_LOOP: // middle of cut loop
		//progess = Callback::ROOT; // could be any node's cut loop 
		return 0;				// just wait for CUT_LAST
      default:
		return 0;				// should never happen(?)
	  }
	  Callback *callback = mip->getCallback();
	  if(callback == 0) return 2; // should never happen
	  if(callback->primalSolution.size() == 0) return 2; // should never happen
	  int status =
	  CPXgetcallbacknodex(mip->env,cbdata,wherefrom,
						  &callback->primalSolution[0],0,
						  callback->primalSolution.size()-1);
	  qolAssert(status==0,mip->getErrorMessage(status));
	  *useraction_p = CPX_CALLBACK_DEFAULT;
	  if(callback->callback(progress) == Callback::ABORT)
		*useraction_p = CPX_CALLBACK_FAIL;
	}catch(...){
	  return 1;					// failed;
	}
	return 0;					// OK
  }
  void CplexFormulation::setCallback(Callback *cb)
  {
	callback = cb;
	callback->setSolver(this);
	callback->primalSolution.resize(nVar(),0.0);	// make sure we have memory for solution
#   if CPX_VERSION < 12030000
	lastErrorCode=CPXsetusercutcallbackfunc(env,&CplexGenericCallback,this);
#   else	
	// from about CPLEXverison 12.3 onwards there are two different cutbacks
	lastErrorCode=CPXsetlazyconstraintcallbackfunc(
					env,&CplexGenericCallback,this);
	qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));	
	lastErrorCode=CPXsetusercutcallbackfunc(env,&CplexGenericCallback,this);
#   endif
	qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));	
  }

  void CplexFormulation::addCut(const Row &row,bool isGlobal)
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
	if( ! isGlobal){
	  lastErrorCode=CPXcutcallbackaddlocal(env,callbackData,callbackWherefrom,
										   val.size(),row.rhs,row.op,&ind[0],&val[0]);
	  qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));
	  return;
	}
	int usecut;
	switch(callbackWherefrom){
	case CPX_CALLBACK_MIP_CUT_FEAS:
	case CPX_CALLBACK_MIP_CUT_UNBD:
	  usecut = CPX_USECUT_FORCE; break;
	case CPX_CALLBACK_MIP_CUT_LOOP: // won't happen
	  usecut = CPX_USECUT_FILTER; break; // allow CPLEX to arbitrarily discard
	default:
	  usecut = CPX_USECUT_PURGE; // CPLEX may purge in future iterations
	}	  
	lastErrorCode=CPXcutcallbackadd(env,callbackData,callbackWherefrom,
									val.size(),row.rhs,row.op,&ind[0],&val[0],
									usecut);
	qolAssert(lastErrorCode==0,getErrorMessage(lastErrorCode));	
  }



  
} // end namespace qol

#endif
