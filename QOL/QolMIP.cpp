#include "QolMIP.h"
#include <cstdio>
#include <map>

namespace qol {
  //namespace MIP {
  /*Variable::Variable (MIP *m) : qol::Variable(m->addVar()), model(m) {}
  Variable::Variable(MIP *m,qol::Variable v) : qol::Variable(v), model(m) {}

  Variable &Variable::operator <=(double ub) 
  { model->setVarUB(*this,ub); return *this; }
  Variable &Variable::operator >=(double lb)
  { model->setVarLB(*this,lb); return *this; }	
  Variable &Variable::operator ==(double fix) ///< fix variable (set lower & upper bnd)
  { model->setVarUB(*this,fix); model->setVarLB(*this,fix); return *this;}
  Variable &Variable::operator <=(int ub)
  { model->setVarUB(*this,ub); return *this; }
  Variable &Variable::operator >=(int lb)
  { model->setVarLB(*this,lb); return *this; }	
  Variable &Variable::operator ==(int fix) ///< fix variable (set lower & upper bnd)
  { model->setVarUB(*this,fix); model->setVarLB(*this,fix); return *this;}

  void Variable::setObjCoeff(double cost) 
  { model->setObjCoeff(*this,cost); }

  double Variable::getObjCoeff() const { return model->getObjCoeff(*this); }

  void Variable::setName(const char *fmt,...)
  { // like printf
  char buffer[1028];
  va_list args;
  va_start(args,fmt);
  vsprintf(buffer,fmt,args);
  va_end(args);
  model->setVarName (*this, std::string (buffer));
  }

  void Variable::setName(const std::string &name) {
  model->setVarName (*this, name);
  }

  const std::string Variable::getName() const {
  return (model->getVarName (*this));
  }

  void Variable::setType(char type) { model->setVarType(*this,type); }

  char Variable::getType() const { return model->getVarType(*this); }

  /// make variable integer (default binary)
  void Variable::setInteger() 
  { setType(qol::Variable::INTEGER); }

  /// makes binary & changes bounds to 0/1
  void Variable::setBinary() 
  { 0.0 <= *this; *this <= 1.0; setType(qol::Variable::BINARY); }

  /// make variable continuous (default binary)
  void Variable::setContinuous() { setType(qol::Variable::CONTINUOUS); }

  // solution query methods
  /// solution value of variable (0 if unsolved)
  double Variable::primal() const 
  { return model->getPrimal(*this); }

  /// reduced cost of variable (0 if unsolved)
  double Variable::dual() const 
  { return model->getDual(*this); }

  */

    int Row::setCoeffVec(std::vector<int> &idx,DblVec &val) const
    {
	std::vector<CoeffIdxPair> tmp(lhs); // copy so we can modify
	std::sort(tmp.begin(),tmp.end());
	idx.clear(); val.clear();
	idx.reserve(lhs.size()); val.reserve(lhs.size());
	for(size_t i=0;i<tmp.size();++i)
	    if(idx.empty() || idx.back() != tmp[i].second){
		idx.push_back(tmp[i].second);
		val.push_back(tmp[i].first);
	    } else
		val.back() += tmp[i].first;
	/* -- this should be marginally more efficient but has a bug
	std::map<int,size_t> map; // index to placement
	for(Index i=0;i<lhs.size();++i){
	    size_t &mappedI=map[lhs[i].second];
	    if(i==0) mappedI=map.size()-1;
	}
	idx.resize(map.size());
	for(auto it=map.begin();it!=map.end();++it)
	    idx[it->second]=it->first; // column for each element of map
	val.resize(map.size());
	val=0;
	for(Index i=0;i<lhs.size();++i)
	val[map[lhs[i].second]] += lhs[i].first;
	*/
	return int(idx.size());
    }

  //ConstraintMIP MIP::getConstr(Index i) {return ConstraintMIP(this,i);}
    
  ConstraintMIP::ConstraintMIP(MIP *m,qol::Constraint c) : qol::Constraint(c),model(m) {}

  void ConstraintMIP::setName(const char *fmt,...) { ///< like printf
    char buffer[1028];
    va_list args;
    va_start(args,fmt);
    vsprintf(buffer,fmt,args);
    va_end(args);
    model->setConstrName (*this, std::string (buffer));
  }

  void ConstraintMIP::setName(const std::string &name) {
    model->setConstrName(*this,name); }

    std::string ConstraintMIP::getName() const {
	//std::cout << "ConstraintMIP::getName()..., model = " << model << std::endl;
	//std::cout << "\tgetConstrName() = " << &model->getConstrName << std::endl;
	//std::string name = model->getConstrName(*this);
	//std::cout << "Got name " << name << std::endl;
	return model->getConstrName(*this);
  }

  double ConstraintMIP::getRHS() const {
    return model->getRHS(*this); }
  void ConstraintMIP::setRHS(double rhs)  {
	  return model->setRHS(*this,rhs); }

  char ConstraintMIP::getSense() const {
    return model->getSense(*this); }

  
  // solution query methods
  /*  
  double ConstraintMIP::getPrimal() const  /// activity level (row * solution_vec)
  { return model->getPrimal(*this); }
  double ConstraintMIP::getDual() const // dual variable value (row price)
  { return model->getDual(*this); }
  */
  //}
}

