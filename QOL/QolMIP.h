#ifndef __QOL_ABSTRACTFORMULATION_H__
#define __QOL_ABSTRACTFORMULATION_H__

//#include "AligmentAllocator.h" // no longer used
#include "QolUtil.h"
#include "QolData.h"
#include <utility> // for definition of pair
#include <algorithm>
#include <memory> // std::auto_ptr
#include <assert.h>
#include <cstdarg>
#include <cmath>

#include <iostream>
#include <sstream>

/// Wrapper code to set up an MILP with an arbitrary solver using the OSI
/// interface. 
/// Basic usage:
/// \code
/// using namespace qol
/// Formulation m;
/// Variable x(m),y(m),z(m);
/// x >= 0; y >= 0; z >= 0; // set lower bounds
/// x <= 1; y <= 2; z <= 1; // set upper bounds
/// x.setObjCoeff(2); y.setObjCoeff(-1.0);
/// m.addConstraint( x + 2*y == 2.5 )
/// Expression sum;
/// for(i=0;i<3;++i) sum += m.getVar(i)
/// m.addConstraint( sum <= 2 );
/// m.loadProblem();
/// \endcode
/// This will solve a problem of the form
/// \f[ \min 2 x - y \f]
/// Subject to
/// \f$  x+2\times y = 2.5\f$ \f[ \sum_{var\in\{x,y,z\}} var \le 2 \f]
/* the following doxygen markup only works in LaTeX not in html
\f{eqnarray*}{
x+2y &=& 2.5\\ \sum_{var\in\{x,y,z\}} var &\le& 2 
\f}
*/

namespace qol {

  /// Constraint defines the index of a constraint (or a dual variable)
  /// For a complete definition of the coefficients we need a Row object.
  class MIP;

  class ConstraintMIP : public qol::Constraint {
    MIP *model;
  public:
    //const int idx; ///< numbered from 0 in Model
    ConstraintMIP() : model(0) {}
    /// Construct a Constraint from a generic constraint
    ConstraintMIP(MIP *m,qol::Constraint c) ;//: qol::Constraint(c),model(m) {}
    ConstraintMIP(const ConstraintMIP &c) : qol::Constraint(c),model(c.model) {}

    void setName(const char *fmt,...); ///< like printf
    void setName(const std::string &name);
    std::string getName() const;
    double getRHS() const;
    void setRHS(double val);	// allow changing of RHS
	char getSense() const;
	const MIP *getModel() const {return model;}
	MIP *getModel() {return model;}  

	// solution query methods via MIPSolver object:
	//         mipSolver->getPrimal(constraint)
	// as model may not have a solver attached so can't provide solution info
    //double getPrimal() const; /// activity level (row * solution_vec)
    //{ return model->getPrimal(*this); }
    //double getDual() const; // dual variable value (row price)
    //{ return model->getDual(*this); }
  };

  class CoeffIdxPair { //: public std::pair<double,int> {
  public:
    double first; Index second;
    CoeffIdxPair() : first(0),second(0) {}
	CoeffIdxPair(double coef,Index idx) : first(coef),second(idx) {}
	CoeffIdxPair operator -() const { return CoeffIdxPair(-first,second);}  
    bool operator < (const CoeffIdxPair &other) const
    { return (second < other.second) ||
    (second==other.second && fabs(first) > fabs(other.first) ); }
  };

  /// Row is a helper class only. Not intended to be used directly.
  struct Row {
    std::vector<CoeffIdxPair> lhs;
    double rhs;
    char op;  ///< 'L', 'E', or 'G'
    Row() : rhs(0.0),op('\0') {lhs.clear();}
    Row(const std::vector<CoeffIdxPair> &LHS,char OP,double RHS) :
	  lhs(LHS), rhs(RHS), op(OP) {}
	/// create a row from coefficient index pairs on both left & right hand
	/// side. Either of these may be null for the general case
	Row(const std::vector<CoeffIdxPair> *LHS,
		const std::vector<CoeffIdxPair> *RHS,
		char OP,double RHSconst) : rhs(RHSconst),op(OP) {
	  size_t size=0;
	  if(LHS) size += LHS->size();
	  if(RHS) size += RHS->size();
	  lhs.resize(size);
	  std::vector<CoeffIdxPair>::iterator lhsIt = lhs.begin();
	  if(LHS)
		lhsIt = std::copy(LHS->begin(),LHS->end(),lhsIt);
	  if(RHS){
		for(std::vector<CoeffIdxPair>::const_iterator it=RHS->begin();
			it != RHS->end();++it,++lhsIt)
		  *lhsIt = -(*it);
	  }
	}
    Row(ConstraintRow c) {
      lhs.clear();
      for (size_t i=0; i<c.size(); ++i) {
        lhs.push_back (CoeffIdxPair (c[i].first, c[i].second));
      }
      op = c.sense;
      rhs = c.rhs;
    }
	// convert to a (cleaned-up) list of indices & values (adds duplicate
	// variable entries but doesn't guarantee sorted variable indices)
	  int setCoeffVec(std::vector<int> &idx,DblVec &val) const;
  };
	inline std::ostream &operator<<(std::ostream &os,const Row &row)
	{   for(size_t i=0;i<row.lhs.size();++i){
			os << " + " << row.lhs[i].first << " x" << row.lhs[i].second;
		}
		switch(row.op){
		case 'L': os << " <= "; break;
		case 'E': os << " == "; break;
		case 'G': os << " >= "; break;
		default: os << " ?? ";
		}
		os << row.rhs;
		return os;
	}

  /// An expression is the sum of coefficients * variables
  class Expression {
	double constant;			// constant part of expression
	std::vector<CoeffIdxPair> terms;
  public:
    Expression(double val=0.0) : constant(val) {}
    Expression(const qol::Variable &v): constant(0) { terms.push_back(CoeffIdxPair(1.0,int(v))); }
    Expression(double coef,const qol::Variable &v) : constant(0)
    { terms.push_back(CoeffIdxPair(coef,int(v))); }
	Expression(const CoeffIdxPair &c) : constant(0) {terms.push_back(c);}
    Expression(const Expression &e) : constant(e.constant), terms(e.terms) {}
    Expression(size_t size) : constant(0.0){
    	reserve(size);
    }
    void reserve(size_t size) { 
    	terms.reserve(size); 
    }
	Expression &operator =(const Expression &e) {
	  constant = e.constant;
	  terms.reserve(std::max(e.size(),terms.capacity()));
	  terms.resize(e.size());
	  std::copy(e.terms.begin(),e.terms.end(),terms.begin());
	  return *this;
	}
	Expression &operator =(const qol::Variable &v) {
	  constant = 0.0;
	  terms.clear();
	  terms.push_back(CoeffIdxPair(1.0,int(v)));
	  return *this;
	}
	Expression &operator =(const CoeffIdxPair &c) {
	  constant = 0.0;
	  terms.clear();
	  terms.push_back(c);
	  return *this;
	}
	Expression &operator =(double c) {
	  constant = c;
	  terms.clear();
	  return *this;
	}
    Expression &operator +=(double c){
      constant += c; return *this; 
    }
    Expression &operator +=(const qol::Variable &v){
      terms.push_back(CoeffIdxPair(1.0,int(v))); return *this; 
    }
    Expression &operator +=(const CoeffIdxPair &cip){
      if(cip.first != 0.0) terms.push_back(cip); return *this; 
    }
    Expression &operator +=(const Expression &e){
	  constant += e.constant;
      reserve(size() + e.size());
      std::copy(e.terms.begin(),e.terms.end(),back_inserter(terms));
      return *this;
    }
    Expression &operator -=(int c){
	  constant -= c; return *this;
    }
    Expression &operator -=(double c){
	  constant -= c; return *this;
    }
    Expression &operator -=(const qol::Variable &v){
      terms.push_back(CoeffIdxPair(-1.0,int(v))); return *this; 
    }
	Expression &operator -=(const CoeffIdxPair &cip){
	  if(cip.first != 0.0){
		terms.push_back(cip);
		terms.back().first = -cip.first;
	  }
	  return *this; 
    }
    Expression &operator -=(const Expression &e){
	  constant -= e.constant;
      terms.reserve(terms.size() + e.terms.size());
      for(Index i=0;i<e.terms.size();++i)
        terms.push_back(CoeffIdxPair(-e.terms[i].first,e.terms[i].second));
      return *this;
    }
    Expression &operator *=(double c) {
	  constant *= c;
	  if(c == 0.0){
		terms.clear();
	  }else
		for(Index i=0;i<terms.size();++i) terms[i].first*=c;
      return *this;
    }
    Expression &operator /=(double c) {
	  constant /= c;
      for(Index i=0;i<terms.size();++i) terms[i].first/=c;
      return *this;
    }
    Expression operator +(const Expression &e) const
		  { Expression exp(size()+e.size());
			exp=*this; exp+=e; return exp; }
    Expression operator -(const Expression &e) const
    { Expression exp(size()+e.size()); exp=*this; exp-=e; return exp; }
    Row operator <=(double rhs)	{ return Row(terms,'L',rhs-constant); }
    Row operator ==(double rhs) { return Row(terms,'E',rhs-constant); }
    Row operator >=(double rhs) { return Row(terms,'G',rhs-constant); }
    Row operator <=(Variable v)	{ Expression e;e+=v;return Row(&terms,&e.terms,'L',-constant); }
    Row operator ==(Variable v) { Expression e;e+=v;return Row(&terms,&e.terms,'E',-constant); }
    Row operator >=(Variable v) { Expression e;e+=v;return Row(&terms,&e.terms,'G',-constant); }
    Row operator <=(const CoeffIdxPair &v)	{ Expression e;e+=v;return Row(&terms,&e.terms,'L',-constant); }
    Row operator ==(const CoeffIdxPair &v) { Expression e;e+=v;return Row(&terms,&e.terms,'E',-constant); }
    Row operator >=(const CoeffIdxPair v) { Expression e;e+=v;return Row(&terms,&e.terms,'G',-constant); }
    Row operator <=(const Expression &e) {
	  return Row(&terms,&e.terms,'L',e.constant - constant);
	}
    Row operator ==(const Expression &e) {
	  return Row(&terms,&e.terms,'E',e.constant - constant);
	}
    Row operator >=(const Expression &e) {
	  return Row(&terms,&e.terms,'G',e.constant - constant);
	}
	/// constant part of the expression
	double constValue() const {return constant;}
    bool isEmpty(){return (terms.size()==0);}
	size_t size() const {return terms.size(); }
	CoeffIdxPair operator[](size_t i) const { return terms[i]; }
	const std::vector<CoeffIdxPair> &getTerms() const { return terms; }
  };  

  // define all sensible combinations of adding, subtracting
  // variables/expressions or multiplying them by a constant
  inline CoeffIdxPair operator *(int c,const qol::Variable &v)
  { return CoeffIdxPair(c,v); }	// avoid compiler casting v to int
  inline CoeffIdxPair operator *(double c,const qol::Variable &v)
  { return CoeffIdxPair(c,v); }
  inline CoeffIdxPair operator *(const qol::Variable &v,double c)
  { return CoeffIdxPair(c,v); }
  inline CoeffIdxPair operator *(const qol::Variable &v,int c)
  { return CoeffIdxPair(c,v); }	// avoid compiler casting v to int
  inline Expression operator *(double c,const Expression &e)
  { Expression exp(e); exp *= c;   return exp; }
  inline CoeffIdxPair operator /(const qol::Variable &v,double c)
  { return CoeffIdxPair(1.0/c,v); }
  inline Expression operator /(const Expression &exp,double c)
  { Expression e(exp); e /= c;  return e;	}
  
  inline Expression operator +(const double val,const Variable &var)
  { Expression e(1,var); e+=val; return e;}
  inline Expression operator +(const Variable &var,const double val)
  { Expression e(1,var); e+=val; return e;}
  inline Expression operator +(const int val,const Variable &var)
  { Expression e(1,var); e+=(double)val; return e;}
  inline Expression operator +(const Variable &var,const int val)
  { Expression e(1,var); e+=(double)val; return e;}
  inline Expression operator +(const double val,const CoeffIdxPair &c)
  { Expression e(c.first,c.second); e += val; return e;}
  inline Expression operator +(const CoeffIdxPair &c,const double val)
  { Expression e(c.first,c.second); e += val; return e;}
  inline Expression operator +(const Expression &e,const double val)
  { if(val==0.0) return e; Expression exp(e); exp+=val; return exp;}
  inline Expression operator +(const double val,const Expression &e)
  { if(val==0.0) return e; Expression exp(e); exp+=val; return exp;}
  inline Expression operator +(const Variable &v,const CoeffIdxPair &c)
  { Expression exp(size_t(2)); exp+=v; exp+=c; return exp; }
  inline Expression operator +(const CoeffIdxPair &c,const Variable &v)
  { Expression exp(size_t(2)); exp+=v; exp+=c; return exp; }
  inline Expression operator +(const CoeffIdxPair &c1,const CoeffIdxPair &c2)
  { Expression exp((size_t)2); exp+=c1; exp+=c2; return exp; }
  inline Expression operator +(const CoeffIdxPair &c,const Expression &e)
  { Expression exp(e); exp+=c; return exp; }
  inline Expression operator +(const Expression &e,const CoeffIdxPair &c)
  { Expression exp(e); exp+=c; return exp; }
  inline Expression operator +(const qol::Variable &v,const qol::Variable &w)
  { Expression exp((size_t)2); exp+=v; exp += w; return exp; }
  inline Expression operator +(const qol::Variable &v,const Expression &e)
  { Expression exp(e); exp += v; return exp; }
  inline Expression operator +(const Expression &e,const qol::Variable &v)
  { Expression exp(e); exp += v; return exp; }
  inline CoeffIdxPair operator -(const qol::Variable &v)
  { return CoeffIdxPair(-1.0,v); }
  inline Expression operator -(const double val,const Variable &var)
  { Expression e(-1,var); e+=val; return e;}
  inline Expression operator -(const Variable &var,const double val)
  { Expression e(1,var); e-=val; return e;}
  inline Expression operator -(const int val,const Variable &var)
  { Expression e(-1,var); e+=(double)val; return e;}
  inline Expression operator -(const Variable &var,const int val)
  { Expression e(1,var); e-=(double)val; return e;}
  inline Expression operator -(const double val,const CoeffIdxPair &c)
  { Expression e(val); e -= c; return e;}
  inline Expression operator -(const CoeffIdxPair &c,const double val)
  { Expression e(c); e -= val; return e;}
  inline Expression operator -(const Expression &e,const double val)
  { Expression exp(e); exp-=val; return exp;}
  inline Expression operator -(const double val,const Expression &e)
  { Expression exp(val); exp-=e; return exp;}
  inline Expression operator -(const Expression &e,const int val)
  { Expression exp(e); exp-=val; return exp;}
  inline Expression operator -(const int val,const Expression &e)
  { Expression exp(e.size()); exp += val; exp-=e; return exp;}
  inline Expression operator -(const qol::Variable &v,const Expression &e)
  { Expression exp(e.size()+1); exp=v; exp -= e; return exp; }
  inline Expression operator -(const Expression &e,const qol::Variable &v)
  { Expression exp(e); exp -= v; return exp; }
  inline Expression operator -(const qol::Variable &v,const qol::Variable &w)
  { Expression exp((size_t)2); exp+=v; exp -= w; return exp; }
  inline Expression operator -(const qol::Variable &v,const CoeffIdxPair &c)
  { Expression exp((size_t)2); exp=v; exp -= c; return exp; }
  inline Expression operator -(const CoeffIdxPair &e,const qol::Variable &v)
  { Expression exp((size_t)2); exp+=e; exp -= v; return exp; }
  inline Expression operator -(const qol::CoeffIdxPair &v,
							   const qol::CoeffIdxPair &w)
  { Expression exp((size_t)2); exp+=v; exp -= w; return exp; }
  inline Expression operator -(const qol::CoeffIdxPair &v,const Expression &e)
  { Expression exp(e.size()+1); exp=v; exp -= e; return exp; }
  inline Expression operator -(const Expression &e,const qol::CoeffIdxPair &v)
  { Expression exp(e.size()+1); exp=e; exp -= v; return exp; }
	

  /// constraints of the form  constant <= \sum coeff*variable
  inline Row operator <=(double lhs,const Expression &e) {
	return Row(0,&e.getTerms(),'L',e.constValue() - lhs);
  }
  inline Row operator ==(double lhs,const Expression &e) {
	return Row(0,&e.getTerms(),'E',e.constValue() - lhs);
  }
  inline Row operator >=(double lhs,const Expression &e) {
	return Row(0,&e.getTerms(),'G',e.constValue() - lhs);
  }
  /// constraints of the form  coeff*variable <= (\sum) coeff*variable
  inline Row operator <=(const CoeffIdxPair &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(lhs);
	  return Row(&LHS,&e.getTerms(),'L',e.constValue());
  }
  inline Row operator ==(const CoeffIdxPair &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(lhs);
	  return Row(&LHS,&e.getTerms(),'E',e.constValue());
  }
  inline Row operator >=(const CoeffIdxPair &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(lhs);
	  return Row(&LHS,&e.getTerms(),'G',e.constValue());
  }
  inline Row operator <=(const Variable &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(CoeffIdxPair(1.0,lhs));
	  return Row(&LHS,&e.getTerms(),'L',e.constValue());
  }
  inline Row operator ==(const Variable &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(CoeffIdxPair(1.0,lhs));
	  return Row(&LHS,&e.getTerms(),'E',e.constValue());
  }
  inline Row operator >=(const Variable &lhs,const Expression &e) {
	  std::vector<CoeffIdxPair> LHS; LHS.push_back(CoeffIdxPair(1.0,lhs));
	  return Row(&LHS,&e.getTerms(),'G',e.constValue());
  }
  inline Row operator <=(const Variable &lhs,const CoeffIdxPair &e) {
	  return lhs - e <= 0;
  }
  inline Row operator ==(const Variable &lhs,const CoeffIdxPair &e) {
	  return lhs-e == 0;
  }
  inline Row operator >=(const Variable &lhs,const CoeffIdxPair &e) {
	  return lhs-e >= 0;
  }
  inline Row operator <=(const CoeffIdxPair &lhs,const CoeffIdxPair &e) {
	  return lhs - e <= 0;
  }
  inline Row operator ==(const CoeffIdxPair &lhs,const CoeffIdxPair &e) {
	  return lhs-e == 0;
  }
  inline Row operator >=(const CoeffIdxPair &lhs,const CoeffIdxPair &e) {
	  return lhs-e >= 0;
  }
  inline Row operator <=(const CoeffIdxPair &lhs,const Variable &e) {
	  return lhs - e <= 0;
  }
  inline Row operator ==(const CoeffIdxPair &lhs,const Variable &e) {
	  return lhs-e == 0;
  }
  inline Row operator >=(const CoeffIdxPair &lhs,const Variable &e) {
	  return lhs-e >= 0;
  }
	
  /// constraints of the form  variable <= variable
  inline Row operator <=(Variable v,Variable w) {
	  Expression e(size_t(2)); e+=v; e-=w;
	  return Row(e.getTerms(),'L',0.0);
  }
  /* Could add the following - but syntax can be confusing: does v==w check if
	 the two are the same variable (same index) or create a constraint???
  inline Row operator ==(Variable v,Variable w) {
	  Expression e((size_t)2); e+=v; e-=w;
	  return Row(e.getTerms(),'E',0.0);
  }   */
  inline Row operator >=(Variable v,Variable w) {
	  Expression e(size_t(2)); e+=v; e-=w;
	  return Row(e.getTerms(),'G',0.0);
  }



  /// Abstract MIP Formulation.
  /// Note Formulation is the most basic implementation of qol::Formulation
  /// It doesn't solve anything just stores the data.
  /// In practice typically want to use things like CplexFormulation or GrbFormulation
  /// to solve with particular solvers (CPLEX, Gurobi)

  class MIP {
  public:
    // these members are public to make it easier for
    // solvers to dump it into their solution structure
    std::vector<double> collb,colub; // lower & upper bounds
    std::vector<char> coltype; // C, I, B
    std::vector<std::string> colname; // variable name
    std::vector<std::string> rowname; // variable name
    std::vector<double> cost; // objective coefficient
    //matrix
    //matrix[0] will contain a vector of the contraints where variable with
    //index 0 is used, together with the coefficient for the variable for each constraint.
    //if matrix[0] is [(1.0, 0), (2.0, 1)], means the coefficient 
    //of the variable for the first constraint is 1.0, while this coefficient
    //is 2.0 for the second contraint
    std::vector<std::vector<CoeffIdxPair> >matrix; // row indices for each column
    std::vector<char> sense;
    std::vector<double> rhs;
    //if alligment is needed, then use something like
    //std::vector<double, AlignmentAllocator<double, 16> > rhs;
    Index varSize;

    MIP () {}
    virtual ~MIP(){}// = 0;

    virtual void setVarSize (Index s) { varSize = s; }
    virtual Index getVarSize () const { return varSize; }

    virtual Index nVar() const = 0;
	virtual Index nConstr() const = 0;
	/// when resizing all new variables default to binary variables  
    virtual void setNumVar(Index n) = 0; // resize to given size
    virtual qol::Variable getVar(Index i) const {return qol::Variable(i); }
	virtual qol::ConstraintMIP getConstr(Index i) {return qol::ConstraintMIP(this,i);}

    /// <summary>
    /// adds a new variable to the problem (defaults to binary)
    /// </summary>
    /// <return>
    /// Index of the variable
    /// </return>
    virtual qol::Variable addVar(){
      const Index n=nVar();
      setNumVar(n+1);
      return getVar(n);
    }

    /// <summary>
    /// adds a new variable to the problem and initialises it
    /// </summary>
    /// <param name="lb"> lower bound of the variable </param>
    /// <param name="ub"> upper bound of the variable </param>
    /// <param name="obj"> objective coefficient of the variable </param>
    /// <param name="type"> type of the variable </param>
    /// <param name="name"> name of the variable </param>
    /// <return>
    /// Index of the variable
    /// </return>
    virtual qol::Variable addVar (double lb, double ub, double obj, char type, std::string name ="") {
      const Index n=nVar();
      setNumVar(n+1);
      qol::Variable var = getVar(n);
      setVarLB (var, lb);
      setVarUB (var, ub);
      setObjCoeff (var, obj);
      setVarType (var, type);
      setVarName (var, name);
      return var;
    }
    /// <summary>
    /// initialises a variable
    /// </summary>
    /// <param name="lb"> lower bound of the variable </param>
    /// <param name="ub"> upper bound of the variable </param>
    /// <param name="obj"> objective coefficient of the variable </param>
    /// <param name="type"> type of the variable </param>
    /// <param name="name"> name of the variable </param>
    /// <return>
    /// Index of the variable
    /// </return>
    virtual void setVar (qol::Variable var, double lb, double ub, double obj, char type, std::string name ="") {
      setVarLB (var, lb);
      setVarUB (var, ub);
      setObjCoeff (var, obj);
      setVarType (var, type);
	  if(name != "")
		  setVarName (var, name);
    }

    /// <summary>
    /// adds an array of variables
    /// </summary>
    /// <param name="len"> number of variable </param>
    /// <return>
    /// Array of variables
    /// </return>
    virtual std::vector<qol::Variable> addVarVec(Index len) {
      const Index n=nVar(); 
      setNumVar(n+len);
      std::vector<qol::Variable> vec (len);
      for(int i=0;i<(int)len;++i)
        vec[i] = getVar(n+i);
      return vec;
    }

    //virtual Index nSubprob() const = 0;// { return 1; }

    // variable modification

    virtual void setObjCoeff(qol::Variable v,double obj) = 0;
    virtual double getObjCoeff(qol::Variable v) const = 0;

    /// <summary>
    /// set the upper bound of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    /// <param name="ub">upper bound of the variable</param>
    virtual void setVarUB(Variable v,double ub) = 0;

    /// <summary>
    /// get the upper bound of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    virtual double getVarUB(Variable v) const = 0;

    /// <summary>
    /// set the lower bound of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    /// <param name="lb">lower bound of the variable</param>
    virtual void setVarLB(Variable v,double lb) = 0;

    /// <summary>
    /// get the lower bound of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    virtual double getVarLB(Variable v) const = 0;

    /// <summary>
    /// set the type of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    /// <param name="type">type of the variable</param>
    virtual void setVarType(Variable v,char type) = 0;

    /// <summary>
    /// get the type of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    virtual char getVarType(Variable v) const = 0;

    /// <summary>
    /// set the name of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    /// <param name="name">name of the variable</param>
    virtual void setVarName(Variable v,const std::string &name) = 0;

    /// <summary>
    /// get the name of a variable
    /// </summary>
    /// <param name="v">index of the variable</param>
    virtual const std::string getVarName (Variable v) const = 0;

    /// constraint modification
    virtual ConstraintMIP addConstraint(const Row &row) = 0;
	/// allow .add() since it's obvious we are adding a constraint
	/// also makes for easier modification from CPLEX/Concert code
    virtual ConstraintMIP add(const Row &row) { return addConstraint(row); }
    virtual void setSense(Constraint cnstr,char sense) = 0;
    virtual char getSense(const Constraint cnstr) const = 0;
    virtual void setRHS(Constraint cnstr,double rhs) = 0;
    virtual double getRHS(Constraint cnstr) const = 0 ;
    virtual void setConstrName (Constraint cnstr,const std::string &name) = 0;
	  virtual std::string getConstrName(Constraint cnstr) const = 0;//{std::cerr << "No getConstrName()\n"; return "";}
    /// write formulation out to file
    virtual void writeLP(const char *filename) = 0;
    /// write formulation out to file
    void writeLP(const std::string filename) { writeLP(filename.c_str()); } // convenience method for string names
    virtual double getTolerance () {return eps;}
  };
}
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
