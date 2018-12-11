/// This header file provides the basic concepts for defining a problem
#ifndef __QOL_DATA__
#define __QOL_DATA__

#include "QolUtil.h"
#include <list>
#include <vector>
#include <string>

namespace qol {
  /// Variable is really just an index
  class Variable {
  protected:
    Index index;
  public:
#     if __cplusplus >= 201103L
	  enum : char {CONTINUOUS='C',BINARY='B',INTEGER='I'};
#     else	  
	  enum {CONTINUOUS='C',BINARY='B',INTEGER='I'};
#     endif	  
	  //static const char CONTINUOUS='C',BINARY='B',INTEGER='I';
    Variable(Index i=InvalidIndex) : index(i) {}
    operator Index() const { return index; }
    operator int() const { return int(index); } // so variable can be an index
    bool operator == (const Variable &v) const { return (int(v)==int(*this));}
    bool operator != (const Variable &v) const { return (int(v)!=int(*this));}
	bool operator < (const Variable &v) const { return (index < v.index); } // mainly for use of variables in maps
	bool operator < (const int i) const { return (int(*this)<i); }
    bool operator > (const int i) const { return (int(*this)>i); }
    bool operator < (const Index i) const { return (Index(*this)<i); }
    bool operator > (const Index i) const { return (Index(*this)>i); }
    bool valid() const {return index!=InvalidIndex; }
  };

  inline std::ostream& operator << (std::ostream& os, const Variable &v) {
    return os<<int(v);
  }

  class VarVec : public std::vector<qol::Variable> {
  public:
    VarVec(Index s=0,qol::Variable v=0) : std::vector<qol::Variable>(s,v) {}

    /// assign every value to v
    qol::Variable operator=(const qol::Variable v) {std::fill(begin(),end(),v); return v;}

    /// allow assignment of straight vectors
    VarVec &operator =(const std::vector<qol::Variable> &v){
      *static_cast<std::vector<qol::Variable> *>(this) = v;
      return *this;
    }
  };

  /// ExtendedVariable extends the Variable class to include the name of the variable
  /// and the primal and dual values. This is used in the solutionFull class
  class ExtendedVariable:public Variable {
  protected:
    std::string name;
    double primal, dual;
  public:
    ExtendedVariable(Index i) : Variable(i), primal(0.0), dual(0.0) {}
    ExtendedVariable(Index i, std::string n, double p, double d) : Variable(i), name(n), primal(p), dual(d) {}
    friend std::ostream& operator << (std::ostream& os, const ExtendedVariable &v);
  };
  inline std::ostream& operator << (std::ostream& os, const ExtendedVariable &v) {
    return os << v.index << " " << v.name << " " << v.primal << " " << v.dual;
  }
  // Class that extends an array of Extended variables. Used by solutionFull
  class ExtendedVarVec : public std::vector<qol::ExtendedVariable> {
  public:
    ExtendedVarVec(Index s=0,qol::ExtendedVariable v=0) : std::vector<qol::ExtendedVariable>(s,v) {}

    /// assign every value to v
    qol::ExtendedVariable operator=(const qol::ExtendedVariable v) {std::fill(begin(),end(),v); return v;}

    /// allow assignment of straight vectors
    ExtendedVarVec &operator =(const std::vector<qol::ExtendedVariable> &v){
      *static_cast<std::vector<qol::ExtendedVariable> *>(this) = v;
      return *this;
    }
  };

  class Constraint {
    Index index;
  public:
    static const char LE='L',EQ='E',GE='G';
    Constraint(Index i=InvalidIndex) : index(i) {}
    operator Index() const { return index; }
    operator int() const { return int(index); }
  };

  typedef std::pair<double,Variable> CoeffVarPair;

  class ConstraintRow : public Constraint, public std::vector<CoeffVarPair> {
  public:
    double rhs;
    char sense; // one of LE, EQ or GE
  };
}

#endif
