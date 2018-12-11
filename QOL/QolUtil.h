/// Utility functions and type definitions
#ifndef __QOL_UTIL__
#define __QOL_UTIL__

#include <limits>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <exception>
#include <sstream>
#include "QolRandom.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define  omp_get_max_threads() 1
#endif

namespace qol {
  ///////////////////////////////////////////////////////////////////////////
  //Types and variables used in the library

  ///Problem status
	enum Status {HEURISTIC=1,OPTIMAL=3,INFEASIBLE=4,UNBOUNDED=5,FAILED=6,ABORTED=7};

  ///Formulation type
  enum FormulationType {QOLFORM = 1, COIN = 2, GUROBI = 3, CPLEX = 4};

  ///Model type
  enum ModelType {QOLMODEL = 1, BENDERS = 2, LAGRANGIAN = 3, VOLUME = 4};
  
  //Algorithm type
  enum AlgorithmType {QOLOPTIMISER = 1, QOLPSO = 2, QOLACO = 3};

  //Solution type
  enum SolutionType {PRIMAL = 1, FULL =2};

  //Heuristic type
  // NeighbourhoodSearch - need a primal feasible solution
  // Repair - need a primal, not necessarily feasible solution
  // Constructive - just need some weights (costs) for each variable
  enum HeuristicType {REPAIR = 1, CONSTRUCTIVE = 2, NEIGHBOURHOODSEARH = 3}; 

  typedef size_t Index;

  ///////////////////////////////////////////////////////////////////////////
  ///Constants
  const double eps = std::numeric_limits <double>::epsilon();
  /// double::infinty is too big. but there is no solver standare
  /// Gurobi 5.5 uses 1e100
  /// CPLEX 12.5 uses 1e20
  const double inf = 1e100;//std::numeric_limits<double>::infinity();
  const Index InvalidIndex = std::numeric_limits<Index>::max();

#ifdef _GLIBCXX_USE_NOEXCEPT
#  define QOL_EXCEPT_THROW _GLIBCXX_USE_NOEXCEPT
#else // should work under Visual C++ (?)
#  define QOL_EXCEPT_THROW throw()
#endif  
  /// Exception Class to use when things go really wrong
  class Exception : public std::exception {
	std::string message;
  public:
    Exception(const std::string &msg) : message(msg) {}
	virtual ~Exception() QOL_EXCEPT_THROW {}
	virtual const char *what() const QOL_EXCEPT_THROW
	{ return message.c_str(); }
  };
# ifdef DEBUG
#   define qolAssert(condition,message) {if( !(condition) ){std::stringstream msgbuf; msgbuf << message;throw Exception(msgbuf.str());}}
# else
#  define qolAssert(condition,message)
# endif  
  ///////////////////////////////////////////////////////////////////////////
  /** Integer vector with  some extra convenience methods */
  class IntVec : public std::vector<int> {
  public:
    IntVec(size_t s=0,int v=0) : std::vector<int>(s,v) {}

    /// fill with value v
    int operator=(const int v) { std::fill(begin(),end(),v); return v; }

    /// allow assignment of straight vectors
    IntVec &operator =(const std::vector<int> &v){
      *static_cast<std::vector<int> *>(this) = v;
      return *this;
    }

    /// vector addition
    IntVec &operator += (const IntVec v) {
      const_iterator vi = v.begin();
      for(iterator ti=begin();ti!=end()&&vi!=v.end();++ti,++vi)*ti += *vi;
      return *this;
    }

    /// vector scalar multiplication
    IntVec &operator *=(int v) {
      for(iterator ti=begin();ti!=end();++ti)*ti *= v;
      return *this;
    }

    int min() const { ///< minimum value in vector
      return *std::min_element(begin(),end());
    }			    

    int max() const { ///< maximum value in vector
      return *std::max_element(begin(),end());
    }

    int sum() const { ///< sum of values in vector
      return std::accumulate(begin(),end(),0);
    }
  };

  ///////////////////////////////////////////////////////////////////////////
  /** Floating point vector with some extra convenience methods */
  class DblVec : public std::vector<double> {
  public:
    DblVec(Index s=0,double v=0) : std::vector<double>(s,v) {}

    /// assign every value to v
    double operator=(const double v) {std::fill(begin(),end(),v); return v;}

    /// allow assignment of straight vectors
    DblVec &operator =(const std::vector<double> &v){
      *static_cast<std::vector<double> *>(this) = v;
      return *this;
    }

    /// vector addition
    DblVec &operator += (const DblVec v) {
      const_iterator vi = v.begin();
      for(iterator ti=begin();ti!=end()&&vi!=v.end();++ti,++vi)*ti += *vi;
      return *this;
    }

    /// vector scalar multiplication
    DblVec &operator *=(double v) {
      for(iterator ti=begin();ti!=end();++ti)*ti *= v;
      return *this;
    }

    double min() const { ///< minimum value in vector
      return *std::min_element(begin(),end());
    }

    double max() const { ///< maximum value in vector
      return *std::max_element(begin(),end());
    }

    double sum() const { ///< sum of values in vector
      return std::accumulate(begin(),end(),0.0);
    }
  };

  typedef std::pair<int,double> IdxVal;
  /** SparseVec defines a sparse integer vector. Defined in terms of indices
  and the value to be taken at that index
  */
  class SparseVec: public std::vector<IdxVal> {
  public:
    void add(int index,double value) { push_back(IdxVal(index,value)); }
  };
} //namespace qol
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
