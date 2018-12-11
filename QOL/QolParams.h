#ifndef __QOL_PARAM__
#define __QOL_PARAM__

#include <string>
#include <map>

namespace qol {
  enum ParamType {
    //Generic params (CPU, threads,...)
	  CPUTIME           = 01, // for cplex this changes the clock type
    WALLTIME          = 02,	  // for cplex this changes the clock type
    CPU_NUMBER        = 03,
	/// VERBOSITY for CPLEX & Gurobi = screen indicator off/on (0/1 only)
    VERBOSITY         = 04, // 0: LOW, 1: MED, 2: HIGH 
    THREADS           = 05, // passed to CPLEX & Gurobi
	/// Number of nodes in branch & bound to print for Gurobi/CPLEX
    PRINTFREQ         = 06, // node logging frequency

    //Generic optim. params
    EPS               = 21,
    ABSGAP            = 22,		// absolute optimality gap
    RELGAP            = 23,		// relative optimality gap
    EPINT             = 24,		
    TIMELIMIT         = 25,
    ITLIMIT           = 26,    	// iteration limi
    //
    ITERATIONS        = 51,

    //
    HEURISTIC_FREQ    = 71,   //How often we run the heuristic

    /// PSO
    SUBGRADFACTOR     = 101,
    VELOCITY          = 102,
    GLOBALFACTOR      = 103,
    PERTURBFACTOR     = 104,
    NPARTICLES        = 105,
    LBCHECKFREQ       = 106,
    SUBGRADFMIN       = 107,
    SUBGRADFMULT      = 108,

    /// LAGRANGIAN
    LAGITER           = 301,
    LAGTIMER          = 302,    //Time in seconds for each execution of the lagrangian

    //
    UNKNOWN           = 999999
  };

  class Param {
  private:
    ParamType ptype;
    double defaultval, minval, maxval, val;
    std::string pname;

  public:

    Param (ParamType t, std::string n, double dfault, double min, 
           double max, double v) : ptype(t), defaultval (dfault), 
           minval (min), maxval (max), val (v), pname (n) {}

    Param &operator = (const double v) {
      val = v;
      return *this;
    }

    Param &operator = (const Param p) {
      ptype = p.ptype;
      defaultval = p.defaultval;
      minval = p.minval;
      maxval = p.maxval;
      val = p.val;
      pname = p.pname;
      return *this;
    }

    inline std::string name() const {
      return pname;
    }

    inline ParamType type() const {
      return ptype;
    }

    inline double min() const { ///< minimum value in vector
      return minval;
    }
    inline double max() const { ///< maximum value in vector
      return maxval;
    }
    inline double value() const {
      return val;
    }
    inline void setval (double _v) { val = _v; }
    inline double defaultVal() const {
      return defaultval;
    }
  };

  const double InvalidParam = -999999;

  class Parameters : public std::map <ParamType, Param*> {
  private:
    std::map<ParamType, Param*>::iterator it;

  public:
    Parameters () { clear(); }
    ~Parameters () { clear(); }

    void addParam (Param *p) {
      insert (std::pair<ParamType, Param*>(p->type(), p));
    }

    void addParam (ParamType t, std::string n, double dfault, double min, 
                   double max, double v) {
      Param* p = new Param (t, n, dfault, min, max, v);
      addParam (p);
    }

    void setParamVal (ParamType type, double v) {
      it = find (type);
      if (it==end()) 
		addParam (type, "", 0.0, -inf, inf, v);
	  else
		it->second->setval (v);
    }

    inline double getParamValue (ParamType type) {
      it = find(type);
      if (it==end()) return InvalidParam;
      return it->second->value();
    }

    inline Param* getParam (ParamType type) {
      it = find(type);
      if (it==end()) return NULL;
      return it->second;
    }
  };
} //namespace

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
