// unit test: try to more or less systematically test all of the functionality
// of the QOL interface (or at least the MIP formulation)
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
#include "QolMIP.h"
#include "QolColFormulation.h"
#include "CplexFormulation.h"
#include "GurobiFormulation.h"
#include "CoinFormulation.h"

//using namespace qol;

/// \code
/// using namespace qol
/// MIP::Formulation m;
/// MIP::Variable x(m),y(m),z(m);
/// x >= 0; y >= 0; z >= 0; // set lower bounds
/// x <= 1; y <= 2; z <= 1; // set upper bounds
/// x.setObjCoeff(2); y.setObjCoeff(-1.0);
/// m.addRow( x + 2*y == 2.5 )
/// Expression sum;
/// for(i=0;i<3;++i) sum += m.getVar(i)
/// m.addRow( sum <= 2 );
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

int main(int argc, char* argv[])
{
    qol::FormulationType solver = qol::CPLEX;
    qol::MIPSolver *mip=0;
    std::cout << "Unit testing of QOL interface to MIP solvers\n";
    for(int a=1;a<argc;++a){
	if( argv[a][0] == '-'){
	    switch(argv[a][1]){
#               ifdef _USE_CPLEX_		
		case 'c': case 'C':
		    solver = qol::CPLEX;
		    mip = new qol::CplexFormulation;
		    std::cout << "Testing CplexFormulation class\n";
		    break;
#               endif
#               ifdef _USE_GUROBI_		
		case 'g': case 'G':
		    solver = qol::GUROBI;
		    mip = new qol::GurobiFormulation;
		    std::cout << "Testin GurobiFormulation class\n";
		    break;
#               endif
#               ifdef _USE_COIN__		
		case 'o': case 'O':
		    solver = qol::COIN;
		    mip = new qol::CoinFormulation;
		    std::cout << "Testin GurobiFormulation class\n";
		    break;
#               endif
		case 'h':
		    std::cout << "Usage: choose solver to test. Available options:\n";
#		    ifdef _USE_CPLEX_
		    std::cout << "  -c CPLEX\n";
#                   endif

#		    ifdef _USE_GUROBI_
		    std::cout  << "  -g GUROBI\n";
#                   endif

#		    ifdef _USE_COIN_
		    std::cout << "  -o COIN\n";
#                   endif
		    std::cout << "  -h This help message\n";		    
		    return 0;
	    }
	}
    }
		
    if(mip==0) {
	std::cerr << "ERROR no solver specified. Use -h to see options\n";
	return 2;
    }
    

    //Declaration of variables
    qol::Variable x,y,z;
    std::cout << "Testing addVar() ...\n";
    //Example of how to completely define two variables in the constructor
    x = mip->addVar (0.0, 5.0, 0.5, qol::Variable::CONTINUOUS, "x");
    y = mip->addVar (0, 2, 1, qol::Variable::INTEGER, "y");
    //"Empty" variable which information can be completed later
    z = mip->addVar ();
    std::cout << "Testing setVarType() ...\n";
    mip->setVarType (y, qol::Variable::INTEGER);
    mip->setVarType (z, qol::Variable::CONTINUOUS);
    std::cout << "Testing setVarLB() ...\n";
    mip->setVarLB (z, 0.0);
    std::cout << "Testing setVarUB() ...\n";
    mip->setVarUB (z, 1.0);
    std::cout << "Testing setVarName() ...\n";
    mip->setVarName (z, "z");
    std::cout << "Testing setObjCoeff() ...\n";
    mip->setObjCoeff (z, 1.0);

    //Add two constrains with the two different options provided
    std::cout << "Testing addConstraint() ...\n";
    qol::ConstraintMIP ceq=mip->addConstraint( x - 2.0*y == 2.5 );
    ceq.setName("eq%.1f",2.5);
    //Using an expression
    qol::Expression sum;
    for(int i=0;i<3;++i) 
	sum += mip->getVar(i);
    qol::ConstraintMIP ccomp = mip->addConstraint(-x + y +3 <= 2*y - z + 3.0*z + 11);
    ccomp.setName("comp");    
    qol::Constraint csum = mip->addConstraint( sum <= 7 );
    std::cout << "Testing constraint.setName() ... \n";
    mip->setConstrName(ccomp,"complex");
    mip->setConstrName(csum,"sum");
    std::cout << "Testing set VERBOSITY=high\n";
    qol::Parameters params;
    params.setParamVal(qol::VERBOSITY,1);
    mip->setParameters(params);
    std::cout << "Testing solveRelaxed\n";
    try{
	mip->solveRelaxed();
	std::cout << "Objective = " << mip->getObjective() << std::endl;
    }catch(qol::Exception err){
	std::cout << "\tFailed: " << err.what() << std::endl;
    }
    std::cout << "Testing solveExact\n";
    mip->solveExact();
    std::cout << "Objective bound = " << mip->getObjectiveBound() << std::endl;
    std::cout << "Testing of behaviour after problem loaded into solver ------------\n";
    std::cout << "Testing getConstrName()\n";
    std::string name = mip->getConstrName(ceq);
    name += mip->getConstrName(csum);
    name += mip->getConstrName(ccomp);
    std::cout << "Testing constraint.setName() ... \n";
    mip->setConstrName(ccomp,"complex");
    mip->setConstrName(csum,"sum");
    std::cout << "Testing getPrimal(constr) && getRHS(constr)\n";
    std::cout << "\t" << mip->getConstrName(ceq) << " = "
	      << mip->getPrimal(ceq) << " == " << mip->getRHS(ceq) <<std::endl;
    std::cout << "\t" << mip->getConstrName(csum) << "= "
	      << mip->getPrimal(csum) << " <= " << mip->getRHS(csum) <<std::endl;
    std::cout << "\t" << mip->getConstrName(ccomp) << " = "
	      << mip->getPrimal(ccomp) << " <= " << mip->getRHS(ccomp) <<std::endl;
    std::cout << "Testing getPrimal(var)\n";
    double tmp=0;
    tmp += mip->getPrimal(x);
    tmp += mip->getPrimal(y);
    tmp += mip->getPrimal(z);
    //Write the LP file
    std::cout << "Testing writeLP()\n";
    mip->writeLP("test.lp");
    std::cout << "wrote test.lp\n";
    std::cout << "Done...all tests passed\n";
    return 0;
}

