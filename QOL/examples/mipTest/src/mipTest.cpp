// mipTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include "QolMIP.h"
#include "QolColFormulation.h"

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
    // QolColFormulation is not intended to solve a model. It is only used
    // to store variables and constraints
    qol::QolColFormulation *m = new qol::QolColFormulation();
    //Other Formulations provided (that provide a solve() method)
    // qol::CpxFormulation *m = new qol::CpxFormulation ();
    // qol::GurobiFormulation *m = new qol::GurobiFormulation();

    //Declaration of variables
    qol::Variable x,y,z;
    //Example of how to completely define two variables in the constructor
    x = m->addVar (0.0, 1.0, 0.5, qol::Variable::CONTINUOUS, "x");
    y = m->addVar (0, 2, 1, qol::Variable::INTEGER, "y");
    //"Empty" variable which information can be completed later
    z = m->addVar ();
    m->setVarType (y, qol::Variable::INTEGER);
    m->setVarLB (z, 0.0);
    m->setVarUB (z, 1.0);
    m->setVarType (z, qol::Variable::CONTINUOUS);
    m->setVarName (z, "z");
    m->setObjCoeff (z, 1.0);

    //Add two constrains with the two different options provided
    m->addConstraint( x - 2.0*y == 2.5 ).setName("eq%.1f",2.5);
    //Using an expression
    qol::Expression sum;
    for(int i=0;i<3;++i) {
      sum += m->getVar(i);
    }
	qol::Constraint c;
    m->addConstraint( sum <= 2 ).setName("sum");
	m->addConstraint(-x + y +3 <= 2*y - z + 3.0*z + 11).setName("complex");
    //Write the LP file
    m->writeLP("test.lp");
	std::cout << "wrote test.lp\n";
	return 0;
}

