#ifndef __VOLUME_MODEL_H__
#define __VOLUME_MODEL_H__

#include "QolModel.h"

namespace qol {
  class VolumeModel: public Model {
  public:
    VolumeModel (): Model(){}
    virtual ~VolumeModel (){}


    virtual ModelType type() const {return VOLUME;}

    // virtual Index nSubprob() const ;
    virtual Index nConstr() const ; // explicit/linking constraints

    virtual qol::ConstraintMIP /*ConstraintRow*/ getConstraintRow(Constraint c);

    virtual Status solveRelaxed();
    virtual Status solveExact();

    virtual void setVarBounds(Variable var,double lb,double ub); // fix/unfix

    /// setting variable cost may be sensible as part of Lagrangian
    /// scheme or maybe to do some perturbation in a Benders scheme.
    virtual void setVarCost(Variable var,double obj, Index problem); // objective coefficient
  };
}

#endif
