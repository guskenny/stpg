#ifndef __QOLMIP_INSTANCE_H__
#define __QOLMIP_INSTANCE_H__

#include "QolMIP.h"

namespace qol {

  /// <summary>
  /// MIPInstance is a qol::MIP also without solver, just for
  /// storing data.
  /// </summary>
  class MIPInstance: public MIP {
  public:
    MIPInstance():MIP(){};

    virtual ~MIPInstance(){};

    /// <summary>
    /// Writes a LP file. This is a generic format.
    /// <param name="filename">Name of the output file.</param>
    /// </summary>
    virtual void writeLP (const char* filename) {
      const Index numVars = nVar();
      const Index numConstr = nConstr();

      colname.resize(numVars,"xxx");
      rowname.resize(numConstr,"xxx");
      // do something here about missing names
      FILE *fp= fopen(filename,"w");
      if(!fp) return;
      fprintf(fp, "\\Problem name: %s\n\n", "Formulation");
      std::vector<std::vector<CoeffIdxPair> > row(nConstr());  // col index for each row
      fprintf(fp, "Minimize\n");
      
      for(Index i=0;i<numVars;++i) {
		  if( fabs(cost[i]) >= eps) {
			  fprintf(fp," %c ",(cost[i] > 0 ? '+' : '-'));
			  if( fabs(fabs(cost[i])-1.0) > eps )
				  fprintf(fp,"%g ",fabs(cost[i]));
			  fprintf(fp,"%s",colname[i].c_str());
			  if( i % 20 == 19) fprintf(fp,"\n");
		  }
		  const Index matrixSize = matrix[i].size();
		  for(Index j=0;j<matrixSize;++j){
			  if(fabs(matrix[i][j].first) < eps) continue;
			  row[matrix[i][j].second].push_back(
				  CoeffIdxPair(matrix[i][j].first,i));
        }
      }
      fprintf(fp,"\n");

      fprintf(fp, "Subject To\n");
      std::vector<Index> idx(nVar(),0);
      for(Index c=0;c<numConstr;++c){
		  if(rowname[c] == "")
			  fprintf(fp,"C#%lu:",c);
		  else
			  fprintf(fp,"%s:",rowname[c].c_str());
        int cnt=0;
        const Index rowSize = row[c].size();
        for(Index j=0;j<rowSize;++j){
          Index i = row[c][j].second;
          const double coeff = matrix[i][idx[i]].first;
          ++idx[i];
          fprintf(fp," %c ",(coeff > 0 ? '+' : '-'));
          if( fabs(fabs(coeff)-1.0) > eps )
            fprintf(fp,"%g ",fabs(coeff));
          fprintf(fp,"%s",colname[i].c_str());
          if( ++cnt % 10 == 9) fprintf(fp,"\n");
        }
        switch(sense[c]){
                    case Constraint::LE: fprintf(fp," <= "); break;
                    case Constraint::EQ: fprintf(fp," = "); break;
                    case Constraint::GE: fprintf(fp," >= "); break;
                    default: fprintf(fp," ?? ");
        }
        fprintf(fp,"%g\n",rhs[c]);
      }

      fprintf(fp, "Bounds\n");
      for(Index i=0;i<numVars;++i){
        if( colub[i] < inf){
          if(collb[i] > -inf)
            fprintf(fp,"%g <= ",collb[i]);
          fprintf(fp,"%s <= %g\n",colname[i].c_str(),colub[i]);
        }else 
          if(collb[i] > -inf)
            fprintf(fp,"%s >= %g\n",colname[i].c_str(),collb[i]);
      }

      fprintf(fp,"Binaries\n");
      int cnt = 0;
      for(Index i=0;i<numVars;++i) {
        if(coltype[i] == Variable::BINARY){
          fprintf(fp,"%s ",colname[i].c_str());
          if(++cnt%20==19) fprintf(fp,"\n");
        }
      }

      fprintf(fp,"\nIntegers\n");
      cnt = 0;
      for(Index i=0;i<numVars;++i) {
        if(coltype[i] == Variable::INTEGER){
          fprintf(fp,"%s ",colname[i].c_str());
          if(++cnt%20==19) fprintf(fp,"\n");
        }
      }

      fprintf(fp,"\nGenerals\n");
      cnt = 0;
      for(Index i=0;i<numVars;++i) {
        if(coltype[i] == Variable::CONTINUOUS){
          fprintf(fp,"%s ",colname[i].c_str());
          if(++cnt%20==19) fprintf(fp,"\n");
        }
      }

      fprintf(fp,"\nEnd\n");
      fclose(fp);
    }
  };
}

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
* Local Variables:
* tab-width: 4
* eval: (c-set-style "stroustrup")
* End:
*/
