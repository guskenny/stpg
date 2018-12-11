/// Useful class for logging
#ifndef __QOL_LOG__
#define __QOL_LOG__

#include <cstring>
#include <cstdio>

#define FLERR __FILE__,__LINE__

namespace qol {
  class Log {
  protected:
    FILE * screen, * logfile;
    int verbosity;
  public:
    Log (int _v, std::string logfilename, std::string screenfilename) : verbosity (_v) {
      if (logfilename != "")
        logfile = fopen(logfilename.c_str(),"w");
      if (screenfilename!="")
        screen = fopen (screenfilename.c_str(), "w");
      else
        screen = stdout;
    }

    ~Log () {
      if (logfile) fclose(logfile);
      if (screen && screen != stdout) fclose(screen);
    }

    inline void error(const char *str) {
      if (screen) {
        fprintf(screen,"ERROR: %s (%s:%d)\n",str,FLERR);
        if (screen!=stdout) fclose(screen);
      }
      if (verbosity && logfile) {
        fprintf(logfile,"ERROR: %s (%s:%d)\n", str,FLERR);
        fclose (logfile);
      }
      exit (1);
    }

    inline void warning(const char *str) {
      if (screen) fprintf(screen,"WARNING: %s (%s:%d)\n",str,FLERR);
      if (verbosity && logfile) fprintf(logfile,"WARNING: %s (%s:%d)\n",
        str,FLERR);
    }

    inline void message(const char *str) {
      if (screen) fprintf(screen,"%s \n",str);
      if (verbosity >1 && logfile) fprintf(logfile,"%s\n",str);
    }

  };

} //namespace qol
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
* Local Variables:
* tab-width: 4
* eval: (c-set-style "stroustrup")
* End:
*/
