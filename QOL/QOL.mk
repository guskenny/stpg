# change any of the below to avoid using some of the libraries
ifndef USE_CPLEX
USE_CPLEX=1
endif
ifndef USE_GUROBI
USE_GUROBI=1
endif
ifndef USE_COIN
USE_COIN=0
endif

CXXFLAGS+=-std=c++11

## find CPLEX, Gurobi, boost, etc
ifdef BOOST_LIBRARYDIR
  CPPFLAGS+=-I$(BOOST_LIBRARYDIR)/..
  QOL_LIBS+=-L$(BOOST_LIBRARYDIR)
endif


################## FIND CPLEX DIRECTORY ############################
ifeq ($(shell which cplex),"cplex not found")
   CPPFLAGS+=-U_USE_CPLEX
   USE_CPLEX:=0
else
   CPLEX_ROOT=$(patsubst %/cplex,%/../../..,$(shell which cplex))
endif
ifeq ($(USE_CPLEX),1)
  ifdef CPLEX_ROOT
    ifneq ($(wildcard $(CPLEX_ROOT)/lib/*/static_pic),)
	QOL_LIBS+=-L$(wildcard $(CPLEX_ROOT)/lib/*/static_pic)
    else 
        # version 12.5 and above:
        QOL_LIBS+=-L$(wildcard $(CPLEX_ROOT)/cplex/lib/*/static_pic)
    endif
    QOL_LIBS+= -lcplex 
    CPPFLAGS+=-I$(CPLEX_ROOT)/include  -I$(CPLEX_ROOT)/cplex/include -D _USE_CPLEX_
  else 
    $(warning  "Cannot find CPLEX directory")
    USE_CPLEX:=0
  endif
  QOL_SRC+=\
  CplexFormulation.cpp \

endif
################## FIND GUROBI DIRECTORY ############################
ifdef GUROBI_HOME
   CPPFLAGS+=-D_USE_GUROBI_ -I$(GUROBI_HOME)/include
   # set GUROBI_VER=65 (for version 6.5.1 or whatever
   GUROBI_VER=$(patsubst $(GUROBI_HOME)/lib/libgurobi%.so,%,$(wildcard $(GUROBI_HOME)/lib/libgurobi*.so))
   QOL_LIBS+= -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi$(GUROBI_VER)
else
   USE_GUROBI:=0
endif
################## FIND COIN DIRECTORY ############################
ifeq ($(USE_COIN),1)
  ifneq ($(wildcard /mos/software/gcc4.6.1/lib/libOsiCbc*),)
    QOL_LIBS+= -L/mos/software/gcc4.6.1/lib/ \
	-lOsi -lOsiCbc -lOsiClp -lCbc -lCgl -lClp -lCoinUtils -lz
    CPPFLAGS+=-I/mos/software/gcc4.6.1/include -D_USE_COIN_
  else
    $(warning "Cannot find COIN")
    USE_COIN:=0
  endif
  QOL_SRC+=CoinFormulation.cpp
endif

