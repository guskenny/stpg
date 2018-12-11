all: qolTest grbTest cpxTest

include ../../Makefile
CPPFLAGS+= -I$(CPLEX_ROOT)/concert/include
CPPFLAGS+=-DIL_STD  -Wno-sign-compare
QOL_LIBS:= -L$(wildcard $(CPLEX_ROOT)/concert/lib/*/static_pic) \
	-lconcert -lilocplex $(QOL_LIBS) -lgurobi_c++ -lgurobi$(GUROBI_VER)
$(QOL_LIB): $(wildcard $(QOL_DIR)/*.h $(QOL_DIR)/*.cpp)
	cd $(QOL_DIR); make

qolTest: qolTest.o $(QOL_LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(QOL_LIBS)

grbTest:

cpxTest:
