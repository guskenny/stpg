all: mipTest unitTest bendersatsp

include ../../../Makefile
$(QOL_LIB): $(wildcard $(QOL_DIR)/*.h $(QOL_DIR)/*.cpp)
	cd $(QOL_DIR); make
mipTest: mipTest.o $(QOL_LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(QOL_LIBS)
unitTest: unitTest.o $(QOL_LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(QOL_LIBS)
bendersatsp: bendersatsp.o $(QOL_LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(QOL_LIBS)
