LimitedResourceSched:

include ../../../Makefile
$(QOL_LIB): $(wildcard $(QOL_DIR)/*.h $(QOL_DIR)/*.cpp)
	cd $(QOL_DIR); make
SRC:=$(filter-out %td%fx.h,$(wildcard *.cpp))
OBJ:=$(SRC:.cpp=.o)
debug: LimitedResourceSched
debug: NAME:=_debug

release: LimitedResourceSched
release: NAME:=_release

LimitedResourceSched: $(OBJ) $(wildcard *.h) $(QOL_LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $(OBJ) $(QOL_LIBS)
	'cp' -f LimitedResourceSched ../bin/limitedresched$(NAME).exe
