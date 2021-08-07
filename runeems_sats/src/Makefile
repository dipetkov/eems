
EIGEN_INC = /usr/local/opt/eigen/include/eigen3
BOOST_LIB = /usr/local/lib
BOOST_INC = /usr/local/include

EXE = runeems_sats
OBJ = runeems_sats.o eems.o util.o mcmc.o draw.o habitat.o graph.o

CXXFLAGS = -std=c++14 -I${BOOST_INC} -I${EIGEN_INC} -O3 -DNDEBUG
LDFLAGS = \
	-lboost_system \
	-lboost_program_options \
	-lboost_filesystem

all:
	@echo "Use 'make linux' on a Linux machine, or 'make darwin' on a Mac."

linux: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ -L${BOOST_LIB} -Wl,-R${BOOST_LIB} $(LDFLAGS)

darwin: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ -L${BOOST_LIB} $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE) *~
