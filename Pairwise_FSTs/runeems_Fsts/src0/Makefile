
EIGEN_INC = /opt/local/include/eigen3
BOOST_LIB = /opt/local/lib
BOOST_INC = /opt/local/include

EXE = runeems_Fsts_v0
OBJ = runeems_Fsts_v0.o eems.o util.o mcmc.o draw.o habitat.o graph.o

CXXFLAGS = -I${BOOST_INC} -I${EIGEN_INC} -O3 -DNDEBUG
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
