
EIGEN_INC = /opt/local/include/eigen3
BOOST_LIB = /opt/local/lib
BOOST_INC = /opt/local/include

EXE = runeems_sats
OBJ = runeems_sats.o eems.o util.o mcmc.o draw.o habitat.o graph.o

CXXFLAGS = -I${BOOST_INC} -I${EIGEN_INC} -O3 -DNDEBUG
LDFLAGS = -L${BOOST_LIB} \
	-lboost_system \
	-lboost_program_options \
	-lboost_filesystem

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE) *~
