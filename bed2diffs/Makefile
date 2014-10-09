
.PHONY: all

PLINKIO = /usr/local/src/plinkio
EIGEN = /usr/local/src/eigen-6b38706d90a9

EXE = bed2diffs
OBJ = bed2diffs.o data.o

CXXFLAGS = -I${PLINKIO}/include -I${EIGEN} -Wall -Wextra -O2 -m64 -DEIGEN_NO_DEBUG
LDFLAGS = -L${PLINKIO}/lib -lplinkio

all: $(EXE)
$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJ) $(EXE)
