
.PHONY: all

PLINKIO = /usr/local
EXE = bed2diffs_v1
OBJ = bed2diffs_v1.o data.o

CXXFLAGS = -I${PLINKIO}/include -O3 -Wall -Werror -fopenmp
LDFLAGS = -L${PLINKIO}/lib -lplinkio

all: $(EXE)
$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJ) $(EXE) *~
