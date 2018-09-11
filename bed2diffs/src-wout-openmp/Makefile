
PLINKIO = /usr/local
EXE = bed2diffs_v1
OBJ = bed2diffs_v1.o data.o

CXXFLAGS = -I${PLINKIO} -O3 -Wall -Werror
LDFLAGS = -lplinkio


all:
	@echo "Use 'make linux' on a Linux machine, or 'make darwin' on a Mac."

linux: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ -L${PLINKIO}/lib -Wl,-R${PLINKIO}/lib $(LDFLAGS)

darwin: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ -L${PLINKIO}/lib $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE) *~
