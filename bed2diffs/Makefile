
.PHONY: all

PLINKIO = /usr/local/src/plinkio
EIGEN = /usr/local/src/eigen-6b38706d90a9

all: bed2diffs
OBJ = bed2diffs.o data.o

CXXFLAGS = -I${PLINKIO}/include -I${EIGEN}
LDFLAGS = -lplinkio -L${PLINKIO}/lib

bed2diffs: CXXFLAGS += -Wall -Wextra -O2 -m64 -DEIGEN_NO_DEBUG
bed2diffs: $(OBJ)
	$(CXX) $(CXXFLAGS) -o bed2diffs $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) bed2diffs
