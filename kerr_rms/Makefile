CPP=g++
LIBS=-lgsl -lgslcblas
LPATH=-L/opt/homebrew/opt/gsl/lib
IPATH=-I/opt/homebrew/opt/gsl/include
SRC=kerr_rms.cpp
FILES=$(addprefix ./src/,$(SRC))
BINDIR=./bin/
OUT=$(BINDIR)kerrrms

build: $(FILES)
	$(CPP) -o $(OUT) $(FILES) $(LIBS) $(IPATH) $(LPATH)
