CPP=g++
LIBS=-lgsl -lgslcblas
LPATH=-L/opt/homebrew/opt/gsl/lib
IPATH=-I/opt/homebrew/opt/gsl/include
SRC=main.cpp kerr.cpp numeric.cpp kdisk.cpp kd_rms.cpp kd_flux.cpp Noviflux_PageThorne74.cpp
FILES=$(addprefix ./src/,$(SRC))
BINDIR=./bin/
OUT=$(BINDIR)kdfirewall

build: $(FILES)
	$(CPP) -o $(OUT) $(FILES) $(LIBS) $(IPATH) $(LPATH)
