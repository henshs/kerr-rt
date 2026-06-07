CPP=g++
LIBS=
SRC=main.cpp kdimage.cpp
SRCDIR=./src/
BINDIR=./bin/

FILES=$(addprefix $(SRCDIR),$(SRC))
OUT=$(BINDIR)kdimage

build: $(FILES)
	$(CPP) -o $(OUT) $(FILES) $(LIBS)
	
