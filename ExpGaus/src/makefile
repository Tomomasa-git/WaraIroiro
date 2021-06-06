CC = gcc 
CXX = g++ 

CFLAGS  = -O2

BINDIR = ../bin
LIBDIR = ../lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

#TARGET1=     SkewnessMod
#OBJS1=       SkewnessMod.o
#TARGET1 = Chebyshev   Setting
#OBJS1   = Chebyshev.o Setting.o
#TARGET1 = Test_2DGr   Setting
#OBJS1   = Test_2DGr.o Setting.o
#TARGET1=TestRandExpGaus Setting
#OBJS1=TestRandExpGaus.o Setting.o
#TARGET1 = TransientExp Setting
#OBJS1   = TransientExp.o Setting.o
#TARGET1 = TransDecayExp   Setting
#OBJS1   = TransDecayExp.o Setting.o
TARGET1 = RandTransExp   Setting
OBJS1   = RandTransExp.o Setting.o

all: $(TARGET1)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*
