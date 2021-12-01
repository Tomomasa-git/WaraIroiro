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

#TARGET1=     HypChart   Setting
#OBJS1=       HypChart.o Setting.o
TARGET1=     HypChart_MThesis   Setting
OBJS1=       HypChart_MThesis.o Setting.o

all: $(TARGET1)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*
