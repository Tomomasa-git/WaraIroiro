CC=gcc
CXX=g++

CFLAGS  = -O2

BINDIR = ../bin
LIBDIR = ../lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET = main   HyperTritonLifetime   Setting
OBJS   = main.o HyperTritonLifetime.o Setting.o

.PHONY: all
all: $(TARGET)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)

$(TARGET): $(patsubst %,$(LIBDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(TARGET).d $(TARGET).o $(TARGET)

