CC = gcc 
CXX = g++ 

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit #Minuitを使いたい時は"-lMinuit"を追加する。<--ココ重要。
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1  = CalcElossFunc 
OBJS1    = CalcElossFunc.o 

all: $(TARGET1) 
#     $(TARGET2) 
#     $(TARGET3) 
#     $(TARGET4) 
#     $(TARGET5)  
#     $(TARGET6)  \
#     $(TARGET7)  \
#     $(TARGET8)  \
#     $(TARGET9)  \
#     $(TARGET10) \
#     $(TARGET11) \
#     $(TARGET12) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
##$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
##	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#####
##$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
##	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#####
##$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
##	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
###
##$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
##	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET8): $(patsubst %,$(LIBDIR)/%,$(OBJS8))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET9): $(patsubst %,$(LIBDIR)/%,$(OBJS9))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET10): $(patsubst %,$(LIBDIR)/%,$(OBJS10))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET11): $(patsubst %,$(LIBDIR)/%,$(OBJS11))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)
#
#$(TARGET12): $(patsubst %,$(LIBDIR)/%,$(OBJS12))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

