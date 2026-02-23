.PHONY: run all clean

SHELL = /bin/bash
MODULES := main diagnostic have cubes cubes-border

###########################################################################
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/

SRCS := $(MODULES:%=%.cpp)
OBJS := $(SRCS:%.cpp=%.obj)

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
###############################################################################

GXX = g++-14

EXTRA_DEP = Makefile

# Several PRs exist like PR124207 for wrong bounds warnings.
WARN  = -W -Wall -Wno-stringop-overflow
DEBUG = -g0

HOST_FLAGS += $(WARN) -O3 $(DEBUG) -fopenmp $(FLAGS)

HOST_CXXFLAGS = -std=c++17 -fno-exceptions $(HOST_FLAGS) $(CXXFLAGS)

%.obj: %.cpp
%.obj: %.cpp $(EXTRA_DEP) $(DEPDIR)/%.d | $(DEPDIR)
	$(GXX) -c $< $(HOST_CXXFLAGS) -o $@ $(I_GMP) $(DEPFLAGS)

all: dorun

dorun: $(OBJS) $(EXTRA_DEP)
	$(GXX) $(HOST_CXXFLAGS) $(OBJS) -o $@ -Wl,--gc-sections -lm

NICE ?= 20

run: dorun
	time nice -$(NICE) ./$< $(ARGS)
#	convert out-1.ppm out-1.png
#	convert out-2.ppm out-2.png

.PHONY: progress rotor alloc

progress:
	g++ test-progress.cpp -o prog.x -Wall -O3 -fopenmp $(FLAGS)
	./prog.x

rotor:
	g++ test-rotation.cpp -o rotor.x -Wall -O3 -fopenmp $(FLAGS)
	./rotor.x

canon:
	$(GXX) test-canonical.cpp -o canon.x -Wall -O3 -fopenmp $(FLAGS)
	./canon.x

ALLOC_CPP = test-allocator.cpp diagnostic.cpp
alloc:
	$(GXX) $(ALLOC_CPP)  -o alloc.x -Wall -O3 -fopenmp $(FLAGS)
	./alloc.x

clean:
	rm -f -- $(wildcard *.[iso] *.ii *.obj *.exe *.x *.x.* dorun)
	rm -f -- $(wildcard *.lst *.lss *.out *.map)
	rm -rf -- $(wildcard .deps)


demo.x: demo.cpp Makefile
	g++ -Wall -O3 -g2 -std=c++11 $< -o $@ -fopenmp

xxx: demo.x
	time ./demo.x $(ARGS)
