.PHONY: run all clean

SHELL = /bin/bash
MODULES := main # cluster diagnostic options

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

HOST_FLAGS += -W -Wall -Wno-frame-address -save-temps -fverbose-asm -ffunction-sections -O3 -g0 -fopenmp $(FLAGS)

HOST_CXXFLAGS = -std=c++17 -fno-exceptions $(HOST_FLAGS) $(CXXFLAGS)

%.obj: %.cpp
%.obj: %.cpp $(EXTRA_DEP) $(DEPDIR)/%.d | $(DEPDIR)
	$(GXX) -c $< $(HOST_CXXFLAGS) -o $@ $(I_GMP) $(DEPFLAGS)

all: dorun

dorun: $(OBJS) $(EXTRA_DEP)
	$(GXX) $(HOST_CXXFLAGS) $(OBJS) -o $@ -Wl,--gc-sections -lm

run: dorun
	@echo "OMP_SCHEDULE=$(OMP_SCHEDULE)"
	time nice -20 ./$< $(ARGS)
#	convert out-1.ppm out-1.png
#	convert out-2.ppm out-2.png

clean:
	rm -f -- $(wildcard *.[iso] *.ii *.obj *.exe *.x *.x.* dorun)
	rm -f -- $(wildcard *.lst *.lss *.out *.map)
	rm -rf -- $(wildcard .deps)


demo.x: demo.cpp Makefile
	g++ -Wall -O3 -g2 -std=c++11 $< -o $@ -fopenmp

xxx: demo.x
	time ./demo.x $(ARGS)
