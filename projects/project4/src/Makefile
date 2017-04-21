ADVFLAGS = -std=c++11
OPTIM = -O3
LINK_TARGET = moldyn.exe
OBJS =                    \
  lib/atom.o              \
  lib/io.o                \
  lib/lennardjones.o      \
  lib/statisticssampler.o \
  lib/system.o            \
  lib/unitconverter.o     \
  lib/velocityverlet.o    \
  lib/main.o

REBUILDABLES = $(OBJS) $(LINK_TARGET)

all: $(LINK_TARGET)

$(LINK_TARGET) : $(OBJS) lib/vec3.o
	g++ -o $@ $^


lib/%.o: %.cpp
	g++ $(ADVFLAGS) -o $@ -c $<


lib/vec3.o: math/vec3.cpp
	g++ $(ADVFLAGS) -o lib/vec3.o -c math/vec3.cpp

clean:
	rm -f $(REBUILDABLES)
