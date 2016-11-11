OBJS = SoG4Tubs.o SoHelicalTrack.o SoEllipsoid.o SoArrow.o
CFILES = $(OBJS:.so=.C)
CXX=/usr/local/bin/g++
CXXFLAGS = -fPIC -shared -O -c $(INCLUDES)

all: $(OBJS)

gl: all
clean:
	rm -f $(OBJS) *~
