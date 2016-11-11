# make sure, env variable LD_LIBRARY_PATH = /usr/local/lib

#ML=$(shell [ -e Makefile.locale ]  && echo 1)
#ifeq ($(ML),1)
#	include Makefile.locale
#endif

OBJS = Detector.o DetElement.o Hit.o Plane.o Track.o StateVector.o FullVector.o any2str.o getopt.o getopt1.o
GLOBJS = GLDetector.o GLDetElement.o GLHit.o GLTrack.o GLStateVector.o GLFullVector.o Plane.o any2str.o getopt.o getopt1.o

CFILES := $(OBJS:.o=.cc)

OBJS+= getopt.o getopt1.o
CFILES += getopt.c getopt1.c

ifneq ($(DEBUG),0)
	LDFLAGS += -pg -g
	CXXFLAGS += -pg -g
endif
CXXFLAGS += -DETCDIR=\"$(DESTDIR)/etc/\" $(SOLFLAG) $(INCLUDES) $(SETOPTIMIZE) -DCOLORS $(SETDEBUG) $(SETVERBOSE)
gl : override CXXFLAGS += -fPIC -I$(INVENTOR)include -I/usr/local/include -DHEPVIS # -DGL
gl : CFILES += SoArrow.cc SoEllipsoid.cc

GLDetector.o: Detector.cc
	$(CXX) $(CXXFLAGS) -c -o GLDetector.o Detector.cc

GLDetElement.o: DetElement.cc
	$(CXX) $(CXXFLAGS) -c -o GLDetElement.o DetElement.cc

GLHit.o: Hit.cc
	$(CXX) $(CXXFLAGS) -c -oGLHit.o Hit.cc

GLTrack.o: Track.cc
	$(CXX) $(CXXFLAGS) -c -o GLTrack.o Track.cc

GLStateVector.o: StateVector.cc
	$(CXX) $(CXXFLAGS) -c -o GLStateVector.o StateVector.cc

GLFullVector.o: FullVector.cc
	$(CXX) $(CXXFLAGS) -c -o GLFullVector.o FullVector.cc
	

all: nogl gl
nogl: $(OBJS)
gl: $(GLOBJS)
	
clean:
	cd hepvis && $(MAKE) $@
	rm -f $(OBJS) *~
mrproper: clean
	rm -f .\#* .*.swp core .ignore_starter.m
ctags:
	ctags --no-warn --c++ *.cc *.h
	echo "!_TAG_FILE_SORTED	0" > /tmp/tags
	cat tags |  sed -e 's/^\w*:://' >> /tmp/tags
	cp /tmp/tags .

FORCE:
