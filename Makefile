# make sure, env variable LD_LIBRARY_PATH = /usr/local/lib

# yes. make -e 'DEBUG=0' turns off debugging.


export VERBOSE=5
export DEBUG=1
export OPTIMIZE=1


export INVENTOR=
export CLHEP=

export DESTDIR=.
export GL_LDFLAGS = -lXm -lXt -lXi -lXext -lX11 -lGL -lGLU -lm -ldl # -limage
export GL_LDFLAGS += -lInventor -lInventorXt

export VERSION=$(shell ./getversion)

export HEPVIS_BASE= ./hepvis
export HEPVIS_INCLUDE = -I$(HEPVIS_BASE)
export HEPVIS_SO=$(HEPVIS_BASE)/
export HEPVIS_O=$(HEPVIS_BASE)/

export INCLUDES = -I/usr/X11R6/include $(HEPVIS_INCLUDE) \
           -I$(CLHEP)/include\
		   -I../interface -I../hepvis/
export LIBS=-L$(CLHEP)/lib -lCLHEP\

export SETVERBOSE=-DVERBOSE=$(VERBOSE)
	
ifeq ($(DEBUG),0)
	SETDEBUG=
else 
	SETDEBUG=-DDEBUG
endif

ifeq ($(OPTIMIZE),0)
	SETOPTIMIZE=
else 
	SETOPTIMIZE=-O
endif
export SETDEBUG
export SETOPTIMIZE

export CXX=/usr/bin/g++
# CC=/usr/local/bin/gcc
# CXX=/opt/SUNWspro/bin/cc

export SYSOS=$(shell uname | tr A-Z a-z )

ifeq ($(SYSOS),sunos)
	export CXX=/opt/SUNWspro/bin/CC
else
	export CXX=/usr/bin/g++
endif

ifeq ($(SYSOS),sunos)
	export SOLFLAG = -DSOLARIS
endif

#ML=$(shell [ -e Makefile.locale ]  && echo 1)
#ifeq ($(ML),1)
#	include Makefile.locale
#endif

# do we have to pass -limage to the linker?
CHECK_IMAGE=$(shell [ -e $(INVENTOR)/lib/libimage.so ] && echo 1)
ifeq ($(CHECK_IMAGE),1)
	GL_LDFLAGS += -limage
endif
	

	
CXXFLAGS += -DETCDIR=\"$(DESTDIR)/etc/\" $(SOLFLAG) $(INCLUDES) $(SETOPTIMIZE) \
	-DCOLORS $(SETDEBUG) $(SETVERBOSE)
LDFLAGS += $(LIBS)

NN=$(shell ls /afs/hephy.at/scratch/w/wwalten/data/ | sort | tail -1)
NN2= $(shell echo $(NN)+1 | bc)
NC=$(shell echo $(NN2) | wc -c )
NC2= $(shell echo $(NC)-1 | bc)
NP=$(shell echo 00$(NN2) | cut -b $(NC2)-)

BINARY = traxx


GLOBJS = src/GLDetector.o src/GLDetElement.o src/GLHit.o src/Plane.o \
       src/GLStateVector.o src/GLFullVector.o src/any2str.o src/getopt.o \
	   src/getopt1.o src/GLTrack.o

GLOBJS+= test/gltraxx.o
OBJS = src/Detector.o src/DetElement.o src/Hit.o src/Plane.o src/Track.o \
       src/StateVector.o src/FullVector.o src/any2str.o src/getopt.o \
	   src/getopt1.o

OBJS+= test/traxx.o
	
ifneq ($(DEBUG),0)
	override LDFLAGS += -pg -g
	override CXXFLAGS += -pg -g
endif
gl : override CXXFLAGS += -fPIC -I$(INVENTOR)include -I/usr/local/include -DHEPVIS # -DGL
gl : override BINARY = gltraxx
gl : override LDFLAGS += -DHEPVIS -L. -L$(HEPVIS_BASE)/lib/ -lhepvis -fPIC -L/usr/X11R6/lib -L$(INVENTOR)lib $(GL_LDFLAGS)
gl : override GLOBJS += $(HEPVIS_SO)SoHelicalTrack.o $(HEPVIS_SO)SoG4Tubs.o $(HEPVIS_SO)SoArrow.o $(HEPVIS_SO)SoEllipsoid.o
gl : CFILES += SoArrow.cc SoEllipsoid.cc

nogl : override CXXFLAGS += -Wall
5test : override CXXFLAGS += -Wall -g

nogl-nodebug : override CXXFLAGS += -Wall

all: nogl gl
	
test2: FORCE
	@echo $(CHECK_IMAGE)

vim_nogl:
	@echo '----< make nogl >----'
	make nogl
nogl:
	cd src && $(MAKE) $@
	cd test && $(MAKE) $@
	$(CXX) $(LDFLAGS) -o $(BINARY) $(OBJS)
	@echo 
	@echo "----< Made non-OpenGL target. >----"

nogl-nodebug:
	$(CXX) $(LDFLAGS) -o $(BINARY) $(OBJS)
	@echo 
	@echo "----< Made non-OpenGL target with debugging symbols. >----"
	
gl:
	cd src && $(MAKE) $@
	cd test && $(MAKE) $@
	cd hepvis && $(MAKE) $@
	$(CXX) $(LDFLAGS) -o $(BINARY) $(GLOBJS)
	@echo 
	@echo "----< Made OpenGL target. >----"
#	@gl-restart
	
clean:
	rm -f $(OBJS) $(GLOBJS) *~ traxx gltraxx
mrproper: clean
	rm -f .\#* .*.swp core .ignore_starter.m
doc: FORCE
	doxygen
	cd latex; make ps
ctags:
	ctags --no-warn --c++ *.cc *.h
	echo "!_TAG_FILE_SORTED	0" > /tmp/tags
	cat tags |  sed -e 's/^\w*:://' >> /tmp/tags
	cp /tmp/tags .
install:
	mkdir -p $(DESTDIR)/usr/bin
	mkdir -p $(DESTDIR)/etc/traxx
	mkdir -p $(DESTDIR)/usr/man/man3
	test -f gltraxx && install gltraxx $(DESTDIR)/usr/bin
	test -f traxx && install traxx $(DESTDIR)/usr/bin
	test -d man && install man/man3/* $(DESTDIR)/usr/man/man3/
	test -f etc/detector.conf && cp etc/detector.conf $(DESTDIR)/etc/traxx

test:
	@echo $(DEF_MAKE)
	
getversion: FORCE
	g++ -o getversion test/getversion.cpp
stat_clean:
	rm -f *riemann? *global? *kalman?
tar: getversion mrproper
	export VERSION=$(shell ./getversion)
	tar --exclude=.gdb_history --exclude=.gdb_init --exclude=latex --exclude=man --exclude=tmp -czvf ../traxx-$(VERSION).tar.gz .
publish_tar: getversion mrproper
	export VERSION=$(shell ./getversion)
	tar --exclude=latex --exclude=man --exclude=tmp -czvf /afs/hephy.at/user/w/wwalten/www/traxx-$(VERSION).tar.gz .
	

reference: FORCE
	@traxx -a2 -m2 -f etc/forward-det.old.conf -o reference.forward
	@traxx -a2 -m2 -f etc/barrel-det.conf -o reference.barrel

check: check_barrel check_forward check2_barrel check2_forward

check_forward:
	@traxx -p -v0 -a2 -m2 -f etc/forward-det.old.conf -o new.for
	@diff new.for reference.forward
	@echo "Message: [Makefile] forward check  I done"
	@mv new.for /tmp/

check_barrel:
	@traxx -v0 -a2 -m2 -f etc/barrel-det.conf -o new.bar
	@diff new.bar reference.barrel
	@echo "Message: [Makefile] barrel  check  I done"
	@mv new.bar /tmp/

check2_barrel:
	@traxx -v0 -a2 -m2 -f etc/barrel-det.conf -R 0 -o seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -G 0 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -K 0 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -R 1 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -G 1 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -K 1 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -R 2 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -G 2 -O | grep ^SV >> seg.bar
	@traxx -a2 -m2 -f etc/barrel-det.conf -G 3 -O | grep ^SV >> seg.bar
	@diff seg.bar reference.barrel
	@echo "Message: [Makefile] barrel  check II done"
	@mv seg.bar /tmp/

check2_forward:
	@traxx  -p -v0 -a2 -m2 -f etc/forward-det.old.conf -R 0 -o seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -G 0 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -K 0 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -R 1 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -G 1 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -K 1 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -R 2 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -G 2 -O | grep ^SV >> seg.for
	@traxx  -p -a2 -m2 -f etc/forward-det.old.conf -G 3 -O | grep ^SV >> seg.for
	@diff seg.for reference.forward
	@echo "Message: [Makefile] forward check II done"
	@mv seg.for /tmp/
	
run:
	run.pl

publish:
	mkdir $(SCRATCH)/data/$(NP)
	cp data/* $(SCRATCH)/data/$(NP)

prop_scram:
	cp *.cc ~/OTraxx/src/
	cp *.h  ~/OTraxx/interface/
	cp *.cpp ~/OTraxx/test/
gprof_rieman1:
	cp profile profile.old
	./traxx -p -n 5000 -R1 -f etc/barrel.conf -T
	gprof traxx > profile
	
gprof_kalman1:
	cp profile.k1 profile.k1.old
	./traxx -p -n 5000 -K1 -f etc/barrel.conf -T
	gprof traxx > profile.k1

killrun:
	killall run.pl
	killall traxx

makefile_template:
	cat Makefile | grep -v export\ CELL | sed -e 's/export INVENTOR=
	
	
distclean:	clean check_clean
	rm -f data/*
	rm -f core gmon.out profile*
	rm -f tmp/[^CVS]*

check_clean:
	rm -f new.bar new.for seg.bar seg.for
FORCE:
